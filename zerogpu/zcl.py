"""
zcl.py --- Zero 计算引擎的 OpenCL 后端（AMD RX 570 / 任意 OpenCL 设备）。

设计约束（本机实测，见 bench/）：
  * 本机是 Intel Mac + AMD RX 570(8GB, 32CU)。**没有 CUDA，没有 MPS，没有 MLX。**
  * 朴素 OpenCL matmul 只有 6.5 GFLOP/s，而 NumPy BLAS 有 67 GFLOP/s
      => 稠密线性代数**必须**交给 CPU/numpy，GPU 只在特定形态上赢。
  * GPU 真正赢的形态（实测）：
      - 位并行字操作 (uint64 批量)      : 6.0x
      - 大量独立智能体的分支整数步进     : 33x
  * 所以本模块只提供 GPU 擅长的原语；稠密谱分解等仍走 numpy/scipy。

用法:
    from zcl import Engine
    e = Engine()                      # 自动选 GPU，失败回落 CPU
    buf = e.to_device(arr)
    e.run('canon_rot', n, buf_in, buf_out, np.int32(T), np.uint64(mask))
    out = e.from_device(buf_out, n, np.uint64)
"""
from __future__ import annotations

import os
import time
import numpy as np

try:
    import pyopencl as cl
    import pyopencl.array as cl_array  # noqa: F401
    HAVE_CL = True
except Exception:  # pragma: no cover
    cl = None
    HAVE_CL = False


# 标量一律用 numpy 标量传入（pyopencl 直接支持 np.int32/np.uint64/...），
# 不需要也不应依赖 pyopencl.cltypes 里的 C 类型常量。


class Engine:
    """一个 OpenCL 上下文 + 内核缓存 + 计时器。"""

    def __init__(self, prefer: str = "gpu", platform_idx: int = 0, device_idx: int | None = None,
                 verbose: bool = False):
        if not HAVE_CL:
            raise RuntimeError("pyopencl 不可用：请用系统 python3 (它装了 pyopencl)")
        self.verbose = verbose
        P = cl.get_platforms()[platform_idx]
        devs = P.get_devices()
        if device_idx is not None:
            dev = devs[device_idx]
        else:
            dev = None
            if prefer == "gpu":
                for d in devs:
                    if d.type & cl.device_type.GPU:
                        dev = d
                        break
            if dev is None:
                dev = devs[0]
        self.device = dev
        self.ctx = cl.Context([dev])
        self.queue = cl.CommandQueue(self.ctx)
        self._progs: dict[str, object] = {}
        self._kernels: dict[tuple[str, str], object] = {}
        self._build_opts: str | None = None
        # 探测 fp64
        self.has_fp64 = "cl_khr_fp64" in (dev.extensions or "")
        # 注意：不要传 "-cl-khr-fp64"——那不是一个合法的构建选项，AMD 编译器会直接报错。
        self._build_opts = "-cl-mad-enable"

    # ------------------------------------------------------------- 设备信息
    @property
    def name(self) -> str:
        return f"{self.device.name.strip()}"

    def info(self) -> dict:
        d = self.device
        return {
            "name": d.name.strip(),
            "type": cl.device_type.to_string(d.type),
            "cus": d.max_compute_units,
            "global_mem_gb": round(d.global_mem_size / 2**30, 2),
            "max_alloc_gb": round(d.max_mem_alloc_size / 2**30, 2),
            "clock_mhz": d.max_clock_frequency,
            "fp64": self.has_fp64,
            "max_wg": d.max_work_group_size,
        }

    # --------------------------------------------------------------- 内核
    def program(self, source: str, key: str | None = None):
        key = key or source
        if key not in self._progs:
            t0 = time.time()
            self._progs[key] = cl.Program(self.ctx, source).build(options=self._build_opts)
            if self.verbose:
                print(f"[zcl] built program '{key[:24]}...' in {time.time()-t0:.2f}s")
        return self._progs[key]

    def kernel(self, key: str, source: str, name: str):
        """取（并缓存）内核对象——避免 RepeatedKernelRetrieval 的开销。"""
        ck = (key, name)
        if ck not in self._kernels:
            prg = self.program(source, key)
            self._kernels[ck] = cl.Kernel(prg, name)
        return self._kernels[ck]

    # -------------------------------------------------------------- 缓冲
    def buffer(self, nbytes: int, readonly: bool = False):
        mf = cl.mem_flags
        flags = mf.READ_ONLY if readonly else (mf.READ_WRITE | mf.ALLOC_HOST_PTR)
        return cl.Buffer(self.ctx, flags, int(nbytes))

    def to_device(self, arr: np.ndarray, readonly: bool = True):
        arr = np.ascontiguousarray(arr)
        mf = cl.mem_flags
        flags = (mf.READ_ONLY if readonly else mf.READ_WRITE) | mf.COPY_HOST_PTR
        return cl.Buffer(self.ctx, flags, hostbuf=arr)

    def empty(self, n: int, dtype) -> "cl.Buffer":
        return cl.Buffer(self.ctx, cl.mem_flags.READ_WRITE, int(n) * np.dtype(dtype).itemsize)

    def from_device(self, buf, n: int, dtype) -> np.ndarray:
        out = np.empty(int(n), dtype=dtype)
        cl.enqueue_copy(self.queue, out, buf).wait()
        return out

    # ---------------------------------------------------------------- 执行
    def run(self, kernel, global_size, *args, local_size=None, wait=True,
            arg_dtypes=None):
        """跑一个内核。args 里的 numpy 标量会自动转成合适的 CL 类型。"""
        k = kernel
        if isinstance(kernel, tuple):
            k = self.kernel(*kernel)
        conv = []
        for i, a in enumerate(args):
            if isinstance(a, np.generic):
                conv.append(a)
            elif arg_dtypes and i < len(arg_dtypes) and arg_dtypes[i] is not None:
                conv.append(arg_dtypes[i](a))
            else:
                conv.append(a)
        ev = k(self.queue, global_size, local_size, *conv)
        if wait:
            ev.wait()
        return ev

    def finish(self):
        self.queue.finish()

    def timeit(self, fn, repeat: int = 5, warmup: int = 2):
        for _ in range(warmup):
            fn()
        self.finish()
        t0 = time.time()
        for _ in range(repeat):
            fn()
        self.finish()
        return (time.time() - t0) / repeat


# ------------------------------------------------------------------ 公共内核
COMMON = r"""
// 循环左移一个 T 位字（位 0..T-1 有效），并取最小值扫描
__kernel void canon_rot(const __global ulong* in, __global ulong* out,
                        const int T, const ulong mask) {
    size_t i = get_global_id(0);
    ulong w = in[i] & mask;
    ulong best = w;
    ulong r = w;
    const int s = T - 1;
    for (int j = 1; j < T; ++j) {
        r = ((r << 1) | (r >> s)) & mask;
        best = min(best, r);
    }
    out[i] = best;
}

// 循环相邻对换：把位置 i 与 (i+1)%T 的两位交换（仅当它们不同时才有变化）
__kernel void swap_at(const __global ulong* in, __global ulong* out,
                      const int T, const ulong mask, const int i) {
    size_t g = get_global_id(0);
    ulong w = in[g] & mask;
    int j = (i + 1) % T;
    ulong bi = (w >> i) & 1UL;
    ulong bj = (w >> j) & 1UL;
    if (bi != bj) {
        w ^= (1UL << i) | (1UL << j);
    }
    out[g] = w;
}

// 对每个字，把所有 T 个循环相邻对换的结果写出来（用于建图）
__kernel void all_swaps(const __global ulong* in, __global ulong* out,
                        const int T, const ulong mask) {
    size_t g = get_global_id(0);
    ulong w = in[g] & mask;
    for (int i = 0; i < T; ++i) {
        int j = (i + 1) % T;
        ulong bi = (w >> i) & 1UL;
        ulong bj = (w >> j) & 1UL;
        ulong r = w;
        if (bi != bj) r ^= (1UL << i) | (1UL << j);
        out[g * T + i] = r;
    }
}

// popcount（零和检查用）
__kernel void popcnt_k(const __global ulong* in, __global int* out, const ulong mask) {
    size_t i = get_global_id(0);
    out[i] = (int)popcount(in[i] & mask);
}
"""
