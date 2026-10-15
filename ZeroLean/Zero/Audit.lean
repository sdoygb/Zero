/-
  ZeroLean 自动审计。

  遍历环境中所有 `Zero.*` 定理，逐条打印其依赖公理（`LEDGER|定理|公理1,公理2`）。
  CI（scripts/check.py）解析这些行 —— 于是**新增定理自动纳入审计**，作者不必记得写
  `#print axioms`（早期版本靠手写，结果漏掉了故意植入的违规样例）。
-/

import Zero

open Lean

run_cmd do
  let env ← getEnv
  let mut count := 0
  for (name, info) in env.constants.toList do
    if (`Zero).isPrefixOf name then
      match info with
      | .thmInfo _ =>
        let axs ← collectAxioms name
        let s := String.intercalate "," (axs.toList.map (fun a => a.toString))
        IO.println s!"LEDGER|{name}|{s}"
        count := count + 1
      | _ => pure ()
  IO.println s!"LEDGER-COUNT|{count}"
