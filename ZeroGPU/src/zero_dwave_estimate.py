#!/usr/bin/env python3
"""Zb073: finite-level BCS d-wave pairing strength on a square-lattice torus.

Estimator (linear response, Zb063 eq. 0.2):
    lam_d(mu) = (1/Nk) * sum_k g_k^2 / (2 |eps_k - mu|),   g_k = cos kx - cos ky
    V_c = 1 / lam_d        (Thouless criterion for the d-wave channel)

Compared against the OSC validation 4x4 Hubbard ED table (n, P_d) --
see _staging/OSC_进展.md and the OSC dataset validation/ directory.
"""
import json, math, sys

def spectrum(N, t=1.0):
    ks = [2*math.pi*i/N for i in range(N)]
    return [(-2*t*(math.cos(kx)+math.cos(ky)), kx, ky) for kx in ks for ky in ks]

def mu_for_filling(eigs, n):
    # occupation (2/Nk) sum f(eps-mu) = n ; T -> 0 step
    lo, hi = min(eigs)-1.0, max(eigs)+1.0
    for _ in range(200):
        mid = 0.5*(lo+hi)
        occ = 2.0*sum(1 for e in eigs if e < mid)/len(eigs)
        if occ < n: lo = mid
        else: hi = mid
    return 0.5*(lo+hi)

def lam_d(N, n, t=1.0, eta=1e-3):
    sp = spectrum(N, t)
    eigs = [e for e,_,_ in sp]
    mu = mu_for_filling(eigs, n)
    s = 0.0
    for e,kx,ky in sp:
        g = math.cos(kx)-math.cos(ky)
        s += g*g/(2.0*max(abs(e-mu), eta))
    return s/len(sp), mu

if __name__ == "__main__":
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 4
    ED = [(0.5,0.048076),(0.625,0.060122),(0.75,0.009831),(0.875,0.019219),(1.0,0.006197)]
    rows = []
    for n, pd in ED:
        lam, mu = lam_d(N, n)
        rows.append(dict(n=n, lam_d=lam, Vc=1.0/lam, mu=mu, ed_Pd=pd))
    # Spearman
    def rank(xs):
        order = sorted(range(len(xs)), key=lambda i: xs[i])
        r = [0]*len(xs)
        for pos,i in enumerate(order): r[i] = pos+1
        return r
    a = rank([r["lam_d"] for r in rows]); b = rank([r["ed_Pd"] for r in rows])
    d2 = sum((x-y)**2 for x,y in zip(a,b)); nn = len(rows)
    rho = 1 - 6*d2/(nn*(nn*nn-1))
    out = dict(cluster="%dx%d"%(N,N), rows=rows, spearman_lam_vs_ED=rho)
    json.dump(out, open("/Users/oygb/Downloads/Zero/ZeroGPU/results/zb073_dwave_estimate_N%d.json"%N, "w"), indent=1)
    print("cluster N=%d  spearman(lam_d, ED P_d) = %.3f" % (N, rho))
    print("%6s %10s %10s %10s %10s" % ("n","mu","lam_d","V_c","ED P_d"))
    for r in rows:
        print("%6.3f %10.4f %10.4f %10.4f %10.5f" % (r["n"], r["mu"], r["lam_d"], r["Vc"], r["ed_Pd"]))
