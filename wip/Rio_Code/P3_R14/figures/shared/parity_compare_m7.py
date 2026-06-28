"""Compare parity_source_m7.npz and parity_port_m7.npz."""
import os
import sys

import numpy as np


SRC = os.path.join(os.path.dirname(__file__), "parity_source_m7.npz")
PORT = r"D:/LRHWork/Model_RH/PyBaMM-port-work/wip/parity_port_m7.npz"


def main():
    if not os.path.exists(SRC):
        print(f"missing {SRC}")
        return 1
    if not os.path.exists(PORT):
        print(f"missing {PORT}")
        return 1
    s = np.load(SRC)
    q = np.load(PORT)
    t_s, t_q = s["t"], q["t"]
    V_s, V_q = s["V"], q["V"]
    cE_s, cE_q = s["c_e_avg"], q["c_e_avg"]
    cEC_s, cEC_q = s["c_EC_avg"], q["c_EC_avg"]

    n = min(len(t_s), len(t_q))

    print("== Time grids ==")
    print(f"  source: shape={t_s.shape}, [{t_s[0]:.2f}, {t_s[-1]:.2f}]s")
    print(f"  port:   shape={t_q.shape}, [{t_q[0]:.2f}, {t_q[-1]:.2f}]s")

    print("== Voltage [V] ==")
    print(f"  source: V[0]={V_s[0]:.4f}, V[-1]={V_s[-1]:.4f}, drop={V_s[0]-V_s[-1]:+.4f}")
    print(f"  port:   V[0]={V_q[0]:.4f}, V[-1]={V_q[-1]:.4f}, drop={V_q[0]-V_q[-1]:+.4f}")
    dV = V_q[:n] - V_s[:n]
    print(f"  delta:  mean={dV.mean():+.4f}, max|d|={np.max(np.abs(dV)):+.4f}")

    print("== X-averaged c_e [mol/m3] ==")
    print(f"  source: c_e[0]={cE_s[0]:.2f}, c_e[-1]={cE_s[-1]:.2f}")
    print(f"  port:   c_e[0]={cE_q[0]:.2f}, c_e[-1]={cE_q[-1]:.2f}")
    dCe = cE_q[:n] - cE_s[:n]
    print(f"  delta:  mean={dCe.mean():+.3f}, max|d|={np.max(np.abs(dCe)):+.3f}")

    print("== X-averaged c_EC [mol/m3] ==")
    print(f"  source: c_EC[0]={cEC_s[0]:.4f}, c_EC[-1]={cEC_s[-1]:.4f}, "
          f"drop={cEC_s[0]-cEC_s[-1]:+.4f}")
    print(f"  port:   c_EC[0]={cEC_q[0]:.4f}, c_EC[-1]={cEC_q[-1]:.4f}, "
          f"drop={cEC_q[0]-cEC_q[-1]:+.4f}")
    dCEC = cEC_q[:n] - cEC_s[:n]
    print(f"  delta:  mean={dCEC.mean():+.4f}, max|d|={np.max(np.abs(dCEC)):+.4f}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
