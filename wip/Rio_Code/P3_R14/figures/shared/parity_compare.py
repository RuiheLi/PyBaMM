"""Compare parity_source.npz and parity_port.npz."""
import os
import sys

import numpy as np


SRC = os.path.join(os.path.dirname(__file__), "parity_source.npz")
PORT = r"D:/LRHWork/Model_RH/PyBaMM-port-work/wip/parity_port.npz"


def main():
    if not os.path.exists(SRC):
        print(f"missing {SRC}")
        return 1
    if not os.path.exists(PORT):
        print(f"missing {PORT}")
        return 1
    s = np.load(SRC)
    q = np.load(PORT)
    t_s = s["t"]
    t_q = q["t"]
    V_s = s["V"]
    V_q = q["V"]
    cE_s = s["c_e_avg"]
    cE_q = q["c_e_avg"]
    cEC_s = s["c_EC_avg"]
    cEC_q = q["c_EC_avg"]

    print("== Time grids ==")
    print(f"  source: shape={t_s.shape}, [{t_s[0]:.2f}, {t_s[-1]:.2f}]s")
    print(f"  port:   shape={t_q.shape}, [{t_q[0]:.2f}, {t_q[-1]:.2f}]s")

    n = min(len(t_s), len(t_q))

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
    print(f"  source: c_EC[0]={cEC_s[0]:.2f}, c_EC[-1]={cEC_s[-1]:.2f}")
    print(f"  port:   c_EC[0]={cEC_q[0]:.2f}, c_EC[-1]={cEC_q[-1]:.2f}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
