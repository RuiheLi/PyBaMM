"""Compare parity_source_m9b.npz and parity_port_m9b.npz."""
import os
import sys

import numpy as np


SRC = os.path.join(os.path.dirname(__file__), "parity_source_m9b.npz")
PORT = r"D:/LRHWork/Model_RH/PyBaMM-port-work/wip/parity_port_m9b.npz"


def main():
    if not os.path.exists(SRC):
        print(f"missing {SRC}")
        return 1
    if not os.path.exists(PORT):
        print(f"missing {PORT}")
        return 1
    s = np.load(SRC)
    q = np.load(PORT)

    print("== Time / Voltage envelope ==")
    print(f"  source: t_end={float(s['t'][-1]):.2f}s, "
          f"V_min={float(s['V'].min()):.4f}, V_max={float(s['V'].max()):.4f}, "
          f"V[0]={float(s['V'][0]):.4f}")
    print(f"  port:   t_end={float(q['t'][-1]):.2f}s, "
          f"V_min={float(q['V'].min()):.4f}, V_max={float(q['V'].max()):.4f}, "
          f"V[0]={float(q['V'][0]):.4f}")

    print("== X-averaged c_e [mol/m3] ==")
    print(f"  source: [{float(s['c_e_avg'].min()):.2f}, "
          f"{float(s['c_e_avg'].max()):.2f}], range="
          f"{float(s['c_e_avg'].max()-s['c_e_avg'].min()):.2f}")
    print(f"  port:   [{float(q['c_e_avg'].min()):.2f}, "
          f"{float(q['c_e_avg'].max()):.2f}], range="
          f"{float(q['c_e_avg'].max()-q['c_e_avg'].min()):.2f}")
    print(f"  delta min={float(q['c_e_avg'].min()-s['c_e_avg'].min()):+.2f}, "
          f"delta max={float(q['c_e_avg'].max()-s['c_e_avg'].max()):+.2f}")

    print("== X-averaged c_EC [mol/m3] ==")
    print(f"  source: [{float(s['c_EC_avg'].min()):.2f}, "
          f"{float(s['c_EC_avg'].max()):.2f}], range="
          f"{float(s['c_EC_avg'].max()-s['c_EC_avg'].min()):.2f}")
    print(f"  port:   [{float(q['c_EC_avg'].min()):.2f}, "
          f"{float(q['c_EC_avg'].max()):.2f}], range="
          f"{float(q['c_EC_avg'].max()-q['c_EC_avg'].min()):.2f}")
    print(f"  delta min={float(q['c_EC_avg'].min()-s['c_EC_avg'].min()):+.2f}, "
          f"delta max={float(q['c_EC_avg'].max()-s['c_EC_avg'].max()):+.2f}")

    print("== SEI lithium loss [mol] ==")
    print(f"  source: Q_sei[-1]={float(s['Q_sei'][-1]):.4e}")
    print(f"  port:   Q_sei[-1]={float(q['Q_sei'][-1]):.4e}")
    print(f"  ratio (source/port): {float(s['Q_sei'][-1]/q['Q_sei'][-1]):.2f}x")
    return 0


if __name__ == "__main__":
    sys.exit(main())
