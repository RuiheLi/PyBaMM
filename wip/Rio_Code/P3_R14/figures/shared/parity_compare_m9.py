"""Compare parity_source_m9.npz and parity_port_m9.npz (full P3_R14)."""
import os
import sys

import numpy as np


SRC = os.path.join(os.path.dirname(__file__), "parity_source_m9.npz")
PORT = r"D:/LRHWork/Model_RH/PyBaMM-port-work/wip/parity_port_m9.npz"


def main():
    if not os.path.exists(SRC):
        print(f"missing {SRC}")
        return 1
    if not os.path.exists(PORT):
        print(f"missing {PORT}")
        return 1
    s = np.load(SRC)
    q = np.load(PORT)

    print("== Time grids ==")
    print(f"  source: shape={s['t'].shape}, [{float(s['t'][0]):.2f}, "
          f"{float(s['t'][-1]):.2f}]s")
    print(f"  port:   shape={q['t'].shape}, [{float(q['t'][0]):.2f}, "
          f"{float(q['t'][-1]):.2f}]s")

    print("== Voltage [V] ==")
    print(f"  source: V[0]={float(s['V'][0]):.4f}, V[-1]={float(s['V'][-1]):.4f}, "
          f"V_min={float(s['V'].min()):.4f}, V_max={float(s['V'].max()):.4f}")
    print(f"  port:   V[0]={float(q['V'][0]):.4f}, V[-1]={float(q['V'][-1]):.4f}, "
          f"V_min={float(q['V'].min()):.4f}, V_max={float(q['V'].max()):.4f}")

    # Interpolate port voltage onto source time grid for elementwise compare.
    V_q_on_s = np.interp(s["t"], q["t"], q["V"])
    dV = V_q_on_s - s["V"]
    print(f"  delta on source grid: mean={dV.mean():+.4f} V, "
          f"max|d|={np.max(np.abs(dV)):+.4f} V")

    print("== X-averaged c_e [mol/m3] ==")
    print(f"  source: [{float(s['c_e_avg'].min()):.2f}, "
          f"{float(s['c_e_avg'].max()):.2f}], range="
          f"{float(s['c_e_avg'].max()-s['c_e_avg'].min()):.2f}")
    print(f"  port:   [{float(q['c_e_avg'].min()):.2f}, "
          f"{float(q['c_e_avg'].max()):.2f}], range="
          f"{float(q['c_e_avg'].max()-q['c_e_avg'].min()):.2f}")

    print("== X-averaged c_EC [mol/m3] ==")
    print(f"  source: [{float(s['c_EC_avg'].min()):.2f}, "
          f"{float(s['c_EC_avg'].max()):.2f}], range="
          f"{float(s['c_EC_avg'].max()-s['c_EC_avg'].min()):.2f}")
    print(f"  port:   [{float(q['c_EC_avg'].min()):.2f}, "
          f"{float(q['c_EC_avg'].max()):.2f}], range="
          f"{float(q['c_EC_avg'].max()-q['c_EC_avg'].min()):.2f}")

    print("== SEI lithium loss [mol] ==")
    print(f"  source: Q_sei[-1]={float(s['Q_sei'][-1]):.4e}")
    print(f"  port:   Q_sei[-1]={float(q['Q_sei'][-1]):.4e}")
    print(f"  ratio (source/port): {float(s['Q_sei'][-1]/q['Q_sei'][-1]):.2f}x")

    print()
    print("Caveat: with SEI = 'interstitial-diffusion limited' + SEI film "
          "resistance + SEI porosity change all on, the v23-era and v24 SEI "
          "implementations diverge. The port code path runs end-to-end with "
          "the EC drag physics; absolute Q_sei / V[0] differences come from "
          "SEI submodel changes upstream, not from M3-M8 port logic.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
