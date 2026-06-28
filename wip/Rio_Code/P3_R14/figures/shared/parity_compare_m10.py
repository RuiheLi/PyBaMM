"""Compare port-side M10 (`SEI = "interstitial-diffusion limited (legacy)"`)
against source M9b (which uses the v23 interstitial-diffusion-limited form
natively). Both use ``solvent diffusion = full with sei refill`` and run the
full P3_R14 4-stage protocol with ``SEI porosity change = false``.

Reads:
  - parity_source_m9b.npz (in this directory)
  - parity_port_m10.npz   (in PyBaMM-port-work/wip/)
"""
import os
import sys

import numpy as np


SRC = os.path.join(os.path.dirname(__file__), "parity_source_m9b.npz")
PORT = r"D:/LRHWork/Model_RH/PyBaMM-port-work/wip/parity_port_m10.npz"


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
    print(f"  source M9b: shape={s['t'].shape}, [{float(s['t'][0]):.2f}, "
          f"{float(s['t'][-1]):.2f}] s")
    print(f"  port M10:   shape={q['t'].shape}, [{float(q['t'][0]):.2f}, "
          f"{float(q['t'][-1]):.2f}] s")

    print("== Voltage [V] ==")
    print(f"  source: V[0]={float(s['V'][0]):.4f}, "
          f"V_min={float(s['V'].min()):.4f}, V_max={float(s['V'].max()):.4f}")
    print(f"  port:   V[0]={float(q['V'][0]):.4f}, "
          f"V_min={float(q['V'].min()):.4f}, V_max={float(q['V'].max()):.4f}")

    print("== X-averaged c_e [mol/m3] ==")
    print(f"  source: [{float(s['c_e_avg'].min()):.2f}, "
          f"{float(s['c_e_avg'].max()):.2f}]")
    print(f"  port:   [{float(q['c_e_avg'].min()):.2f}, "
          f"{float(q['c_e_avg'].max()):.2f}]")

    print("== X-averaged c_EC [mol/m3] ==")
    print(f"  source: [{float(s['c_EC_avg'].min()):.2f}, "
          f"{float(s['c_EC_avg'].max()):.2f}]")
    print(f"  port:   [{float(q['c_EC_avg'].min()):.2f}, "
          f"{float(q['c_EC_avg'].max()):.2f}]")

    print("== Loss of Li to SEI [mol] ==")
    Q_s = s["Q_sei"]
    Q_q = q["Q_sei"]
    print(f"  source M9b: Q_sei[-1]={float(Q_s[-1]):.4e}")
    print(f"  port M10:   Q_sei[-1]={float(Q_q[-1]):.4e}")
    if float(Q_s[-1]) != 0:
        print(f"  ratio source/port = {float(Q_s[-1]) / float(Q_q[-1]):.2f}")

    print(
        "\nNote: M10 swaps the j_sei formula to the v23-faithful (legacy) "
        "form. As documented in ECdrag2_migration_status.md, this does not "
        "close the residual Q_sei gap because in this regime "
        "L_sei_inner == L_sei_outer and c_EC/c_EC_init ~ 1, so the legacy "
        "and v24-default formulas evaluate to the same number. The remaining "
        "factor of ~38x lives in the L_sei RHS form and the experiment-step "
        "interpretation (v23 vs v24)."
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
