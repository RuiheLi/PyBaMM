"""Apply one-off fixes to three P3_R14 notebooks (category-C bugs)."""
import json
import re

ROOT = r"D:/LRHWork/Model_RH/PyBaMM-ECdrag2/wip/Rio_Code/P3_R14"


def patch_rate_performance_get_sol(nb):
    """Cell 11: align Rate_Dis_All with Cap_* / Vol_* lengths before plotting."""
    old = """for i in range(len(Rate_Dis_All)):
    Cap_HD.append(Get_cap(Sol_HD[i]))
    Cap_SD.append(Get_cap(Sol_SD[i]))
    Cap_LD.append(Get_cap(Sol_LD[i]))
    Vol_HD.append(Get_vol(Sol_HD[i]))
    Vol_SD.append(Get_vol(Sol_SD[i]))
    Vol_LD.append(Get_vol(Sol_LD[i]))
fig, axs = plt.subplots(2,1, figsize=(30/2.53,12/2.54),tight_layout=True)
axs[0].plot(Rate_Dis_All,Cap_HD,"-",marker = 'o',label="HD")
axs[0].plot(Rate_Dis_All,Cap_LD,"--",marker = '*',label="LD")
axs[0].plot(Rate_Dis_All,Cap_SD,"-.",marker = '^',label="SD")
axs[1].plot(Rate_Dis_All,Vol_HD,"-",marker = 'o',label="HD")
axs[1].plot(Rate_Dis_All,Vol_LD,"--",marker = '*',label="LD")
axs[1].plot(Rate_Dis_All,Vol_SD,"-.",marker = '^',label="SD")"""
    new = """for i in range(len(Rate_Dis_All)):
    Cap_HD.append(Get_cap(Sol_HD[i]))
    Cap_SD.append(Get_cap(Sol_SD[i]))
    Cap_LD.append(Get_cap(Sol_LD[i]))
    Vol_HD.append(Get_vol(Sol_HD[i]))
    Vol_SD.append(Get_vol(Sol_SD[i]))
    Vol_LD.append(Get_vol(Sol_LD[i]))
# Guard: partial sweeps can leave len(Cap_*) < len(Rate_Dis_All) and break plt.plot
_rates = np.asarray(Rate_Dis_All, dtype=float)
def _nplt(r, y):
    y = np.asarray(y, dtype=float)
    n = min(len(r), len(y))
    return r[:n], y[:n]
fig, axs = plt.subplots(2,1, figsize=(30/2.53,12/2.54),tight_layout=True)
r0, y0 = _nplt(_rates, Cap_HD); axs[0].plot(r0, y0,"-",marker = 'o',label="HD")
r0, y0 = _nplt(_rates, Cap_LD); axs[0].plot(r0, y0,"--",marker = '*',label="LD")
r0, y0 = _nplt(_rates, Cap_SD); axs[0].plot(r0, y0,"-.",marker = '^',label="SD")
r0, y0 = _nplt(_rates, Vol_HD); axs[1].plot(r0, y0,"-",marker = 'o',label="HD")
r0, y0 = _nplt(_rates, Vol_LD); axs[1].plot(r0, y0,"--",marker = '*',label="LD")
r0, y0 = _nplt(_rates, Vol_SD); axs[1].plot(r0, y0,"-.",marker = '^',label="SD")"""
    for c in nb["cells"]:
        if c.get("cell_type") != "code":
            continue
        src = "".join(c.get("source", []))
        if old not in src:
            continue
        c["source"] = [l + "\n" for l in src.replace(old, new).split("\n")]
        if c["source"] and c["source"][-1] == "\n":
            c["source"][-1] = ""
        return True
    return False


def patch_debug_xitilde(nb):
    old = """Save_Fig = True; Crate_index = -1
# SD/LDx/HDx datasets were not generated in this debug notebook (it only
# builds DD_Xitilde_Fun_Crate / DD_Xitilde_3_Crate). To keep the plotting
# code self-contained, fall back to the Xitilde-function dataset for all
# three slots when the legacy names aren't defined.
try:
    SD_Crate, DD_LDx_Crate, DD_HDx_Crate
except NameError:
    SD_Crate = DD_LDx_Crate = DD_HDx_Crate = DD_Xitilde_Fun_Crate
Plot_Concentration_1_Crate(SD_Crate,DD_LDx_Crate,DD_HDx_Crate,"""
    new = """Save_Fig = True; Crate_index = -1
# This notebook only builds DD_Xitilde_Fun_Crate / DD_Xitilde_3_Crate.
# Map them onto the plotting API (SD / LDx / HDx) — avoids fragile
# ``try: SD_Crate, ...`` tuple probes that behave badly in batch runs.
SD_Crate = DD_Xitilde_Fun_Crate
DD_LDx_Crate = DD_Xitilde_Fun_Crate
DD_HDx_Crate = DD_Xitilde_3_Crate
Plot_Concentration_1_Crate(SD_Crate,DD_LDx_Crate,DD_HDx_Crate,"""
    for c in nb["cells"]:
        if c.get("cell_type") != "code":
            continue
        src = "".join(c.get("source", []))
        if old not in src:
            continue
        c["source"] = [l + "\n" for l in src.replace(old, new).split("\n")]
        if c["source"] and c["source"][-1] == "\n":
            c["source"][-1] = ""
        return True
    return False


def clear_error_outputs_cell_containing(nb, needle: str):
    for c in nb["cells"]:
        if c.get("cell_type") != "code":
            continue
        src_needle = needle in "".join(c.get("source", []))
        if not src_needle:
            continue
        outs = []
        for o in c.get("outputs", []):
            if o.get("output_type") == "error":
                continue
            outs.append(o)
        c["outputs"] = outs
        c["execution_count"] = None
        return True
    return False


def main():
    p1 = f"{ROOT}/Rate_Performance_Get_Sol.ipynb"
    with open(p1, encoding="utf-8") as f:
        nb = json.load(f)
    assert patch_rate_performance_get_sol(nb), "Rate_Performance_Get_Sol patch failed"
    with open(p1, "w", encoding="utf-8", newline="\n") as f:
        json.dump(nb, f, ensure_ascii=False, indent=2)
        f.write("\n")
    print("patched", p1)

    p2 = f"{ROOT}/Debug_Xitilde.ipynb"
    with open(p2, encoding="utf-8") as f:
        nb = json.load(f)
    assert patch_debug_xitilde(nb), "Debug_Xitilde patch failed"
    with open(p2, "w", encoding="utf-8", newline="\n") as f:
        json.dump(nb, f, ensure_ascii=False, indent=2)
        f.write("\n")
    print("patched", p2)

    p3 = f"{ROOT}/1EC1DMC_Rate_Performance.ipynb"
    with open(p3, encoding="utf-8") as f:
        nb = json.load(f)
    # Source already has Return_sol/Save_sol; drop stale error blobs only.
    hit = clear_error_outputs_cell_containing(nb, "DD_LDx_Crate = Scan_Crate_Paper")
    if hit:
        with open(p3, "w", encoding="utf-8", newline="\n") as f:
            json.dump(nb, f, ensure_ascii=False, indent=2)
            f.write("\n")
        print("cleared stale error output in", p3)
    else:
        print("no stale error block to clear in", p3)


if __name__ == "__main__":
    main()
