"""Patch Debug_Xitilde.ipynb: replace fragile FindClose with _nearest_time_index."""
import json

P = r"D:/LRHWork/Model_RH/PyBaMM-ECdrag2/wip/Rio_Code/P3_R14/Debug_Xitilde.ipynb"

OLD_HEADER = '''def FindClose(time_DD, t_select):
    """Nearest index for t_select in list time_DD."""
    for i in range(len(time_DD)):
        if abs(time_DD[i] - t_select) < 1e-3:
            break
    return i

# compare conccnetration
def Plot_Concentration_1_Crate('''

NEW_HEADER = '''def _nearest_time_index(time_DD, t_select):
    """Index of element in time_dd closest to t_select (robust vs float noise)."""
    t = np.asarray(time_DD, dtype=float)
    return int(np.argmin(np.abs(t - float(t_select))))


# compare conccnetration
def Plot_Concentration_1_Crate('''


def main():
    with open(P, encoding="utf-8") as f:
        nb = json.load(f)
    c18 = nb["cells"][18]
    s = "".join(c18.get("source", []))
    if OLD_HEADER not in s:
        raise SystemExit("cell 18 header not found")
    s = s.replace(OLD_HEADER, NEW_HEADER).replace("FindClose(", "_nearest_time_index(")
    c18["source"] = [ln + "\n" for ln in s.split("\n")]
    if c18["source"] and c18["source"][-1] == "\n":
        c18["source"][-1] = ""

    c19 = nb["cells"][19]
    s9 = "".join(c19.get("source", []))
    c19["source"] = [
        ln + "\n" for ln in s9.replace("FindClose(", "_nearest_time_index(").split("\n")
    ]
    if c19["source"] and c19["source"][-1] == "\n":
        c19["source"][-1] = ""

    with open(P, "w", encoding="utf-8", newline="\n") as f:
        json.dump(nb, f, ensure_ascii=False, indent=2)
        f.write("\n")
    print("OK:", P)


if __name__ == "__main__":
    main()
