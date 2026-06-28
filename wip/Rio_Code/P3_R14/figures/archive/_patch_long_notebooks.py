"""Patch long-run notebooks: solver max_num_steps; Reproduce_sol_seg mesh + env rate."""
import json

ROOT = r"D:/LRHWork/Model_RH/PyBaMM-ECdrag2/wip/Rio_Code/P3_R14"


def patch_reproduce_sol_seg():
    p = f"{ROOT}/Reproduce_sol_seg.ipynb"
    with open(p, encoding="utf-8") as f:
        nb = json.load(f)

    for c in nb["cells"]:
        if c.get("cell_type") != "code":
            continue
        s = "".join(c.get("source", []))
        orig = s

        if '"Mesh list":[ [20,10,20,100,20]' in s:
            s = s.replace(
                '"Mesh list":[ [20,10,20,100,20], ],   # Simon uses 30',
                '"Mesh list":[ [20,10,20,80,20], ],   # was 100; 80 matches Xitilde notebook & avoids extreme stiffness at high C (see cell 5 comments)',
            )

        if "Rate_Dis = 1.9" in s and "P3_REPRO_SOL_SEG_RATE" not in s:
            s = s.replace(
                "Sol_DD = []; Rate_Dis = 1.9",
                'Sol_DD = []; Rate_Dis = float(os.environ.get("P3_REPRO_SOL_SEG_RATE", "1.9"))  # set e.g. 1.5 to debug faster',
            )

        if "solver = pybamm.CasadiSolver(return_solution_if_failed_early=True)," in s:
            s = s.replace(
                "solver = pybamm.CasadiSolver(return_solution_if_failed_early=True),",
                "solver = pybamm.CasadiSolver(return_solution_if_failed_early=True, extra_options_setup={\"max_num_steps\": 200000}),",
            )

        if s != orig:
            c["source"] = [line + "\n" for line in s.split("\n")]
            if c["source"] and c["source"][-1] == "\n":
                c["source"][-1] = ""

    with open(p, "w", encoding="utf-8", newline="\n") as f:
        json.dump(nb, f, ensure_ascii=False, indent=2)
        f.write("\n")
    print("patched", p)


def patch_onecyc_2c():
    p = f"{ROOT}/1EC1DMC_OneCycAge_2C.ipynb"
    with open(p, encoding="utf-8") as f:
        nb = json.load(f)

    for c in nb["cells"]:
        if c.get("cell_type") != "code":
            continue
        s = "".join(c.get("source", []))
        orig = s
        if "solver = pybamm.CasadiSolver(return_solution_if_failed_early=True)," in s:
            s = s.replace(
                "solver = pybamm.CasadiSolver(return_solution_if_failed_early=True),",
                "solver = pybamm.CasadiSolver(return_solution_if_failed_early=True, extra_options_setup={\"max_num_steps\": 200000}),",
            )

        if s != orig:
            c["source"] = [line + "\n" for line in s.split("\n")]
            if c["source"] and c["source"][-1] == "\n":
                c["source"][-1] = ""

    with open(p, "w", encoding="utf-8", newline="\n") as f:
        json.dump(nb, f, ensure_ascii=False, indent=2)
        f.write("\n")
    print("patched", p)


if __name__ == "__main__":
    patch_reproduce_sol_seg()
    patch_onecyc_2c()
