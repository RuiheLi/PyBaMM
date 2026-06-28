"""Inject progress timing into Reproduce_sol_seg.ipynb and 1EC1DMC_OneCycAge_2C.ipynb."""
import json

ROOT = r"D:/LRHWork/Model_RH/PyBaMM-ECdrag2/wip/Rio_Code/P3_R14"


def _set_source(cell, text: str) -> None:
    cell["source"] = [line + "\n" for line in text.split("\n")]
    if cell["source"] and cell["source"][-1] == "\n":
        cell["source"][-1] = ""


def patch_reproduce():
    p = f"{ROOT}/Reproduce_sol_seg.ipynb"
    with open(p, encoding="utf-8") as f:
        nb = json.load(f)
    old = (
        'Sol_DD = []; Rate_Dis = float(os.environ.get("P3_REPRO_SOL_SEG_RATE", "1.9"))  # set e.g. 1.5 to debug faster\n'
        "for i in range(len(Para_DD)):   # len(Para_DD)\n"
        "    Sol_DD.append(RunOne(Para_DD[i],Rate_Dis=Rate_Dis)) \n"
        '    para = Para_DD[i]\n'
        '    print("D_e,EC=",para["Lithium ion EC cross diffusivity [m2.s-1]"])\n'
    )
    new = (
        "import time\n"
        "from datetime import datetime\n"
        'Sol_DD = []; Rate_Dis = float(os.environ.get("P3_REPRO_SOL_SEG_RATE", "1.9"))  # set e.g. 1.5 to debug faster\n'
        "_n = len(Para_DD)\n"
        "_t_batch = time.perf_counter()\n"
        "for i in range(_n):\n"
        '    print(datetime.now(), f"Progress [{i+1}/{_n}] RunOne START", flush=True)\n'
        "    _t0 = time.perf_counter()\n"
        "    Sol_DD.append(RunOne(Para_DD[i],Rate_Dis=Rate_Dis))\n"
        "    _wall = time.perf_counter() - _t0\n"
        "    _batch = time.perf_counter() - _t_batch\n"
        '    para = Para_DD[i]\n'
        '    print("D_e,EC=", para["Lithium ion EC cross diffusivity [m2.s-1]"])\n'
        '    print(datetime.now(), f"Progress [{i+1}/{_n}] wall_s={_wall:.1f} batch_s={_batch:.1f}", flush=True)\n'
    )
    for c in nb["cells"]:
        if c.get("cell_type") != "code":
            continue
        s = "".join(c.get("source", []))
        if old in s:
            s = s.replace(old, new)
            _set_source(c, s)
            break
    else:
        raise SystemExit("reproduce patch: block not found")
    with open(p, "w", encoding="utf-8", newline="\n") as f:
        json.dump(nb, f, ensure_ascii=False, indent=2)
        f.write("\n")
    print("patched Reproduce_sol_seg.ipynb")


def patch_onecyc():
    p = f"{ROOT}/1EC1DMC_OneCycAge_2C.ipynb"
    with open(p, encoding="utf-8") as f:
        nb = json.load(f)
    old = (
        "sol_DD_HDx = Run_OneCycleAge_Dict(1.0, Para_DD[1],Path_pack)\n"
        "sol_DD_LDx = Run_OneCycleAge_Dict(0.0, Para_DD[0],Path_pack)\n"
        "sol_SD     = Run_OneCycleAge_Dict(2.0, Para_SD[0],Path_pack)\n"
    )
    new = (
        "import time\n"
        "from datetime import datetime\n"
        "_jobs = [\n"
        '    ("HDx", 1.0, Para_DD[1]),\n'
        '    ("LDx", 0.0, Para_DD[0]),\n'
        '    ("SD", 2.0, Para_SD[0]),\n'
        "]\n"
        "_t_batch = time.perf_counter()\n"
        "sol_DD_HDx = sol_DD_LDx = sol_SD = None\n"
        "for _j, (_name, _coef, _para) in enumerate(_jobs):\n"
        '    print(datetime.now(), f"OneCyc [{_j+1}/3] {_name} START", flush=True)\n'
        "    _t0 = time.perf_counter()\n"
        "    _sol = Run_OneCycleAge_Dict(_coef, _para, Path_pack)\n"
        "    _wall = time.perf_counter() - _t0\n"
        "    _batch = time.perf_counter() - _t_batch\n"
        '    if _name == "HDx":\n'
        "        sol_DD_HDx = _sol\n"
        '    elif _name == "LDx":\n'
        "        sol_DD_LDx = _sol\n"
        "    else:\n"
        "        sol_SD = _sol\n"
        "    _rem = (_batch / (_j + 1)) * (3 - _j - 1)\n"
        '    print(datetime.now(), f"OneCyc [{_j+1}/3] {_name} DONE wall_s={_wall:.1f} ETA_remain~{_rem/60:.1f}min", flush=True)\n'
    )
    for c in nb["cells"]:
        if c.get("cell_type") != "code":
            continue
        s = "".join(c.get("source", []))
        if old in s:
            s = s.replace(old, new)
            _set_source(c, s)
            break
    else:
        raise SystemExit("onecyc patch: block not found")
    with open(p, "w", encoding="utf-8", newline="\n") as f:
        json.dump(nb, f, ensure_ascii=False, indent=2)
        f.write("\n")
    print("patched 1EC1DMC_OneCycAge_2C.ipynb")


if __name__ == "__main__":
    patch_reproduce()
    patch_onecyc()
