"""Patch Compare_parameter_value.ipynb for current PyBaMM API."""
import json
from pathlib import Path

path = Path(__file__).resolve().parent / "Compare_parameter_value.ipynb"
nb = json.loads(path.read_text(encoding="utf-8"))


def set_cell(i: int, src: str):
    lines = src.splitlines(keepends=True)
    if not lines:
        nb["cells"][i]["source"] = []
        return
    nb["cells"][i]["source"] = lines


# Cell 1
set_cell(
    1,
    """Para_Ecker2015 = pb.ParameterValues("Ecker2015")
# Ren2018 不在当前仓库参数集中；使用 ORegan2022（LGM50）作为替代
Para_Ren2018 = pb.ParameterValues("ORegan2022")
Para_Chen2020 = pb.ParameterValues("Chen2020")
""",
)

# Cell 2
set_cell(
    2,
    """str_path_0 = os.path.abspath(os.path.join(pb.__path__[0], ".."))
BasicPath = os.path.join(str_path_0, "wip", "Rio_Code", "Compare_parameter_value_output")
Target = "InputData" + os.sep
os.makedirs(os.path.join(BasicPath, Target), exist_ok=True)
""",
)

# Cell 3: drop first duplicate graphite_LGM50_diffusivity_Chen2020 block
c3_lines = "".join(nb["cells"][3]["source"]).splitlines()
idxs = [i for i, l in enumerate(c3_lines) if l.startswith("def graphite_LGM50_diffusivity_Chen2020")]
if len(idxs) >= 2:
    # keep from second def onward, but keep header comment
    header = []
    for i, line in enumerate(c3_lines):
        if line.startswith("def graphite_LGM50_diffusivity_Chen2020"):
            break
        header.append(line)
    # skip first duplicate: from idxs[0] to line before "# Compare electrode" or second def
    # Simpler: take lines[0:idxs[0]] + lines[idxs[1]:]  drops first function body including first def
    merged = c3_lines[: idxs[0]] + c3_lines[idxs[1] :]
    set_cell(3, "\n".join(merged) + "\n")

# Cell 4: remove extra arg; optional rename
c4 = "".join(nb["cells"][4]["source"])
c4 = c4.replace(
    "graphite_LGM50_electrolyte_exchange_current_density_Chen2020  (\n"
    "        1000, NegSOC_i*33133, 298.15, Para_Chen2020_coupled  ) )",
    "graphite_LGM50_electrolyte_exchange_current_density_Chen2020(\n"
    "        1000, NegSOC_i*33133, 298.15))",
)
set_cell(4, c4)

# plt.savefig paths: BasicPath + Target -> join (cell 9, 10, 11 use BasicPath + Target)
for idx in (9, 10, 11):
    s = "".join(nb["cells"][idx]["source"])
    if "BasicPath + Target" in s or "BasicPath +  Target" in s:
        s = s.replace("BasicPath +  Target+", "os.path.join(BasicPath, Target) + ")
        s = s.replace("BasicPath + Target+", "os.path.join(BasicPath, Target) + ")
        s = s.replace(
            'plt.savefig(BasicPath + Target+"',
            'plt.savefig(os.path.join(BasicPath, Target) + "',
        )
        set_cell(idx, s)

path.write_text(json.dumps(nb, ensure_ascii=False, indent=1), encoding="utf-8")
print("Wrote", path)
