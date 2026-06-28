## Long-run resume (after disconnect)

**Before starting a long run**, check Task Manager / `Get-Process python` for duplicate `debug_long_sims.py reproduce_seg`; two instances double CPU use and may fight over outputs.

`debug_long_sims.py` tries to enable **line buffering** on stdout/stderr before importing PyBaMM, so `Tee-Object`/`>>` usually see `_log()` lines promptly without `PYTHONUNBUFFERED`. You may still use `-u` or `PYTHONUNBUFFERED` if anything remains sticky.

Logs are appended by redirecting stderr/stdout. For **maximum** zeal (optional), force unbuffered Python:

```powershell
cd D:\LRHWork\Model_RH\PyBaMM-ECdrag2\wip\Rio_Code\P3_R14
$env:PYTHONUNBUFFERED = "1"
python debug_long_sims.py reproduce_seg --rate 1.9 --mesh-r-n 80 2>&1 | Tee-Object -FilePath run_seg_19_resume.log -Append
```

Or:

```powershell
python -u debug_long_sims.py reproduce_seg --rate 1.9 --mesh-r-n 80 *>> run_seg_19.log
```

Legacy append-only variant:

```powershell
cd D:\LRHWork\Model_RH\PyBaMM-ECdrag2\wip\Rio_Code\P3_R14

# solvent segregation notebook target (~ may take 25–45+ min)
python debug_long_sims.py reproduce_seg --rate 1.9 --mesh-r-n 80 *>> run_seg_19.log

# full aging triple (heavy)
python debug_long_sims.py onecyc_full *>> run_onecyc_full.log
python debug_long_sims.py onecyc_full --light *>> run_onecyc_light.log
```

`debug_long_sims.py` prints timestamps and ETA for `reproduce_seg_full` / `onecyc_full`.
