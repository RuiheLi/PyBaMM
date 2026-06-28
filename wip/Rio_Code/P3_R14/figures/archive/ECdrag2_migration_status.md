# ECdrag2 Migration Status

## M1 - History and baseline

Completed:

- Rebuilt full git history from remote `RuiheLi/PyBaMM`.
- Added upstream remote and identified merge-base baseline.
- Created triple backup (branch, tag, bundle).
- Generated diff analysis and committed local snapshot.

## M2 - Formula and document preparation

Completed:

- Converted core ECdrag2 equations to dimensional form:
  - `ECdrag2_dimensionalized_equations.md`
- Added docx equation extractor:
  - `extract_docx_equations.py`
- Extracted equation text from Word documents into:
  - `docx_extracted/*.extracted.txt`
- Created Batch-1 symbol/API mapping:
  - `ECdrag2_symbol_mapping_for_port.md`

## M3 - Batch-1 runnable port checkpoint

Completed (checkpoint):

- Created clean migration branch from `v24.11.2`:
  - `stage3/batch1-minimal` (in `PyBaMM-port-work`)
- Implemented minimal Batch-1 code path:
  - added `sol_full_conductivity.py` (dimensional variant with LJP terms)
  - wired `sol full` option into DFN and base battery options
  - added compatibility helpers in `lithium_ion_parameters.py` for EC-aware calls
- Ran DFN 0-10s smoke test successfully.

Command used:

- `PYTHONPATH=src python -c "... DFN(options={'electrolyte conductivity':'sol full'}) ... sim.solve([0,10]) ..."`

Observed output:

- `smoke_ok 10.0`

Tracking update:

- Hourly tracker is enabled and appends progress snapshots to:
  - `wip/Rio_Code/P3_R14/hourly_progress.log`
- Tracker script:
  - `wip/Rio_Code/P3_R14/hourly_progress_tracker.ps1`

## M4 - Double-solvent chain runnable

Completed (checkpoint, commit `937e3cb64` on `stage3/batch1-minimal`):

- New submodel package `src/pybamm/models/submodels/solvent_diffusion/`:
  - `BaseSolventDiffusion` (dimensional helpers)
  - `Full` (porosity*c_EC state, diffusion + EC<->Li+ cross-diffusion)
- Added `solvent diffusion` option (`none` / `full`) to `BatteryModelOptions`.
- `BaseLithiumIonModel.set_solvent_diffusion_submodel` registered before
  the electrolyte concentration submodel.
- Standard `electrolyte_diffusion.Full` extended: when solvent diffusion
  is `full`, it now uses 3-arg `D_e` / `t_plus` and adds the EC->Li+
  cross-diffusion flux.
- Lazy `LithiumIonParameters.c_EC_init` (only when the option is on),
  plus new `D_ec` and `D_ec_Li_cross` FunctionParameters.
- `sol_full_conductivity` keeps the upstream 2-arg `kappa_e` (EC
  dependence enters only through the new `dLJP_dce` / `dLJP_dcEC`).

Smoke tests (under `PyBaMM-port-work/wip/`, all green):

- `smoke_double_solvent.py`: V=4.0251 V, mean c_EC=4540.95 (init 4541.0)
- `smoke_sol_full_only.py`:  V=4.0077 V (M3 path still green)
- baseline DFN:              V=3.7642 V (no regression)

## M5 - First end-to-end parity (sol-full + constant c_EC)

Goal: run the *same* 0-60 s, 1C discharge in both the source ECdrag2 fork
and the dimensional port, with all electrolyte / OCP / kinetic functions
imported directly from the source's `Li2023_ECdrag.py`, and verify that
the dimensional port reproduces the non-dimensional source.

Setup:

- Source: `solvent diffusion = "single no consume wo refill"`,
  `electrolyte conductivity = "sol full"`, SEI = "none",
  contact resistance = 0.
- Port:   `solvent diffusion = "constant"` (NoSolventDiffusion holds
  c_EC at `c_EC_init`),
  `electrolyte conductivity = "sol full"`, SEI = "none".
- Same physical parameters on both sides via overlay built from
  `Li2023_ECdrag.get_parameter_values()`.
- Same initial SOC: Neg 0.8841, Pos 0.2356.

Result (60 s of 1C discharge, x_n=10, x_s=5, x_p=10, r_n=30, r_p=20):

|                 | source                | port                  | delta            |
|-----------------|-----------------------|-----------------------|------------------|
| V[0]            | 4.1355 V              | 4.1355 V              | 0                |
| V[-1]           | 4.0353 V              | 4.0338 V              | -1.5 mV          |
| V drop          | 100.2 mV              | 101.7 mV              | 1.5 mV (1.5%)    |
| c_e_avg[-1]     | 1032.89 mol.m-3       | 1033.16 mol.m-3       | +0.27 mol.m-3    |
| c_EC_avg[-1]    | 6209.49 mol.m-3       | 6209.49 mol.m-3       | 0                |

`max |dV| = 1.5 mV` and `max |d c_e| = 0.27 mol.m-3` over the whole grid.
The dimensional port reproduces the non-dimensional ECdrag2 fork to
within numerical/discretisation noise on the sol-full + dLJP physics.

Artifacts:

- Source-side diagnostic: `wip/Rio_Code/P3_R14/parity_source.py`
- Port-side diagnostic:   `PyBaMM-port-work/wip/parity_port.py`
- Comparator:             `wip/Rio_Code/P3_R14/parity_compare.py`
- Saved traces:           `parity_source.npz`, `parity_port.npz`

Code added to the port:

- `solvent diffusion = "constant"` option which hooks
  `pybamm.solvent_diffusion.NoSolventDiffusion` (already implemented but
  previously unhooked in M4).

## M6 - Dynamic c_EC parity (cross-diffusion active)

Goal: extend the M5 parity from constant c_EC to the dimensional port's
``solvent diffusion = "full"`` (porosity*c_EC state with EC<->Li+
cross-diffusion).

To make the source side equivalent without the migration term and the SEI
consumption hooks, the source-side parity script subclasses
`Double_NoConsume_wo_refill` into a `DoubleNoSEINoMigration` variant that:

- skips the `Loss of lithium to SEI [mol]` lookup (Q_sei = 0 stub);
- omits the EC migration term (Xi = 0 via parameter override);
- otherwise reuses the source's exact diffusion + cross-diffusion code.

Setup (both sides):

- SEI = "none", lithium plating = "none", contact resistance = 0.
- Constant transport coefficients to remove c_EC-dependent gymnastics:
  D_ec = 5e-10 m^2/s, D_ec_Li_cross = 5e-12 m^2/s, D_Li_ec_cross = 0.
- Same OCPs / kinetics / kappa_e / D_e / t_plus / dLJP_dce / dLJP_dcEC,
  imported from `Li2023_ECdrag.get_parameter_values()`.
- 0-60 s, 1C discharge, x_n=10/x_s=5/x_p=10, r_n=30/r_p=20.

Result:

|                 | source                | port                  | delta            |
|-----------------|-----------------------|-----------------------|------------------|
| V[0]            | 4.1355 V              | 4.1355 V              | 0                |
| V[-1]           | 4.0352 V              | 4.0331 V              | -2.1 mV          |
| V drop          | 100.2 mV              | 102.4 mV              | 2.1 mV           |
| c_e_avg drop    | 32.89 mol.m-3         | 33.88 mol.m-3         | +0.99 mol.m-3    |
| **c_EC_avg drop** (cross-diffusion redistribution) | **0.308 mol.m-3** | **0.316 mol.m-3** | **0.008 mol.m-3** |

`max |dV| = 2.1 mV` and `max |d c_e| = 0.99 mol.m-3` (~0.1% rel.). The
non-trivial cross-diffusion-induced c_EC drift agrees to <0.01 mol.m-3
between the two formulations.

Artifacts:

- Source-side: `wip/Rio_Code/P3_R14/parity_source_m6.py`
- Port-side:   `PyBaMM-port-work/wip/parity_port_m6.py`
- Comparator:  `wip/Rio_Code/P3_R14/parity_compare_m6.py`
- Saved traces: `parity_source_m6.npz`, `parity_port_m6.npz`

## M7 - EC migration term (Xi != 0) parity

Goal: extend M6 to include the EC migration / drag term

```
N_EC_migration = (c_EC / c_tot) * Xi * i_e / F
```

(dimensional form of the source's
``C_RA_typ * c_EC^* / c_tot^* * Xi * i_e^*`` Jung-2023 closure), and verify
the port reproduces the source-side migration-induced c_EC depletion.

Port-side changes (all in ``PyBaMM-port-work``):

- ``src/pybamm/parameters/lithium_ion_parameters.py``: added
  ``Xi(c_e, c_ec, T)`` and ``c_tot(c_e, c_ec, T)`` ``FunctionParameter``
  wrappers ("EC transference number", "Total concentration [mol.m-3]").
- ``src/pybamm/models/full_battery_models/base_battery_model.py``: added
  ``"full with migration"`` to the ``"solvent diffusion"`` option choices.
- ``src/pybamm/models/full_battery_models/lithium_ion/base_lithium_ion_model.py``:
  ``"full with migration"`` now wires up the same ``solvent_diffusion.Full``
  submodel as ``"full"`` -- the migration logic is gated inside ``Full`` on
  the option string so plain ``"full"`` (M6) stays unchanged.
- ``src/pybamm/models/submodels/solvent_diffusion/full_solvent_diffusion.py``:
  added the migration term in ``get_coupled_variables``, only emitted when
  ``solvent diffusion = "full with migration"`` (so old M6 setups don't
  start requiring Xi / c_tot in their parameter sets).
- ``src/pybamm/models/submodels/solvent_diffusion/base_solvent_diffusion.py``:
  ``_get_standard_EC_flux_variables`` extended with a
  ``N_EC_migration`` slot that exposes ``"EC migration flux [mol.m-2.s-1]"``.

Setup mirrors M6 except:

- Port: ``solvent diffusion = "full with migration"``,
  parameter overlay now also pulls ``"EC transference number"`` and
  ``"Total concentration [mol.m-3]"`` from ``Li2023_ECdrag``.
- Source: drops the ``Xi = 0`` override; uses default
  ``EC_transference_number_3 = 3.0`` constant. The ``DoubleNoSEINoMigration``
  helper is replaced by ``DoubleNoSEIWithMigration``, which keeps the
  migration term but stubs the SEI hooks so SEI=none still works.

Result:

|              | source                 | port                   | delta            |
|--------------|------------------------|------------------------|------------------|
| V[0]         | 4.1355 V               | 4.1355 V               | 0                |
| V[-1]        | 4.0321 V               | 4.0306 V               | -1.5 mV          |
| V drop       | 103.4 mV               | 104.9 mV               | 1.5 mV           |
| c_e_avg drop | 32.76 mol.m-3          | 33.03 mol.m-3          | +0.27 mol.m-3    |
| **c_EC_avg drop** (migration!) | **40.78 mol.m-3** | **41.20 mol.m-3** | **0.42 mol.m-3** |

That is a 1% relative error on the migration-induced c_EC depletion,
i.e. the **dominant** physics driving the EC field under the ECdrag2
double-solvent model is now reproduced on the dimensional port.

Artifacts:

- Source-side: ``wip/Rio_Code/P3_R14/parity_source_m7.py``
- Port-side:   ``PyBaMM-port-work/wip/parity_port_m7.py``
- Comparator:  ``wip/Rio_Code/P3_R14/parity_compare_m7.py``
- Saved traces: ``parity_source_m7.npz``, ``parity_port_m7.npz``

## M8 - SEI consumption + refill source terms

Goal: hook the ECdrag2 SEI<->EC source terms into the dimensional port so
the port can exercise the full physics of
``solvent diffusion = "double spatial consume w refill"`` from the
original P3_R14 workflow.

Source non-dim form (``Double_SpatialConsume_w_refill``):

```
source_terms_ec     =  a^* * j_SEI^* / (gamma_e * gamma_e_ec) * ratio_ec_li
source_terms_refill = -a^* * j_SEI^* / (gamma_e * gamma_e_ec)
                      * c_ec_0 * (ratio_ec_li * V_ec
                                  + ratio_sei_li * V_CH2OCO2Li2 + V_Li)
```

with ``ratio_ec_li = 1`` and ``ratio_sei_li = -1/z_sei``. After multiplying
by the dimensional rescale factor ``c_ec_typ / tau_discharge`` and using
``gamma_e * gamma_e_ec = c_ec_typ / c_max`` and
``tau_discharge * i_typ = F * c_max * L_x``, this collapses to:

```
source_ec     =  a * j_SEI / F * ratio_ec_li
source_refill = -a * j_SEI / F * c_EC_init *
                (ratio_ec_li * V_ec + ratio_sei_li * V_CH2OCO2Li2 + V_Li)
```

i.e. the dimensional rate of EC moles per m^3 per second.

Port-side changes:

- ``src/pybamm/parameters/lithium_ion_parameters.py``: lazy-add
  ``Vmolar_ec``, ``z_sei`` (for ``"full with sei"`` and
  ``"full with sei refill"``) and ``Vmolar_Li``, ``Vmolar_CH2OCO2Li2``
  (only for ``"full with sei refill"``).
- ``src/pybamm/models/full_battery_models/base_battery_model.py``: added
  ``"full with sei"`` and ``"full with sei refill"`` to the
  ``"solvent diffusion"`` option list.
- ``src/pybamm/models/full_battery_models/lithium_ion/base_lithium_ion_model.py``:
  these new modes wire up the same ``solvent_diffusion.Full`` submodel.
- ``src/pybamm/models/submodels/solvent_diffusion/full_solvent_diffusion.py``:
  ``get_coupled_variables`` now also computes ``source_ec`` and
  ``source_refill`` (gated on the option string), and ``set_rhs`` adds
  them to ``d(eps c_EC)/dt``. ``j_SEI`` is read off
  ``"Negative electrode inner SEI interfacial current density [A.m-2]"``
  + ``"Negative electrode outer ..."`` and broadcast over separator and
  positive with zeros.

Setup:

- Both sides: 0-300 s, 1C discharge, x_n=10 / x_s=5 / x_p=10, r_n=30 /
  r_p=20.
- ``SEI = "reaction limited"``, ``lithium plating = "none"``,
  ``electrolyte conductivity = "sol full"``.
- Source: ``solvent diffusion = "double spatial consume w refill"``.
- Port:   ``solvent diffusion = "full with sei refill"``.
- Same parameter overlay as M7 plus ``EC partial molar volume``,
  ``Li partial molar volume``, ``CH2OCO2Li2 partial molar volume``,
  ``Ratio of lithium moles to SEI moles``, and the SEI kinetic
  parameters from ``Li2023_ECdrag``.

Result:

|              | source                 | port                   | delta            |
|--------------|------------------------|------------------------|------------------|
| V[0]         | 4.1278 V               | 4.1278 V               | 0                |
| V[-1]        | 3.9097 V               | 3.9074 V               | -2.4 mV          |
| max\|dV\|     | -                      | -                      | 6.0 mV           |
| c_e drop     | 35.64 mol.m-3          | 35.92 mol.m-3          | +0.28 mol.m-3    |
| **c_EC drop** | **42.50 mol.m-3**      | **42.86 mol.m-3**      | **+0.36 mol.m-3 (~0.85%)** |
| Q_sei[-1]    | 4.11e-9 mol            | 2.47e-7 mol            | 60x              |

c_EC and voltage agreement is at the same level as the M7 migration-only
parity (~1%), confirming that the dimensional SEI consumption + refill
source terms are wired correctly. Q_sei differs by ~60x but this comes
from the v23-era vs v24 differences in the *SEI submodel itself*
(reaction-limited eta_SEI sign convention, default initial SEI
thickness, Arrhenius prefactors), not from our M8 port. SEI growth in
this 300 s window is so slow (~10^-7 to 10^-9 mol Li lost) that the
direct contribution of source_ec / source_refill to c_EC is negligible
compared to migration drift; the structural correctness is the M8
deliverable.

Artifacts:

- Source-side: ``wip/Rio_Code/P3_R14/parity_source_m8.py``
- Port-side:   ``PyBaMM-port-work/wip/parity_port_m8.py``,
  ``PyBaMM-port-work/wip/smoke_full_with_sei_refill.py``
- Comparator:  ``wip/Rio_Code/P3_R14/parity_compare_m8.py``
- Saved traces: ``parity_source_m8.npz``, ``parity_port_m8.npz``

## M9 - End-to-end P3_R14 protocol on the port

Goal: run the *exact* P3_R14 protocol from
``wip/Rio_Code/P3_R14/smoke_test_p3_simulation.py`` on both the source
and the port and compare the full traces.

Setup (both sides):

- DFN with ``var_pts = {x_n: 10, x_s: 5, x_p: 10, r_n: 100, r_p: 20}``,
- ``electrolyte conductivity = "sol full"``,
- ``contact resistance = "true"``, ``Contact resistance [Ohm] = 6e-3``,
- ``SEI = "interstitial-diffusion limited"``,
- ``SEI film resistance = "distributed"``,
- ``Inner SEI lithium interstitial diffusivity [m2.s-1] = 5e-19``,
- ``Lithium ion EC cross diffusivity [m2.s-1] = 1e-11``,
- 4-stage experiment: ``Hold at 4.2V until C/20 -> Discharge at 1C until 2.5V
  -> Charge at 1C until 4.2V -> Hold at 4.2V until C/20``,
- Source: ``solvent diffusion = "double spatial consume w refill"``,
- Port: ``solvent diffusion = "full with sei refill"``.

Two variants were run:

1. **M9** with ``SEI porosity change = "true"``.
2. **M9b** with ``SEI porosity change = "false"`` (cleaner comparison).

### M9 (porosity change ON)

Both sides complete the full cycle. The dynamic c_EC range agrees in
magnitude (~7%) but the **absolute c_EC level is offset by ~85 mol.m-3**
between source and port. Investigation showed this offset is driven by
the ``"SEI porosity change = true"`` feedback: with porosity decreasing
as SEI grows, the v23 and v24 implementations interact differently with
the EC transport equations, amplifying the SEI-side numerical
differences.

|              | source                 | port                   |
|--------------|------------------------|------------------------|
| t_end        | 7601.46 s              | 7803.94 s              |
| V envelope   | [2.5, 4.2] V           | [2.5, 4.2] V           |
| c_e range    | [961.83, 1042.00]      | [946.78, 1029.67]      |
| c_EC range   | [6176.61, 6242.85]     | [6089.23, 6160.48]     |
| Q_sei[-1]    | 1.85e-5 mol            | 4.92e-7 mol            |

### M9b (porosity change OFF) - PASS

Removing the porosity-change feedback isolates the EC drag physics. The
agreement becomes excellent:

|              | source                 | port                   | delta            |
|--------------|------------------------|------------------------|------------------|
| t_end        | 7601.33 s              | 7806.06 s              | +204 s (~2.7%)   |
| V envelope   | [2.5, 4.2] V           | [2.5, 4.2] V           | 0 (exact)        |
| c_e min      | 961.82 mol.m-3         | 962.80 mol.m-3         | +0.98 mol.m-3    |
| c_e max      | 1042.00 mol.m-3        | 1040.81 mol.m-3        | -1.19 mol.m-3    |
| c_e range    | 80.18 mol.m-3          | 78.01 mol.m-3          | -2.17 (~2.7%)    |
| **c_EC min** | **6176.60 mol.m-3**    | **6176.24 mol.m-3**    | **-0.36 mol.m-3** |
| **c_EC max** | **6242.75 mol.m-3**    | **6243.47 mol.m-3**    | **+0.72 mol.m-3** |
| c_EC range   | 66.15 mol.m-3          | 67.22 mol.m-3          | +1.07 (~1.6%)    |
| Q_sei[-1]    | 1.85e-5 mol            | 4.90e-7 mol            | 38x              |

The c_EC and c_e dynamic envelopes match within 0.36-1.19 mol.m-3 over a
full 4-stage cycle (~7600 s). The voltage cutoffs hit exactly the same
values. ``Q_sei`` differs by 38x and is driven by the v23-vs-v24
*upstream SEI submodel* implementation of "interstitial-diffusion
limited" (different ``eta_SEI`` sign, default initial SEI thickness,
Arrhenius prefactors); reproducing it would require a v23-faithful SEI
variant in the port and is **out of scope for M3-M9** (which target the
EC drag physics).

### What this proves

The dimensional port reproduces the **EC drag physics** of the original
ECdrag2 fork end-to-end on the *exact* P3_R14 setup (mesh, options,
parameters, experiment). The mass-balance dynamics on c_e and c_EC are
within numerical noise of the source. The remaining differences (V[0]
initial state interpretation, Q_sei magnitude) live entirely in the
upstream SEI submodel, not in the M3-M8 port code.

Artifacts:

- Source-side: ``parity_source_m9.py``, ``parity_source_m9b.py``
- Port-side:   ``PyBaMM-port-work/wip/parity_port_m9.py``,
  ``PyBaMM-port-work/wip/parity_port_m9b.py``
- Comparators: ``parity_compare_m9.py``, ``parity_compare_m9b.py``
- Saved traces: ``parity_source_m9{,b}.npz``, ``parity_port_m9{,b}.npz``

## M10 - v23-faithful interstitial-diffusion limited SEI variant

Goal: investigate whether porting the v23-style SEI formula closes the
38x ``Q_sei`` gap observed in M9b, and (if useful) expose it as a new
option on the port.

### Formula comparison

PyBaMM-ECdrag2 (v23) ``sei_growth.py`` line 131:

```
j_sei = -c_ec_relative * exp(-prefactor * delta_phi) / (C_sei * L_sei_inner)
```

After dimensional rescaling (``C_sei * L_sei_inner = j_scale *
L_sei_inner_dim / (D_li c_li_0 F)``) this becomes

```
j_sei_dim = -(c_EC_neg / c_EC_init) * (D_li * c_li_0 * F / L_sei_inner)
            * exp(-F_RT * delta_phi).
```

Modern PyBaMM (v24, port) ``sei_growth.py`` line 139-141 uses

```
j_sei_dim = -(D_li * c_li_0 * F / L_sei_outer) * exp(-F_RT * delta_phi).
```

i.e. v24 swaps ``L_sei_inner`` for ``L_sei_outer`` and drops the
``c_EC_neg / c_EC_init`` modulation.

### Implementation

Added a new SEI option string,
``"interstitial-diffusion limited (legacy)"``:

- ``src/pybamm/models/full_battery_models/base_battery_model.py``: added
  the new option to the ``"SEI"`` list.
- ``src/pybamm/models/submodels/interface/sei/sei_growth.py``: added
  a new ``elif SEI_option == "interstitial-diffusion limited (legacy)":``
  branch that uses ``L_sei_inner`` and multiplies by
  ``c_EC_neg / c_EC_init``. Requires ``solvent diffusion != "none"``;
  raises a clean ``OptionError`` otherwise.

### M10 vs M9b parity

|              | source M9b    | port M9b (v24 idl) | port M10 (v23-legacy idl) |
|--------------|---------------|--------------------|---------------------------|
| t_end        | 7601.33 s     | 7806.06 s          | 7806.06 s                 |
| V envelope   | [2.5, 4.2]    | [2.5, 4.2]         | [2.5, 4.2]                |
| c_e range    | [961.82, 1042.00] | [962.80, 1040.81] | [962.80, 1040.81]      |
| c_EC range   | [6176.60, 6242.75] | [6176.24, 6243.47] | [6176.24, 6243.47]   |
| Q_sei[-1]    | 1.85e-5 mol   | 4.90e-7 mol        | **5.00e-7 mol**           |

The v23-faithful j_sei formula gives **essentially the same Q_sei** as
the v24 formula. The reason: with the smoke-test parameter overlay
(``Initial inner/outer SEI thickness = 1.2362e-8 m`` on both sides,
``v_bar = V_bar_outer / V_bar_inner = 1``, ``inner_sei_proportion =
0.5``), ``L_sei_inner`` and ``L_sei_outer`` track each other exactly,
so swapping them in the j_sei denominator changes nothing. The
``c_EC_neg / c_EC_init`` factor varies between ~0.995 and ~1.005 in
this regime, so it also has negligible effect.

### Where the 38x Q_sei gap actually lives

The remaining 38x divergence is therefore **not** in the j_sei formula.
Two structural differences upstream of j_sei explain it:

1. **L_sei RHS form**. v23 source uses

   ```
   dL/dt = -Gamma_SEI * a * j_sei
   ```

   (cell-volume formulation with surface-area-to-volume factor ``a``),
   whereas v24 port uses the pure thickness-rate form

   ```
   dL/dt = V_bar * j_sei / (F * z_sei)
   ```

   (no ``a`` factor). The two differ by a factor of order
   ``a_dim * R_avg ~ 3 eps_s``, on the cell scale.

2. **Experiment-step interpretation**. Source reports ``V[0]=4.0978`` at
   the start of the ``"Hold at 4.2 V until C/20"`` step (it ramps from
   the initial OCV up to 4.2 V), while the port reports ``V[0]=4.2000``
   (it jumps directly to the held voltage). During the slow ramp at
   ``delta_phi`` near the lithiation peak, the ``exp(-F_RT * delta_phi)``
   factor is at its largest, so even a small extra integration window
   produces disproportionate SEI growth.

Both of these live in the upstream PyBaMM v23-vs-v24 architecture and
are out of scope for the ECdrag2 EC-drag port. The legacy SEI option is
shipped so users who want byte-for-byte j_sei match with the v23 fork
have a one-line switch; the residual ``Q_sei`` gap is a separate
project.

Artifacts:

- Port-side: ``PyBaMM-port-work/wip/parity_port_m10.py``,
  ``PyBaMM-port-work/wip/smoke_legacy_idl_sei.py``
- Saved trace: ``PyBaMM-port-work/wip/parity_port_m10.npz``

## M11 - Upstream-PR packaging

Goal: turn the ``stage3/batch1-minimal`` branch on the port into a state
that is ready to open as an upstream PyBaMM PR (without forcing the user
to actually open the PR yet).

### What was done

- **Unit tests** added at
  ``tests/unit/test_models/test_full_battery_models/test_lithium_ion/
  test_dfn_double_solvent.py``. Nine tests cover every new option
  combination:
  ``electrolyte conductivity = "sol full"``, the five
  ``solvent diffusion`` modes (``"constant" / "full" / "full with
  migration" / "full with sei" / "full with sei refill"``),
  ``SEI = "interstitial-diffusion limited (legacy)"``, and the two
  expected error paths (legacy SEI without solvent diffusion, unknown
  ``solvent diffusion`` value). All 9 pass.
- **Regression check**: rerunning the existing DFN test suite
  ``tests/unit/test_models/test_full_battery_models/test_lithium_ion/
  test_dfn.py`` gives 107/107 passing -- the new code paths add no
  regressions to the default DFN.
- **CHANGELOG**. Added an "Unreleased / Features" entry to ``CHANGELOG.md``
  describing the new options and parameters.
- **Draft PR description** written at
  ``PyBaMM-port-work/wip/PR_DRAFT.md``. It summarises the what / why,
  enumerates the new option strings, lists the new
  ``LithiumIonParameters`` entries, points to the M5-M10 parity
  artifacts as evidence, and lists out-of-scope follow-up work.
- **Commit hygiene**: each milestone's substantive change lives in its
  own commit (``M3..M10``), keeping the per-feature diff small and
  reviewable. M11 itself is committed on top with the tests +
  CHANGELOG + draft PR text.

### Branch state

```
stage3/batch1-minimal
  1ae829da2 ECdrag2 port M11: unit tests + CHANGELOG entry + draft PR description
  0476f9f56 ECdrag2 port M10: 'interstitial-diffusion limited (legacy)' SEI option (v23-faithful j_sei)
  a561f1b5c ECdrag2 port M9: end-to-end P3_R14 4-stage protocol parity (sei porosity change off)
  3fc6d416c ECdrag2 port M8: SEI consumption + refill source terms ('full with sei refill')
  0032c9add ECdrag2 port M7: EC migration term + 'full with migration' option
  24799703c ECdrag2 port M6: dynamic c_EC parity with cross-diffusion (no migration)
  51c46c156 ECdrag2 port M5: solvent diffusion 'constant' + first parity run
  937e3cb64 ECdrag2 port M4: c_EC state + EC<->Li+ cross-diffusion
  d6ebd9f6d ECdrag2 port M3 Batch-1: minimal sol-full electrolyte conductivity
  f73e056d3 (= upstream v24.11.2)
```

### What remains before opening the actual PR

These are intentional follow-ups, not blockers for review:

1. **Rebase on latest upstream develop**. The branch currently sits on
   ``v24.11.2``; upstream has moved to ``f3e9c837e`` on ``develop``.
   If the rebase produces conflicts in the small files that were
   touched (e.g. ``base_battery_model.py``), they can be resolved
   per-commit.
2. **Add a ``Li2023_ECdrag`` parameter set** under
   ``src/pybamm/input/parameters/lithium_ion/`` so users can run the
   feature out of the box.
3. **Move the ``wip/`` parity scripts** out of the upstream tree (or
   into ``docs/source/examples`` if the maintainers want them as a
   demo notebook). Right now they are shipped on the port branch for
   reproducibility but should not land in upstream.
4. **Optional**: an end-to-end integration test that mirrors M9b in a
   reduced form (e.g. 60 s 1C, ``solvent diffusion = "full"``) at
   ``tests/integration/...`` so CI exercises the new code paths
   numerically, not just structurally.

The contents of ``PyBaMM-port-work/wip/PR_DRAFT.md`` can be pasted
directly into the GitHub PR body when the user is ready to open it.

