# ECdrag2 Port Mapping (Batch-1)

## Scope

This batch maps two core modules for dimensional migration:

- `pybamm/models/submodels/electrolyte_conductivity/sol_full_conductivity.py`
- `pybamm/models/submodels/electrolyte_diffusion/full_diffusion.py`

Target reference for modern PyBaMM dimensional style:

- `src/pybamm/models/submodels/electrolyte_conductivity/full_conductivity.py`
- `src/pybamm/models/submodels/electrolyte_diffusion/full_diffusion.py`
- `src/pybamm/parameters/lithium_ion_parameters.py`

## A) Variable and API mapping

### A1. Conductivity module (`sol_full_conductivity`)

- Old state variable:
  - `phi_e` (dimensionless electrolyte potential)
- New state variable:
  - `Phi_e` with key `"Electrolyte potential [V]"`

- Old concentration keys:
  - `"Electrolyte concentration"` (dimensionless)
  - `"EC concentration"` (dimensionless)
- New concentration keys:
  - `"Electrolyte concentration [mol.m-3]"`
  - new EC key should be dimensional if retained, e.g. `"EC concentration [mol.m-3]"`

- Old constitutive call:
  - `param.kappa_e(c_e, c_EC, T)` returns dimensionless conductivity
- New constitutive call:
  - should return dimensional conductivity directly, e.g.
  - `param.kappa_e(c_e, c_EC, T)` in `S.m-1`

- Old current variable:
  - `"Electrolyte current density"` (dimensionless)
- New current variable:
  - `"Electrolyte current density [A.m-2]"`

### A2. Diffusion module (`full_diffusion`)

- Old conserved state:
  - `"Porosity times concentration"` (dimensionless `eps*c_e`)
- New conserved state:
  - `"Porosity times concentration [mol.m-3]"`

- Old electrolyte flux:
  - `"Li+ flux"` (dimensionless)
- New electrolyte flux:
  - `"Electrolyte flux [mol.m-2.s-1]"`

- Old source term:
  - `"Sum of electrolyte reaction source terms"` (dimensionless current source)
  - converted through `/gamma_e`
- New source term:
  - `"Sum of electrolyte reaction source terms [A.m-3]"`
  - converted to molar source through `/F`

## B) Scale-factor removals and replacements

### B1. Remove non-dimensional prefactors from PDE-level equations

Replace these legacy factors in governing equations:

- remove `gamma_e / C_e` in conductivity current law
- remove `C_e / gamma_e` in migration flux term
- remove `/C_e` in divergence terms in RHS
- remove explicit `c_e_typ/potential_scale` and `c_ec_typ/potential_scale` multipliers in LJP gradient terms

Reason: modern PyBaMM equations are dimensional, so these are already absorbed by symbol units/scales.

### B2. Keep only physical dimensional coefficients

Use directly:

- `kappa_e_dim(c_e, c_EC, T)` in `S.m-1`
- `D_e_dim(c_e, c_EC, T)` in `m2.s-1`
- `D_Li_ec_cross_dim(c_e, c_EC, T)` in `m2.s-1`
- `t_plus(c_e, c_EC, T)` dimensionless
- `dLJP_dce`, `dLJP_dcEC` in `V/(mol.m-3)`

## C) Equation-level migration templates

### C1. Electrolyte current (conductivity)

Old implemented form:

- `i_e = (kappa_e * tor * gamma_e/C_e) * (...)`

Target dimensional form:

- `i_e_dim = tor * kappa_e_dim(c_e, c_EC, T) * ( -grad(Phi_e) + dLJP_dce * grad(c_e) + dLJP_dcEC * grad(c_EC) )`

with:

- `i_e_dim` in `A.m-2`
- `Phi_e` in `V`
- `c_e`, `c_EC` in `mol.m-3`

### C2. Li+ flux decomposition (full diffusion)

Old:

- `N_e_diffusion = -tor * D_e * grad(c_e)`
- `N_e_migration = (C_e/gamma_e) * t_plus * i_e`
- `N_cross = -gamma_e_ec * tau_diffusion_e/tau_cross * tor * D_Li_ec_cross * grad(c_EC)`

Target dimensional:

- `N_e_diff = -tor * D_e_dim * grad(c_e)`
- `N_e_mig = t_plus * i_e_dim / F`
- `N_e_cross = -tor * D_Li_ec_cross_dim * grad(c_EC)`
- `N_e = N_e_diff + N_e_cross + N_e_mig + c_e * v_box`

### C3. Li+ conservation RHS

Old:

- `d(eps*c_e)/dt = -div(N_e)/C_e + source_terms - c_e*div_v + source_refill`

Target dimensional:

- `d(eps*c_e)/dt = -div(N_e_dim) + sum_s_a_j/F - c_e*div_v + source_refill_dim`

where `sum_s_a_j` uses `[A.m-3]`.

## D) Word-equation consistency anchors

From extracted `Round_240126_Non_dimensionalization - new`:

- `[EQ0001]`: Li+ flux includes diffusion + cross + migration
- `[EQ0002]`: EC flux includes diffusion + cross + migration with `(c_EC/c_T)*Xi`
- `[EQ0008]`: `C_RA_typ` definition used for EC migration scaling

These anchors are consistent with code and with the dimensional rewrites above.

## E) Implementation checklist for next commit

1. Add dimensional EC parameter functions in migrated parameter class:
   - `kappa_e(c_e, c_EC, T)`
   - `D_e(c_e, c_EC, T)`
   - `D_Li_ec_cross(c_e, c_EC, T)`
   - `dLJP_dce(c_e, c_EC, T)`, `dLJP_dcEC(c_e, c_EC, T)`
2. Implement `sol_full` current law using dimensional variables/units.
3. Patch `full_diffusion` to dimensional flux/source forms (`/F` source conversion).
4. Keep BC migration form consistent with dimensional upstream pattern:
   - `-(1 - t_plus)/(tor*D_e*F) * i_boundary_cc`
5. Run 0-10s DFN smoke with double-solvent options.

