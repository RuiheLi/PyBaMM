# ECdrag2 Double-Solvent Model: Dimensionalized Equations (from code)

## Purpose

This note converts the current ECdrag2 implementation from its internal non-dimensional form to dimensional equations, directly from source code behavior in:

- `pybamm/models/submodels/electrolyte_conductivity/sol_full_conductivity.py`
- `pybamm/models/submodels/electrolyte_diffusion/full_diffusion.py`
- `pybamm/models/submodels/solvent_diffusion/Double_SpatialConsume_wo_refill.py`
- `pybamm/parameters/lithium_ion_parameters.py`

The objective is to provide a dimensional equation set that can be ported into the modern dimensional PyBaMM framework.

## 1) Key non-dimensional definitions used by current code

From `lithium_ion_parameters.py`:

- Concentrations:
  - `c_e = c_e_dim / c_e_typ`
  - `c_EC = c_EC_dim / c_ec_typ`
- Potential:
  - `phi_e = (Phi_e_dim + U_ref) / potential_scale`
  - `potential_scale = R T_ref / F`
- Flux scales:
  - Li+ flux scale: `N_e,scale = D_e_typ c_e_typ / L_x`
  - EC flux scale: `N_EC,scale = D_ec_typ c_ec_typ / L_x`
- Conductivity scaling:
  - `kappa_e = kappa_e_dim / kappa_scale`
  - `kappa_scale = F^2 D_e_typ c_e_typ / (R T_ref)`

The code also defines:

- `gamma_e / C_e = F c_e_typ D_e_typ / (i_typ L_x)`
- `gamma_e_ec * tau_diffusion_e / tau_cross * D_Li_ec_cross`
  collapses dimensionally to `D_Li_ec_cross_dim`.

These two identities are what make the final converted equations compact.

## 2) Dimensional electrolyte current equation (sol_full_conductivity)

### 2.1 Code form (non-dimensional)

The implemented expression is:

`i_e = (kappa_e * tor * gamma_e / C_e) * [ -grad(phi_e) + grad(c_e) * dLJP_dce * c_e_typ/potential_scale + grad(c_EC) * dLJP_dcEC * c_ec_typ/potential_scale ]`

where `dLJP_dce` and `dLJP_dcEC` are "Measured dLJP" function parameters.

### 2.2 Dimensionalized form

After substituting scaling relations, the dimensional electrolyte current density is:

`i_e_dim = tor * kappa_e_dim(c_e_dim, c_EC_dim, T_dim) * [ -∇Phi_e_dim + (∂U_LJP/∂c_e) ∇c_e_dim + (∂U_LJP/∂c_EC) ∇c_EC_dim ]`

with:

- `i_e_dim` in A.m^-2
- `kappa_e_dim` in S.m^-1
- `Phi_e_dim` in V
- `c_e_dim`, `c_EC_dim` in mol.m^-3
- `∂U_LJP/∂c` in V / (mol.m^-3)

So this is physically a conductivity-weighted driving force containing:

1. Ohmic migration term (`-∇Phi_e_dim`)
2. Li+ concentration-driven LJP term
3. EC concentration-driven LJP term

## 3) Dimensional Li+ flux equation (electrolyte_diffusion/full_diffusion)

### 3.1 Code decomposition

The model splits Li+ flux into:

- `N_e_diffusion = - tor * D_e * grad(c_e)`
- `N_e_migration = (C_e/gamma_e) * t_plus * i_e`
- `N_cross_diffusion = - gamma_e_ec * tau_diffusion_e/tau_cross * tor * D_Li_ec_cross * grad(c_EC)`
- `N_e_convection = C_e * c_e * v_box`

and `N_e = N_e_diffusion + N_cross_diffusion + N_e_migration + N_e_convection`.

### 3.2 Dimensionalized form

The dimensional Li+ flux becomes:

`N_e_dim = - tor * D_e_dim(c_e_dim, c_EC_dim, T_dim) ∇c_e_dim`
`          - tor * D_Li_ec_cross_dim(c_e_dim, c_EC_dim, T_dim) ∇c_EC_dim`
`          + t_plus(c_e_dim, c_EC_dim, T_dim) * i_e_dim / F`
`          + N_conv,dim`

where:

- `N_e_dim` in mol.m^-2.s^-1
- `D_e_dim`, `D_Li_ec_cross_dim` in m^2.s^-1

`N_conv,dim` is the convection contribution associated with the chosen volume-averaged velocity variable in this model branch.

## 4) Dimensional EC flux equation (solvent_diffusion/Double_SpatialConsume_wo_refill)

### 4.1 Code decomposition

The EC flux is written as:

- `N_EC_diffusion = - tor * D_ec * grad(c_EC)`
- `N_cross_diffusion = - (tau_ec/tau_cross/gamma_e_ec) * tor * D_ec_Li_cross * grad(c_e)`
- `N_EC_migration = C_RA_typ * (c_EC/c_tot) * Xi * i_e`

and `N_EC = N_EC_diffusion + N_cross_diffusion + N_EC_migration`.

### 4.2 Dimensionalized form

The first two terms reduce directly to:

`N_EC,dim = - tor * D_ec_dim(c_e_dim, c_EC_dim, T_dim) ∇c_EC_dim`
`           - tor * D_ec_Li_cross_dim(c_e_dim, c_EC_dim, T_dim) ∇c_e_dim`
`           + N_EC,mig,dim`

For migration, code scaling gives:

`N_EC,mig,dim = Xi(c_e_dim, c_EC_dim, T_dim) * [ c_EC_dim / c_tot_dim(c_e_dim, c_EC_dim, T_dim) ] * i_e_dim / F`

This has consistent units mol.m^-2.s^-1.

## 5) Dimensional conservation PDEs

From the coded RHS structures (`eps*c` state variables):

### 5.1 Li+ balance

`∂(epsilon * c_e_dim)/∂t + ∇·N_e_dim = S_Li,dim + S_refill,dim - c_e_dim * (∇·v_box_dim)`

(`v_box` term appears explicitly in the current code branch.)

### 5.2 EC balance (double spatial consume wo refill branch)

`∂(epsilon * c_EC_dim)/∂t + ∇·N_EC,dim = S_EC,SEI,dim + S_EC,refill,dim`

In this specific `wo_refill` variant, refill is set to zero in the implemented branch.

## 6) Implication for migration to modern dimensional PyBaMM

To port into modern PyBaMM, keep these equations directly in dimensional form and remove the legacy dependencies on:

- `gamma_e`, `C_e`, `tau_*` prefactors inside governing equations
- explicit `potential_scale` and concentration scaling factors in PDE terms

Those legacy factors should only appear (if needed) in helper post-processing, not in the physical PDE definitions.

## 7) Word equation cross-check (docx extraction enabled)

I added a local extractor script:

- `wip/Rio_Code/P3_R14/extract_docx_equations.py`

It reads `word/document.xml` and extracts:

- paragraph text
- OMML equation text (`m:t`) as numbered `[EQxxxx]` entries

Generated extracted files:

- `wip/Rio_Code/P3_R14/docx_extracted/key equations.extracted.txt`
- `wip/Rio_Code/P3_R14/docx_extracted/Round_240126_Non_dimensionalization - new.extracted.txt`
- `wip/Rio_Code/P3_R14/docx_extracted/Run model of P3 in PyBaMM.extracted.txt`

## 8) Cross-check results against manuscript equations

From `Round_240126_Non_dimensionalization - new` extraction:

- Li+ flux form explicitly contains diffusion, cross-diffusion and migration:
  - `[EQ0001] Ne,k* = ... - De,EC* d(c_EC*)/dx* + t+0/F* ie,k*`
- EC flux form explicitly contains cross term and migration:
  - `[EQ0002] NEC,k* = ... - DEC,e* d(ce*)/dx* + (cEC*/cT*) Xi /F * ie,k*`
- Cross diffusivity scaling was intentionally revised:
  - `[EQ0003]..[EQ0005]`
- `C_RA,typ` definition is documented:
  - `[EQ0008] CRA,typ = Lx*I* / (F*cT,typ*DEC,typ*)` (symbol text spacing follows extraction)

These extracted forms are consistent with the code-level dimensionalization in Sections 2-4:

- conductivity equation includes two measured LJP gradient terms;
- Li+ and EC fluxes each include their own cross-diffusion coupling;
- EC migration term is proportional to `(c_EC/c_tot) * Xi * i_e / F`.

## 9) Remaining gap for paper-grade alignment

The extractor captures equation text tokens but cannot fully preserve Word visual structure (fractions, superscripts, stacked operators) in perfect typeset fidelity.

For final publication-grade reconciliation, the next step is:

1. Use these extracted equation IDs as anchors.
2. Manually verify each anchor against the original rendered Word equation.
3. Finalize one canonical dimensional equation set and then map to modern PyBaMM code symbols.

