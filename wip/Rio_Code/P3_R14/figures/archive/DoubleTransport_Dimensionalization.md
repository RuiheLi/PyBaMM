# Double-transport model dimensional formulation (from current ECdrag2 code)

## Purpose

This note rewrites the core equations in the current ECdrag2 implementation from code-level nondimensional form into dimensional form.

Code sources used:

- `pybamm/models/submodels/electrolyte_conductivity/sol_full_conductivity.py`
- `pybamm/models/submodels/electrolyte_diffusion/full_diffusion.py`
- `pybamm/models/submodels/solvent_diffusion/Double_SpatialConsume_wo_refill.py`
- `pybamm/parameters/lithium_ion_parameters.py`

## 1) Nondimensional variables and scales used by current code

From `LithiumIonParameters`:

- \(x = L_x \tilde{x}\)
- \(t = \tau \tilde{t}\), with default \(\tau = \tau_{\mathrm{discharge}}\)
- \(c_e = c_{e,\mathrm{typ}} \tilde{c}_e\)
- \(c_{EC} = c_{EC,\mathrm{typ}} \tilde{c}_{EC}\)
- \(\phi_e = \phi_{\mathrm{scale}} \tilde{\phi}_e\), with \(\phi_{\mathrm{scale}} = RT_{ref}/F\)
- \(i_e = i_{\mathrm{typ}} \tilde{i}_e\)
- \(N_{Li} = \dfrac{D_{e,\mathrm{typ}} c_{e,\mathrm{typ}}}{L_x}\tilde{N}_{Li}\)
- \(N_{EC} = \dfrac{D_{EC,\mathrm{typ}} c_{EC,\mathrm{typ}}}{L_x}\tilde{N}_{EC}\)

Important grouped parameters:

- \(C_e = \tau_{\mathrm{diff},e}/\tau\), \(\tau_{\mathrm{diff},e}=L_x^2/D_{e,\mathrm{typ}}\)
- \(\gamma_e = (\tau_{\mathrm{discharge}}/\tau)\, c_{e,\mathrm{typ}}/c_{max}\)
- \(\gamma_{e,EC} = c_{EC,\mathrm{typ}}/c_{e,\mathrm{typ}}\)
- \(\tau_{EC}=L_x^2/D_{EC,\mathrm{typ}}\)
- \(\tau_{\mathrm{cross}}=L_x^2/D_{EC,Li,\mathrm{cross,typ}}\)
- \(C_{RA,\mathrm{typ}} = \dfrac{L_x i_{\mathrm{typ}}}{F D_{EC,\mathrm{typ}} c_{\mathrm{tot,typ}}}\)

Useful identity from code definitions:

\[
\frac{\gamma_e}{C_e}
=
\frac{F c_{e,\mathrm{typ}} D_{e,\mathrm{typ}}}{i_{\mathrm{typ}} L_x}
\]

## 2) Electrolyte current equation (`sol_full_conductivity.py`)

### 2.1 Nondimensional code form

\[
\tilde{i}_e
=
\tilde{\kappa}_e\,\mathrm{tor}\,\frac{\gamma_e}{C_e}
\left(
-\nabla_{\tilde{x}}\tilde{\phi}_e
+
\nabla_{\tilde{x}}\tilde{c}_e \,
\frac{\partial U_{LJP}}{\partial c_e}\,
\frac{c_{e,\mathrm{typ}}}{\phi_{\mathrm{scale}}}
+
\nabla_{\tilde{x}}\tilde{c}_{EC} \,
\frac{\partial U_{LJP}}{\partial c_{EC}}\,
\frac{c_{EC,\mathrm{typ}}}{\phi_{\mathrm{scale}}}
\right)
\]

where:

- `dLJP_dce` = \(\partial U_{LJP}/\partial c_e\)
- `dLJP_dcEC` = \(\partial U_{LJP}/\partial c_{EC}\)

### 2.2 Dimensional form

After substituting all scales:

\[
i_e
=
\kappa_e(c_e,c_{EC},T)\,\mathrm{tor}
\left(
-\nabla \Phi_e
+
\frac{\partial U_{LJP}}{\partial c_e}\nabla c_e
+
\frac{\partial U_{LJP}}{\partial c_{EC}}\nabla c_{EC}
\right)
\]

This is the dimensional electrolyte Ohm/Stefan-Maxwell form with measured liquid-junction-potential gradients replacing classical \( \chi RT/Fc \) closure.

## 3) Lithium-ion electrolyte flux (`electrolyte_diffusion/full_diffusion.py`)

### 3.1 Nondimensional code form

\[
\tilde{N}_{Li}
=
\underbrace{-\mathrm{tor}\,\tilde{D}_e\nabla_{\tilde{x}}\tilde{c}_e}_{\text{diffusion}}
+
\underbrace{\frac{C_e}{\gamma_e}\,t_+\,\tilde{i}_e}_{\text{migration}}
+
\underbrace{\tilde{N}_{Li,EC\ cross}}_{\text{EC cross-diffusion}}
+
\underbrace{C_e \tilde{c}_e \tilde{v}_{box}}_{\text{convection}}
\]

Double-solvent mode cross term:

\[
\tilde{N}_{Li,EC\ cross}
=
-\gamma_{e,EC}\frac{\tau_{\mathrm{diff},e}}{\tau_{\mathrm{cross}}}
\mathrm{tor}\,\tilde{D}_{Li,EC,\mathrm{cross}}
\nabla_{\tilde{x}}\tilde{c}_{EC}
\]

### 3.2 Dimensional form

Dimensionalized term-by-term:

\[
N_{Li}
=
-\mathrm{tor}\,D_e(c_e,c_{EC},T)\nabla c_e
+
t_+(c_e,c_{EC},T)\frac{i_e}{F}
+
N_{Li,EC\ cross}
+
c_e\,v_{box}
\]

\[
N_{Li,EC\ cross}
=
-\mathrm{tor}\,D_{Li,EC,\mathrm{cross}}(c_e,c_{EC},T)\nabla c_{EC}
\]

So the code corresponds to a dimensional Stefan-Maxwell-like flux with explicit cross-diffusion by solvent concentration gradient.

## 4) EC flux (`solvent_diffusion/Double_SpatialConsume_wo_refill.py`)

### 4.1 Nondimensional code form

\[
\tilde{N}_{EC}
=
\underbrace{-\mathrm{tor}\,\tilde{D}_{EC}\nabla_{\tilde{x}}\tilde{c}_{EC}}_{\text{EC diffusion}}
+
\underbrace{-\frac{\tau_{EC}}{\tau_{\mathrm{cross}}}\frac{1}{\gamma_{e,EC}}
\mathrm{tor}\,\tilde{D}_{EC,Li,\mathrm{cross}}\nabla_{\tilde{x}}\tilde{c}_{e}}_{\text{Li gradient cross-diffusion}}
+
\underbrace{C_{RA,\mathrm{typ}}
\frac{\tilde{c}_{EC}}{\tilde{c}_{tot}}
\Xi\,\tilde{i}_e}_{\text{migration}}
\]

### 4.2 Dimensional form

\[
N_{EC}
=
-\mathrm{tor}\,D_{EC}(c_e,c_{EC},T)\nabla c_{EC}
-\mathrm{tor}\,D_{EC,Li,\mathrm{cross}}(c_e,c_{EC},T)\nabla c_e
+
\Xi(c_e,c_{EC},T)\frac{c_{EC}}{c_{tot}}\,\frac{i_e}{F}
\]

This is a two-species coupled transport structure (EC + Li+) with both cross-diffusion and migration.

## 5) PDE forms in dimensional notation

Using porosity-weighted concentrations:

\[
\frac{\partial (\epsilon c_e)}{\partial t}
=
-\nabla\cdot N_{Li}
+ S_{Li}
- c_e\,\nabla\cdot v_{box}
+ S_{refill}
\]

\[
\frac{\partial (\epsilon c_{EC})}{\partial t}
=
-\nabla\cdot N_{EC}
+ S_{EC}
\]

where source terms are driven mainly by SEI side reaction terms in your implementation (`j_inner`, `j_outer`, stoichiometric ratios, and refill option switches).

## 6) Practical migration guidance (to dimensional-only modern PyBaMM)

For a dimensional-only target, keep the physics but drop explicit nondimensional groups in equation assembly:

1. Replace `param.gamma_e / param.C_e` prefactor in current equation by direct dimensional form shown above.
2. Replace `param.C_e / param.gamma_e * i_e` in Li migration term by `i_e / F`.
3. Replace EC migration prefactor `C_RA_typ * (...) * i_e` by \((c_{EC}/c_{tot})\Xi\,i_e/F\).
4. Keep `D_e`, `kappa_e`, `D_EC`, cross-diffusion and LJP-derivative closures as dimensional function parameters of \((c_e, c_{EC}, T)\).

## 7) Note about paper Word files

The current toolchain in this session cannot directly parse `.docx` binary files.
If you export the key Word equations to PDF or plain text, this note can be aligned line-by-line against the manuscript equations and updated with exact equation numbering.
