# 1D Semiconductor Heterojunction Device Simulator

A 1D numerical semiconductor device simulator written in Python to model electrostatic potential, electric fields, charge density, and carrier concentrations ($n, p$) across multi-material heterojunctions under applied bias voltages.

The simulator solves the **coupled non-linear Poisson equation** using **Finite Difference Discretization (FDM)** and the **Newton-Raphson Method**.

---

## ðŸ“‹ Features

- **Multi-Material Support:** Pre-configured material database for common semiconductors:
  - Silicon ($\text{Si}$)
  - Germanium ($\text{Ge}$)
  - Gallium Arsenide ($\text{GaAs}$)
  - Aluminum Gallium Arsenide ($\text{AlGaAs}$)
  - Indium Phosphide ($\text{InP}$)
  - Gallium Indium Arsenide ($\text{GaInAs}$)
- **Heterostructure Interfaces:** Handles continuous flux boundary conditions and material parameter discontinuities ($\varepsilon_r, N_c, N_v, E_g$) at material boundaries.
- **Robust Numerical Engine:**
  - Finite Difference spatial grid ($500$ nodes per material layer by default).
  - Newton-Raphson non-linear solver with dynamic tridiagonal Jacobian matrix construction.
  - Convergence controlled by thermal voltage threshold ($\Delta V < V_t$).
- **Comprehensive Visualizations:** Automatically generates Matplotlib figures for:
  - Electrostatic Potential Profile ($V$)
  - Electric Field Distribution ($E$)
  - Linear Carrier Concentration ($n, p$)
  - Logarithmic Carrier Concentration ($\log_{10} n, \log_{10} p$)
  - Net Charge Density ($\rho$)

---

## ðŸ”¬ Physics & Numerical Methods

### 1. Poisson Equation
The core physical model is governed by the 1D non-linear Poisson equation:

$$\frac{d}{dx} \left( \varepsilon(x) \frac{d V(x)}{dx} \right) = -q \left( N_d^+(x) - N_a^-(x) + p(x, V) - n(x, V) \right)$$

Where the carrier concentrations follow Maxwell-Boltzmann statistics:

$$n(x) = N_c \cdot \exp\left( \frac{q V(x) + E_f}{k_B T} \right)$$

$$p(x) = N_v \cdot \exp\left( \frac{-q V(x) - E_g - E_f}{k_B T} \right)$$

### 2. Discretization & Solver
Using a 3-point central finite difference scheme, the continuous spatial domain is discretized into $N$ nodes. The resulting non-linear system $F(V) = 0$ is iteratively solved using the Newton-Raphson updates:

$$J(V^{(k)}) \cdot \Delta V^{(k)} = -F(V^{(k)})$$

$$V^{(k+1)} = V^{(k)} + \Delta V^{(k)}$$

Where $J$ is the tridiagonal Jacobian matrix built dynamically during each iteration.

---

## ðŸ› ï¸ Installation & Setup

### Prerequisites
Make sure you have Python 3.8+ installed along with the required libraries:

```bash
pip install numpy matplotlib
```

## ðŸš€ Usage

Run the main simulator script:

```bash
python main.py
```

### Interactive Prompts
Upon execution, the terminal will prompt for simulation parameters:

1. **Initial Voltage ($V_{\text{initial}}$):** Applied bias at the left terminal (e.g., `0.0`).
2. **Final Voltage ($V_{\text{final}}$):** Applied bias at the right terminal (e.g., `0.7`).
3. **Number of Materials:** Total distinct layers in the heterostructure (e.g., `2` for a $p\text{-GaAs} / n\text{-AlGaAs}$ junction).
4. **Per Material Inputs:**
   - Material Name (`Si`, `Ge`, `GaAs`, `AlGaAs`, `InP`, or `GaInAs`)
   - Layer length in meters (e.g., `500e-9`)
   - Acceptor Concentration $N_a$ in $\text{cm}^{-3}$ (e.g., `1e17`)
   - Donor Concentration $N_d$ in $\text{cm}^{-3}$ (e.g., `0` or `1e16`)

---

## ðŸ“‚ Repository Structure

```text
â”œâ”€â”€ main.py              # Main simulation code (solver, loops, plotting)
â”œâ”€â”€ README.md            # Project documentation
â””â”€â”€ LICENSE              # Project License (e.g., MIT)
```

---

## ðŸ“Š Output Visualizations

The simulator generates five distinct figures upon convergence:

1. **Potential Profile ($V$ vs $x$):** Shows electrostatic potential across layers and junction offset.
2. **Electric Field ($E$ vs $x$):** Peak field located at the heterojunction depletion region.
3. **Carrier Concentrations (Linear):** Plot of $n(x)$ and $p(x)$ across the device.
4. **Carrier Concentrations (Log Scale):** Clear view of minority and majority carrier profiles across several orders of magnitude.
5. **Charge Density ($\rho$ vs $x$):** Space-charge region behavior at heterojunction interfaces.

---

## ðŸ”® Future Improvements

- [ ] Implement Thomas Algorithm (TDMA) for $O(N)$ tridiagonal matrix solving instead of full matrix inversion `np.linalg.inv()`.
- [ ] Incorporate Fermi-Dirac integrals for highly degenerate doping regions.
- [ ] Add transient time-dependent drift-diffusion solver for dynamic I-V characterization.
- [ ] Support custom material input JSON/YAML configuration files instead of interactive CLI prompts.

---

## ðŸ‘¤ Author

**Bipul Biswas**  
M.Tech Student, Electrical Engineering  
IIT Hyderabad  
- **Email:** ee24mtech12003@iith.ac.in  
