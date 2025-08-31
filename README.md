# 2D Incompressible LES Solver – Lid-Driven Cavity

This project implements a **2D incompressible Navier–Stokes solver** in C for the **lid-driven cavity problem**, a classical benchmark in Computational Fluid Dynamics (CFD).  
It uses **Large Eddy Simulation (LES)** with the **Smagorinsky model** to capture turbulent effects at high Reynolds numbers.

---

## ✨ Features
- **Finite Difference Method** on a uniform 2D grid.
- **Projection Method** for enforcing incompressibility:
  - Predictor step: explicit advection + variable viscosity diffusion.
  - Pressure correction: iterative Poisson solve via **Successive Over-Relaxation (SOR)**.
- **Turbulence Modeling**:
  - **Smagorinsky LES** model for eddy viscosity.
- **Configurable Parameters**:
  - Grid resolution (Nx, Ny).
  - Reynolds number.
  - Lid velocity.
  - Time step size.
  - Smagorinsky constant.
- **Data Output**:
  - Periodically saves `CSV` files containing velocity and pressure fields for visualization.

---

## 📂 Project Structure
- **2d-cavity.c** – main C source code.
- **fields_xxxxxx.csv** – simulation snapshots (u, v, p) saved during runtime.

---

## ⚙️ How It Works
1. **Boundary Conditions**: 
   - No-slip on walls.
   - Moving lid at the top boundary.
2. **LES Eddy Viscosity**: 
   - Computed using the Smagorinsky model.
3. **Intermediate Velocities**: 
   - Calculated with explicit advection + diffusion.
4. **Pressure Poisson Equation**: 
   - Solved with SOR for enforcing incompressibility.
5. **Projection Step**: 
   - Velocity field corrected with pressure gradient.
6. **Output**: 
   - Velocity (u, v) and pressure (p) stored for post-processing.

---

## 🚀 Compilation & Run
Compile with GCC:
```bash
gcc -O2 -o cavity 2d-cavity.c -lm
