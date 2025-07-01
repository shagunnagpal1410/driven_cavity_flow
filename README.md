# 💨 2D Lid-Driven Cavity Flow Simulation using Finite Pointset Method (FPM)

This project implements a **2D incompressible Navier–Stokes solver** for the classic **lid-driven cavity flow** problem using the **Finite Pointset Method (FPM)**. It simulates velocity evolution over time using a meshfree approach based on local MLS approximations and Gaussian weights.

---

## ⚙️ Running the Code

You can run the code directly using:

```
Ctrl + Shift + B
```

> Make sure your `tasks.json` is properly set up and includes compilation and execution commands for `driven_cavity_flow.cpp`.

---

## 📦 Output

- Generates a separate CSV file for each time step:
  ```
  Velocity1.csv
  Velocity2.csv
  ...
  ```
- Each file contains:
  - `X`, `Y` → coordinates
  - `u`, `v` → velocity components

---

## 📐 Model Setup

- Domain: Unit square \([0,1] \times [0,1]\)
- Top lid moves with a **parabolic velocity profile**
- All other walls are **stationary (no-slip)**
- Gaussian kernel used for spatial interpolation
- MLS-based approximations for derivatives

---

## 🔄 Method Overview

1. Discretizes the 2D domain into a **point cloud**
2. Applies **Chorin projection method**:
   - Computes intermediate velocities (\(u^*, v^*\))
   - Solves pressure Poisson equation
   - Corrects velocity to ensure incompressibility
3. Applies Dirichlet boundary conditions
4. Writes velocity field to CSV after each time step

---

## 🧪 Applications

- Educational demonstration of meshfree methods (FPM)
- Benchmarking against standard CFD solvers
- Research in meshless incompressible flow modeling

---

## 📁 Dependencies

- [Eigen](https://eigen.tuxfamily.org/) for sparse matrix operations

> Make sure Eigen is accessible to your build system (your `tasks.json` should include the correct include path).

