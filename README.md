# Design Optimization and Probabilistic Analyses in Structural Dynamics

This repository contains the MATLAB implementations and documentation for the project "Design Optimization and Probabilistic Analyses in Structural Dynamics." The work focuses on optimizing a beam profile under specific constraints using deterministic and probabilistic methods. The repository is structured into three main tasks, each housed in a separate folder, following the project outline and problem requirements.

## Repository Structure

```
Design-Optimization-and-Probabilistic-Analyses-in-Structural-Dynamics/
├── Aufgabe 1 - Deterministic Optimization/
├── Aufgabe 2 - Probabilistic Analysis/
└── Aufgabe 3 - Robustness Optimization/
```

### Task Overview

1. **Deterministic Optimization**  
   The goal is to optimize the deflection of a beam under specific loading and geometric constraints. MATLAB's `fmincon` function is utilized to adjust beam parameters while adhering to boundary conditions. The solution involves calculating derivatives analytically and validating them numerically.

2. **Probabilistic Analysis**  
   Statistical methods, including Monte Carlo Simulation and First Order Second Moment (FOSM) method, are applied to evaluate the variability in beam deflection due to stochastic properties of the elasticity modulus. Distribution fitting and validation using Kolmogorov-Smirnov tests are part of this task.

3. **Robustness Optimization**  
   Combines deterministic and probabilistic approaches to minimize both the deflection and its standard deviation simultaneously. A multi-objective optimization strategy is employed to achieve robust designs under stochastic conditions.

## Project Highlights

- **Deterministic Optimization:**
  - Adjusts beam geometry to minimize center deflection.
  - Analytical derivative calculations validated against finite difference methods.

- **Probabilistic Analysis:**
  - Evaluates deflection variability using statistical techniques.
  - Validates elasticity modulus distributions and compares simulation outcomes.

- **Robustness Optimization:**
  - Simultaneously optimizes mean deflection and its variability.
  - Integrates deterministic optimization with probabilistic insights.

## How to Navigate the Repository

Each folder contains:
- MATLAB scripts and functions required for the task.
- A detailed report summarizing the methodologies and results.
- Plots and visualizations of key findings.

### Folder Details

1. **`Aufgabe 1 - Deterministic Optimization/`**
   - Contains MATLAB code for deterministic optimization.
   - Includes analytical and numerical validation of derivatives.
   - Reference and optimized designs are compared.

2. **`Aufgabe 2 - Probabilistic Analysis/`**
   - Implements Monte Carlo simulations and FOSM.
   - Tests various probability distributions using Kolmogorov-Smirnov.
   - Outputs statistical summaries and plots.

3. **`Aufgabe 3 - Robustness Optimization/`**
   - Integrates deterministic and probabilistic approaches.
   - Multi-objective optimization with weighting factors.
   - Final results compared to deterministic and probabilistic outcomes.

## Prerequisites


- Optimization Toolbox for `fmincon`.

---

For further questions or contributions, feel free to contact me via email: [dogukaan.dere@outlook.de](mailto:dogukaan.dere@outlook.de) or [dogukaan.dere@tuhh.de](mailto:dogukaan.dere@tuhh.de).
