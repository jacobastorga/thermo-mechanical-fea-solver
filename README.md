# Thermo-Mechanical Finite Difference & FEA Solver

[![Python](https://img.shields.io/badge/python-3.8%20%7C%203.9-blue.svg)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![FEA/Thermal](https://img.shields.io/badge/simulation-thermo--elastic-orange.svg)]()
[![Tests](https://img.shields.io/badge/tests-passing-brightgreen.svg)]()

A 2D finite difference and thermo-elastic stress solver for electronic semiconductor packages, heat sinks, and structural components subjected to internal volumetric heat dissipation.

---

## Technical Capabilities

1. **2D Steady-State Heat Conduction:** Solves Fourier's law with internal heat generation $k 
abla^2 T + Q = 0$ using Successive Over-Relaxation (SOR).
2. **Convective Boundary Conditions:** Implements mixed Dirichlet and Robin boundaries $-k rac{\partial T}{\partial n} = h (T_s - T_\infty)$.
3. **Thermo-Elastic Stress Tensor:** Calculates thermal strain $\epsilon_{	ext{th}} = lpha \Delta T$ and Von Mises equivalent stress under plane strain constraints.

---

## Quickstart

```bash
git clone https://github.com/jacobastorga/thermo-mechanical-fea-solver.git
cd thermo-mechanical-fea-solver
pip install -r requirements.txt
pytest tests/ -v
python examples/simulate_package.py
```

---

## Author
**Jacob Astorga**  
Mechanical & Systems Engineer  
Portfolio: [github.com/jacobastorga](https://github.com/jacobastorga)
