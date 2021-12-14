import os
import sys
import numpy as np

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from src.thermal_fdm import SteadyStateThermalSolver2D
from src.thermo_elastic_stress import ThermoElasticStressEvaluator

def main():
    print("Simulating 2D Silicon Die Thermal & Stress Profile in Semiconductor Package...")
    solver = SteadyStateThermalSolver2D(nx=31, ny=31, lx=0.015, ly=0.015, k_thermal=148.0)
    # Heat generation: 5W over chip area -> Q = 2.2e7 W/m^3
    t_field = solver.solve(top_temp=350.15, bottom_temp=310.15, left_temp=325.15, right_temp=325.15, internal_heat_gen=2e6)
    
    stress_eval = ThermoElasticStressEvaluator(youngs_modulus=130e9, poissons_ratio=0.28, cte_alpha=2.6e-6, t_reference=298.15)
    stress_res = stress_eval.calculate_thermal_stress_field(t_field)

    print(f"Max Package Temperature: {np.max(t_field) - 273.15:.2f} C")
    print(f"Min Package Temperature: {np.min(t_field) - 273.15:.2f} C")
    print(f"Max Equivalent Von Mises Stress: {stress_res['max_von_mises_mpa']:.2f} MPa")

if __name__ == "__main__":
    main()
