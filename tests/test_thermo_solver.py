import unittest
import numpy as np
from src.thermal_fdm import SteadyStateThermalSolver2D
from src.convection_boundary import RobinConvectionBoundary
from src.thermo_elastic_stress import ThermoElasticStressEvaluator

class TestThermoMechanicalSolver(unittest.TestCase):

    def test_linear_heat_conduction(self):
        solver = SteadyStateThermalSolver2D(nx=21, ny=21, lx=0.02, ly=0.02, k_thermal=150.0)
        t_field = solver.solve(top_temp=373.15, bottom_temp=293.15, left_temp=333.15, right_temp=333.15)
        self.assertEqual(t_field.shape, (21, 21))
        # Center should be bounded between bottom and top temp
        center_t = t_field[10, 10]
        self.assertGreater(center_t, 293.15)
        self.assertLess(center_t, 373.15)

    def test_convective_heat_flux(self):
        conv = RobinConvectionBoundary(h_coeff=50.0, t_ambient=295.0)
        flux = conv.compute_boundary_flux(surface_temp=345.0)
        self.assertEqual(flux, 50.0 * (345.0 - 295.0))
        nu = conv.compute_nusselt_number(characteristic_len=0.05, k_fluid=0.025)
        self.assertEqual(nu, 100.0)

    def test_thermal_stress_evaluation(self):
        evaluator = ThermoElasticStressEvaluator(youngs_modulus=100e9, poissons_ratio=0.25, cte_alpha=1e-5, t_reference=300.0)
        temp_field = np.ones((5, 5)) * 350.0  # +50 K delta
        res = evaluator.calculate_thermal_stress_field(temp_field)
        self.assertGreater(res["max_von_mises_mpa"], 0.0)
        self.assertAlmostEqual(res["max_von_mises_mpa"], 100.0, places=1)

if __name__ == "__main__":
    unittest.main()
