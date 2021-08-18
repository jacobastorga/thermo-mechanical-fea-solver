import numpy as np
from typing import Dict, Any

class ThermoElasticStressEvaluator:
    """Evaluates thermal expansion strain epsilon_th = alpha * delta_T and thermo-mechanical stress."""

    def __init__(self, youngs_modulus: float = 130e9, poissons_ratio: float = 0.28,
                 cte_alpha: float = 2.6e-6, t_reference: float = 298.15):
        self.e = youngs_modulus  # Pa (e.g. Silicon: 130 GPa)
        self.nu = poissons_ratio
        self.alpha = cte_alpha   # 1/K (e.g. Silicon: 2.6e-6 /K)
        self.t_ref = t_reference

    def calculate_thermal_stress_field(self, temp_field: np.ndarray) -> Dict[str, np.ndarray]:
        """Calculates thermal strain and constrained thermal stress sigma_th = -E * alpha * delta_T / (1 - 2*nu)."""
        delta_t = temp_field - self.t_ref
        thermal_strain = self.alpha * delta_t

        # Plane strain thermal stress tensor components (biaxial constraint)
        factor = (self.e * self.alpha * delta_t) / (1.0 - 2.0 * self.nu)
        sigma_xx = -factor
        sigma_yy = -factor
        
        # Von Mises equivalent stress under biaxial thermal load
        von_mises = np.abs(factor)

        return {
            "delta_t": delta_t,
            "thermal_strain": thermal_strain,
            "sigma_xx_pa": sigma_xx,
            "sigma_yy_pa": sigma_yy,
            "von_mises_pa": von_mises,
            "max_von_mises_mpa": float(np.max(von_mises) / 1e6)
        }
