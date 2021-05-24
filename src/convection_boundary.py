import numpy as np

class RobinConvectionBoundary:
    """Applies Robin boundary condition -k*(dT/dn) = h*(T_surface - T_ambient)."""

    def __init__(self, h_coeff: float = 25.0, t_ambient: float = 298.15):
        self.h_coeff = h_coeff
        self.t_ambient = t_ambient

    def compute_boundary_flux(self, surface_temp: float) -> float:
        """Calculates outward convective heat flux q = h * (T_s - T_inf)."""
        return float(self.h_coeff * (surface_temp - self.t_ambient))

    def compute_nusselt_number(self, characteristic_len: float, k_fluid: float = 0.026) -> float:
        """Calculates dimensionless Nusselt number Nu = h * L / k."""
        return float((self.h_coeff * characteristic_len) / k_fluid)
