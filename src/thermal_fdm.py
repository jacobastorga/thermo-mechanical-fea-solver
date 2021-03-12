import numpy as np
from typing import Tuple, Dict, Any

class SteadyStateThermalSolver2D:
    """Solves 2D steady-state heat equation k*(d^2T/dx^2 + d^2T/dy^2) + Q = 0."""

    def __init__(self, nx: int, ny: int, lx: float, ly: float, k_thermal: float = 150.0):
        self.nx = nx
        self.ny = ny
        self.lx = lx
        self.ly = ly
        self.dx = lx / (nx - 1)
        self.dy = ly / (ny - 1)
        self.k_thermal = k_thermal

    def solve(self, top_temp: float, bottom_temp: float, left_temp: float, right_temp: float,
              internal_heat_gen: float = 0.0, max_iter: int = 5000, tol: float = 1e-5) -> np.ndarray:
        """Solves steady-state temperature field via successive over-relaxation (SOR)."""
        t = np.ones((self.ny, self.nx), dtype=np.float64) * ((top_temp + bottom_temp) / 2.0)
        t[-1, :] = top_temp
        t[0, :] = bottom_temp
        t[:, 0] = left_temp
        t[:, -1] = right_temp

        dx2 = self.dx ** 2
        dy2 = self.dy ** 2
        factor = 2.0 * (dx2 + dy2)
        q_term = (internal_heat_gen / self.k_thermal) * dx2 * dy2

        omega = 1.7  # Relaxation parameter
        for _ in range(max_iter):
            t_old = t.copy()
            for j in range(1, self.ny - 1):
                for i in range(1, self.nx - 1):
                    rhs = dy2 * (t[j, i+1] + t[j, i-1]) + dx2 * (t[j+1, i] + t[j-1, i]) + q_term
                    t_new = rhs / factor
                    t[j, i] = (1.0 - omega) * t[j, i] + omega * t_new

            diff = np.max(np.abs(t - t_old))
            if diff < tol:
                break

        return t
