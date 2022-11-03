//
// Created by mblomquist on 10/13/22.
//

#include "LevelSet.h"
#include <stdlib.h>

LevelSet::LevelSet(Grid2d &newGrid, std::vector<double> &phi_0) {

    grid = newGrid;
    phi.resize(phi_0.size());
    phi = phi_0;

}

void LevelSet::setGrid(Grid2d &newGrid) {
    grid = newGrid;
}

void LevelSet::setPhi(std::vector<double> &phi_0) {
    phi.resize(phi_0.size());
    phi = phi_0;
}

void LevelSet::getGrid(Grid2d &grid_o) {
    grid_o = grid;
}

void LevelSet::getPhi(std::vector<double> &phi_o) {
    phi_o = phi;
}

void LevelSet::advance_sl(std::vector<double> &vel_u, std::vector<double> &vel_v, double dt) {

    SL_method sl_solver;
    sl_solver.set_grid(grid);
    sl_solver.set_velocity(vel_u, vel_v);
    sl_solver.update_sol(phi, dt);

}

void
LevelSet::reinitialize(const std::vector<double> &phi_0, std::vector<double> &phi_n, std::vector<double> &phi_np1) {

    for (int i = 0; i < grid.get_N(); i++) {
        for (int j = 0; j < grid.get_M(); j++) {

            double dx = grid.get_dx();
            double dy = grid.get_dy();
            double dpx_m, dpx_p, dpy_m, dpy_p;

            double p_0 = phi_0[grid.n_from_ij(i, j)];
            double p_n = phi_n[grid.n_from_ij(i, j)];

            if (abs(p_n) < 1.e-1 * dx) {
                phi_np1[grid.n_from_ij(i, j)] = 0.;
            } else {

                double p_l = phi_n[grid.n_from_ij(i - 1, j)];
                double p_r = phi_n[grid.n_from_ij(i + 1, j)];
                double p_b = phi_n[grid.n_from_ij(i, j - 1)];
                double p_t = phi_n[grid.n_from_ij(i, j + 1)];

                double Sp_0 = sgn(p_0);

                if (p_n * p_l < 0.) dpx_m = p_n / abs(p_0);
                else dpx_m = (p_n - p_l) / dx;

                if (p_r * p_n < 0.) dpx_p = -p_n / abs(p_0);
                else dpx_p = (p_r - p_n) / dx;

                if (p_n * p_b < 0.) dpy_m = p_n / abs(p_0);
                else dpy_m = (p_n - p_b) / dy;

                if (p_t * p_n < 0.) dpy_p = -p_n / abs(p_0);
                else dpy_p = (p_t - p_n) / dy;

                if (i == 0) {
                    dpx_m = dpx_p;
                }
                if (i == grid.get_N()) {
                    dpx_p = dpx_m;
                }
                if (j == 0) {
                    dpy_m = dpy_p;
                }
                if (j == grid.get_M()) {
                    dpy_p = dpy_m;
                }

                double dt = 0.5 * dx;

                // Godunov scheme
                if (p_0 > 0.) {
                    if (dpx_p > 0.) dpx_p = 0.;
                    if (dpx_m < 0.) dpx_m = 0.;
                    if (dpy_p > 0.) dpy_p = 0.;
                    if (dpy_m < 0.) dpy_m = 0.;
                } else {
                    if (dpx_p < 0.) dpx_p = 0.;
                    if (dpx_m > 0.) dpx_m = 0.;
                    if (dpy_p < 0.) dpy_p = 0.;
                    if (dpy_m > 0.) dpy_m = 0.;
                }

                p_n = p_n -
                      dt * Sp_0 * (sqroot(max(dpx_p * dpx_p, dpx_m * dpx_m) + max(dpy_p * dpy_p, dpy_m * dpy_m)) - 1.);

                if (p_n * p_0 < 0.)
                    phi_np1[grid.n_from_ij(i, j)] = -p_n;
                else
                    phi_np1[grid.n_from_ij(i, j)] = p_n;

            }
        }
    }
}

void LevelSet::perturb_ls(std::vector<double> &phi_p, double eps) {

    for (int i = 0; i < grid.get_N(); i++){
        for (int j = 0; j < grid.get_M(); j++){

            double& p = phi_p[grid.n_from_ij(i,j)];

            if (abs(p) > eps){
                if (p > 0.)
                    p += 0.01 * (rand() % 10);
                if (p < 0.)
                    p -= 0.01 * (rand() % 10);
            }
        }
    }
}