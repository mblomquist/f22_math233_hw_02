//
// Created by mblomquist on 10/13/22.
//

#include "LevelSet.h"

LevelSet::LevelSet(Grid2d & newGrid, std::vector<double> & phi_0) {

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

void LevelSet::advance_sl(std::vector<double> & vel_u, std::vector<double> & vel_v, double dt) {

    SL_method sl_solver;
    sl_solver.set_grid(grid);
    sl_solver.set_velocity(vel_u, vel_v);
    sl_solver.update_sol(phi, dt);

}

void LevelSet::reinitialize(const std::vector<double> &phi_0, std::vector<double> &phi_n, std::vector<double> &phi_np1) {

    for(int i = 0; i < grid.get_N(); i++){
        for(int j = 0; j < grid.get_M(); j++){

            double dx = grid.get_dx();
            double dy = grid.get_dy();
            double dpx_m, dpx_p, dpy_m, dpy_p;

            double p_0 = phi_0[grid.n_from_ij(i,j)];
            double p_n = phi_n[grid.n_from_ij(i,j)];

            double Sp_0 = sgn(p_0);

            // compute one sided derivatives
            dpx_m = (phi_n[grid.n_from_ij(i,j)] - phi_n[grid.n_from_ij(i-1,j)]) / dx;
            dpx_p = (phi_n[grid.n_from_ij(i+1,j)] - phi_n[grid.n_from_ij(i,j)]) / dx;
            dpy_m = (phi_n[grid.n_from_ij(i,j)] - phi_n[grid.n_from_ij(i,j-1)]) / dy;
            dpy_p = (phi_n[grid.n_from_ij(i,j+1)] - phi_n[grid.n_from_ij(i,j)]) / dy;

            if (i == 0){
                dpx_m = 0.;
            }
            if (i == grid.get_N()){
                dpx_p = 0.;
            }
            if (j == 0){
                dpy_m = 0.;
            }
            if (j == grid.get_M()){
                dpy_p = 0.;
            }

            double dt = min(abs(dpx_p), abs(dpx_m));
            dt = min(dt, abs(dpy_p));
            dt = min(dt, abs(dpy_m));
            dt = 0.2 * dt;

            // Godunov scheme
            if (Sp_0 > 0.) {
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

            p_n = p_n - dt * Sp_0 * (sqroot(max(dpx_p*dpx_p,dpx_m*dpx_m) + max(dpy_p*dpy_p,dpy_m*dpy_m)) - 1.);

            if (p_n * p_0 < 0.)
                phi_np1[grid.n_from_ij(i,j)] = -p_n;
            else
                phi_np1[grid.n_from_ij(i,j)] = p_n;
        }
    }
}