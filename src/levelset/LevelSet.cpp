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

void LevelSet::reinitialize(std::vector<double> &phi_0, std::vector<double> phi_n, std::vector<double> &phi_np1) {

}