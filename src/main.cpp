#include <iostream>
#include "grid/Grid2d.h"
#include "advection/SL_method.h"
#include "tools/math_tools.h"
#include "math.h"

double cf_phi(double x, double y){
    return pow(pow((x - 0.25),2)+ pow(y,2), 0.5) - 0.2;
}

void cf_vel0(double x, double y, double & vel_u, double & vel_v){
    vel_u = -y;
    vel_v = x;
}

int main() {
    std::cout << "Hello, F22 MATH 233 Homework 2!" << std::endl;

    // Problem 1 - Semi-Lagrangian Method for the advection equation

    double xmin = -1.0;
    double xmax = 1.0;
    double ymin = -1.0;
    double ymax = 1.0;

    int N = 100;
    int M = N;

    Grid2d newGrid(N, M, xmin, xmax, ymin, ymax);

    double dx = newGrid.get_dx();

    newGrid.print_VTK_format("../results/f22_hw02_problem1_0.vtk");

    std::vector<double> vel_u, vel_v, phi;
    vel_u.resize(N*M);
    vel_v.resize(N*M);
    phi.resize(N*M);

    for (int i = 0; i < N; i++){
        for (int j = 0; j < M; j++){
            cf_vel0(newGrid.x_from_n(newGrid.n_from_ij(i,j)),
                    newGrid.y_from_n(newGrid.n_from_ij(i,j)),
                    vel_u[newGrid.n_from_ij(i,j)],
                    vel_v[newGrid.n_from_ij(i,j)]);

            phi[newGrid.n_from_ij(i,j)] = cf_phi(newGrid.x_from_n(newGrid.n_from_ij(i,j)),
                                                 newGrid.y_from_n(newGrid.n_from_ij(i,j)));
        }
    }

    newGrid.print_VTK_format(vel_u, "vel_u", "../results/f22_hw02_problem1_0.vtk");
    newGrid.print_VTK_format(vel_v, "vel_v", "../results/f22_hw02_problem1_0.vtk");
    newGrid.print_VTK_format(phi, "phi", "../results/f22_hw02_problem1_0.vtk");

    double cfl = 0.5;
    double dt = dx/cfl;

    SL_method sl_solver;
    sl_solver.set_grid(newGrid);
    sl_solver.set_velocity(vel_u, vel_v);

    int steps = 314;
    int pr, mprint;

    mprint = 10;
    pr = 0;

    for (int i = 0; i < steps; i++){
        sl_solver.update_sol(phi, dt);

        if (i % mprint == 0){
            pr++;

            std::string fileName = "../results/f22_hw02_problem1_";
            std::string fileType = ".vtk";

            fileName = fileName + std::to_string(pr);

            fileName.append(fileType);

            std::cout << fileName << std::endl;

            newGrid.print_VTK_format(fileName);
            newGrid.print_VTK_format(vel_u, "vel_u", fileName);
            newGrid.print_VTK_format(vel_v, "vel_v", fileName);
            newGrid.print_VTK_format(phi, "phi", fileName);
        }
    }

    // Problem 2 - Reinitialization Equation


    // Problem 3 - Level Set Method


    // Problem 4 - Extra Credit


    return 0;
}
