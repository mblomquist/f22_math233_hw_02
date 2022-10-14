#include <iostream>
#include "grid/Grid2d.h"
#include "advection/SL_method.h"
#include "tools/math_tools.h"
#include "levelset/LevelSet.h"
#include "math.h"

#define PI 3.1415926535897932384626433832795028841971

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

    std::vector<double> vel_u, vel_v, phi, phi_0;
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

    phi_0 = phi;

    newGrid.print_VTK_format(vel_u, "vel_u", "../results/f22_hw02_problem1_0.vtk");
    newGrid.print_VTK_format(vel_v, "vel_v", "../results/f22_hw02_problem1_0.vtk");
    newGrid.print_VTK_format(phi, "phi", "../results/f22_hw02_problem1_0.vtk");

    double ratio = 0.1;
    double dt = dx/ratio;

    SL_method sl_solver;
    sl_solver.set_grid(newGrid);
    sl_solver.set_velocity(vel_u, vel_v);

    double t_final = 2.*PI;
    int steps = floor(t_final / dt);
    int pr, mprint;

    mprint = 10;
    pr = 0;

    for (int i = 0; i < steps + 2; i++){

        if (i*dt > t_final){
            double temp_t = (i-1)*dt;
            dt = t_final - (i-1)*dt;
            std::cout << "Time: " << temp_t + dt << std::endl;
        } else {
            std::cout << "Time: " << i*dt << std::endl;
        }

        sl_solver.update_sol(phi, dt);

        if (i % mprint == 0){
            pr++;

            std::string fileName = "../results/f22_hw02_problem1_";
            std::string fileType = ".vtk";

            fileName = fileName + std::to_string(pr);

            fileName.append(fileType);

            //std::cout << fileName << std::endl;

            newGrid.print_VTK_format(fileName);
            newGrid.print_VTK_format(vel_u, "vel_u", fileName);
            newGrid.print_VTK_format(vel_v, "vel_v", fileName);
            newGrid.print_VTK_format(phi, "phi", fileName);
        }
    }


    // Problem 2 - Reinitialization Equation
    double l1, l2, linf;

    l1   = norm_l1(phi, phi_0);
    l2   = norm_l2(phi, phi_0);
    linf = norm_linf(phi, phi_0);

    std::cout << "\nFor ratio: " << ratio << "\n" << std::endl;
    std::cout << "L1 Norm: " << l1 << std::endl;
    std::cout << "L2 Norm: " << l2 << std::endl;
    std::cout << "Linf Norm: " << linf << std::endl;

    // Problem 3 - Level Set Method
    LevelSet ls_phi(newGrid, phi_0);
    std::vector<double> p2_phi;
    p2_phi.resize(phi_0.size());
    ls_phi.getPhi(p2_phi);

    std::string fileName = "../results/f22_hw02_problem2_";
    std::string fileType = ".vtk";

    fileName = fileName + std::to_string(0);
    fileName.append(fileType);

    newGrid.print_VTK_format(fileName);
    newGrid.print_VTK_format(vel_u, "vel_u", fileName);
    newGrid.print_VTK_format(vel_v, "vel_v", fileName);
    newGrid.print_VTK_format(p2_phi, "phi", fileName);

    for (int i = 0; i < 100; i++){
        ls_phi.advance_sl(vel_u, vel_v, .1);
    }

    fileName = "../results/f22_hw02_problem2_";

    fileName = fileName + std::to_string(1);
    fileName.append(fileType);

    ls_phi.getPhi(p2_phi);

    newGrid.print_VTK_format(fileName);
    newGrid.print_VTK_format(vel_u, "vel_u", fileName);
    newGrid.print_VTK_format(vel_v, "vel_v", fileName);
    newGrid.print_VTK_format(p2_phi, "phi", fileName);


    // Problem 4 - Extra Credit


    return 0;
}
