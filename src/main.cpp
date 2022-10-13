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

    int N = 10;
    int M = 10;

    Grid2d newGrid(N, M, xmin, xmax, ymin, ymax);

    newGrid.print_VTK_format("../results/f22_hw02_problem1.vtk");

    std::vector<double> vel_u, vel_v;
    vel_u.resize(N*M);
    vel_v.resize(N*M);

    for (int i = 0; i < N; i++){
        for (int j = 0; j < M; j++){
            cf_vel0(newGrid.x_from_n(newGrid.n_from_ij(i,j)),
                    newGrid.y_from_n(newGrid.n_from_ij(i,j)),
                    vel_u[newGrid.n_from_ij(i,j)],
                    vel_v[newGrid.n_from_ij(i,j)]);
        }
    }

    newGrid.print_VTK_format(vel_u, "vel_u", "../results/f22_hw02_problem1.vtk");
    newGrid.print_VTK_format(vel_v, "vel_v", "../results/f22_hw02_problem1.vtk");

    // Problem 2 - Reinitialization Equation


    // Problem 3 - Level Set Method


    // Problem 4 - Extra Credit


    return 0;
}
