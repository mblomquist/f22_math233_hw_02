#include <iostream>
#include "grid/Grid2d.h"
#include "tools/math_tools.h"
#include "levelset/LevelSet.h"
#include "math.h"

#define PI 3.1415926535897932384626

double cf_phi(double x, double y){
    return pow(pow((x - 0.25),2)+ pow(y,2), 0.5) - 0.2;
//    return pow(pow((x),2)+ pow(y,2), 0.5) - 0.5; // this is for problem 2
}

void cf_vel0(double x, double y, double & vel_u, double & vel_v){
    vel_u = -y;
    vel_v = x;
}

void measure_error(std::vector<double> phi, std::vector<double> phi_0){

    double l1, l2, linf;

    l1   = norm_l1_ls(phi, phi_0, 0.02);
    l2   = norm_l2_ls(phi, phi_0, 0.02);
    linf = norm_linf_ls(phi, phi_0, 0.02);

    std::cout << "L1 Norm: " << l1 << std::endl;
    std::cout << "L2 Norm: " << l2 << std::endl;
    std::cout << "Linf Norm: " << linf << std::endl;
}

void print_results_vtk(Grid2d newGrid, std::string fileName, double cfl, int itr,
                       std::vector<double> phi, std::vector<double> vel_u, std::vector<double> vel_v){
    // build fileName
    fileName = fileName + std::to_string(cfl);
    fileName.append("_");
    fileName = fileName + std::to_string(itr);
    fileName.append(".vtk");

    // print to file
    newGrid.print_VTK_format(fileName);
    newGrid.print_VTK_format(vel_u, "vel_u", fileName);
    newGrid.print_VTK_format(vel_v, "vel_v", fileName);
    newGrid.print_VTK_format(phi, "phi", fileName);

}

int main() {
    std::cout << "\n\nProblem 1!\n\n" << std::endl;

    // -------------------------------------- //
    // Problem 1 - Semi-Lagrangian Method for the advection equation
    // -------------------------------------- //
    double xmin = -1.0;
    double xmax = 1.0;
    double ymin = -1.0;
    double ymax = 1.0;

    int N = 64;
    int M = N;

    double ratio = 0.1;
    double t_final = 2.*PI;

    std::vector<double> vel_u, vel_v, phi, phi_0;
    vel_u.resize(N*M);
    vel_v.resize(N*M);
    phi.resize(N*M);
    phi_0.resize(N*M);

    Grid2d newGrid(N, M, xmin, xmax, ymin, ymax);

    double dx = newGrid.get_dx();
    double dt = dx/ratio;

    // initialize vel_u, vel_v, phi_0
    for (int i = 0; i < N; i++){
        for (int j = 0; j < M; j++){
            cf_vel0(newGrid.x_from_n(newGrid.n_from_ij(i,j)),
                    newGrid.y_from_n(newGrid.n_from_ij(i,j)),
                    vel_u[newGrid.n_from_ij(i,j)],
                    vel_v[newGrid.n_from_ij(i,j)]);

            phi_0[newGrid.n_from_ij(i,j)] = cf_phi(newGrid.x_from_n(newGrid.n_from_ij(i,j)),
                                                 newGrid.y_from_n(newGrid.n_from_ij(i,j)));
        }
    }

    print_results_vtk(newGrid, "../results/f22_hw02_problem1_", (1./ratio), 0,
                      phi_0, vel_u, vel_v);

    LevelSet ls_phi(newGrid, phi_0);

    // compute steps to t_final
    int steps = floor(t_final / dt);
    int iprint = 0;
    int mprint = 1; //floor(steps * ratio) + 1;

    for (int i = 1; i <= steps+1; i++){

        if (i*dt > t_final){
            double temp_t = (i-1)*dt;
            dt = t_final - (i-1)*dt;
            std::cout << "Percent Complete: " << (temp_t + dt) / t_final * 100. << " %," << std::endl;
        } else {
            std::cout << "Percent Complete: " << i*dt/t_final * 100. << " %," << std::endl;
        }

        ls_phi.advance_sl(vel_u, vel_v, dt);

        if (i % mprint == 0){
            iprint++;
            ls_phi.getPhi(phi);
            print_results_vtk(newGrid, "../results/f22_hw02_problem1_", (1./ratio), iprint,
                              phi, vel_u, vel_v);
        }
    }

    // update phi
    ls_phi.getPhi(phi);

    // measure the error
    std::cout << "\nProblem 1 - Measured Error for CFL " << (1./ratio) << std::endl;
    std::cout << "dt: " << dx/ratio << std::endl;
    std::cout << "dx: " << dx << std::endl;
    measure_error(phi, phi_0);

    // -------------------------------------- //
    // Problem 2 - Reinitialization Equation
    // -------------------------------------- //
    std::cout << "\n\nProblem 2!\n\n" << std::endl;

    std::vector<double> phi_p, phi_pi, phi_pf;
    phi_p.resize(N*M);
    phi_pi.resize(N*M);
    phi_pf.resize(N*M);
    phi_p = phi_0;

    print_results_vtk(newGrid, "../results/f22_hw02_problem2_", (1./ratio), 0,
                      phi_p, vel_u, vel_v);

    ls_phi.perturb_ls(phi_p, 1.e-1);
    ls_phi.setPhi(phi_p);
    phi_pi = phi_p;

    print_results_vtk(newGrid, "../results/f22_hw02_problem2_", (1./ratio), 1,
                      phi_p, vel_u, vel_v);

    int mitr = 100;

    for (int i = 2; i < mitr; i++){
        ls_phi.reinitialize(phi_p, phi_pi, phi_pf);
        print_results_vtk(newGrid, "../results/f22_hw02_problem2_", (1./ratio), i,
                          phi_pf, vel_u, vel_v);
        phi_pi = phi_pf;
    }



    // -------------------------------------- //
    // Problem 3 - Level Set Method
    // -------------------------------------- //

// Switch for Problem 3
#if 1
    std::cout << "\n\nProblem 3!\n\n" << std::endl;

    std::vector<double> phi2, phi_t;
    phi2.resize(N*M);
    phi_t.resize(N*M);

    ls_phi.setPhi(phi_0);
    ls_phi.getPhi(phi);

    print_results_vtk(newGrid, "../results/f22_hw02_problem3_", (1./ratio), 0,
                      phi, vel_u, vel_v);

    dt = dx/ratio;
    iprint = 0;

    for (int i = 1; i <= steps+1; i++){

        if (i*dt > t_final){
            double temp_t = (i-1)*dt;
            dt = t_final - (i-1)*dt;
            std::cout << "Percent Complete: " << (temp_t + dt) / t_final * 100. << " %," << std::endl;
        } else {
            std::cout << "Percent Complete: " << i*dt/t_final * 100. << " %," << std::endl;
        }

        ls_phi.advance_sl(vel_u, vel_v, dt);
        ls_phi.getPhi(phi);

        phi_t = phi;

        double err_p3 = 1.;
        int maxitrs = 10;
        int itrs = 0;

        while (err_p3 > 1.e-8 && itrs < maxitrs){
            ls_phi.reinitialize(phi_t, phi, phi2);
            err_p3 = norm_linf_ls(phi, phi2, 0.02);
            phi = phi2;
            itrs++;
        }

        ls_phi.setPhi(phi2);

        if (i % mprint == 0){
            iprint++;
            ls_phi.getPhi(phi2);
            print_results_vtk(newGrid, "../results/f22_hw02_problem3_", (1./ratio), iprint,
                              phi2, vel_u, vel_v);
        }
    }

    // measure the error
    std::cout << "\nProblem 3 - Measured Error for CFL " << (1./ratio) << std::endl;
    measure_error(phi, phi_0);
#endif

    // -------------------------------------- //
    // Problem 4 - Extra Credit
    // -------------------------------------- //
//    std::cout << "\n\nProblem 4!\n\n" << std::endl;


    return 0;
}
