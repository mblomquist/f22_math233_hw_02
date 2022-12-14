//
// Created by Matt Blomquist on 10/6/22.
//

#ifndef F22_MATH233_HW_02_GRID2D_H
#define F22_MATH233_HW_02_GRID2D_H

#include <iostream>
#include <fstream>
#include <vector>

class Grid2d {
    // Objective: creates a 2D grid
private:
    double dx;  // spacing in x
    double dy;  // spacing in y
    long N;     // number of nodes in x
    long M;     // number of nodes in y
    double xmin;
    double xmax;
    double ymin;
    double ymax;
public:
    double get_dx();    // return dx to user
    double get_dy();    // return dy to user

    double get_M(){return M;}
    double get_N(){return N;}

    double get_xmin();
    double get_ymin();
    double get_xmax();
    double get_ymax();

    // Constructor to create object
    Grid2d();
    // Constructor to initialize values
    Grid2d(long NN, long MM, double xlo, double xhi, double ylo, double yhi);
    // Return coord (i,j) from grid index n
    long i_from_n(long n);
    long j_from_n(long n);
    // Return grid index n from coord (i, j)
    long n_from_ij(long i, long j);
    // Return position of grid index n
    double x_from_n(long n);
    double y_from_n(long n);

    // output file in VTK format
    void print_VTK_format(std::string output_file);
    void print_VTK_format(std::vector<double> &F, std::string data_name,
                          std::string file_name);
};

#endif //F22_MATH233_HW_02_GRID2D_H
