//
// Created by Matt Blomquist on 10/6/22.
//

#ifndef F22_MATH233_HW_02_MATH_TOOLS_H
#define F22_MATH233_HW_02_MATH_TOOLS_H


#include "../grid/Grid2d.h"
#include <vector>


double bilinear_interpolation(Grid2d & grid,std::vector<double> & func,double x, double y);

double quadratic_interpolation(Grid2d & grid,std::vector<double> & func,double x, double y);

double minmod(double x, double y);

double norm_l1(std::vector<double> & x, std::vector<double> & y);
double norm_l2(std::vector<double> & x, std::vector<double> & y);
double norm_linf(std::vector<double> & x, std::vector<double> & y);

double norm_l1_ls(std::vector<double> & x, std::vector<double> & y, double eps);
double norm_l2_ls(std::vector<double> & x, std::vector<double> & y, double eps);
double norm_linf_ls(std::vector<double> & x, std::vector<double> & y, double eps);

double sgn(double x);

double max(double a, double b);
double min(double a, double b);

double sqroot(double a);

#endif //F22_MATH233_HW_02_MATH_TOOLS_H
