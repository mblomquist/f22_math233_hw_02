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

#endif //F22_MATH233_HW_02_MATH_TOOLS_H
