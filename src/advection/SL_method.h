//
// Created by Matt Blomquist on 10/6/22.
//

#ifndef F22_MATH233_HW_02_SL_METHOD_H
#define F22_MATH233_HW_02_SL_METHOD_H


#include "../grid/Grid2d.h"
#include <vector>

// Semi-Langrangian Method
class SL_method {
private:
    Grid2d sl_grid;
    std::vector<double> sol;
    std::vector<double> vel_u;
    std::vector<double> vel_v;
    void find_trajectory(int n, double & x_d, double & y_d, double dt);

public:
    void set_grid(Grid2d & new_grid){sl_grid = new_grid;} // set grid
    std::vector<double> get_sol(){ return sol; }        // access solution
    void set_velocity(std::vector<double> & vel_u0, std::vector<double> & vel_v0);


};


#endif //F22_MATH233_HW_02_SL_METHOD_H
