//
// Created by mblomquist on 10/13/22.
//

#ifndef F22_MATH233_HW_02_LEVELSET_H
#define F22_MATH233_HW_02_LEVELSET_H

#include <vector>
#include "../grid/Grid2d.h"
#include "../advection/SL_method.h"
#include "../tools/math_tools.h"

class LevelSet {

private:
    Grid2d grid;
    std::vector<double> phi;

public:
    LevelSet();
    LevelSet(Grid2d & newGrid, std::vector<double> & phi_0);

    void setGrid(Grid2d & newGrid);
    void setPhi(std::vector<double> & phi_0);

    void getGrid(Grid2d & grid_o);
    void getPhi(std::vector<double> & phi_o);

    void advance_sl(std::vector<double> & vel_u, std::vector<double> & vel_v, double dt);

    void reinitialize(std::vector<double> & phi_0, std::vector<double> phi_n, std::vector<double> & phi_np1);

};


#endif //F22_MATH233_HW_02_LEVELSET_H
