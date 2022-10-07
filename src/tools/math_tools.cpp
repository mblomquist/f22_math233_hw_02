//
// Created by Matt Blomquist on 10/6/22.
//

#include "math_tools.h"
#include <cmath>

// suppose we have grid in 2D [xmin,xmax] x [ymin,ymax]
// find cell in which (x,y) belongs
// find weighted avg of values (?)

double bilinear_interpolation(Grid2d & grid,std::vector<double> & func,double x, double y){
    //  get xmin and get ymin
    // if (x,y) are outside domain, throw warning/error message
    if (x < grid.get_xmin() || x > grid.get_xmax() || y < grid.get_ymin() || y > grid.get_ymax() ) {
        throw std::invalid_argument("ERROR: (x,y) is outside of grid");
    }

    double phi;
    double dx = grid.get_dx();
    double dy = grid.get_dy();

    std::cout <<"dx: " << dx << " dy: " << dy << std::endl;

    // get i and j to find which cell C x belongs to
    //  i = floor((x - xmin)/dx)
    //  j = floor((y-ymin)/dx)
    int i = floor( (x - grid.get_xmin()) / dx);
    int j = floor( (y - grid.get_ymin()) / dy);

    std::cout << "i: " << i << " j: " << j << std::endl;

    double x_i = grid.get_xmin() + i * dx;
    double y_j = grid.get_ymin() + j * dy;
    double x_ip1 = x_i + dx;
    double y_jp1 = y_j + dy;

    // Use quadratic interpolation (formula in Lab 3) to get value at x
    // (i.e. think weighted avg)
    // (i, j), (i + 1, j), (i, j + 1), (i + 1, j + 1) are the corners of the cell C
    phi  = func[grid.n_from_ij(i    ,   j  )]  * ( x_ip1 - x   ) * ( y_jp1 - y   ) / (dx*dy) ;
    phi += func[grid.n_from_ij(i+1,   j  )]  * ( x     - x_i ) * ( y_jp1 - y   ) / (dx*dy) ;
    phi += func[grid.n_from_ij(i    , j+1)]  * ( x_ip1 - x   ) * ( y     - y_j ) / (dx*dy) ;
    phi += func[grid.n_from_ij(i+1, j+1)]  * ( x     - x_i ) * ( y     - y_j ) / (dx*dy) ;

    return phi;
}

double minmod(double x, double y){
    if ((x*y) >= 0.)
        return 0.;
    else {
        if (abs(x) < abs(y))
            return x;
        else
            return y;
    }
}

double quadratic_interpolation(Grid2d & grid,std::vector<double> & func,double x, double y){

    double phi = 0.;

    // check if point is outside the domain
    if (x < grid.get_xmin())
        x = grid.get_xmin();
    if (x > grid.get_xmax())
        x = grid.get_xmax();
    if (y < grid.get_ymin())
        y = grid.get_ymin();
    if (y > grid.get_ymax())
        y = grid.get_ymax();

    double dx = grid.get_dx();
    double dy = grid.get_dy();

    int i = floor( (x - grid.get_xmin()) / dx);
    int j = floor( (y - grid.get_ymin()) / dy);

    double x_i = grid.get_xmin() + i * dx;
    double y_j = grid.get_ymin() + j * dy;
    double x_ip1 = x_i + dx;
    double y_jp1 = y_j + dy;

    double phi_xx_l = (func[grid.n_from_ij(i-1,j)] - 2.*func[grid.n_from_ij(i,j)] + func[grid.n_from_ij(i+1,j)])/dx/dx;
    double phi_xx_r = (func[grid.n_from_ij(i,j)] - 2.*func[grid.n_from_ij(i+1,j)] + func[grid.n_from_ij(i+2,j)])/dx/dx;
    double phi_yy_t = (func[grid.n_from_ij(i,j-1)] - 2.*func[grid.n_from_ij(i,j)] + func[grid.n_from_ij(i,j-1)])/dy/dy;
    double phi_yy_b = (func[grid.n_from_ij(i,j)] - 2.*func[grid.n_from_ij(i,j+1)] + func[grid.n_from_ij(i,j+2)])/dy/dy;

    phi = func[grid.n_from_ij(i,j)] * (x_ip1-x)*(y_jp1-y)/dx/dy
        + func[grid.n_from_ij(i,j+1)] * (x_ip1-x)*(y-y_j)/dx/dy
        + func[grid.n_from_ij(i+1,j)] * (x-x_i)*(y_jp1-y)/dx/dy
        + func[grid.n_from_ij(i+1,j+1)] * (x-x_i)*(y-y_j)/dx/dy
        - (x-x_i)*(x_ip1-x)/2 * minmod(phi_xx_l, phi_xx_r)
        - (y-y_j)*(y_jp1-y)/2 * minmod(phi_yy_t, phi_yy_b);

    return phi;
}