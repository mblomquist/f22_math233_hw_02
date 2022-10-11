//
// Created by Matt Blomquist on 10/6/22.
//

#include "SL_method.h"

void SL_method::find_trajectory(int n, double &x_d, double &y_d, double dt) {
    // RK1 Euler Method
    double x_0 = sl_grid.x_from_n(n);
    double y_0 = sl_grid.y_from_n(n);
    x_d = x_0 - dt * vel_u[n];
    y_d = y_0 - dt * vel_v[n];
}

void SL_method::find_departure_2nd(int n, double & x_d, double & y_d, double dt) {
    // RK2 Mid-Point Rule
    double x_0 = sl_grid.x_from_n(n);
    double y_0 = sl_grid.y_from_n(n);

    double x_star = x_0 - 0.5 * dt * vel_u[n];
    double y_star = y_0 - 0.5 * dt * vel_v[n];

    // interpolate vel_u, vel_v at (x_star, y_star)
    double vel_u_star = quadratic_interpolation(sl_grid, vel_u, x_star, y_star);
    double vel_v_star = quadratic_interpolation(sl_grid, vel_v, x_star, y_star);

    x_d = x_0 - dt * vel_u_star;
    y_d = y_0 - dt * vel_v_star;
}

void SL_method::set_velocity(std::vector<double> &vel_u0, std::vector<double> &vel_v0) {
    vel_u.resize(vel_v0.size());
    vel_v.resize(vel_v0.size());

    vel_u = vel_u0;
    vel_v = vel_v0;
}