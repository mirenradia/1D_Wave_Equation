#ifndef CONSTANTS_H
#define CONSTANTS_H

const double c {1.0}; //c is the wavespeed
const double dx {0.01}; //spatial step size
const double dt {0.01}; //temporal step size
const double alpha {c*dt/dx}; //Courant factor
const double xMax {5.0}; //max value of x
const double xMin {-5.0}; //min value of x
const int nxSteps {static_cast<int>((xMax-xMin)/dx)}; //width of spatial grid
const double tMax {12.0}; //max value of t
const int ntSteps {static_cast<int>(tMax/dt)}; //number of time steps
const double a {1.0}; //initial Gaussian width in A*q-xexp(-(x-x_0)^2/(2*a^2))
const double A {1.0}; //initial Gaussian amplitude
const double x_0 {0}; //initial Gaussian centre

#endif
