#ifndef _FITTING_ROOT_H_
#define _FITTING_ROOT_H_

// This file includes methods taken from the following files (which are no more updated):
// linear_fit.h
// plane_fit.h
// circular_fit.h
// helix_fit.h

#include <vector>

// #include <Eigen/Dense>

#include "point.h"

#include "TMatrixD.h"

#include "Rtypes.h"

#define def_lft SVD

// Linear fitting ////////////////////////////////////////////////////////////////////////////////

enum fit_type { // Here are different types of fitting
    SVD, // The default type
    X, // Linear regression for y = y(x), z = z(x); BE CAREFUL WITH THESE 3 METHODS:
    Y, // -||-                                        THEY DON'T ALWAYS WORK!
    Z, // -||-                                        BAD CASE - LINE ORTHOGONAL TO PARAMETER AXIS
    AVG, // Three previous averaged with weights (see cpp-file); normally works ~ as SVD
    OPT, // Chooses the one with smallest chi^2 from 5 previous after adding each point; works slow
    CIRC, // Circle fitting 
    HEL // Helix fitting
};

const char * fit_type_to_string(const fit_type t);
const bool is_linear(const fit_type t);

void least_squares(const std::vector <double> x, const std::vector <double> y,
                   double& k, double& b, double& eta, bool& use);
void svd_fit_3(const std::vector <point> points, point& dir); // dir points from begin to end
void reg_fit(const std::vector <point> points, const unsigned type, point& dir); // type: 0 for AVG, 1 for X, 2 for Y, 3 for Z, else AVG
void linear_fit(const std::vector <point> points, point& dir, const fit_type lft);



// Plane fitting /////////////////////////////////////////////////////////////////////////////////

// A * x + B * y + C * z + D = 0
// A ^ 2 + B ^ 2 + C ^ 2 = 1

void plane_fit(std::vector <point> points, double& A, double& B, double& C, double& D);



// Circle fitting ////////////////////////////////////////////////////////////////////////////////

// First three methods minimize the error of F(x) = 0
// Second one is the most stable!! Two others can even give the negative radius!!

void circle_fit_2d(const std::vector <double> x, const std::vector <double> y, double& X, double& Y, double& R);
void circle_fit_2d_(const std::vector <double> x, const std::vector <double> y, double& X, double& Y, double& R);
void circle_fit_2d__(const std::vector <double> x, const std::vector <double> y, double& X, double& Y, double& R);

void circle_fit_2d_root(const std::vector <double> x, const std::vector <double> y, double& X, double& Y, double& R);
void myfcn_root(Int_t &, Double_t *, Double_t &f, Double_t *par, Int_t);

// The last method solves nonlinear least squares problem using Gauss-Newton algorithm (10 steps by default)
// The solution of previous problem is used as a starting vector
void circle_fit_2d___(const std::vector <double> x, const std::vector <double> y, double& X, double& Y, double& R);

// Here we first fit points by plane using least squares (see plane_fit.h) and then use circle_fit_2d___
void circle_fit_3d(const std::vector <point>& points, point& mid, point& normal, double& r);



// Helix fitting /////////////////////////////////////////////////////////////////////////////////

void helix_fit(std::vector <point> points, point& mid, const point normal, double& r, double& step, bool& right, point& on);

void stdout_hide(int &fd, fpos_t &pos);
void stdout_unhide(const int &fd, const fpos_t &pos);

#endif

