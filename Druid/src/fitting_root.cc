// This file was based on fitting.cpp
// Here all the linear algebra is replaced by ROOT functions (instead of Eigen package)

#include <unistd.h>
#include <vector>
#include <sstream>
#include <iostream>
#include <fstream>
#include <cmath>

#include "TDecompSVD.h"
#include "TMatrixDEigen.h"
#include "TVirtualFitter.h"

#include "point.h"
#include "fitting_root.h"

// #include <Eigen/Dense>

// using namespace Eigen;
using namespace std;

#define def_lft SVD // Default linear fitting method

vector <double> x_root, y_root;

// Linear fitting ////////////////////////////////////////////////////////////////////////////////

const char * fit_type_to_string(const fit_type t) { // Use this function for printing the fit_type of the segment
    switch(t){
        case SVD:
            return "SVD";
        case X:
            return "X";
        case Y:
            return "Y";
        case Z:
            return "Z";
        case AVG:
            return "AVG";
        case OPT:
            return "OPT";
        case CIRC:
            return "CIRC";
        case HEL:
            return "HEL";
        default:
            return "Error!";
         
      }
 }  

const bool is_linear(const fit_type ft) { // Returns 'true' if the fitting method is of linear type
    if ( ( ft == SVD ) || // If one chose the linear fitting
             ( ft == AVG ) ||
             ( ft == OPT ) ||
             ( ft == X ) ||
             ( ft == Y ) ||
             ( ft == Z ) ) {
        return true;
    }
    else {
        return false;
    }
}

void least_squares(const std::vector <double> x, const std::vector <double> y, // Performs usual least squares procedure for y = k * x + b line
                   double& k, double& b, double& eta, bool& use){ //eta == 0(bad)..1(good)
    double SUMx = 0, SUMy = 0, SUMxy = 0, SUMxx = 0, AVGx = 0, AVGy = 0, SUMerr_num = 0, SUMerr_den = 0;
    int n = x.size();
    double nd = (double) n;
    for (int i = 0; i < n; i++){
        SUMx += x[i];
        SUMy += y[i];
        SUMxy += x[i] * y[i];
        SUMxx += x[i] * x[i];
        }
    AVGx = SUMx / nd;
    AVGy = SUMy    / nd;
    k = (nd * SUMxy - SUMx * SUMy) / (nd * SUMxx - SUMx * SUMx);

    b = AVGy - k * AVGx;

    for (int i=0; i<n; i++){
        SUMerr_num += (y[i] - (k * x[i] + b)) * (y[i] - (k * x[i] + b));
        SUMerr_den += (y[i] - AVGy) * (y[i] - AVGy);
    }
    
    //eta = SUMerr_num  / (nd - 2.);
    eta = SUMerr_num / SUMerr_den;

    const double delta = 1e-10;

    if (fabs(SUMerr_den) < delta) {
        if (fabs(SUMerr_num) < delta) {
            eta = 1; // Most likely, points already were on the line
        }
        else {
            eta = 2345;
            use = 0;
        }
    }
    
    if ( (n * SUMxx - SUMx * SUMx) == 0 ){ // Here we check if the denominator of k is fishily small
        eta = 2345;
        use = 0;
    }

    //if (eta > 1 && use == 1) {
    //    cout << "HOHOHO!!!" << endl; // Hope, this will never appear on the screen :)
    //}
}

void svd_fit_3(const vector <point> points, point& dir) { // Fits the std::vector of points in 3D by a straight line with help of SVD-decomposition. The orthogonal squared distance (ODR) is minimized

    // Hiding ROOT's output:
        //freopen("/dev/null", "a", stdout);
    int fd;
    fpos_t pos;
    stdout_hide(fd, pos);
    
    point centroid = average(points);
    unsigned n = points.size(); double nd = (double) n;
    TMatrixD M_(n,3);
    // MatrixXd M(n, 3);
    // stringstream str;
    double x, y, z; 
    
    for (unsigned i = 0; i < n; i++) {
        M_(i, 0) = points[i].getX() - centroid.getX();
        M_(i, 1) = points[i].getY() - centroid.getY();
        M_(i, 2) = points[i].getZ() - centroid.getZ();
        // M(i, 0) = points[i].getX() - centroid.getX();
        // M(i, 1) = points[i].getY() - centroid.getY();
        // M(i, 2) = points[i].getZ() - centroid.getZ();
    }
    
    // JacobiSVD<MatrixXd> svd(M, ComputeThinV);
    TMatrixD v_ = TDecompSVD(M_).GetV();
    

    // str.precision(17);
    // str << svd.matrixV().col(0).row(0) << " " << svd.matrixV().col(0).row(1) << " " << svd.matrixV().col(0).row(2);
    // str >> x >> y >> z;

    // dir = point(x, y, z);
    dir = point(v_(0,0), v_(1,0), v_(2,0));    

    if (angle(dir, (points.back() - points.front())) > m_pi / 2) {
        dir = - dir;
    }

    dir = dir.normalized();
    
    // Fixing the output back
        // freopen("/dev/stdout", "a", stdout);
    stdout_unhide(fd, pos);
}

void reg_fit(const std::vector <point> points, const unsigned type, point& dir) {
    dir = point(0,0,0);
    unsigned n = points.size();
    point dir1, dir2, dir3;
    double k_x_y, k_x_z, b_x_y, b_x_z, eta_x_y, eta_x_z; // ?_<parameter>_?
    double k_y_z, k_y_x, b_y_z, b_y_x, eta_y_z, eta_y_x;
    double k_z_x, k_z_y, b_z_x, b_z_y, eta_z_x, eta_z_y;
    double eta1, eta2, eta3;
    double w1 = 0, w2 = 0, w3 = 0, w_den = 0;
    bool use1 = 1, use2 = 1, use3 = 1;
    vector <double> X, Y, Z;
    for (unsigned i = 0; i < n; i++) {
        X.push_back(points[i].getX());
        Y.push_back(points[i].getY());
        Z.push_back(points[i].getZ());
    }
    point fl = points.back() - points.front();
    
    least_squares(X, Y, k_x_y, b_x_y, eta_x_y, use1); // x is a parameter
    least_squares(X, Z, k_x_z, b_x_z, eta_x_z, use1);
    dir1 = point(1, k_x_y, k_x_z).normalized();
    eta1 = (eta_x_y + eta_x_z) / 2.;
    if (angle(fl, dir1) > m_pi/2) {dir1 = - dir1;}
    if (type == 1) {
        if (use1 != 0) {
            dir = dir1;
            //err = eta1;
        }
        else {
            dir = point(0,0,0);
            cout << "Error in simple linear regression.\n";
            //err = 1234;
        }
        return;
    }
    
    least_squares(Y, Z, k_y_z, b_y_z, eta_y_z, use2); // y is a parameter
    least_squares(Y, X, k_y_x, b_y_x, eta_y_x, use2);
    dir2 = point(k_y_x, 1, k_y_z).normalized();
    eta2 = (eta_y_z + eta_y_x) / 2.;
    if (angle(fl, dir2) > m_pi/2) {dir2 = - dir2;}
    if (type == 2) {
        if (use2 != 0) {
            dir = dir2;
            //err = eta2;
        }
        else {
            dir = point(0,0,0);
            cout << "Error in simple linear regression.\n";
            //err = 1234;
        }
        return;
    }

    least_squares(Z, X, k_z_x, b_z_x, eta_z_x, use3); // z is a parameter
    least_squares(Z, Y, k_z_y, b_z_y, eta_z_y, use3);
    dir3 = point(k_z_x, k_z_y, 1).normalized();
    eta3 = (eta_z_x + eta_z_y) / 2.;
    if (angle(fl, dir3) > m_pi/2) {dir3 = - dir3;}
    if (type == 3) {
        if (use3 != 0) {
            dir = dir3;
            //err = eta3;
        }
        else {
            dir = point(0,0,0);
            cout << "Error in simple linear regression.\n";
            //err = 1234;
        }
        return;
    }
    
    // Some special cases:
    if ((use1 == 1) && (use2 == 1) && (use3 == 1)){
        if ((eta1 == 0) && (eta2 != 0) && (eta3 != 0)){
            dir = dir1;
            //err = 0;
        }
        else if ((eta1 != 0) && (eta2 == 0) && (eta3 != 0)){
            dir = dir2;
            //err = 0;
        }
        else if ((eta1 != 0) && (eta2 != 0) && (eta3 == 0)){
            dir = dir3;
            //err = 0;
        }
        else if ((eta1 == 0) && (eta2 == 0) && (eta3 != 0)){
            dir = 0.5*dir1 + 0.5*dir2;
            //err = 0;
        }
        else if ((eta1 != 0) && (eta2 == 0) && (eta3 == 0)){
            dir = 0.5*dir2 + 0.5*dir3;
            //err = 0;
        }
        else if ((eta1 == 0) && (eta2 != 0) && (eta3 == 0)){
            dir = 0.5*dir3 + 0.5*dir1;
            //err = 0;
        }
        else if ((eta1 == 0) && (eta2 == 0) && (eta3 == 0)){
            dir = (dir1 + dir2 + dir3) * (1./3.);
            //err = 0;
        }
        else {
            w_den = (1./eta1) + (1./eta2) + (1./eta3);
            w1 = (1./eta1) / w_den;
            w2 = (1./eta2) / w_den;
            w3 = (1./eta3) / w_den;
            dir = dir1 * w1 + dir2 * w2 + dir3 * w3;
            //err = ((eta1 + eta2 + eta3) / 3.);
        }
    }
    else if ((use1 == 1) && (use2 == 1) && (use3 == 0)){
        if ((eta1 == 0) && (eta2 != 0)){
            dir = dir1;
            //err = 0;
        }
        if ((eta1 != 0) && (eta2 == 0)){
            dir = dir2;
            //err = 0;
        }
        if ((eta1 == 0) && (eta2 == 0)){
            dir = 0.5 * dir1 + 0.5 * dir2;
            //err = 0;
        }
        else {
            w_den = (1./eta1) + (1./eta2);
            w1 = (1./eta1) / w_den;
            w2 = (1./eta2) / w_den;
            dir = dir1 * w1 + dir2 * w2;
            //err = ((eta1 + eta2) / 2.);
        }
    }
    else if ((use1 == 0) && (use2 == 1) && (use3 == 1)){
        if ((eta2 == 0) && (eta3 != 0)){
            dir = dir2;
            //err = 0;
        }
        if ((eta2 != 0) && (eta3 == 0)){
            dir = dir3;
            //err = 0;
        }
        if ((eta2 == 0) && (eta3 == 0)){
            dir = 0.5 * dir2 + 0.5 * dir3;
            //err = 0;
        }
        else {
            w_den = (1./eta2)+ (1./eta3);
            w2 = (1./eta2) / w_den;
            w3 = (1./eta3) / w_den;
            dir = dir2 * w2 + dir3 * w3;
            //err = ((eta2 + eta3) / 2.);
        }
    }
    else if ((use1 == 1) && (use2 == 0) && (use3 == 1)){
        if ((eta3 == 0) && (eta1 != 0)){
            dir = dir3;
            //err = 0;
        }
        if ((eta3 != 0) && (eta1 == 0)){
            dir = dir1;
            //err = 0;
        }
        if ((eta3 == 0) && (eta1 == 0)){
            dir = 0.5 * dir3 + 0.5 * dir1;
            //err = 0;
        }
        else {
            w_den = (1./eta3) + (1./eta1);
            w3 = (1./eta3) / w_den;
            w1 = (1./eta1) / w_den;
            dir = dir3 * w3 + dir1 * w1;
            //err = ((eta3 + eta1) / 2.);
        }
    }
    else if ((use1 == 1) && (use2 == 0) && (use3 == 0)){
        dir = dir1;
        //err = (eta1);
    }
    else if ((use1 == 0) && (use2 == 1) && (use3 == 0)){
        dir = dir2;
        //err = (eta2);
    }
    else if ((use1 == 0) && (use2 == 0) && (use3 == 1)){
        dir = dir3;
        //err = (eta3);
    }
    dir = dir.normalized();
}

void linear_fit(const std::vector <point> points, point& dir, const fit_type lft) {
    if ( lft == SVD ) {
        svd_fit_3(points, dir);
    }
    else if ( lft == X ) {
        reg_fit(points, 1, dir);
    }
    else if ( lft == Y ) {
        reg_fit(points, 2, dir);
    }
    else if ( lft == Z ) {
        reg_fit(points, 3, dir);
    }
    else if ( lft == AVG ) {
        reg_fit(points, 0, dir);
    }
    else if ( lft == OPT ) {
        point dir_n;
        double err_n, err_c;
        point centroid = average(points);
        
        svd_fit_3(points, dir);
        err_c = chi_2_lin(points, centroid, dir);

        for (unsigned i = 0; i <= 3; i++) {
            reg_fit(points, i, dir_n);
            err_n = chi_2_lin(points, centroid, dir_n);
            if (err_n < err_c) {
                dir = dir_n;
            }
        }
    }

    if ( angle(dir, points[points.size() - 1] - points[0]) > m_pi / 2. ) {
        dir = -dir;
    }
    
    dir = dir.normalized();
}



// Circle fitting ////////////////////////////////////////////////////////////////////////////////

void circle_fit_3d(const std::vector <point>& points, point& mid, point& n3, double& R) {
    vector <point> ppoints; // Points projected on the plane..
    point p;
    double A, B, C, D;
    plane_fit(points, A, B, C, D); //.. on this plane :)
    point in = project_onto_plane(point(0,0,0), A, B, C, D);
    //point in(-D*A, -D*B, -D*C);
    for (int i = 0; i < (int) points.size(); i++) {
        p = project_onto_plane(points[i], A, B, C, D);
        //p.p_print();
        ppoints.push_back(p);
    }

    // Now we want to form a orthonormal basis with n3(A,B,C) and n1&n2 in the plane
    // To do that we need ONE vector which is not equal to n3
    n3 = point(A, B, C).normalized();
    point n1, n2, not_n3;
    if ( (fabs(A-B) < 1e-10)  && (fabs(A-C) < 1e-10) ) { // This can happen only if A == B == C == (+/-) 1/sqrt(3)
        not_n3 = point(-A, B, C); // In both cases it is not n3
    }
    else { // Otherwise point(B,C,A) != point(A,B,C)
        not_n3 = point(B, C, A); // If some of them are different, it is not n3
    }
    n1 = vector_p(not_n3, n3);
    n2 = vector_p(n3, n1);
    n1 = n1.normalized();
    n2 = n2.normalized();
    /*n1.p_print();
    n2.p_print();
    vector_p(n1,n2).p_print();
    n3.p_print();*/
    
    vector <double> x, y; // Now we fill find the coordinates of points in (n1,n2) basis

    for (int i = 0; i < (int) ppoints.size(); i++) {
        x.push_back(scalar_p(ppoints[i],n1));
        y.push_back(scalar_p(ppoints[i],n2));
        //cout << "x.push_back(" << x[i] << ");\ny.push_back(" << y[i] << ");\n";
        //cout << x[i] << " " << y[i] << " " << endl;
    }

    double X, Y; // Coordinates of the center of the circle in the (n1,n2) basis
    double xm, ym, zm;
    
    circle_fit_2d_root(x, y, X, Y, R);
    double Z = scalar_p(n3, in);

    xm = X * n1.getX() + Y * n2.getX() + Z * n3.getX();
    ym = X * n1.getY() + Y * n2.getY() + Z * n3.getY();
    zm = X * n1.getZ() + Y * n2.getZ() + Z * n3.getZ();
    mid = point(xm, ym, zm);
    //cout << zm << endl;
}

void circle_fit_2d(const std::vector <double> x, const std::vector <double> y, double& X, double& Y, double& R) {
    // Hiding ROOT's output:
        //freopen("/dev/null", "a", stdout);
    int fd;
    fpos_t pos;
    stdout_hide(fd, pos);

    // MatrixXd Bt(3, x.size());
    // MatrixXd d(x.size(), 1);
    // MatrixXd B(x.size(), 3);
    
    TMatrixD Bt_(3, x.size());
    TMatrixD d_(x.size(), 1);
    TMatrixD B_(x.size(), 3);
    
    

    for (int i = 0; i < (int) x.size(); i++) {
        // Bt(0,i) = x[i];
        // Bt(1,i) = y[i];
        // Bt(2,i) = 1;
        // d(i, 0) = x[i] * x[i] + y[i] * y[i];
        
        Bt_(0,i) = x[i];
        Bt_(1,i) = y[i];
        Bt_(2,i) = 1;
        d_(i, 0) = x[i] * x[i] + y[i] * y[i];
    }
    
    // B = Bt.transpose();
    B_.Transpose(Bt_);
    
    // MatrixXd z = (Bt * B).inverse() * Bt * d;
    TMatrixD z_ = (Bt_ * B_).Invert() * Bt_ * d_; //Mult( Mult( TMatrixD( Mult( Bt_, B_ ) ), Bt_ ), d_ );
    // X = (double) z(0,0) / 2.;
    // Y = (double) z(1,0) / 2.;
    // R = (double) sqrt(z(2,0) + z(0,0) * z(0,0) / 4. + z(1,0) * z (1,0) / 4.);
    
    X = (double) z_(0,0) / 2.;
    Y = (double) z_(1,0) / 2.;
    R = (double) sqrt(z_(2,0) + z_(0,0) * z_(0,0) / 4. + z_(1,0) * z_(1,0) / 4.);
    
    // Fixing the output back
    // freopen("/dev/stdout", "a", stdout);
    stdout_unhide(fd, pos);
    
    return;
}

void circle_fit_2d_(const std::vector <double> x, const std::vector <double> y, double& X, double& Y, double& R) {
    // TODO: check that points are not on one line (that includes >2 points check)
    double cx = 0, cy = 0, Suuu = 0, Svvv = 0, Suvv = 0, Svuu = 0, Suu = 0, Svv = 0, Suv = 0, r;
    double uc, vc, xc, yc;
    unsigned n = x.size(); double nd = (double) n;
    vector <double> u, v;
    for (unsigned i = 0; i < n; i++) {
        cx += x[i];
        cy += y[i];
    }
    cx /= nd;
    cy /= nd;
    for (unsigned i = 0; i < n; i++) {
        u.push_back( (x[i] - cx) );
        v.push_back( (y[i] - cy) );
    }
    for (unsigned i = 0; i < n; i++) {
        Suuu += ( u[i] * u[i] * u[i] );
        Svvv += ( v[i] * v[i] * v[i] );
        Suvv += ( u[i] * v[i] * v[i] );
        Svuu += ( v[i] * u[i] * u[i] );
        Suu += ( u[i] * u[i] );
        Suv += ( u[i] * v[i] );
        Svv += ( v[i] * v[i] );
    }
    uc = (Suuu * Svv + Suvv * Svv - Svvv * Suv - Svuu * Suv) / (Suu * Svv - Suv * Suv) /2.;
    vc = (Svvv * Suu + Svuu * Suu - Suuu * Suv - Suvv * Suv) / (Suu * Svv - Suv * Suv) /2.;
    //cout << (Suu * Svv - Suv * Suv) << endl;
    xc = uc + cx;
    yc = vc + cy;
    r = sqrt( uc * uc + vc * vc + ( Suu + Svv ) / nd );

    X = xc;
    Y = yc;
    R = r;
}

void circle_fit_2d__(const std::vector <double> x, const std::vector <double> y, double& X, double& Y, double& R) {
    double a, b, r = 0;
    double A, B, C, D, E;
    double Sx = 0, Sxx = 0, Sy = 0, Syy = 0, Sxy = 0, Sxxx = 0, Syyy = 0, Syxx = 0, Sxyy = 0;
    unsigned n = x.size(); double nd = (double) n;
    for (unsigned i = 0; i < n; i++) {
        Sx += x[i];
        Sxx += ( x[i] * x[i] );
        Sxxx += ( x[i] * x[i] * x[i] );
        Sy += y[i];
        Syy += ( y[i] * y[i] );
        Syyy += ( y[i] * y[i] * y[i] );
        Sxy += ( x[i] * y[i] );
        Sxyy += ( x[i] * y[i] * y[i] );
        Syxx += ( y[i] * x[i] * x[i] );
    }

    A = nd * Sxx - Sx * Sx;
    B = nd * Sxy - Sx * Sy;
    C = nd * Syy - Sy * Sy;
    D = 0.5 * ( nd * Sxyy - Sx * Syy + nd * Sxxx - Sx * Sxx );
    E = 0.5 * ( nd * Syxx - Sy * Sxx + nd * Syyy - Sy * Syy );

    a = (D * C - B * E) / (A * C - B * B);
    b = (A * E - B * D) / (A * C - B * B);

    for (unsigned i = 0; i < n; i++) {
        r += sqrt ( ( x[i] - a ) * ( x[i] - a ) + ( y[i] - b ) * ( y[i] - b ));
    }
    r /= nd;

    X = a;
    Y = b;
    R = r;
}

void circle_fit_2d___(const std::vector <double> x, const std::vector <double> y, double& X, double& Y, double& R) {
    // Hiding ROOT's output:
        //freopen("/dev/null", "a", stdout);
    int fd;
    fpos_t pos;
    stdout_hide(fd, pos);

     if ( x.size() < 3 ) {return;}
    // Here we will use the Gauss-Newton method for solving Nonlinear Least Squares problem
    const unsigned STEPS = 10; // Number of iterations
    // First we calculate a starting vector using one of previous algorithms
    double a, b, r;
    unsigned n = x.size();
    circle_fit_2d(x, y, a, b, r);
    if ( x.size() == 3 ) {
        X = a;
        Y = b;
        R = r;
        stdout_unhide(fd, pos);
        return;
    }
    //r =  (point(x[0], y[0], 0) - point(x[x.size()-1], y[y.size()-1], 0)).norm();
    //a = (x[0] + x[x.size() - 1]) / 2.;
    //b = (y[0] + y[y.size() - 1]) / 2.;

     // MatrixXd u(3, 1);
     // MatrixXd tmp(3, 1);
     // u << a, b, r;
     // MatrixXd J(n, 3), Jt(3, n), JtJm1(3, 3);
     // MatrixXd f(n, 1);
    
    TMatrixD u_(3, 1);
    TMatrixD tmp_(3, 1);
    u_(0,0) = a; u_(1,0) = b; u_(2,0) = r;
    TMatrixD J_(n, 3), Jt_(3, n), JtJm1_(3, 3);
    TMatrixD f_(n, 1);


    for (unsigned i = 0; i < STEPS; i++) { // These are the iterations of the method
        // On each step we construct the Jacobian
         // for (unsigned j = 0; j < n; j++) {
             // J(j, 0) = ( u(0,0) - x[j] ) / sqrt( (u(0,0) - x[j]) * (u(0,0) - x[j]) + (u(1,0) - y[j]) * (u(1,0) - y[j]) );
             // J(j, 1) = ( u(1,0) - y[j] ) / sqrt( (u(0,0) - x[j]) * (u(0,0) - x[j]) + (u(1,0) - y[j]) * (u(1,0) - y[j]) );
             // J(j, 2) = - 1;
             // f(j, 0) = sqrt( (u(0,0) - x[j]) * (u(0,0) - x[j]) + (u(1,0) - y[j]) * (u(1,0) - y[j]) ) - u(2,0);
         // }
        
         // Jt = J.transpose();
         // JtJm1 = (Jt * J).inverse();

         // tmp = u - JtJm1 * Jt * f; // And perform the step!
         // u = tmp;
        
        
        for (unsigned j = 0; j < n; j++) {
            J_(j, 0) = ( u_(0,0) - x[j] ) / sqrt( (u_(0,0) - x[j]) * (u_(0,0) - x[j]) + (u_(1,0) - y[j]) * (u_(1,0) - y[j]) );
            J_(j, 1) = ( u_(1,0) - y[j] ) / sqrt( (u_(0,0) - x[j]) * (u_(0,0) - x[j]) + (u_(1,0) - y[j]) * (u_(1,0) - y[j]) );
            J_(j, 2) = - 1;
            f_(j, 0) = sqrt( (u_(0,0) - x[j]) * (u_(0,0) - x[j]) + (u_(1,0) - y[j]) * (u_(1,0) - y[j]) ) - u_(2,0);
        }
        
        Jt_.Transpose(J_);
        JtJm1_ = (Jt_ * J_).Invert();

        tmp_ = u_ - JtJm1_ * Jt_ * f_; // And perform the step!
        u_ = tmp_;
    }
    
     // X = u(0,0);
     // Y = u(1,0);
     // R = u(2,0);
    
    X = u_(0,0);
    Y = u_(1,0);
    R = u_(2,0);
        
   
    if ( R < 0 ) {
        cout << "AAAAAAAAAAAAAAAAAAAAAAAAAAAAA!!!! " << R << endl;
        for (unsigned i= 0; i < x.size(); i++) {
            //cout << "x.push_back(" << x[i] << ");\ny.push_back(" << y[i] << ");\n";
            cout << x[i] << " " << y[i] << " " << endl;
        }
        cout << endl;
    }
    
    // Fixing the output back
        // freopen("/dev/stdout", "a", stdout);
    stdout_unhide(fd, pos);
 }
 
 void circle_fit_2d_root(const std::vector <double> x, const std::vector <double> y, double& X, double& Y, double& R) {
    if ( x.size() < 3 ) {return;}

    circle_fit_2d_(x, y, X, Y, R);
      
    if ( x.size() == 3 ) {return;}
    
    x_root = x;
    y_root = y;
    
    // Hiding ROOT's output:
    int fd;
    fpos_t pos;
    stdout_hide(fd, pos); 

    
    //TVirtualFitter::SetDefaultFitter("Minuit");  //default is Minuit
 
    TVirtualFitter *fitter = TVirtualFitter::Fitter(0, 3);
    

    fitter->SetFCN(myfcn_root);
    

    
    const double err = 1e-5;

    fitter->SetParameter(0, "x0",   X, err, 0, 0);
    fitter->SetParameter(1, "y0",   Y, err, 0, 0);
    fitter->SetParameter(2, "R",    R, err, 0, 0);

    
    Double_t arglist[1] = {0};
    //gMinuit->SetPrintLevel(-1);
    
    
    fitter -> ExecuteCommand("MIGRAD", arglist, 0);
    
    X = fitter->GetParameter(0);
    Y = fitter->GetParameter(1);
    R = fitter->GetParameter(2);
    
    // Fixing the output back
    stdout_unhide(fd, pos);
    //cout << "Smilee\n";
 }
 
 void myfcn_root(Int_t &, Double_t *, Double_t &f, Double_t *par, Int_t) {
    //minimisation function computing the sum of squares of residuals
    unsigned n = x_root.size();
    f = 0;
    double u, v, dr;
    for (unsigned i = 0; i < n; i++) {
        u = x_root[i] - par[0];
        v = y_root[i] - par[1];
        dr = par[2] - sqrt(u * u + v * v);
        f += (dr * dr);
    }
}



// Plane fitting /////////////////////////////////////////////////////////////////////////////////

void plane_fit(vector <point> points, double& A, double& B, double& C, double& D) {
    // Hiding ROOT's output:
        //freopen("/dev/null", "a", stdout);
    int fd;
    fpos_t pos;
    stdout_hide(fd, pos);
    
    int n = (int) points.size(); double nd = (double) points.size(); double tmp;
    point centroid(0,0,0);
    // stringstream str;
    // Matrix3d M = Matrix3d::Zero(3,3);
    TMatrixD M_(3,3);
    M_.Zero();

    for (int i = 0; i < n; i++) {
        centroid += points[i];
    }
    centroid *= (1. / nd);
    
    for (int i = 0; i < n; i++) {
        // tmp = M(0,0); tmp += (points[i].getX() - centroid.getX()) * (points[i].getX() - centroid.getX());
        // M(0,0) = tmp;

        // tmp = M(1,1); tmp += (points[i].getY() - centroid.getY()) * (points[i].getY() - centroid.getY());
        // M(1,1) = tmp;
        
        // tmp = M(2,2); tmp += (points[i].getZ() - centroid.getZ()) * (points[i].getZ() - centroid.getZ());
        // M(2,2) = tmp;

        // tmp = M(0,1); tmp += (points[i].getX() - centroid.getX()) * (points[i].getY() - centroid.getY());
        // M(0,1) = tmp; M(1,0) = tmp;
        
        // tmp = M(0,2); tmp += (points[i].getX() - centroid.getX()) * (points[i].getZ() - centroid.getZ());
        // M(0,2) = tmp; M(2,0) = tmp;

        // tmp = M(1,2); tmp += (points[i].getY() - centroid.getY()) * (points[i].getZ() - centroid.getZ());
        // M(1,2) = tmp; M(2,1) = tmp;
        
        
        
        
        
        tmp = M_(0,0); tmp += (points[i].getX() - centroid.getX()) * (points[i].getX() - centroid.getX());
        M_(0,0) = tmp;

        tmp = M_(1,1); tmp += (points[i].getY() - centroid.getY()) * (points[i].getY() - centroid.getY());
        M_(1,1) = tmp;
        
        tmp = M_(2,2); tmp += (points[i].getZ() - centroid.getZ()) * (points[i].getZ() - centroid.getZ());
        M_(2,2) = tmp;

        tmp = M_(0,1); tmp += (points[i].getX() - centroid.getX()) * (points[i].getY() - centroid.getY());
        M_(0,1) = tmp; M_(1,0) = tmp;
        
        tmp = M_(0,2); tmp += (points[i].getX() - centroid.getX()) * (points[i].getZ() - centroid.getZ());
        M_(0,2) = tmp; M_(2,0) = tmp;

        tmp = M_(1,2); tmp += (points[i].getY() - centroid.getY()) * (points[i].getZ() - centroid.getZ());
        M_(1,2) = tmp; M_(2,1) = tmp;
    }

//     cout << M << endl << endl;
    // SelfAdjointEigenSolver <Matrix3d> eigensolver(M);
    TMatrixDEigen eigensolver_(M_);
    
    // str.precision(17);
    // double a, b, c;
    // str << eigensolver.eigenvectors().col(0).row(0) << " " << eigensolver.eigenvectors().col(0).row(1) << " " << eigensolver.eigenvectors().col(0).row(2);
    // str >> a >> b >> c;
    // A = a; B = b; C = c;
    // cout << A<< " " << B << " "<< C << " " << D << endl;
    
    // We need the eigenvector of M_ which corresponds to the minimum eigenvalue
    TMatrixD eigenvalues_ = eigensolver_.GetEigenValues();
    unsigned mini = 0;
    if ( eigenvalues_(1,1) < eigenvalues_(0,0) ) {
        if ( eigenvalues_(2,2) < eigenvalues_(1,1) ) {
            mini = 2;
        }
        else {
            mini = 1;
        }
    }
    else if ( eigenvalues_(2,2) < eigenvalues_(0,0) ) {
        mini = 2;
    }
    
    // Got it!
    TMatrixD eigenvectors_  = eigensolver_.GetEigenVectors();
    
    mini = 2;
    
    A = eigenvectors_(0, mini);
    B = eigenvectors_(1, mini);
    C = eigenvectors_(2, mini);     
   
    
    
    if (A <= 0 && B <= 0 && C <= 0) {
        A = -A; B = -B; C = -C;
    }
    D = - (A * centroid.getX() + B * centroid.getY() + C * centroid.getZ());

    // Fixing the output back
    stdout_unhide(fd, pos);
}



// Helix fitting /////////////////////////////////////////////////////////////////////////////////

void helix_fit(vector <point> points, point& mid, const point n, double& r, double& step, bool& right,  point& on) {
    if ( points.size() < 3 ){
        return;
    }
    point normal = n.normalized();
    // 1 ///// Projecting points onto the plane and circular fitting //////////////////////////////////////
    // First we define a plane orthogonal to the axis:
    double A = normal.getX(), B = normal.getY(), C = normal.getZ(), D = 0;

    vector <point> points_p;

    for (unsigned i = 0; i < points.size(); i++) {
        points_p.push_back( project_onto_plane(points[i], A, B, C, D) );
    }
   

    circle_fit_3d(points_p, mid, normal, r); // Yes, I was lazy...
    // 2 ///// Projecting points onto the cylinder surface and linear fitting ////////////////////////////
    vector <double> z; // Projections of points onto axis
    vector <double> phir; // Projections of points onto circle
    
    z.push_back(0.);
    phir.push_back(0.);
    point yi, yii;
    double tmp;
    for (unsigned i = 1; i < points.size(); i++) {
        z.push_back( scalar_p(project_onto_line(points[i], mid, normal) - project_onto_line(points[0], mid, normal), normal) );
        yi  = project_onto_circle(points[i-1]  , mid, normal, r);
        yii = project_onto_circle(points[i], mid, normal, r);
        tmp =  angle(yi - mid, yii - mid) * r ;
        if ( scalar_p(vector_p(yi - mid, yii - mid), normal) < 0 ) {
            tmp = - tmp;
        }
        phir.push_back(phir.back() + tmp);
    }

    double zav = 0, phirav = 0;
    
    for (unsigned i = 0; i < points.size(); i++) {
        zav += z[i];
        phirav += phir[i];
    }
    
    zav /= ( (double) points.size() );
    phirav /= ( (double) points.size() );

    double k, b, eta;
    bool use;

    least_squares(z, phir, k, b, eta, use);

    if (k > 0) {
        right = true;
        step = 2. * m_pi * r / k;
    }
    else {
        right = false;
        step = - 2. * m_pi * r / k;
    }
    on = euler_rot(points[0], mid, normal, phirav / r) + normal * zav;
}

void stdout_hide(int &fd, fpos_t &pos) {
    fflush(stdout); 
    fgetpos(stdout, &pos);
    fd = dup(fileno(stdout));
    freopen("/dev/null", "a", stdout);
}

void stdout_unhide(const int &fd, const fpos_t &pos) {
    fflush(stdout);
    dup2(fd, fileno(stdout)); 
    close(fd);
    clearerr(stdout);
    fsetpos(stdout, &pos);
}


