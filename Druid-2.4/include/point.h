#ifndef _POINT_H_
#define _POINT_H_

#include <vector>
#include <iostream>

// Settings ///////////////////////////////////////////////////////////////////////////////////////////////////////////

#define m_pi 3.1415926535897932 // :)

#define def_delta 1e-11 // Used in the eq(...) functions

// In all functions the equation of plane in 3D is: A * x + B * y + C * z + D = 0
// With A ^ 2 + B ^ 2 + C ^ 2 = 1
// So, point(A, B, C) is the normal vector to the plane


class point // The class is used both as point and as a vector
{
    double m_x;
    double m_y;
    double m_z;
public:
    point();
    point(const double x, const double y, const double z);

    void p_print(const unsigned pr = std::cout.precision()) const; // Here one can set the precision of the output (17 max)
    void p_fprint(const char * fname, const unsigned pr = std::cout.precision()) const;

    point& operator += (const point& rhs);
    point& operator -= (const point& rhs);

    point& operator *= (double d);
    point& operator /= (double d);

    double getX() const;
    double getY() const;
    double getZ() const;
    double norm() const;
    const point normalized() const;
};

const bool operator == (const point& p1, const point& p2);
const bool operator != (const point& p1, const point& p2);

const point operator - (const point& p);

const point operator + (const point& p1, const point& p2);
const point operator - (const point& p1, const point& p2);

const point operator * (double d, const point& p);
const point operator * (const point& p, double d);

const point operator / (const point& p, double d);

const double angle(const point x, const point y); // In radians!
const double scalar_p(const point x, const point y);

const point vector_p(const point x, const point y);

const double distance_point_line (const point p, const point a, const point n); // Line is expressed as x = a + t * n
const double distance_point_circle (const point p, const point mid, const point normal, const double R);
const double distance_point_plane(const point p, const double A, const double B, const double C, const double D);

const point project_onto_plane(const point q, const double A, const double B, const double C, const double D);
const point project_onto_line(const point q, const point a, const point n);
const point project_onto_circle(const point q, const point mid, const point normal, const double R);
const point project_onto_cylinder(const point q, const point mid, const point normal, const double R);
const point project_onto_helix(const point q, const point mid, const point normal, const double r, const double step, const bool right, const point on);

//void plane_3_points(const point a, const point b, const point c, double& A, double& B, double& C, double& D);

const point euler_rot(const point p, const point n, const double angle); // Angle in radians! Axis from origin
const point euler_rot(const point p, const point a, const point n, const double angle); // Angle in radians! Axis through a

const point tangent_circle_point(const point mid, const point normal, const double R, const point p, const point inner); // The vector points 'outside' the arc (that's why we need the inner point)
const point tangent_helix_point(const point q, const point mid, const point normal, const double r, const double step, const bool right, const point on, const point inner);
const point average(const std::vector <point>& points);

// Chi^2 WITHOUT diving by N or (N-1) or smth else:
const double chi_2_lin(const std::vector <point>& points, const point a, const point n);
const double chi_2_circ(const std::vector <point>& points, const point mid, const point normal, const double R);

const bool eq(const double a, const double b);
const bool eq(const point& a, const point& b);
const bool collinear(const point& a, const point& b);

template <typename T> int sgn(T val) { // Signum function for number types
    return (T(0) < val) - (val < T(0));
}

// To compare two vectors of points, first sort each of them using the command
// std::sort(v.begin(), v.end(), ComparePoints);
// Where ComparePoints is:

const bool ComparePoints(const point& p1, const point& p2);

const point point_helix_parametric(const point n1, const point n2, const point n3, const double r, const double step, const bool right, const double t);

#endif

