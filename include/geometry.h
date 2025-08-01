#ifndef _GEOMETRY_H_
#define _GEOMETRY_H_

#include <vector>

#include "point.h"
// #include "fitting.h"
#include "fitting_root.h"

using namespace std;

// Settings ///////////////////////////////////////////////////////////////////////////////////////////////////////////


// chi^2 defaults:
// The sum of squares is divided by (N - M), where M is:

#define def_figure_chi2 3
#define def_line_chi2 2
#define def_plane_chi2 3
#define def_circle_chi2 2
#define def_cylinder_chi2 3
#define def_helix_chi2 2

// For each class please define following methods (if projecting of point is defined):
// unsigned get_min_fit_points() const; // If different from 2
// virtual const point project(const point q) const;
// And line:
// using figure::project;
// And also the comparison function:
// const boo eq(<your_figure_class> f1, <your_figure_class> f2);

class figure {
public:
    virtual const unsigned get_min_fit_points() const; // The default value is 2
    virtual const point project(const point q) const = 0;
    
    const double chi2(const vector <point>& points) const;
    const vector <point> project(const vector <point>& points) const;
};



class line : public figure {
    point a; // Line goes through this point
    point n; // The direction
public:
    line();
    line(const point sa, const point sn);
    line(const vector <point>& points, fit_type lft = def_lft); // Fitting

    using figure::project;
    const unsigned get_min_fit_points() const;

    const point get_a() const;
    const point get_n() const;

    const vector <point> points_to_draw(const point p1, const point p2, const unsigned n) const;
    const point project(const point q) const;
};

const line line_2_points(const point a1, const point a2);
const bool eq(const line& l1, const line& l2);



class plane : public figure {
    double A, B, C, D;
public:
    plane();
    plane(const point p1, const point p2, const point p3);
    plane(const double cA, const double cB, const double cC, const double cD);
    plane(const point p, const point n); // Through p, orthogonal to n
    plane(const vector <point>& points);

    using figure::project;
    virtual const unsigned get_min_fit_points() const;
    
    const double get_A() const;
    const double get_B() const;
    const double get_C() const;
    const double get_D() const;
    const point  get_n() const;

    const point project(const point q) const;
};

const bool eq(const plane& p1, const plane& p2);



class circle : public figure {
    line axis;
    double r;
public:
    circle();
    circle(const line caxis, const double cr);
    circle(const point cmid, const point cn, const double cr);
    circle(const vector <point>& points); // Fitting

    using figure::project;
    virtual const unsigned get_min_fit_points() const;

    const line get_axis() const;
    const point get_mid() const;
    const point get_n() const;
    const double get_r() const;

    const point tangent(const point q) const;
    const point tangent(const point q, const point in) const;

    const vector <point> points_to_draw(const point p1, const point p2, const unsigned n, const bool small) const; // small == 1  ---  draw the small arc

    const point project(const point q) const;
};

const bool eq(const circle& c1, const circle& c2);



class cylinder : public figure {
    circle cc;
public:
    cylinder();
    cylinder(const circle ccc);
    cylinder(const point cmid, const point cn, const double cr);
    cylinder(const vector <point>& points, const point n); // Fitting

    using figure::project;
    virtual const unsigned get_min_fit_points() const;
    
    const line get_axis() const;
    const circle get_cc() const;
    const point get_mid() const;
    const point get_n() const;
    const double get_r() const;

    const point project(const point q) const;
};

const bool eq(const cylinder& c1, const cylinder& c2);



class helix : public figure {
    cylinder hc;
    double step;
    bool right; // Right-handed screw?
    point on;
public:
    helix();
    helix (const cylinder chc, const double cstep, const bool cright, const point con);
    helix (const circle cc, const double cstep, const bool cright, const point con);
    helix(const point cmid, const point cn, const double cr, const double cstep, const bool cright, const point con);
    helix(const vector <point>& points, const point n); // Fitting
    helix(const point con, const point tangent, const double r, const point n, const bool cright); // 'tangent' taken at point 'on'

    using figure::project;
    virtual const unsigned get_min_fit_points() const;

    const cylinder get_hc() const;
    const line get_axis() const;
    const circle get_cc() const;
    const point get_mid() const;
    const point get_n() const;
    const double get_r() const;
    const double get_step() const;
    const bool get_right() const;
    const point get_on() const;

    const double chi2(const vector <point>& points) const;

    const point intersect_plane(const point q) const; // Intersection with plane going through q and orthogonal to the axis

    const point tangent(const point q) const;
    const point tangent(const point q, const point in) const; // 'in' defines normal direction

    const vector <point> points_to_draw(const point p1, const point p2, const unsigned n) const;
    const point project(const point q) const;
};

const bool eq(const helix& h1, const helix& h2);

#endif

