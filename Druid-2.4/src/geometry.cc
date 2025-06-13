#include <vector>
#include <iostream>
#include <cmath>

#include "geometry.h"
// #include "fitting.h"
#include "fitting_root.h"

using namespace std;



const double figure::chi2(const vector <point>& points) const {
    double chi = 0;
    unsigned s = points.size();
    
    point x, y, z;
    for (unsigned i = 0; i < points.size(); i++) {
        x = points[i];
        y = this -> project(points[i]);
        z = x - y;
        chi += ( z.getX() * z.getX() + z.getY() * z.getY() + z.getZ() * z.getZ() );
    }
    
    double DEN;
    if ( s <= (this -> get_min_fit_points() ) ) {
        DEN = 1;
    }
    else {
        DEN = s - ( (double) (this -> get_min_fit_points() ) );
    }

    chi /= DEN;

    return chi;
}

const vector <point> figure::project(const vector <point>& points) const {
    vector <point> ppoints;
    for (unsigned i = 0; i < points.size(); i++) {
        ppoints.push_back( this->project(points[i]) );
    }
    return ppoints;
}

const unsigned figure::get_min_fit_points() const {
    return def_figure_chi2;
}

const unsigned line::get_min_fit_points() const {
    return def_line_chi2;
}

line::line()
    : a(),
      n()
{}

line::line(const point sa, const point sn) {
    a = sa;
    n = sn.normalized();
}

line::line(const vector <point>& points, fit_type lft) {
    if ( points.size() < 2 ) {
        a = point(0,0,0);
        n = point(0,0,0);
        cout << "Error: not enough points to fit the line.\n";
        return;
    }
    // First we check if there are at least 2 different points
    bool ok = false;
    for (unsigned i = 0; i < points.size(); i++) {
        for (unsigned j = i + 1; j < points.size(); j++) {
            if ( ! eq(points[i], points[j]) ) {
                ok = true;
                break;
            }
        }
    }
    if ( ok == false ) {
        a = point(0,0,0);
        n = point(0,0,0);
        cout << "Error: at least 2 points must be different.\n";
        return;
    }

    linear_fit(points, n, lft);
    a = average(points);
}

const point line::get_a() const {
    return a;
}

const point line::get_n() const {
    return n;
}

const point line::project(const point q) const {
    return project_onto_line(q, a, n);
}

const line line_2_points(const point a1, const point a2) {
    if ( eq(a1, a2) ) {
        cout << "Error: points must be different.\n";
        return line(point(0,0,0), point(0,0,0));
    }
    else {
        return line(a1, (a2 - a1).normalized());
    }
}

const bool eq(const line& l1, const line& l2) {
    point n1 = l1.get_n(), n2 = l2.get_n();
    point a1 = l1.get_a(), a2 = l2.get_a();
    if ( ! (eq(n1, n2) || eq(n1, -n2)) ) {
        return false;
    }

    if ( eq(a1, a2 ) ) {
        return true;
    }

    point a1a2 = (a2 - a1).normalized();

    if ( eq(a1a2, n1) || eq(a1a2, -n1) ) {
        return true;
    }
    else {
        return false;
    }

}



const unsigned plane::get_min_fit_points() const {
    return def_plane_chi2;
}

plane::plane()
    : A(),
      B(),
      C(),
      D()
{}

plane::plane(const point p, const point q, const point s) {
    double px = p.getX(), py = p.getY(), pz = p.getZ();
    double qx = q.getX(), qy = q.getY(), qz = q.getZ();
    double sx = s.getX(), sy = s.getY(), sz = s.getZ();
    double d = px * qy * sz - px * qz * sy + py * qz * sx - py * qx * sz + pz * qx * sy - pz * qy * sx;
    if ( eq(d, 0.) ) {
        A = 0, B = 0, C = 0, D = 0;
        cout << "Error: points must not lie on the same line.\n";
        return;
    }
    
    A = - (1 * qy * sz - 1 * qz * sy + py * qz * 1 - py * 1 * sz + pz * 1 * sy - pz * qy * 1) / d;
    B = - (px * 1 * sz - px * qz * 1 + 1 * qz * sx - 1 * qx * sz + pz * qx * 1 - pz * 1 * sx) / d;
    C = - (px * qy * 1 - px * 1 * sy + py * 1 * sx - py * qx * 1 + 1 * qx * sy - 1 * qy * sx) / d;

    double norm = sqrt(A * A + B * B + C * C);
    A /= norm;
    B /= norm;
    C /= norm;
    D = - A * px - B * py - C * pz;
    
}

plane::plane(const double cA, const double cB, const double cC, const double cD) {
    if (eq(cA,0) && eq(cB,0) && eq(cC,0)) {
        cout << "Error: at least one of A, B, C should be nonzero.\n";
        A = 0, B = 0, C = 0, D = 0;
        return;
    }
    A = cA, B = cB, C = cC, D = cD;
}

plane::plane(const point p, const point n) {
    point nn = n.normalized();
    A = nn.getX();
    B = nn.getY();
    C = nn.getZ();
    D = - A * p.getX() - B * p.getY() - C * p.getZ();
}

plane::plane(const vector <point>& points) {
    unsigned s = points.size();
    if (s < 3) {
        cout << "Error: not enough points to fit the plane.\n";
        return;
    }
    else if (s == 3) {
        plane tmp(points[0], points[1], points[2]);
        A = tmp.get_A(); B = tmp.get_B(); C = tmp.get_C(); D = tmp.get_D();
    }
    else {
        plane_fit(points, A, B, C, D);
    }
}

const double plane::get_A() const {
    return A;
}

const double plane::get_B() const {
    return B;
}

const double plane::get_C() const {
    return C;
}

const double plane::get_D() const {
    return D;
}

const point  plane::get_n() const {
    return point(A, B, C);
}

const point plane::project(const point q) const {
    return project_onto_plane(q, A, B, C, D);
}

const bool eq(const plane& p1, const plane& p2) {
    double A1 = p1.get_A(), B1 = p1.get_B(), C1 = p1.get_C(), D1 = p1.get_D();
    double A2 = p2.get_A(), B2 = p2.get_B(), C2 = p2.get_C(), D2 = p2.get_D();
    if ( (eq(A1, A2) && eq(B1, B2) && eq(C1, C2) && eq(D1, D2)) || (eq(A1, -A2) && eq(B1, -B2) && eq(C1, -C2) && eq(D1, -D2)) ) {
        return true;
    }
    else {
        return false;
    }
}


const unsigned circle::get_min_fit_points() const {
    return def_circle_chi2;
}

circle::circle()
    : axis(),
      r()
{}

circle::circle(const line caxis, const double cr)
    : axis(caxis),
      r(cr)
{}


circle::circle(const point cmid, const point cn, const double cr) {
    axis = line(cmid, cn.normalized());
    r = cr;
}

circle::circle(const vector <point>& points) {
    if ( eq(line(points).chi2(points), 0) ) {
        axis = line(point(0,0,0), point(0,0,0));
        r = 0;
    }
    else {
        point mid, n;
        circle_fit_3d(points, mid, n, r);
        axis = line(mid, n);
    }
}

const line circle::get_axis() const {
    return axis;
}

const point circle::get_mid() const {
    return axis.get_a();
}

const point circle::get_n() const {
    return axis.get_n();
}

const double circle::get_r() const {
    return r;
}

const point circle::project(const point q) const {
    return project_onto_circle(q, axis.get_a(), axis.get_n(), r);
}

const point circle::tangent(const point q) const {
    return tangent_circle_point(axis.get_a(), axis.get_n(), r, q, point(q.getX() + 10, 0, 0));
}

const point circle::tangent(const point q, const point in) const {
    return tangent_circle_point(axis.get_a(), axis.get_n(), r, q, in);
}

const bool eq(const circle& c1, const circle& c2) {
    if ( eq(c1.get_axis(), c2.get_axis()) && eq(c1.get_r(), c2.get_r()) ) {
        return true;
    }
    else {
        return false;
    }
}



const unsigned cylinder::get_min_fit_points() const {
    return def_cylinder_chi2;
}

cylinder::cylinder()
    : cc()
{}

cylinder::cylinder(const circle ccc)
    : cc(ccc)
{}

cylinder::cylinder(const point cmid, const point cn, const double cr) {
    cc = circle(cmid, cn, cr);
}

cylinder::cylinder(const vector <point>& points, const point n) {
    plane pl(average(points), n);
    vector <point> ppoints;
    ppoints = (pl).project(points);
    circle circ(ppoints);
    cc = circle(circ.get_mid(), n, circ.get_r());

}

const line cylinder::get_axis() const {
    return cc.get_axis();
}

const circle cylinder::get_cc() const {
    return cc;
}

const point cylinder::get_mid() const {
    return cc.get_mid();
}

const point cylinder::get_n() const {
    return cc.get_n();
}

const double cylinder::get_r() const {
    return cc.get_r();
}

const point cylinder::project(const point q) const {
    return project_onto_cylinder(q, get_mid(), get_n(), get_r());
}

const bool eq(const cylinder& c1, const cylinder& c2) {
    if ( eq(c1.get_r(), c2.get_r()) && eq(c1.get_axis(), c2.get_axis()) ) {
        return true;
    }
    else {
        return false;
    }
}



const unsigned helix::get_min_fit_points() const {
    return def_helix_chi2;
}

helix::helix()
    : hc(),
      step(),
      right(),
      on()
{}

helix::helix (const cylinder chc, const double cstep, const bool cright, const point con) {
    point pcon = chc.project(con);
    if ( !eq(pcon, con) ) {
        cout << "Warning: 'on' point was projected onto the cylinder.\n";
    }
    else {
        pcon = con;
    }
    hc = chc;
    step = cstep;
    right = cright;
    on = pcon;
}

helix::helix (const circle cc, const double cstep, const bool cright, const point con) {
    hc = cylinder(cc);

    point pcon;
    if ( eq(cstep, 0) ) {
        step = 0;
        pcon = cc.project(con);
        if ( !eq(pcon, con) ) {
            // cout << "Warning: 'on' point was projected on the circle (step == 0).\n";
        }
        else {
            pcon = con;
        }
    }
    else {
        if ( !eq(pcon, con) ) {
            pcon = hc.project(con);
            // cout << "Warning: 'on' point was projected onto the cylinder.\n";
        }
        else {
            pcon = con;
        }
    }
    step = cstep;
    right = cright;
    on = pcon;
}

helix::helix(const point cmid, const point cn, const double cr, const double cstep, const bool cright, const point con) {
    hc = cylinder(cmid, cn, cr);

    point pcon;
    if ( eq(cstep, 0) ) {
        step = 0;
        pcon = project_onto_circle(con, cmid, cn, cr);
        if ( !eq(pcon, con) ) {
            // cout << "Warning: 'on' point was projected on the circle (step == 0).\n";
        }
        else {
            pcon = con;
        }
    }
    else {
        pcon = project_onto_cylinder(con, cmid, cn, cr);
        if ( !eq(pcon, con) ) {
            pcon = hc.project(con);
            // cout << "Warning: 'on' point was projected onto the cylinder.\n";
        }
        else {
            pcon = con;
        }
    }
    step = cstep;
    right = cright;
    on = pcon;
}

helix::helix(const vector <point>& points, const point n) {
    point mid;
    double r;
    plane check(points);
    bool bad = true;
    if ( points.size() < 3 ) {
        cout << "Not enough points to build helix!\n";
        return;
    }

    for (unsigned i = 0; i < points.size(); i++) {
        if ( !eq(check.project(points[i]), points[i]) ) {
            bad = false;
            break;
        }
    }
    
    bad = bad && collinear(n, check.get_n());

    if ( bad == false ) {
        helix_fit(points, mid, n, r, step, right, on);
    }
    else {
        point ptrash;
        circle_fit_3d(points, mid, ptrash, r);
        right = 0;
        step = 0;
        on = project_onto_circle(points[0], mid, n, r);
    }
    hc = cylinder(mid, n, r);
}

helix::helix(const point con, const point t, const double r, const point normal, const bool cright) {
    point n = normal.normalized();
    point mid;

    if ( eq(vector_p(normal, t), point(0,0,0)) ) {
        cout << "Warning: step == 0; normal is sign-sensitive. Treated as magnetic field direction.\n";
        step = 0;
        right = 0;
        mid = con + ( vector_p(t, normal) ).normalized() * r;
        
    }
    else {
        point tt;
        if ( angle(t, n) > (m_pi / 2.) ) {
            tt = - t;
        }
        else {
            tt = t;
        }
        point mover = ( vector_p(tt, n) ).normalized() * r;
        mid = con + mover;
        point vp = vector_p(tt, n);
        if ( (scalar_p(con - mid, vp) > 0) != cright ) {
            mid = con - mover;
        }

        step = fabs(scalar_p(t, n) / scalar_p(t, ( vector_p(n, con - mid) ).normalized())) * 2 * m_pi * r;
        right = cright;
    }
    hc = cylinder(mid, n, r);
    
    on = con;

}

const cylinder helix::get_hc() const {
    return hc;
}

const line helix::get_axis() const {
    return hc.get_axis();
}

const circle helix::get_cc() const {
    return hc.get_cc();
}

const point helix::get_mid() const {
    return hc.get_mid();
}

const point helix::get_n() const {
    return hc.get_n();
}

const double helix::get_r() const {
    return hc.get_r();
}

const double helix::get_step() const {
    return step;    
}

const bool helix::get_right() const {
    return right;
}

const point helix::get_on() const {
    return on;
}

const point helix::intersect_plane(const point q) const {
    return project(get_axis().project(q));
}

const point helix::tangent(const point q) const {
    return tangent_helix_point(q, get_mid(), get_n(), get_r(), get_step(), get_right(), get_on(), point(q.getX() + 10, q.getY(), q.getZ()));
}

const point helix::tangent(const point q, const point in) const {
    return tangent_helix_point(q, get_mid(), get_n(), get_r(), get_step(), get_right(), get_on(), in);
}

const point helix::project(const point q) const {
    return project_onto_helix(q, get_mid(), get_n(), get_r(), step, right, on);
}

const bool eq(const helix& h1, const helix& h2) {
    if ( (!eq(h1.get_hc(), h2.get_hc())) ||
         (!eq(h1.get_step(), h2.get_step())) ||
         (!eq(h1.get_right(), h2.get_right())) ||
         (!eq(h1.get_on(), h2.project(h1.get_on())))) {
        return false;
    }
    return true;
}



const vector <point> helix::points_to_draw(const point p1, const point p2, const unsigned n) const { 
    vector <point> opoints;
    if ( n < 2 ) {
        return opoints;
    }
    point be = get_axis().project(project(get_hc().project(p1)));
    point en = get_axis().project(project(get_hc().project(p2)));
    //point en = get_axis().project(p2);
    if ( eq(step, 0.) ) {
        return get_cc().points_to_draw(p1, p2, n,  true);
    }
    else {
        point ste = ( en - be ) / ( (double) (n - 1) );
        if ( eq((be - en).norm(), 0) ) {
            return circle(be, get_n(), get_r()).points_to_draw(p1, p2, n, true);
        }
        for (unsigned i = 0; i < n; i++) {
            opoints.push_back( project( be + ste * ((double) i) ) );
        }
        return opoints;
    }
}

const vector <point> line::points_to_draw(const point p1, const point p2, const unsigned n) const { 
    vector <point> opoints;
    if ( n < 2 ) {
        return opoints;
    }
    point be = project(p1);
    point en = project(p2);
    point ste = ( en - be ) / ( (double) (n - 1) );
    
    for (unsigned i = 0; i < n; i++) {
        opoints.push_back( be + ste * ((double) i) );
    }
    return opoints;
}

const vector <point> circle::points_to_draw(const point p1, const point p2, const unsigned n, const bool small) const { 
    vector <point> opoints;
    if ( n < 2 ) {
        return opoints;
    }
    point be = project(p1);
    point en = project(p2);

    double an = angle(be - get_mid(), en - get_mid());
    if ( small == false ) {
        an = 2 * m_pi - an;
    }
    point normal;
    if ( scalar_p(vector_p(be - get_mid(), en - get_mid()), get_n())  < 0) {
        normal = -get_n();
    }
    else {
        normal = get_n();
    }
    for (unsigned i = 0; i < n; i++) {
        opoints.push_back( euler_rot(be, get_mid(), normal, an / ((double) (n-1) ) * ((double) i )));
    }
    return opoints;
}

const double helix::chi2(const vector <point>& points) const {
    double A, B, C, D;
    plane_fit(points, A, B, C, D);
    if ( eq(angle(point(A,B,C), get_n()), 0) || eq(angle(point(A,B,C), -get_n()), 0)) {
        return circle(points).chi2(points);
    }

    return figure::chi2(points);
}
