#include "point.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>

using namespace std;

point::point()
    : m_x(),
      m_y(),
      m_z()
      //id()
{}

point::point(const double x, const double y, const double z)
    : m_x(x),
      m_y(y),
      m_z(z)
{}

void point::p_print(const unsigned pr) const
{
    std::streamsize ss = std::cout.precision();
    if (pr <= 17) {
        std::cout.precision(pr);
    }
    else {
        std::cout.precision(17);
    }
    std::cout << m_x << " " << m_y << " " << m_z << '\n';
    std::cout.precision (ss);
}

void point::p_fprint(const char * fname, const unsigned pr) const
{
    ofstream output;
    output.open(fname, ios::app);
    std::streamsize ss = output.precision();
    if (pr <= 17) {
        output.precision(pr);
    }
    else {
        output.precision(17);
    }
    output << m_x << " " << m_y << " " << m_z << '\n';
    output.precision(ss);
    output.close();
}

/*void point::p_id_fprint(const char * fname) const
{
    ofstream output;
    output.open(fname, ios::app);
    std::cout << id << " " << m_x << " " << m_y<< " " << m_z << '\n';
    output.close();
}*/

point& point::operator += (const point& rhs)
{
    m_x += rhs.m_x;
    m_y += rhs.m_y;
    m_z += rhs.m_z;
    return *this;
}

point& point::operator -= (const point& rhs)
{
    return operator += (-rhs);
}

point& point::operator*=( double d )
{
    m_x *= d;
    m_y *= d;
    m_z *= d;
    return *this;
}

point& point::operator/=( double d )
{
    m_x /= d;
    m_y /= d;
    m_z /= d;
    return *this;
}

double point::getX() const
{
    return m_x;
}

double point::getY() const
{
    return m_y;
}

double point::getZ() const
{
    return m_z;
}

/*int point::get_id() const
{
    return id;
}*/

double point::norm() const
{
    return std::sqrt(m_x*m_x + m_y*m_y + m_z*m_z);
}

const point point::normalized() const
{
    if ( eq(point(this->getX(), this->getY(), this->getZ()), point(0,0,0)) ) {
        return point(0,0,0);
    }
    point a = operator * (*this, 1/norm());
    return point(a.getX(), a.getY(), a.getZ());
}

const point operator - (const  point& p)
{
    return point(-p.getX(), -p.getY(), -p.getZ());
}

const point operator + (const point& p1, const point& p2)
{
    point temp(p1);
    temp += p2;
    return point(temp.getX(), temp.getY(), temp.getZ());
}
const point operator - (const point& p1, const point& p2)
{
    point temp(-p2);
    temp += p1;
    return point(temp.getX(), temp.getY(), temp.getZ());
}

const bool operator == (const point& p1, const point& p2)
{
    if ( (p1.getX() == p2.getX()) && (p1.getY() == p2.getY()) && (p1.getZ() == p2.getZ()) )
        return true;
    else
        return false;
}

const bool operator != (const point& p1, const point& p2)
{
    return !(p1 == p2);
}


const point operator*( double d, const point& p )
{
    point temp(p);
    temp *= d;
    return point(temp.getX(), temp.getY(), temp.getZ());
}

const point operator*( const point& p, double d )
{
    return operator * (d, p);
}

const point operator/( const point& p, double d )
{
    point tmp;
    tmp = p * (1. / d);
    return tmp;
}

const double angle(const point x, const point y)
{
    return acos((x.getX()*y.getX() + x.getY()*y.getY() + x.getZ()*y.getZ()) / sqrt((x.getX()*x.getX()+x.getY()*x.getY()+x.getZ()*x.getZ())*(y.getX()*y.getX()+y.getY()*y.getY()+y.getZ()*y.getZ())));
}

const double scalar_p(const point x, const point y)
{
    return x.getX()*y.getX() + x.getY()*y.getY() + x.getZ()*y.getZ();
}

const double distance_point_line (const point p, const point a, const point n) { // See wiki "Distance_from_a_point_to_a_line"
    point q, t1, t2;
    t1 = a - p;
    t2 = n * scalar_p(n, (a - p));
    q = (t1 - t2);
    return q.norm();
}

const point vector_p(const point x, const point y)
{
    double z1, z2, z3;
    z1 = x.getY()*y.getZ() - x.getZ()*y.getY();
    z2 = x.getZ()*y.getX() - x.getX()*y.getZ();
    z3 = x.getX()*y.getY() - x.getY()*y.getX();
    return point(z1, z2, z3);
}

const point project_onto_plane(const point q, const double A, const double B, const double C, const double D) {
    double tmp = - D - A * q.getX() - B * q.getY() - C * q.getZ();
    return (q + tmp * point(A, B, C));
}

const double distance_point_plane(const point p, const double A, const double B, const double C, const double D) {
    if ( fabs(A * A + B * B + C * C) > 1e-10 ) {
        return fabs(A * p.getX() + B * p.getY() + C * p.getZ() + D) / sqrt(A * A + B * B + C * C);
    }
    else {
        return 0;
    }
}

const double distance_point_circle (const point p, const point mid, const point normal, const double R) {
    point pp = project_onto_plane(p, normal.getX(), normal.getY(), normal.getZ(), - scalar_p(mid, normal));
    point pc = mid + (pp - mid).normalized() * R;
    return (p - pc).norm();
}

const point project_onto_line(const point s, const point a, const point n) {
    point b = a + 3 * n;
    point pr;
    pr = a + (b - a) * scalar_p((b - a), (s - a)) / ( ((b - a).norm()) * ((b - a).norm()));
    /*double t = ( (s.getX() - a.getX()) * (b.getX() - a.getX()) + 
                 (s.getY() - a.getY()) * (b.getY() - a.getY()) +
                 (s.getZ() - a.getZ()) * (b.getZ() - a.getZ()) ) /
               ( (b.getX() - a.getX()) * (b.getX() - a.getX()) +
                 (b.getY() - a.getY()) * (b.getY() - a.getY()) +
                 (b.getZ() - a.getZ()) * (b.getZ() - a.getZ()) );
    pr = a + (b - a) * t;*/
    return pr;
}

const point euler_rot(const point p, const point n, const double t) {
    double  x = n.getX(),  y = n.getY(),  z = n.getZ();
    double px = p.getX(), py = p.getY(), pz = p.getZ();
    double qx, qy, qz;
    qx =  (cos(t) + (1 - cos(t)) * x * x) * px +
         ((1 - cos(t)) * x * y - sin(t) * z) * py +
         ((1 - cos(t)) * x * z + sin(t) * y) * pz;
    qy = ((1 - cos(t)) * y * x + sin(t) * z) * px +
          (cos(t) + (1 - cos(t)) * y * y) * py +
         ((1 - cos(t)) * y * z - sin(t) * x) * pz;
    qz = ((1 - cos(t)) * z * x - sin(t) * y) * px +
         ((1 - cos(t)) * z * y + sin(t) * x) * py +
          (cos(t) + (1 - cos(t)) * z * z) * pz;
    return point(qx, qy, qz);
}

const point euler_rot(const point p, const point a, const point n, const double angle) {
    point delta = project_onto_line(point(0,0,0), a, n);
    return ( delta + euler_rot((p - delta), n, angle) );
}

const point project_onto_circle(const point q, const point mid, const point normal, const double R) {
    point pp = project_onto_plane(q, normal.getX(), normal.getY(), normal.getZ(),
                                  - scalar_p(mid, normal));  // The projection of p onto the plane
    pp = mid + (pp - mid).normalized() * R;
    return pp;
}

const point project_onto_cylinder(const point q, const point mid, const point normal, const double R) {
    point p = project_onto_line(q, mid, normal);
    return ( p +(q - p).normalized() * R );
}

/*const point project_onto_helix(const point q, const point mid, const point normal, const double r, const double step, const bool right, const point on) {
    point n = normal.normalized();
    //if ( eq(project_onto_line(q, mid, n), q) ) {
        

        //return ...;
    //}

    double dist = 1e300, distc;
    point outp, outc;
    point qr;
    const unsigned it = 5;
    double itd = (double) it;
    for (unsigned i = 0; i < it; i++ ) {
        double sign = (((double) right) - .5) * 2.;

        //qr = euler_rot(q, mid, n, 2 * m_pi / it * ((double) i));

        point qp = project_onto_cylinder(q, mid, n, r);
        double A, B, C, D; // Plane orthogonal to normal and containing qp
        A = n.getX(), B = n.getY(), C = n.getZ();
        D = - A * qp.getX() - B * qp.getY() - C * qp.getZ();

        point onp = project_onto_plane(on, A, B, C, D);
        
        double moveon;
        double ry = (onp - on).norm();
        modf((onp - on).norm() / step, &moveon);
        point onc = on;

        if ( (on + n * moveon * step - qp).norm() < (on - qp).norm() ) {
            onc = on + n * moveon * step;
        }
        if ( (on - n * moveon * step - qp).norm() < (on - qp).norm() ) {
            onc = on - n * moveon * step;
        }

        double dz = scalar_p(project_onto_line(qp, mid, n) - project_onto_line(onc, mid, n), n);
        double dphi =  angle(project_onto_circle(onc , mid, n, r) - mid, project_onto_circle(qp, mid, n, r) - mid) * r * sign;
        vector <double> dphis;
        dphis.push_back(dphi);
        dphis.push_back(-dphi);
        if ( sign > 0 ) {
            dphis.push_back(dphi + 2*m_pi);
            dphis.push_back(-dphi - 2*m_pi);
        }
        else {
            dphis.push_back(dphi - 2*m_pi);
            dphis.push_back(-dphi + 2*m_pi);
        }
        //if ( scalar_p(vector_p(project_onto_circle(onc , mid, n, r) - mid, project_onto_circle(qp, mid, n, r) - mid), n * sign) < 0 || dz < 0) {
            //dphi = - dphi;
        //}
        /*if ( (project_onto_line(onc, mid, n) - project_onto_line(qp, mid, n)).norm() > step / 2 ) {
            if (sign > 0){
                if (dz < 0) {
                    dphi = -(dphi + 2*m_pi);
                }
                else {
                    dphi = 2*m_pi + dphi;
                }
            }
            else {
                if (dz < 0) {
                    dphi = - dphi + 2*m_pi;
                }
                else {
                    dphi =  -(2*m_pi - dphi);
                }
            }
        }
        
        for (unsigned j = 0; j < dphis.size(); j++) {
            dphi = dphis[j];
            point tmp = project_onto_line(point(dz, dphi, 0), point(0, 0, 0), point(1, 2 * m_pi * r / step * sign, 0));
            double movez = tmp.getX();
            double movephi = tmp.getY();
            outc = euler_rot(onc, mid, normal, movephi / r) + n * movez;
            distc = ( q - outc ).norm();
            if ( distc < dist ) {
                outp = outc;
                dist = distc;
            }
        }
        //distc = (outc - q).norm();
        //if ( distc < dist ) {
        //     outp = outc;
        //}
    //}

    return outp;
}*/

//void plane_3_points(const point a, const point b, const point c, double& A, double& B, double& C, double& D) {
//    if ( eq(angle(b - a, b - c), 0) || eq(angle(b - a, b - c), m_pi) ) {
//        A = 0, B = 0, C = 0, D = 0;
//        return;
//    }
//    
//}

const point average(const vector <point>& points) {
    if (points.size() == 0) {
        return point(0,0,0);
    }
    else {
        point av(0,0,0);
        for (unsigned i = 0; i < points.size(); i++) {
            av += points[i];
        }
        av /= ((double) points.size());
        return av;
    }
}

const double chi_2_lin(const std::vector <point>& points, const point a, const point n) {
    double tmp, chi = 0;
    for (unsigned i = 0; i < points.size(); i++) {
        tmp = distance_point_line(points[i], a, n);
        chi += ( tmp * tmp );
    }
    return chi;
}

const double chi_2_circ(const std::vector <point>& points, const point mid, const point normal, const double R) {
    double tmp, chi = 0;
    for (unsigned i = 0; i < points.size(); i++) {
        tmp = distance_point_circle(points[i], mid, normal, R);
        chi += ( tmp * tmp );
    }
    return chi;
}

const point tangent_circle_point(const point mid, const point normal, const double R, const point p, const point inner) {
    point pp = project_onto_circle(p, mid, normal, R);
    point k = ( vector_p(pp - mid, normal) ).normalized();
    if ( angle(k, (inner - pp)) < m_pi / 2 ) {
        k = -k;
    }
    return k;
}

const point tangent_helix_point(const point q, const point mid, const point normal, const double r, const double step, const bool right, const point on, const point inner) {
    if ( eq(step, 0) ) {
        return tangent_circle_point(mid, normal, r, q, inner);
    }
    
    point qp = project_onto_helix(q, mid, normal, r, step, right, on);
    if ( eq(qp, q) ) {
        qp = q;
    }
    else {
        // cout << "Warning: the point was projected to the helix during calculating tangent direction!\n";
    }
    
    point n = normal.normalized();
    point qa = project_onto_line(qp, mid, n);
    point n1 = tangent_circle_point(qa, n, r, qp, inner); // Component in the plane going through q and orthogonal to axis
    point n2 = n * ( step / 2. / m_pi / r );
    
    if ( angle(n1, n2) < ( m_pi / 2 ) ) {
        n2 = - n2;
    }

    point k = n1 + n2;

    if ( ( scalar_p((qp-qa), vector_p(n1, k)) > 0 ) != right) {
        k = n1 - n2;
    }

    if ( angle(k, qp - inner) > (m_pi / 2.) ) {
        k = -k;
    }
    
    return k;

}

const bool eq(const double a, const double b) {
    if ( fabs(a-b) < def_delta ) {
        return true;
    }
    else {
        return false;
    }
}

const bool eq(const point& a, const point& b) {
    if ( (a - b).norm() < def_delta ) {
        return true;
    }
    else {
        return false;
    }
}

const bool collinear(const point& a, const point& b) {
    if ( ( eq(a, b) ) || ( eq(a, -b) ) ) {
        return true;
    }
    else {
        return false;
    }
}














const point project_onto_helix(const point p, const point mid, const point normal, const double r, const double step, const bool right, const point on) {
    if ( eq(step, 0.) ) {
        return project_onto_circle(p, mid, normal, r);
    }

    // Forming a basis:
    point n3 = normal.normalized();
    point O = project_onto_line(on, mid, n3);
    point P = project_onto_line(p, mid, n3);
    point n1 = (on - O).normalized();
    point n2 = vector_p(n3, n1);
    point pm = p - O;
    double px = scalar_p(n1, pm), py = scalar_p(n2, pm), pz = scalar_p(n3, pm); // Coordinates in new basis
    point Pn(px, py, pz);
    double dist_on_axis = scalar_p(P - O, n3);
    double h = step / 2. / m_pi;
    double t_mid = dist_on_axis / h;
    const unsigned t_num = 5;
    double t_num_d = (double) t_num;
    double t_step = 2. * m_pi / ( t_num_d - 1);
    double sign = (((double) right) - .5) * 2.;

    double dist_best = 1e300, dist_curr;
    point pbest, pcurr;

    double tb, t, tn, f, fp;
    const unsigned maxstep = 5;
    unsigned curr_step = maxstep;

    for (unsigned i = 0; i < t_num; i++) { // Here we try few starting points, because the function has few minima
        tb = t_mid - m_pi + t_step * ((double) i);
        tn = tb;
        do {
            t = tn;
            f = r * ( px * sin(t) - py * sign * cos(t) ) +  h * h * t - pz * h;
            fp = r * ( px * cos(t) + py * sign * sin(t) ) + h * h;
            tn = t - f / fp;
            curr_step--;
        } while ( ( !eq(t, tn) ) && (curr_step > 0) );
        pcurr = point_helix_parametric(point(1,0,0), point(0,1,0), point(0,0,1), r, step, right, tn);
        dist_curr = ( pcurr - Pn ).norm();
        if ( dist_curr < dist_best ) {
            dist_best = dist_curr;
            pbest = pcurr;
        }
        curr_step = maxstep;
    }

    // Now we must move to old coordinates

    return ( O + n1 * pbest.getX() + n2 * pbest.getY() + n3 * pbest.getZ() );    
}

const bool ComparePoints(const point& p1, const point& p2) {
    if ( p1.getX() != p2.getX() ) {
        return ( p1.getX() < p2.getX() );
    }
    if ( p1.getY() != p2.getY() ) {
        return ( p1.getY() < p2.getY() );
    }
    return ( p1.getZ() < p2.getZ() );
}

const point point_helix_parametric(const point n1, const point n2, const point n3, const double r, const double step, const bool right, const double t) {
    double sign = (((double) right) - .5) * 2.;
    return ( n1 * r * cos(t) + n2 * r * sign * sin(t) + n3 * step / 2. / m_pi * t );
}

