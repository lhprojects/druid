//#include <vector>

#include "point.h"
#include "geometry.h"
// #include "fitting.h"
#include "fitting_root.h"
#include "segment3.h"

using namespace std;

// Constructors ///////////////////////////////////////////////////////////////////////////////////////////////////////

segment::segment()
    : m_points(),
      m_line(),
      m_helix(),
      ft(def_fit_type),
      dir_begin(),
      dir_end(),
      m_err()
{}

segment::segment(const vector <point>& in_points, const fit_type ftin, const point n)
    : m_points(in_points),
      m_helix(),
      ft(),
      dir_begin(),
      dir_end(),
      m_err()
{
    if ( ( ftin == HEL ) && ( eq(n, point(0,0,0)) ) ) {
        ft = def_fit_type;
        ADDED(point(0,0,0));
        //m_line = line(in_points);
        return;
    }

    ft = ftin;
    
    /*if ( is_linear(ft) ) {
        m_line = line(in_points, ft);
    }
    else {
        m_line = line(in_points);
    }
    m_err = m_line.chi2(m_points);
    if ( m_err > def_on_line ) {
        if ( ( ft == HEL ) && ( m_points.size() > def_min_helix) ) {
            m_helix = helix(in_points, n);
        }
        else if ( (ft == CIRC) && ( m_points.size() > def_min_helix ) ) {
            m_helix = helix(circle(in_points), 0., false, point(0,0,0));
        }
    }*/
    ADDED(n);
}


      
// Private methods ////////////////////////////////////////////////////////////////////////////////////////////////////

void segment::ADDED(const point n)
{
    unsigned s = m_points.size();

    if ( s == 0 ) {
        set_n(n);
        return;
    }
    else if ( s == 1 ) {
        set_n(n);
        return;
    }
    else if ( s == 2 ) {
        m_line = line(m_points[0], (m_points[1] - m_points[0]));
        dir_begin = m_line.get_n();
        dir_end = m_line.get_n();
        set_n(n);
        return;
    }
    else {
        if ( is_linear(ft) ) {
            m_line = line(m_points, ft);
             m_err = m_line.chi2(m_points);
            dir_begin = m_line.get_n();
            dir_end = m_line.get_n();
            set_n(n);
            return;
        }
        else { // For circular and helix fitting
            m_line = line(m_points, def_lft); // We still use default linear fitting
            m_err = m_line.chi2(m_points); // If its error is small enough

            if ( ( m_err < def_on_line ) ) { //||
                 //( eq(m_helix.get_r(), 0) ) ||
                 //( ( ft == HEL ) && ( m_points.size() < def_min_helix ) ) ||
                 //( ( ft == CIRC) && ( m_points.size() < def_min_circ) ) ) {
                dir_begin = m_line.get_n();
                dir_end = m_line.get_n();
                set_n(n);
                return;
            }
            else if ( ft == CIRC ) {
                m_helix = helix(circle(m_points), 0., false, point(0,0,0));
                m_err = m_helix.get_cc().chi2(m_points);
                dir_begin = m_helix.get_cc().tangent(m_points.front(), 2 * m_points.front() - m_points.back());
                dir_end   = m_helix.get_cc().tangent(m_points.back(), m_points.front());
                return;
            }
            else if ( ft == HEL ) {
                if ( !eq(m_helix.get_step(), 0) ) {
                    m_helix = helix(m_points, n);
                    m_err = m_helix.chi2(m_points);
                    dir_begin = m_helix.tangent(m_points.front(), 2 * m_points.front() - m_points.back());
                    dir_end = m_helix.tangent(m_points.back(), m_points.front());
                }
                else {
                    plane check(m_points);
                    if ( collinear(point(check.get_A(), check.get_B(), check.get_C()), n) ) {
                        m_helix = helix(circle(m_points), 0., false, point(0,0,0));
                        m_err = m_helix.get_cc().chi2(m_points);
                        dir_begin = m_helix.get_cc().tangent(m_points.front(), 2 * m_points.front() - m_points.back());
                        dir_end   = m_helix.get_cc().tangent(m_points.back(), m_points.front());
                    }
                    else {
                        m_helix = helix(m_points, n);
                        m_err = m_helix.chi2(m_points);
                        dir_begin = m_helix.tangent(m_points.front(), 2 * m_points.front() - m_points.back());
                        dir_end = m_helix.tangent(m_points.back(), m_points.front());
                    }
                }
                set_n(n);
                return;
            }
        }
        cout << "This never appears :)\n";
    }
}

void segment::set_n(const point n) {
    m_helix = helix(m_helix.get_mid(), n, m_helix.get_r(), m_helix.get_step(), m_helix.get_right(), m_helix.get_on());
}

// Public methods in alphabetical order (NOT printing, NOT get()) ////////////////////////////////////////////////////

void segment::add_point (const point& p)
{
    if ( (m_points.size() > 0) && (m_points.back() == p) ) {
        return;
    }

    m_points.push_back(p);

    ADDED(m_helix.get_n());
}

void segment::add_points(const std::vector <point>& more_points)
{
    if ( more_points.size() == 0 ) {
        return;
    }

    if ( more_points.size() == 1 ) {
        add_point(more_points[0]);
        return;
    }
    
    for ( unsigned i=0; i < more_points.size(); i++ ) {
        if (m_points.back() != more_points[i]) {
            m_points.push_back(more_points[i]);
        }
    }

    ADDED(m_helix.get_n());
}

const std::vector <point> segment::all() const
{
    return m_points;
}

void segment::clear()
{
    point n = m_helix.get_n();
    m_helix = helix(point(0,0,0), n, 0, 0, false, point(0,0,0));
    m_points.clear();
    m_line = line();
    m_err = 0.;
    dir_begin = point();
    dir_end = point();
}

const point segment::get_dir_begin() const
{
    return dir_begin;
}

const point segment::get_dir_end() const
{
    return dir_end;
}

const point segment::get_dir_lin() const
{
    return m_line.get_n();
}

const double segment::chi2() const {
    return m_err;
}

const point segment::first() const
{
    return m_points.front();
}

const point segment::last() const
{
    return m_points.back();
}

const unsigned segment::size() const
{
    return m_points.size();
}

// get() methods /////////////////////////////////////////////////////////////////////////////////////////////////////

const fit_type segment::get_ft() const {
    return ft;
}

const fit_type segment::get_ft_used() const {
    if ( is_linear(ft) ) {
        return ft;
    }
    if ( eq(m_helix.get_r(), 0) ) {
        return def_lft;
    }
    if ( !eq(m_helix.get_step(), 0) ) {
        return HEL;
    }
    return CIRC;
}

const line segment::get_line() const {
    return m_line;
}

const helix segment::get_helix() const {
    return m_helix;
}

const point segment::get_point(unsigned i) const
{
    return m_points[i];
}



// Printing methods //////////////////////////////////////////////////////////////////////////////////////////////////

void segment::print_p(const unsigned pr) const {
    for (int i = 0; i < (int) m_points.size(); i++) {
        m_points[i].p_print(pr);
    }
}

void segment::fprint_p(const char * fname, const unsigned int pr) const {
    for (int i = 0; i < (int) m_points.size(); i++) {
        m_points[i].p_fprint(fname, pr);
    }
}



void segment::ends_print(const unsigned int pr) const {
    if ( m_points.size() == 0 ) {return;}
    else if ( m_points.size() == 1 ) {
        m_points[0].p_print(pr);
        return;
    }
    else if ( ( m_points.size() == 2 ) || ( def_ends_p == 0) ) {
        m_points.front().p_print(pr);
        m_points.back ().p_print(pr);
        return;
    }
    
    if ( is_linear(ft) ) {
        m_line.project(m_points.front()).p_print(pr);
        m_line.project(m_points.back ()).p_print(pr);
    }
    else if ( ft == CIRC ) {
        m_helix.get_cc().project(m_points.front()).p_print(pr);
        m_helix.get_cc().project(m_points.back ()).p_print(pr);
    }
    else if ( ft == HEL ) {
        m_helix.project(m_points.front()).p_print(pr);
        m_helix.project(m_points.back ()).p_print(pr);
    }
    else {
        cout << "You will never see this :) \n";
    }
    return;
}

void segment::ends_fprint(const char * fname, const unsigned int pr) const {
    if ( m_points.size() == 0 ) {return;}
    else if ( m_points.size() == 1 ) {
        m_points[0].p_fprint(fname, pr);
        return;
    }
    else if ( ( m_points.size() == 2 ) || ( def_ends_p == 0) ) {
        m_points.front().p_fprint(fname, pr);
        m_points.back ().p_fprint(fname, pr);
        return;
    }
    
    if ( is_linear(ft) ) {
        m_line.project(m_points.front()).p_fprint(fname, pr);
        m_line.project(m_points.back ()).p_fprint(fname, pr);
    }
    else if ( ft == CIRC ) {
        m_helix.get_cc().project(m_points.front()).p_fprint(fname, pr);
        m_helix.get_cc().project(m_points.back ()).p_fprint(fname, pr);
    }
    else if ( ft == HEL ) {
        m_helix.project(m_points.front()).p_fprint(fname, pr);
        m_helix.project(m_points.back ()).p_fprint(fname, pr);
    }
    else {
        cout << "You will never see this :) \n";
    }
    return;
}


void segment::print(const unsigned n, const unsigned int pr) const {
    if ( m_points.size() == 0 ) {return;}
    else if ( m_points.size() == 1 ) {
        m_points[0].p_print(pr);
        return;
    }
    else if ( ( m_points.size() == 2 ) || ( ( is_linear(ft) ) && ( def_ends_p == 0 ) ) ) {
        m_points.front().p_print(pr);
        m_points.back ().p_print(pr);
        return;
    }

    if ( is_linear(get_ft_used()) ) {
        m_line.project(m_points.front()).p_print(pr);        
        m_line.project(m_points.back()).p_print(pr);
        return;
    }
    else if ( get_ft_used() == CIRC ) {
        point axis, normal = m_helix.get_n(), mid = m_helix.get_mid();
        point xb, xe, p;
        xb = m_helix.get_cc().project(m_points[0]);
        xb.p_print(pr);
        xb -= mid;
        xe = m_helix.get_cc().project(m_points[m_points.size() - 1]) - mid;
        double a = angle(xb, xe) / ((double) (n - 1));

        if ( (xb - xe).norm() < (euler_rot(xb, normal, a) - xe).norm() ) {
            axis = - normal;
        }
        else {
            axis = normal;
        }

        for (unsigned i = 1; i <= (n - 2); i++) {
            p = euler_rot(xb, axis, a * i);
            (p + mid).p_print(pr);
        }
        (xe + mid).p_print(pr); 
        return;
    }
    else if ( get_ft_used() == HEL ) {
        vector <point> pp;
        pp = m_helix.points_to_draw(m_points.front(), m_points.back(), n);
        for (unsigned i = 0; i < n; i++) {
            pp[i].p_print(pr);
        }
    }
    else {
        cout << "You will never see this :) \n";
    }
}

void segment::fprint(const char * fname, const unsigned n, const unsigned int pr) const {
    if ( m_points.size() == 0 ) {return;}
    else if ( m_points.size() == 1 ) {
        m_points[0].p_fprint(fname, pr);
        return;
    }
    else if ( ( m_points.size() == 2 ) || ( ( is_linear(ft) ) && ( def_ends_p == 0 ) ) ) {
        m_points.front().p_fprint(fname, pr);
        m_points.back ().p_fprint(fname, pr);
        return;
    }

     if ( is_linear(get_ft_used()) ) {
        m_line.project(m_points.front()).p_fprint(fname, pr);        
        m_line.project(m_points.back()).p_fprint(fname, pr);
        return;
    }
    else if ( get_ft_used() == CIRC ) {
        point axis, normal = m_helix.get_n(), mid = m_helix.get_mid();
        point xb, xe, p;
        xb = m_helix.get_cc().project(m_points[0]);
        xb.p_fprint(fname, pr);
        xb -= mid;
        xe = m_helix.get_cc().project(m_points[m_points.size() - 1]) - mid;
        double a = angle(xb, xe) / ((double) (n - 1));

        if ( (xb - xe).norm() < (euler_rot(xb, normal, a) - xe).norm() ) {
            axis = - normal;
        }
        else {
            axis = normal;
        }

        for (unsigned i = 1; i <= (n - 2); i++) {
            p = euler_rot(xb, axis, a * i);
            (p + mid).p_fprint(fname, pr);
        }
        (xe + mid).p_fprint(fname, pr); 
        return;
    }
    else if ( get_ft_used() == HEL ) {
        vector <point> pp;
        pp = m_helix.points_to_draw(m_points.front(), m_points.back(), n);
        for (unsigned i = 0; i < n; i++) {
            pp[i].p_fprint(fname, pr);
        }
        /*point fir = m_helix.get_axis().project(m_points[0]);
        point ste = ( m_helix.get_axis().project(m_points[m_points.size() - 1]) - fir ) / ((double) (n - 1));
        
        for (unsigned i = 0; i < n; i++) {
            m_helix.project( fir + ste * ((double) i) ).p_fprint(fname, pr);
            ( m_helix.get_axis().project(m_points[0]) + ste ).p_fprint(fname, pr);
        }
        return;*/
    }
    else {
        cout << "You will never see this :) \n";
    }
}


const point segment::project(const point p) const {
    if ( is_linear(ft) ) {
        return m_line.project(p);
    }
    else if ( ft == CIRC ) {
        return m_helix.get_cc().project(p);
    }
    else if ( ft == HEL ) {
        return m_helix.project(p);
    }
    else {
        return p;
    }
}
