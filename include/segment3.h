#ifndef _SEGMENT3_H_
#define _SEGMENT3_H_

#include <vector>
#include <iostream>

#include "geometry.h"

using namespace std;



// Settings ///////////////////////////////////////////////////////////////////////////////////////////////////////////

//#define def_ends_p 0 // Choose this one for printing the first and the last points of the segments in the 'ends_print' functions
#define def_ends_p 1 // Choose this one for printing their projections on the fitted line

#define def_on_line 1e-5 // If the error of linear fitting is smaller than this value, we do not perform circular fitting!

#define def_min_circ 3 // The minimum number of points to perform circular fitting

#define def_min_helix 3 // The minimum number of points to perform helix fitting

#define def_fit_type CIRC // If not set in the constructor, see options in "fitting.h". Type here CIRC, HELIX or one of linear fitting types

#define def_helix_n point(0,0,1) // The default direction for the helix axis (if not set in the constructor)

class segment
{
    vector <point> m_points;
    line m_line;
    helix m_helix;
    fit_type ft; // Fitting type we want to use (there are some limitations)
    point dir_begin, dir_end;
    double m_err;
    
    // Private methods ////////////////////////////////////////////////////////////////////////////////////////////////

    void ADDED(const point n); // This procedure runs each time after adding new points and inside constructors
    void set_n(const point n);
public:
    segment();
    segment(const vector <point>& in_points, const fit_type ftin = def_fit_type, const point n = def_helix_n);

    // Public methods in alphabetical order (NOT printing; NOT get()) ////////////////////////////////////////////////

    void add_point(const point& p);
    void add_points(const vector <point>& more_points);

    const vector <point> all() const; // Returns vector of all points
    
    void clear(); // lft and m_helix.hc.cc.n (normal) remain unchanged!!!

    const point get_dir_begin() const;
    const point get_dir_end() const;
    const point get_dir_lin() const;
    const double chi2() const;

    const bool helix_fit_used() const;

    const point first() const;
    
    const point get_centroid() const;
    const fit_type get_ft() const; // Don't forget to convert it to string using fit_type_to_string(const fit_type)
    const fit_type get_ft_used() const;
    
    const point last() const;
    const point project(const point p) const;
       
    const unsigned size() const;

    // get() methods /////////////////////////////////////////////////////////////////////////////////////////////////
    const line get_line() const;
    const helix get_helix() const;
    const point get_point(unsigned i) const;

    // Printing methods ///////////////////////////////////////////////////////////////////////////////////////////////

    void print_p(const unsigned pr = std::cout.precision()) const; // Print all points of the segment
    void fprint_p(const char * fname, const unsigned pr = std::cout.precision()) const; // Adds to the end of the file

    void ends_print (const unsigned pr = std::cout.precision()) const;
    void ends_fprint (const char * fname, const unsigned pr = std::cout.precision()) const;

    void print (const unsigned n, const unsigned pr = std::cout.precision()) const;
    void fprint (const char * fname, const unsigned n, const unsigned pr = std::cout.precision()) const;
};

#endif


