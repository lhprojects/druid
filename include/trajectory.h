#ifndef _TRAJECTORY_H_
#define _TRAJECTORY_H_

#include <vector>

#include "segment3.h"

using namespace std;

// Settings ///////////////////////////////////////////////////////////////////////////////////////////////////////////

// The default value for the cone segmentation IN DEGREES; add_1()
#define def_crit_add 24.

// The default value for the cone-based merging procedure
#define def_crit_merge .1 // So, normally this feature is not used

// The default value for chi^2 in chi^2-based segmentation
#define def_crit_chi2 15. // See add_2() method

#define def_add_1_add_2_decision 6 // See add_points() method
#define def_max_degree_per_hit 35. // See add_3() method

class trajectory
{
    vector <segment> m_segments;
    double add_point_dir_criteria_degree;
    double chi2_criteria;
    fit_type ft;
    point n;
    
    // Different kinds of logic, see add_point() in trajectory.cpp
    const bool add_1(const point& p); // Return 1 if point should be added to the last segment,
    const bool add_2(const point& p); //        0 if a new one should be created
    const bool add_3(const point& p); // Cannot be used as a separate method!!

    void fix_1();
    void fix_2();

public:
    // See default values for constructors in trajectory. /////////////////////////////////////////////////////////////
    trajectory(const fit_type ftin = def_fit_type, const point nin = def_helix_n, double add_crit = def_crit_add, const double crit_chi2 = def_crit_chi2);
    trajectory(const vector <point>& in_points, const fit_type ftin = def_fit_type, const point nin = def_helix_n, const double add_crit = def_crit_add, const double crit_chi2 = def_crit_chi2);

    // Methods in alphabetical order (NOT printing) ///////////////////////////////////////////////////////////////////

    void add_point(const point& p);
    void add_points(const vector <point>& points);
    void clear(); // !!! add_point_dir_criteria_degree, merge_dir_criteria_degree, take_circ and lft stay unchanged !!!

    const point dir_begin(const unsigned i = 0) const; // Both are along
    const point dir_end(const unsigned i) const;   // segment's direction
    const point dir_lin(const unsigned i) const;
    const point dir_end() const;

    const double chi2() const;
    const double chi2(int i) const;

    const point get_point(const unsigned i, const unsigned j) const;
    const segment get_segment(const unsigned int i) const;
    const unsigned n_points() const; // Here each point is taken into account once (endpoints belong to both segments!)
    const unsigned segment_size(const unsigned i) const;
    const unsigned size() const;

    // Printing methods ///////////////////////////////////////////////////////////////////////////////////////////////

    // Prints all points; segments divided by a linebreak
    void print_p(const unsigned pr = std::cout.precision()) const; // Here one can set the precision of the output (17 max)
    void fprint_p(const char * fname, const unsigned pr = std::cout.precision()) const;

    void ends_print(const unsigned pr = std::cout.precision()) const; // See description in segment2.h
    void ends_fprint(const char * fname, const unsigned pr = std::cout.precision()) const;

    // Next 4 functions print first and last points (or their projections onto the fitted line - see segment2.h) for segments fitted by linear fit
    // and n points projected on the circle arc for segments fitted by the circular fit
    void print(const unsigned n, const unsigned pr = std::cout.precision()) const; // Here one can set the precision of the output (17 max)
    void fprint(const char * fname, const unsigned n, const unsigned pr = std::cout.precision()) const;
};

#endif


