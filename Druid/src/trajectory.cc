#include "trajectory.h"
#include <cmath>
#include <iostream>
#include <fstream>
#include <climits> 

using namespace std;



// Constructors ///////////////////////////////////////////////////////////////////////////////////////////////////////

trajectory::trajectory(const fit_type ftin, const point nin, const double add_crit, const double crit_chi2)
    : m_segments(),
      add_point_dir_criteria_degree(add_crit),
      chi2_criteria(crit_chi2),
      ft(ftin),
      n(nin)
{}

trajectory::trajectory(const std::vector <point>& in_points, const fit_type ftin, const point nin, const double add_crit, const double crit_chi2)
    : m_segments(),
      add_point_dir_criteria_degree(add_crit),
      chi2_criteria(crit_chi2),
      ft(ftin),
      n(nin)
{
    for (unsigned i = 0; i < in_points.size(); ++i)
    {
        add_point(in_points[i]);
        //cout << fit_type_to_string(m_segments.back().get_ft_used()) << endl;
    }
}

// Private methods ////////////////////////////////////////////////////////////////////////////////////////////////////
// Used inside add_point() method


const bool trajectory::add_1(const point& p) {
    double add_point_dir_criteria_rad = add_point_dir_criteria_degree * m_pi / 180;
    // Is direction (last point -> new point) very different from the direction of the last segment?

    point x; // 'Direction' of the last segment    
             // Now we decide what to take as x:
    point y; // Last point -> new point
    
    x = m_segments.back().get_dir_end();
    y = p -  m_segments.back().project( m_segments.back().last() );
    if (angle(x,y) < add_point_dir_criteria_rad) { // Add to the last segment
        return true;
    }
    else { // Create a new segment
        return false;
    }

}

const bool trajectory::add_2(const point& p) { // Works good only for short branches
    // Does adding a new point increase chi_2 by a critical value?
    segment tmp = m_segments.back();
    tmp.add_point(p);
    double new_err = tmp.chi2();
    if ( new_err  < chi2_criteria ) {
        return true;
    }
    else {
        return false;
    }
}

/*const bool trajectory::add_3(const point& p) { // He think that trajectory cannot bend too much
    if ( m_segments.size() < 2 ) {
        return true;
    }

    segment tmp = m_segments.back();
    tmp.add_point(p);
    point mid = ( tmp.get_helix().get_axis().project(tmp.first()) + tmp.get_helix().get_axis().project(tmp.last()) );
    point x = ( tmp.get_helix().project(tmp.first()) - mid );
    point y = ( tmp.get_helix().project(tmp.last()) - mid );
    double angle_degree_per_hit = angle(x, y) * 180 / m_pi / ( (double) tmp.size() );
    if ( angle_degree_per_hit <  def_max_degree_per_hit ) {
        return true;
    }
    else {
        return false;
    }
}*/   
/*
void trajectory::fix_1() {
    // Here we try to merge the last segment with one of previous, if they have same direction
    
    if (m_segments.size() < 4){ // 4 is the minimum number of segments for performing this procedure
        return;
    }

    double merge_dir_criteria_rad = def_crit_merge * m_pi / 180;
    unsigned i = 1;
    unsigned c;
    point x, y, z;
    const unsigned DEPTH = 2;
    unsigned segs = m_segments.size();

    // Here are the conventions for x, y, z:
    // x: direction of the last segment at the beginning (lin_dir or tangent to circle at starting point, depending on take_circ value)
    // y: direction of the [segs - 3 - i] segment at the end (lin_dir or tangent to circle at ending point, depending on take_circ value)
    // z: vector pointing from the last point of [segs - 3 - i] segment to the first point of the last segment

    if (m_segments.back().size() < 3 || m_segments.back().size() < take_circ) { // Yes, we have to recalculate it - 
        x = m_segments.back().dir_lin();                                        // Because we've just added a point to the segment     
    }
    else {
        x = m_segments.back().dir_circ_end();
    }
    
    while ( (i <= DEPTH) && ((i + 3) <= segs) ) {
        c = segs - 3 - i; // This counter shows to which segment we compare the direction of the last one

        if (m_segments[c].size() < 3 || m_segments[c].size() < take_circ) {
            y = m_segments[c].dir_lin();
        }
        else {
            y = m_segments[c].dir_circ_end();
        }

        z = m_segments.back().first() - m_segments[c].last();

        if ( ( angle(x,y) < merge_dir_criteria_rad ) && ( angle(z,y) < merge_dir_criteria_rad ) ) {
            m_segments[c].add_points(m_segments.back().all()); // Merging the last segment with the one with similar direction
            for (unsigned j = 0; j < (i + 2); j++){ // Simple check: we can remove the minimum number of 3 segments
                m_segments.pop_back();
            }
            return; // We perform this procedure only once per adding a point to the trajectory
        }
        i++;
    }
}*/

void trajectory::fix_2() {

}

// Methods in alphabetical order (NOT printing) ///////////////////////////////////////////////////////////////////////

void trajectory::add_point( const point& p ) // The most important method!
{
    // Part 1 //////////////////////////////////////////////////////////////////////////////////////////
    // First we consider a couple of special cases:
    if (m_segments.size() == 0){
        vector <point> pp;
        pp.push_back(p);
        segment tmp(pp, ft, n);
        m_segments.push_back(tmp);
        return;
    }

    if (m_segments.size() == 1 && m_segments.back().size() == 1){ // OK
        m_segments[0].add_point(p);
        return;
    }
    
    // Part 2 //////////////////////////////////////////////////////////////////////////////////////////

    // Here we use one of possible criteria to decide if we should
    // add the point to the last segment or begin another one.
    bool decision; // The combined version: add_1() for long segments and add_2() for short

        
    // Option 1 /////////////////////////////////////////////////////////////////////////
    // decision = add_1(p); // Choose this to use cone-based segmentation only
    
    // Option 2 /////////////////////////////////////////////////////////////////////////
    //decision = add_2(p); // Choose this to use chi2-based segmentation only

    // Option 3 /////////////////////////////////////////////////////////////////////////
    // An example of more complicated logic (3 options depending of branch length)
    if ( m_segments.back().size() < def_add_1_add_2_decision ) { //
        decision = add_2(p);
    }
    else if ( m_segments.back().size() < ( def_add_1_add_2_decision + 6) ) {
        decision = add_1(p) || add_2(p);
    }
    else {
        decision = add_1(p);
    }

    // Part 3 //////////////////////////////////////////////////////////////////////////////////////////
    
    if ( decision ) { // See different options above
        // Add to the last segment
        m_segments.back().add_point(p);
    }
    else {
        // Create a new one
        vector <point> pp;
        pp.push_back(m_segments.back().last());
        segment tmp(pp, ft, n);
        tmp.add_point(p);
        m_segments.push_back(tmp);
    }

    // Part 4 //////////////////////////////////////////////////////////////////////////////////////////

    // Now we launch procedures which try to find and remove 'bad' points

    // fix_1(); // This procedure looks if direction of two close segments is ~similar and merges them
    // fix_2(); // This procedure 
    // ...
    
    

}

void trajectory::add_points(const std::vector <point>& points) {
    for (unsigned i = 0; i < points.size(); i++) {
        add_point(points[i]);
    }
}

void trajectory::clear()
{
    m_segments.clear();
}

const point trajectory::dir_begin(const unsigned i) const {
    if(m_segments.size() == 0) {
        return point(0,0,0);
    }
    else {
        return m_segments[i].get_dir_begin();
    }
}

const point trajectory::dir_end() const {
    if(m_segments.size() == 0) {
        return point(0,0,0);
    }
    else {
        return m_segments.back().get_dir_end();
    }
}

const point trajectory::dir_end(const unsigned i) const {
    if(m_segments.size() == 0) {
        return point(0,0,0);
    }
    else {
        return m_segments[i].get_dir_end();
    }
}

const point trajectory::dir_lin(const unsigned i) const {
    if(m_segments.size() == 0) {
        return point(0,0,0);
    }
    else {
        return m_segments[i].get_dir_lin();
    }    
}


const double trajectory::chi2() const {
    double err = 0;
    if ( size() == 0 ) {
        return err;
    }
    
    for (unsigned i = 0; i < m_segments.size(); i++) {
        err += m_segments[i].chi2();
    }
    return err;
}

const double trajectory::chi2(int i) const
{
    return m_segments[i].chi2();
}

const point trajectory::get_point(const unsigned i, const unsigned j) const
{
    return m_segments[i].get_point(j);
}

const segment trajectory::get_segment(const unsigned int i) const
{
    return m_segments[i];
}

const unsigned trajectory::n_points() const
{
    unsigned n = 0;
    for (unsigned i = 0; i < m_segments.size(); i++){
        n += m_segments[i].size();
    }
    n -= (m_segments.size() - 1);
    return n;
}

const unsigned trajectory::segment_size( const unsigned i ) const
{
    return m_segments[i].size();
}

const unsigned trajectory::size() const
{
    return m_segments.size();
}

// Printing methods ///////////////////////////////////////////////////////////////////////////////////////////////////

void trajectory::print_p(const unsigned pr) const
{
    for (unsigned i = 0; i < m_segments.size(); i++) {
        m_segments[i].print_p(pr);
        cout << endl;
    }
}

void trajectory::fprint_p(const char * fname, const unsigned pr) const {
    for (unsigned i = 0; i < m_segments.size(); i++) {
        m_segments[i].fprint_p(fname, pr);
        ofstream output;
        output.open(fname, ios::app);
        output << endl;
        output.close();
    }
}

void trajectory::ends_print(const unsigned pr) const {
    for (unsigned i = 0; i < m_segments.size(); i++) {
        m_segments[i].ends_print(pr);
        cout << endl;
    }
}

void trajectory::ends_fprint(const char * fname, const unsigned pr) const {
    for (unsigned i = 0; i < m_segments.size(); i++) {
        m_segments[i].ends_fprint(fname, pr);
        ofstream output;
        output.open(fname, ios::app);
        output << endl;
        output.close();
    }
}

void trajectory::print(const unsigned n, const unsigned pr) const {
    for (unsigned i = 0; i < m_segments.size(); i++) {
        m_segments[i].print(n, pr);
        cout << endl;
    }
}

void trajectory::fprint(const char * fname, const unsigned n, const unsigned pr) const {
    for (unsigned i = 0; i < m_segments.size(); i++) {
        m_segments[i].fprint(fname, n, pr);
        ofstream output;
        output.open(fname, ios::app);
        output << endl;
        output.close();
    }
}
