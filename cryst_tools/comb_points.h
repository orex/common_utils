/* 
 * File:   comb_points.h
 * Author: kirill
 *
 * Created on December 16, 2013, 9:27 PM
 */

#ifndef COMB_POINTS_H
#define	COMB_POINTS_H

#include <vector>
#include <list>
#include <set>
#include <map>
#include <boost/random.hpp>

class cmb_group
{
public:  
  std::set<int> indexes;
  bool unique_conn;
  double max_dist;
};

typedef std::vector<cmb_group> groups_vc;

typedef std::vector<std::set<int> > index_conn;

typedef std::vector<double> vc_dist;
typedef std::map<double, int> map_dist;

struct cmb_dist
{
  int index_cntr;
  vc_dist   dst_array;
  map_dist  dst_map;
};


class points_clusters
{
protected:
  boost::mt19937 rnd_gen;
  void get_dist_vc_map(int index_cntr, vc_dist &dst_array, map_dist &dst_map);
  void delete_singles(std::vector<int> &data);
protected:
  int  get_possible_connections(index_conn &ic, double tol, 
                                int min_cntr_points = 4);
  int  verify_connections(index_conn &ic, double tol);
  void assign_groups(const index_conn &ic, groups_vc &grp);
public:
  int possible_connections;
  int total_connection;  
  virtual int get_points_size() const = 0;
  virtual double get_distance(int i, int j) const = 0;
public:
  points_clusters();
  void create_groups(groups_vc &vc, double tol, int min_cntr_points);
  void assign_max_dist(groups_vc &vc);
};

#endif	/* COMB_POINTS_H */

