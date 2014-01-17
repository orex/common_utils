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
typedef std::multimap<double, int> map_dist;

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
  double tol_list;
  void get_dist_vc_map(int index_cntr, vc_dist &dst_array, map_dist &dst_map);
  void delete_singles(std::vector<int> &data);
protected:
  int  get_possible_connections(index_conn &ic, double tol, 
                                int min_cntr_points = 4);
  int  verify_connections(index_conn &ic);
  void assign_groups(const index_conn &ic, groups_vc &grp);
  void create_groups_internal(groups_vc &vc, double tol_list_v, int min_cntr_points);  
public:
  int possible_connections;
  int total_connection;  
  virtual int get_points_size() const = 0;
  virtual double get_distance(int i, int j) const = 0;
  virtual bool is_connected(int i, int j) const
  {return get_distance(i, j) <= tol_list; };
public:
  points_clusters();
  virtual void assign_max_dist(groups_vc &gc);
  virtual double min_dist_between_groups(groups_vc &gc); 
};

#endif	/* COMB_POINTS_H */

