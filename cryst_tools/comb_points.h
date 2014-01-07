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
#include <boost/random.hpp>

class cmb_group
{
public:  
  std::set<int> indexes;
  bool unique_conn;
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
  void get_dist_vc_map(int index_cntr, vc_dist &dst_array, map_dist dst_map);
  void delete_singles(std::vector<int> &data);
protected:
  virtual int get_points_size() const = 0;
  virtual double get_distance(int i, int j) const = 0;
  groups_vc split_group(const cmb_group &cg, double tol);
  bool cluster_combinations(groups_vc &vc, double tol);
  int  get_possible_connections(index_conn &ic, double tol, 
                                int min_cntr_points = 4);
  int  verify_connections(index_conn &ic, double tol);
  void create_groups(const index_conn &ic, groups_vc &grp);
public:
  points_clusters();  
};

#endif	/* COMB_POINTS_H */

