/* 
 * File:   comb_points.cpp
 * Author: kirill
 * 
 * Created on December 16, 2013, 9:27 PM
 */

#include "comb_points.h"
#include "others/rnd_utils.h"
#include <boost/range/algorithm/random_shuffle.hpp>
#include <map>

#include <cassert>
#include <algorithm>
#include <ctime>

using namespace std;

points_clusters::points_clusters()
{
  possible_connections = 0;
  total_connection = 0;  
  rnd_gen = create_rnd_gen();
}

void points_clusters::delete_singles(std::vector<int> &data)
{
  sort(data.begin(), data.end());
  int pos_cmb = 0;
  for(int i = 0; i < int(data.size()) - 1; i++)
  {
    if( data[i] == data[i + 1])
    {  
      data[pos_cmb] = data[i];
      pos_cmb++;
    }  
  }
  data.resize(pos_cmb);
}

void points_clusters::get_dist_vc_map(int index_cntr, vc_dist &dst_array, map_dist &dst_map)
{
  dst_map.clear();
  dst_array.resize(get_points_size(), 0.0);
  for(int i = 0; i < get_points_size(); i++)
  {
    double dist = get_distance(index_cntr, i);
    dst_array[i] = dist;
    dst_map[dist] = i;
  }
}

int points_clusters::get_possible_connections(index_conn &ic, double tol,
                                               int min_cntr_points)
{
  vector<cmb_dist> dist_cntr;
  
  dist_cntr.resize(min_cntr_points);
  
  for(int i = 0; i < dist_cntr.size(); i++)
  {
    dist_cntr[i].index_cntr = 
        get_rnd_int_value_in_interval(rnd_gen, 0, get_points_size() - 1);
    get_dist_vc_map(dist_cntr[i].index_cntr, dist_cntr[i].dst_array, dist_cntr[i].dst_map);
  }
  
  int conn_num = 0;
  
  ic.resize(get_points_size());
  for(int i = 0; i < ic.size(); i++)
  {
    std::vector<int> ms;
    ms.clear();
    for(int j = 0; j < dist_cntr.size(); j++)
    {
      double dst = dist_cntr[j].dst_array[i];
      map_dist::const_iterator min_it = dist_cntr[j].dst_map.lower_bound(dst - tol);
      map_dist::const_iterator max_it = dist_cntr[j].dst_map.upper_bound(dst + tol);
      for(map_dist::const_iterator it = min_it; it != max_it; ++it)
      {  
        if( it->second != i )
          ms.push_back(it->second);
      }  
      
      if(j != 0)
        delete_singles(ms);
    }
    ic[i].clear();
    ic[i].insert(ms.begin(), ms.end());
    conn_num += ms.size();
  }
  return conn_num;
}

int points_clusters::verify_connections(index_conn &ic, double tol)
{
  int result = 0;
  
  for(int i = 0; i < ic.size(); i++)
  {
    set<int>::iterator it_next;
    for(set<int>::iterator it = ic[i].begin(); it != ic[i].end(); it = it_next)
    {
      it_next = it;
      ++it_next;
      if(get_distance(i, *it) > tol)
        ic[i].erase(it);
      else
      {  
        //if( (i == 4))
          //cout << "4 connected to " << *it << endl;
        result++;
      }  
      
      
    }  
  }
  
  return result;
}

void points_clusters::assign_groups(const index_conn &ic, groups_vc &grp)
{
  vector<int> group_num;
  map<int, int> group_st; 
  
  group_num.resize(ic.size());
  for(int i = 0; i < group_num.size(); i++)
    group_num[i] = i;
  
  //Set groups num by lowers connection  
  bool changed;
  do
  {  
    changed = false;
    for(int i = 0; i < ic.size(); i++)
    {
      // // Avoid connections down and self connections
      //set<int>::const_iterator b_it = upper_bound(ic[i].begin(), ic[i].end(), i);
      int min_grp_num = group_num[i];
      int max_grp_num = group_num[i];
      for(set<int>::const_iterator it = ic[i].begin(); it != ic[i].end(); ++it)
      {  
        min_grp_num = min(min_grp_num, group_num[*it]);
        max_grp_num = max(max_grp_num, group_num[*it]);
      }
      assert(min_grp_num <= max_grp_num);
      if( min_grp_num <  max_grp_num ) 
      {  
        for(set<int>::const_iterator it = ic[i].begin(); it != ic[i].end(); ++it)
          group_num[*it] = min_grp_num;
        changed = true;
      }  
    }
  }
  while( changed );
  
  //set groups num 0..1..grp_num - 1
  int grp_num = 0;
  for(int i = 0; i < group_num.size(); i++)
  {
    if(group_st.count(group_num[i]) == 0)
    {
      group_st[group_num[i]] = grp_num;
      grp_num++;
    }  
  }
  grp.resize(grp_num);
  
  for(int i = 0; i < group_num.size(); i++)
  {
    int index_group = group_st[group_num[i]];
    grp[index_group].indexes.insert(i);
  }

  for(int i = 0; i < grp.size(); i++)
  {
    int min_con = ic[*grp[i].indexes.begin()].size();
    int max_con = ic[*grp[i].indexes.begin()].size();    
    
    for(set<int>::const_iterator it  = grp[i].indexes.begin(); 
                                 it != grp[i].indexes.end(); ++it)
    {  
      min_con = min(min_con, int(ic[*it].size()) );
      max_con = max(max_con, int(ic[*it].size()) );
    }
    
    assert(min_con <= max_con);
    assert(max_con < grp[i].indexes.size());
    
    grp[i].unique_conn = grp[i].indexes.size() == min_con + 1;
  }
}

void points_clusters::create_groups(groups_vc &vc, double tol, int min_cntr_points)
{
  index_conn ic;
  possible_connections = get_possible_connections(ic, tol, min_cntr_points);
  total_connection = verify_connections(ic, tol);
  assign_groups(ic, vc);
}

void points_clusters::assign_max_dist(groups_vc &vc)
{
  for(int i = 0; i < vc.size(); i++)
  {
    double md = 0;
    for(set<int>::const_iterator it  = vc[i].indexes.begin(); 
                                 it != vc[i].indexes.end(); ++it)
    {
      for(set<int>::const_iterator jt  = vc[i].indexes.begin(); 
                                   jt != vc[i].indexes.end(); ++jt)
        md = max(md, get_distance(*it, *jt));
    }  
    vc[i].max_dist = md;
  }
}

