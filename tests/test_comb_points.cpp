#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE comb_points_test

#include <boost/test/unit_test.hpp>

#include <iostream>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <complex>
#include <vector>

#include "cryst_tools/comb_points.h"

using namespace std;

class points_ND : public points_clusters 
{
protected:
  vector< vector<double> > vc;
public:
  virtual int get_points_size() const;
  virtual double get_distance(int i, int j) const;
public:
  points_ND(int num_of_points, int dimension);
  void create_groups(groups_vc &vc, double tol_list_v, int min_cntr_points)
  { create_groups_internal(vc, tol_list_v, min_cntr_points); };        
};

points_ND::points_ND(int num_of_points, int dimension)
{
  vc.resize(num_of_points);
  for(int i = 0; i < vc.size(); i++)
  {
    vc[i].resize(dimension);
    for(int j = 0; j < dimension; j++)
    {
      vc[i][j] = double(rand()) / double(RAND_MAX);
    }  
  }  
}

int points_ND::get_points_size() const
{
  return vc.size();
}

double points_ND::get_distance(int i, int j) const
{
  double result = 0;
  
  assert(vc[i].size() == vc[j].size());
  
  for(int k = 0; k < vc[i].size(); k++)
    result += (vc[i][k] - vc[j][k]) * (vc[i][k] - vc[j][k]);
  
  result = sqrt(result);
  
  return result;
}

bool groups_equivalent(const groups_vc &g1, const groups_vc &g2)
{
  bool result;
  
  result = g1.size() == g2.size();
  
  if(!result)
    return result;

  for(int i = 0; i < g1.size(); i++)  
  {
    result = g1[i].unique_conn == g2[i].unique_conn;
    if(!result)
      break;
    
    result = g1[i].indexes == g2[i].indexes;
    if(!result)
      break;
    
    result = g1[i].max_dist == g2[i].max_dist;

    if(!result)
      break;
  }
  
  return result;  
}

bool groups_intersect(const points_ND &pt, const cmb_group &g1, const cmb_group &g2, 
                      const double tol)
{
  bool result = false;
  
  for(set<int>::const_iterator it  = g1.indexes.begin(); 
                               it != g1.indexes.end(); ++it)
  {
    for(set<int>::const_iterator jt  = g2.indexes.begin(); 
                                 jt != g2.indexes.end(); ++jt)
    {
      result = pt.get_distance(*it, *jt) < tol;
      if(result)
      {  
        cout << "intersect on " << *it << " - " << *jt;
        break;
      }  
    }
    if(result)
      break;
  }
  return result;
}

bool group_connected(const points_ND &pt, cmb_group g, double tol)
{
  vector<int> cn_ind;
  
  cn_ind.push_back(*g.indexes.begin());
  g.indexes.erase(g.indexes.begin());
  
  bool changed;
  do 
  {  
    changed = false;
    for(int i = cn_ind.size() - 1; i >= 0; i--)
    {
      for(set<int>::iterator it = g.indexes.begin(); it != g.indexes.end(); ++it)
      {
        if(pt.get_distance(*it, cn_ind[i]) < tol)
        {
          cn_ind.push_back(*it);
          g.indexes.erase(it);
          changed = true;
          break;
        }  
      }
      if(changed)
        break;
    }  
  } while( changed );

  return g.indexes.empty();  
}

bool group_uniform(const points_ND &pt, const cmb_group &g, double tol, double &m_dist)
{
  bool result = true;
  m_dist = 0;
  
  for(set<int>::const_iterator it  = g.indexes.begin(); 
                               it != g.indexes.end(); ++it)
  {
    for(set<int>::const_iterator jt  = g.indexes.begin(); 
                                 jt != g.indexes.end(); ++jt)
    {
      double dist = pt.get_distance(*it, *jt);
      m_dist = max(m_dist, dist);
      result = result && (dist < tol);
    }
  }  
  return result;
}


BOOST_AUTO_TEST_SUITE(CombPointsTest)

BOOST_AUTO_TEST_CASE(Test_random_points_comb_3D)
{
  const double tol = 0.1;
  time_t time_rand = time(NULL);
  srand(time_rand);
  //srand(1389185205);
  cout << "Random initialization. Remember if test fails. " << endl;
  cout << "  " << time_rand << endl;
  
  //Generate random points
  const int points[] = {0, 1, 10, 100, 500, 1000, 1500, 2000, 3000, 5000};
  //const int points[] = {100};
  const int points_size = sizeof(points)/sizeof(points[0]);
  
  for(int i = 0; i < points_size; i++)  
  {
    points_ND pt(points[i], 3);
    groups_vc gvc_ref;
    for(int j = 1; j < 7; j++)
    {
      groups_vc gvc;
      pt.create_groups(gvc, tol, j);
      pt.assign_max_dist(gvc);
      cout << "Num reference points: " << j << endl;
      cout << "Raw num of connection: " << pt.possible_connections << endl;
      cout << "Final connections: " << pt.total_connection << endl;
      cout << "Number of groups: " << gvc.size() << endl;
      
      //Compare with previous groups
      if(j > 1)
        BOOST_CHECK(groups_equivalent(gvc, gvc_ref));
      
      gvc_ref = gvc;
    }
    
    BOOST_CHECK_LE(gvc_ref.size(), points[i]);
    
    //Check all indexes exists and exists ones
    multiset<int> ms;
    for(int j = 0; j < gvc_ref.size(); j++)
    {  
      ms.insert(gvc_ref[j].indexes.begin(), gvc_ref[j].indexes.end());
    }
    BOOST_CHECK_EQUAL(ms.size(), points[i]);
    
    for(int j = 0; j < points[i]; j++)
      BOOST_CHECK(ms.count(j) == 1);
    
    //Check groups are not intersect
    for(int j = 0; j < gvc_ref.size(); j++)    
    {
      for(int k = j + 1; k < gvc_ref.size(); k++)
        BOOST_CHECK_EQUAL(groups_intersect(pt, gvc_ref[j], gvc_ref[k], tol), false);
    }  
    //Check groups are connected
    for(int j = 0; j < gvc_ref.size(); j++)    
      BOOST_CHECK(group_connected(pt, gvc_ref[j], tol));
    
    double m_dist;
    //Check uniforming of groups
    for(int j = 0; j < gvc_ref.size(); j++)    
    {  
      BOOST_CHECK_EQUAL(group_uniform(pt, gvc_ref[j], tol, m_dist), gvc_ref[j].unique_conn);
      //Check distances
      BOOST_CHECK_EQUAL(gvc_ref[j].max_dist, m_dist);
    }  
  }  
}

BOOST_AUTO_TEST_SUITE_END()
