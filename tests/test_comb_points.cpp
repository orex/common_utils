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
  virtual int get_points_size() const = 0;
  virtual double get_distance(int i, int j) const = 0;
public:
  points_ND(int num_of_points, int dimension);
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
  double result;
  
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
  }
  
  return result;  
}

bool groups_intersect(const points_ND &pt, const cmb_group &g1, const cmb_group &g2, 
                      const double tol)
{
  bool result = false;
  
  for(int i = 0; i < g1.indexes.size(); i++)
  {
    for(int j = 0; j < g2.indexes.size(); j++)
    {
      result = pt.get_distance(g1[i], g2[i]) < tol;
      if(result)
        break;
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
  
  bool changed = true;
  while( changed )
  {  
    for(int i = cn_ind.size() - 1; i >= 0; i--)
    {
      for(set<int>::iterator it = g.indexes.begin(); it != g.indexes.end(); it++)
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
  }

  return g.indexes.empty();  
}


BOOST_AUTO_TEST_SUITE(CombPointsTest)

BOOST_AUTO_TEST_CASE(Test_random_points_comb_3D)
{
  const double tol = 0.1;
  time_t time_rand = time(NULL);
  srand(time_rand);
  cout << "Random initialization. Remember if test fails. " << endl;
  cout << "  " << time_rand << endl;
  
  //Generate random points
  const int points[] = {1, 10, 100, 500, 1000, 1500, 2000, 3000, 5000, 10000};
  const int points_size = sizeof(points)/sizeof(points[0]);
  
  for(int i = 0; i < points_size; i++)  
  {
    points_ND pt(points[i], 3);
    groups_vc gvc_ref;
    for(int j = 1; j < 7; j++)
    {
      groups_vc gvc;
      pt.create_groups(gvc, tol, j);
      cout << "Raw num of connection: " << pt.possible_connections << endl;
      cout << "Final connections: " << pt.total_connection << endl;
      cout << "Number of groups: " << gvc.size() << endl;
      
      //Compare with previous groups
      if(j > 1)
        BOOST_CHECK(groups_equivalent(gvc, gvc_ref));
      
      gvc_ref = gvc;
    }
    BOOST_CHECK(gvc_ref.size() <= points_size);
    
    //Check all indexes exists and exists ones
    multiset<int> ms;
    for(int j = 0; j < gvc_ref.size(); j++)
      ms.insert(gvc_ref[j].indexes.begin(), gvc_ref[j].indexes.end());
    
    BOOST_CHECK(ms.size() == points_size);
    
    for(int j = 0; j < points_size; j++)
      BOOST_CHECK(ms.count(j) == 1);
    
    //Check groups are not intersect
    for(int j = 0; j < gvc_ref.size(); j++)    
    {
      for(int k = j + 1; k < gvc_ref.size(); k++)
        BOOST_CHECK(!groups_intersect(pt, gvc_ref[j], gvc_ref[k], tol));
    }  
    //Check groups are connected
    for(int j = 0; j < gvc_ref.size(); j++)    
      BOOST_CHECK(group_connected(pt, gvc_ref[j], tol));
    
  }  
}

BOOST_AUTO_TEST_SUITE_END()
