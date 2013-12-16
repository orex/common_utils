/* 
 * File:   obabel_utils.cpp
 * Author: kirill
 * 
 * Created on October 30, 2013, 1:45 PM
 */

#include "obabel_utils.h"

#include <cassert>

OpenBabel::vector3 OpenBabel::get_minimal_distance(const OpenBabel::vector3 &dist, 
                                                   OpenBabel::OBUnitCell * unitcell)
{
  std::vector<OpenBabel::vector3> vr;
  
  vr = get_image_distances(dist, unitcell);
  
  return vr[0];
}


bool vector3_compare(OpenBabel::vector3 v1, OpenBabel::vector3 v2) 
{ 
  return v1.length_2() < v2.length_2(); 
}


std::vector<OpenBabel::vector3> OpenBabel::get_image_distances(const OpenBabel::vector3 &dist, 
                                                              OpenBabel::OBUnitCell * unitcell)
{
  std::vector<OpenBabel::vector3> result;

  matrix3x3 cell = unitcell->GetCellMatrix().transpose();
  matrix3x3 cf   = cell.inverse();
  
  vector3 frac = cf * dist;
  
  result.clear();  
  for(int i = -4; i < 5; i++)
  {
    for(int j = -4; j < 5; j++)
    {
      for(int k = -4; k < 5; k++)
      {
        vector3 frac_test, dist_test;
        frac_test.Set(frac.GetX() + i, frac.GetY() + j, frac.GetZ() + k);
        dist_test = cell * frac_test;
        result.push_back(dist_test);
      }  
    }  
  }  
  
  std::sort(result.begin(), result.end(), vector3_compare);
  
  return result;
}


OpenBabel::vector3 OpenBabel::center_mass(const std::vector<OpenBabel::vector3> &atoms_pos,
                                          OpenBabel::OBUnitCell * unitcell,
                                          const double tol)
{
  assert(atoms_pos.size() > 0);
  
  vector3 central_pos = atoms_pos[0];
  
  vector3 dist_central(0, 0, 0);  
  for(int i = 0; i < atoms_pos.size(); i++)
  {
    vector3 dist = atoms_pos[i] - central_pos;
    dist = get_minimal_distance(dist, unitcell);
    assert(dist.length() < tol);
    dist_central += dist;
  }
  
  dist_central /= double(atoms_pos.size());
  
  vector3 result = central_pos + dist_central;
  
  unitcell->WrapCartesianCoordinate(result);
  
  assert(get_minimal_distance(result - central_pos, unitcell).length() < dist_central.length() + 0.001);
  
  return result;
}
