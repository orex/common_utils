/* 
 * File:   obabel_utils.cpp
 * Author: kirill
 * 
 * Created on October 30, 2013, 1:45 PM
 */

#include "obabel_utils.h"

#include <cassert>

OpenBabel::vector3 OpenBabel::get_minimal_distance(OpenBabel::vector3 dist, 
                                                   OpenBabel::OBUnitCell * unitcell)
{
  vector3 result;

  matrix3x3 cell = unitcell->GetCellMatrix().transpose();
  matrix3x3 cf   = cell.inverse();
  
  vector3 frac = cf * dist;
  
  //cout << cell << endl;
  //cout << cf << endl;          
  
  //cout << dist.GetX() << " " << dist.GetY() << " " << dist.GetZ() << endl;
  //cout << frac.GetX() << " " << frac.GetY() << " " << frac.GetZ() << endl;
  result = dist;
  double res_length2 = dist.length_2();
  for(int i = -2; i < 3; i++)
  {
    for(int j = -2; j < 3; j++)
    {
      for(int k = -2; k < 3; k++)
      {
        vector3 frac_test, dist_test;
        frac_test.Set(frac.GetX() + i, frac.GetY() + j, frac.GetZ() + k);
        dist_test = cell * frac_test;
        double cur_dist2 = dist_test.length_2();
        if( cur_dist2 < res_length2)
        {  
          result = dist_test;
          res_length2 = cur_dist2;
        }  
      }  
    }  
  }  
  
  //cout << result.GetX() << " " << result.GetY() << " " << result.GetZ() << endl;
  
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
