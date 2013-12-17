/* 
 * File:   cryst_tools.cpp
 * Author: kirill
 * 
 * Created on December 16, 2013, 4:00 PM
 */

#include "cryst_tools.h"
#include <set>
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

int Matrix_compare(const Matrix3d &m1, const Matrix3d &m2)
{
  int result = 0;
  for(int i = 0; i < 9; i++)  
  {
    double d = m1(i / 3, i % 3) - m2(i / 3, i % 3);
    if( abs(d) > 1E-6)
    {
      result = int(1.1 * d / abs(d));
      break;
    }  
  }
  
  return result;
}

int Vector_compare(const Vector3d &v1, const Vector3d &v2)
{
  int result = 0;
  for(int i = 0; i < 3; i++)  
  {
    double d = v1(i) - v2(i);
    if( abs(d) > 1E-6)
    {
      result = int(1.1 * d / abs(d));
      break;
    }  
  }
  
  return result;
}

struct Matrix3d_comp_cl 
{
  bool operator() (const Matrix3d &m1, const Matrix3d &m2) const
  {  return Matrix_compare(m1, m2) < 0;}
};

struct Vector3d_comp_cl 
{
  bool operator() (const Vector3d &v1, const Vector3d &v2) const
  {  return Vector_compare(v1, v2) < 0;}
};

std::vector<Eigen::Matrix3d> cryst_tools::get_symmetries(const Eigen::Matrix3d &cell,
                                                           const int range)
{
  std::vector<Eigen::Matrix3d> result;
  set<Eigen::Matrix3d, Matrix3d_comp_cl> real_space_symm;
  
  int mult = 2 * range + 1;
  
  Matrix3d r_cell = cell.inverse();
  
  for(int k = 0; k < pow(mult, 9); k++)
  {
    Matrix3d B;
    int p = k;
    for(int i = 0; i < 9; i++)  
    {
      B(i / 3, i % 3) = p % mult - range;
      p /= mult;
    }
    if( abs(B.determinant()) != 1) continue;
    
    
    Matrix3d A;
    A = cell * B * r_cell;

    if( (A.transpose() - A.inverse()).norm() < 1E-5)
    {  
      if(real_space_symm.count(A) == 0)
      {  
        real_space_symm.insert(A);
        result.push_back(B);
      }  
    }
  }
  
  return result;
}

std::vector<Eigen::Vector3d> cryst_tools::get_shifts(const Eigen::Matrix3d &sym_matrix_B,
                                                      const Eigen::Matrix3d &cell,
                                                      const std::vector<Eigen::Vector3d> &frac_coords,
                                                      double tolerance, bool all_symmetric)
{
  std::vector<Eigen::Vector3d> result;  
  //frac_coords[j] == sym_matrix_B * frac_coords[i] + shift
  
  for  
  
   
  return result;
}



