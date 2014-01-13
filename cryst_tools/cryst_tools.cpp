/* 
 * File:   cryst_tools.cpp
 * Author: kirill
 * 
 * Created on December 16, 2013, 4:00 PM
 */

#include <set>
#include <Eigen/Dense>
#include <eigen2/Eigen/src/Core/MatrixBase.h>
#include "cryst_tools.h"
#include "comb_points.h"
#include "openbabel/math/vector3.h"

using namespace std;
using namespace Eigen;
using namespace cryst_tools;

int Matrix_compare(const Matrix3d &m1, const Matrix3d &m2)
{
  int result = 0;
  for(int i = 0; i < 9; i++)  
  {
    double d = m1(i / 3, i % 3) - m2(i / 3, i % 3);
    if( abs(d) > 1E-4)
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
    if( abs(d) > 1E-4)
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

Eigen::Vector3d cryst_tools::min_frac(const Eigen::Vector3d &frac)
{
  Vector3d result;
  
  for(int i = 0; i < 3; i++)
    result[i] = frac[i] - floor(frac[i] + 0.5);
  
  return result;
}

Eigen::Vector3d cryst_tools::norm_frac(const Eigen::Vector3d &frac)
{
  Vector3d result;

  for(int i = 0; i < 3; i++)
    result[i] = frac[i] - floor(frac[i]);
  
  return result;
}


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
    
    if( (A.transpose() - A.inverse()).norm() < 1E-3)
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

class ps_shifts : public points_clusters 
{
protected:
  vector< Vector3d > vc;
  vector< int > base_index;
  Matrix3d cell;
  double tol_cart;
public:
  virtual int get_points_size() const;
  virtual double get_distance(int i, int j) const;
  virtual bool is_connected(int i, int j) const;
public:
  ps_shifts(const Matrix3d &cell_v);
  void ps_add_shift(int index, const Vector3d &shift);
  Vector3d average_vector(const cmb_group &cbg);
  bool total_shift(const cmb_group &cbg, int poins_num);
  void create_groups(groups_vc &vc, double tol_cart_v);    
};

ps_shifts::ps_shifts(const Matrix3d &cell_v)
{
  base_index.resize(0);
  vc.resize(0);
  cell = cell_v;
}

int ps_shifts::get_points_size() const
{
  return vc.size();
}

double ps_shifts::get_distance(int i, int j) const
{
  return min_frac(vc[i] - vc[j]).norm();
}

bool ps_shifts::is_connected(int i, int j) const
{
  return (cell * min_frac(vc[i] - vc[j])).norm() <= tol_cart;
}

void ps_shifts::create_groups(groups_vc &vc, double tol_cart_v)
{
  tol_cart = tol_cart_v;
  
  Matrix3d r_cell = cell.inverse();
  
  Matrix3d mc = r_cell.transpose() * r_cell;

  double sum = 0;
  for(int i = 0; i < 9; i++)  
  {
    sum += abs(mc(i / 3, i % 3));
  }
  tol_list = sqrt(sum) * tol_cart;  
  create_groups_internal(vc, tol_list, 2);
}

void ps_shifts::ps_add_shift(int index, const Vector3d &shift)
{
  vc.push_back(min_frac(shift));
  base_index.push_back(index);
}

Vector3d ps_shifts::average_vector(const cmb_group &cbg)
{
  Vector3d result;
  
  result.setZero();
  
  for(set<int>::const_iterator it  = cbg.indexes.begin();
                               it != cbg.indexes.end(); ++it)
    result += min_frac(vc[*it] - vc[*cbg.indexes.begin()]);

  result = vc[*cbg.indexes.begin()] + result / double(cbg.indexes.size());
  
  return result;
}
bool ps_shifts::total_shift(const cmb_group &cbg, int poins_num)
{
  bool result;
  
  result = cbg.indexes.size() == poins_num;
  
  #ifndef NDEBUG
  if(result)
  {
    vector<int> vi;
    for(set<int>::const_iterator it  = cbg.indexes.begin();
                                 it != cbg.indexes.end(); ++it)
      vi.push_back(base_index[*it]);
    
    sort(vi.begin(), vi.end());
    assert(vi.size() == poins_num);
    for(int i = 0; i < vi.size(); i++)
      assert(vi[i] == i);
  }  
  
  #endif  
  
  return result;
}

std::vector<Eigen::Vector3d> cryst_tools::get_shifts(const Eigen::Matrix3d &sym_matrix_B,
                                                     const Eigen::Matrix3d &cell,
                                                     const std::vector<Eigen::Vector3d> &frac_coords,
                                                     double tolerance, bool all_symmetric)
{
  typedef list<Vector3d> o_list;
  o_list out_list;
  
  ps_shifts ps(cell);
  
  for(int i = 0; i < frac_coords.size(); i++)  
  {
    for(int j = 0; j < frac_coords.size(); j++)
    {
      Vector3d shift;
      shift = sym_matrix_B * frac_coords[i] - frac_coords[j];
      ps.ps_add_shift(i, shift);      
    }  
  }  
  
  groups_vc gvc;
  
  ps.create_groups(gvc, tolerance);
  ps.assign_max_dist(gvc);
  
  out_list.clear();  
  for(int i = 0; i < gvc.size(); i++)
  {  
    assert(gvc[i].max_dist < tolerance);
    assert(gvc[i].unique_conn);
    if(!all_symmetric)
      out_list.push_back(ps.average_vector(gvc[i]));
    else if( ps.total_shift(gvc[i], frac_coords.size()) )
      out_list.push_back(ps.average_vector(gvc[i]));
  }
  
  vector<Vector3d> result;
  result.resize(out_list.size());
  copy(out_list.begin(), out_list.end(), result.begin());
    
  return result;
}

std::vector<Eigen::Vector3d> cryst_tools::shifts_intersect(const Eigen::Matrix3d &cell,  
                                                           const std::vector<Eigen::Vector3d> &frac_c1,
                                                           const std::vector<Eigen::Vector3d> &frac_c2,
                                                           double tolerance)
{
  std::vector<Eigen::Vector3d> result;
  result.clear();
  
  for(int i = 0; i < frac_c1.size(); i++)
  {
    for(int j = 0; j < frac_c2.size(); j++)
    {
      Vector3d dist = min_frac(frac_c2[j] - frac_c1[i]);
      if((cell * dist).norm() <= tolerance)
        result.push_back(frac_c1[i] + 0.5 * dist);
    }
  }  
  
  return result;
}        


std::vector<Eigen::Affine3d> cryst_tools::get_all_symmetries(const Eigen::Matrix3d &cell,
                                                             const vc_sets &frac_coords,
                                                             const double tol,        
                                                             const int range)
{
  vector<Affine3d> result;
  result.clear();
  
  vector<Matrix3d> symm = get_symmetries(cell, range);
  
  for(int i = 0; i < symm.size(); i++)
  {
    vector<Vector3d> shifts_total;
    for(int j = 0; j < frac_coords.size(); j++)
    {  
      vector<Vector3d> shifts = get_shifts(symm[i], cell, frac_coords[j], tol, true);
      
      if(j == 0)
        shifts_total = shifts;
      else
        shifts_total = shifts_intersect(cell, shifts_total, shifts, tol);
    }

    for(int j = 0; j < shifts_total.size(); j++)
    {
      Affine3d mo;
      
      mo.setIdentity();
      
      mo.linear() = symm[i];
      mo.translation() = shifts_total[j];
      
      result.push_back(mo);
    }
  }
  return result;  
}        


  /*
  //remove translations shifts over cell
  for(o_list::iterator it = out_list.begin(); it != out_list.end(); ++it)
  {
    o_list::iterator jt = it;
    ++jt;
    o_list::iterator jt_next;
    for(; jt != out_list.end(); jt = jt_next)
    {
      jt_next = jt;
      ++jt_next;
      double sum= 0;
      for(int i = 0; i < 3; i++)
      {  
        double d = (*jt)[i] - (*it)[i];
        sum += abs(d - round(d));
      }
      if(sum < 1E-4)
        out_list.erase(jt);
    }
  }
  */  
