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
  virtual void assign_max_dist(groups_vc &gc);  
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

void ps_shifts::assign_max_dist(groups_vc &gc)
{
  for(int i = 0; i < gc.size(); i++)
  {
    double md = 0;
    for(set<int>::const_iterator it  = gc[i].indexes.begin(); 
                                 it != gc[i].indexes.end(); ++it)
    {
      for(set<int>::const_iterator jt  = gc[i].indexes.begin(); 
                                   jt != gc[i].indexes.end(); ++jt)
        md = max(md, (cell * min_frac(vc[*it] - vc[*jt])).norm());
    }  
    gc[i].max_dist = md;
  }  
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
      shift = frac_coords[j] - sym_matrix_B * frac_coords[i];
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
                                                             const std::vector<bool> &all_symm, 
                                                             const double tol,
                                                             const int range)
{
  assert(all_symm.size() == frac_coords.size());
  
  vector<Affine3d> result;
  result.clear();
  
  vector<Matrix3d> symm = get_symmetries(cell, range);
  
  for(int i = 0; i < symm.size(); i++)
  {
    vector<Vector3d> shifts_total;
    for(int j = 0; j < frac_coords.size(); j++)
    {  
      vector<Vector3d> shifts = get_shifts(symm[i], cell, frac_coords[j], tol, all_symm[j]);

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



std::vector<Eigen::Affine3d> cryst_tools::get_all_symmetries(const Eigen::Matrix3d &cell,
                                                             const vc_sets &frac_coords,
                                                             const double tol,        
                                                             const int range)
{
  vector<bool> bc;
  bc.resize(frac_coords.size(), true);
  
  return get_all_symmetries(cell, frac_coords, bc, tol, range);
}        

void cryst_tools::min_dist::set_cell(const Eigen::Matrix3d &cell_v)
{
  box_size.setZero();
  
  cell = cell_v;
  r_cell = cell.inverse();
  
  Matrix3d Q = cell.inverse().transpose() * cell.inverse();
  
  r2_dist_direct = 1E12;
  
  for(int i = -4; i < 5; i++)
  {
    for(int j = -4; j < 5; j++)
    {
      for(int k = -4; k < 5; k++)
      {
        Vector3i vi_ind(i, j, k);
        Vector3d v_ind(i, j, k);
        Vector3d Qv = Q * v_ind;
        double d = v_ind.transpose() * Qv - (abs(Qv[0]) + abs(Qv[1]) + abs(Qv[2]));
        
        if(vi_ind.squaredNorm() != 0)
          r2_dist_direct = min(r2_dist_direct, (cell * v_ind).squaredNorm());
        
        if( d < -tol )
        {  
          for(int l = 0; l < 3; l++)
            box_size[l] = max(box_size[l], abs(vi_ind[l]));
        }  
      }  
    }  
  }
  
  cell_assigned = true;  
}

Eigen::Vector3d cryst_tools::min_dist::operator ()(const Eigen::Vector3d &vd) const
{
  assert(cell_assigned);
  
  Eigen::Vector3d result = cell * min_frac(r_cell * vd);
  
  double sq_norm_min = vd.squaredNorm();
  
  if(sq_norm_min <= r2_dist_direct )
    return result;
  
  for(int i = -box_size[0]; i < box_size[0] + 1; i++)
  {
    for(int j = -box_size[1]; j < box_size[1] + 1; j++)
    {
      for(int k = -box_size[2]; k < box_size[2] + 1; k++)
      {
        if(abs(i) + abs(j) + abs(k) == 0 )
          continue;
        
        Vector3d vc = vd - cell * Vector3d(i, j, k);
        double d = vc.squaredNorm();
        if(d < sq_norm_min )
        {
          sq_norm_min = d;
          result = vc;
        }  
      }
    }
  }  
  
  return result;
}

std::vector<Eigen::Vector3d> cryst_tools::min_dist::get_img_dist(const Eigen::Vector3d &vd, 
                                                                 const Eigen::Vector3i &box) const
{
  assert(cell_assigned);

  std::vector<Eigen::Vector3d> result;
  
  for(int i = -box[0]; i <= box[0]; i++)
  {
    for(int j = -box[1]; j <= box[1]; j++)
    {
      for(int k = -box[2]; k <= box[2]; k++)
      {
        Vector3d vc = vd - cell * Vector3d(i, j, k);
        result.push_back(vc);
      }
    }
  }  
  
  return result;
}        

Eigen::Vector3d cryst_tools::min_dist::average_vector(const std::vector<Eigen::Vector3d> &av) const
{
  assert(cell_assigned);
  
  Vector3d result;
  
  result.setZero();
  
  for(int i = 0; i < av.size(); i++)
    result += (*this)(av[i] - av[0]);

  result = av[0] + result / double(av.size());
  result = (*this)(result);
  
  return result;
  
}

//Ewald sum
void cryst_tools::ewald_sum::set_cell(const Eigen::Matrix3d &cell_v)
{
  cell = cell_v;
  r_cell   = cell.inverse();
  res_cell = 2 * M_PI * (cell.inverse()).transpose();
  volume = fabs(cell.determinant());
}

void cryst_tools::ewald_sum::set_precision(double N, double precision_v)
{
  precision = precision_v;
  
  eta = M_PI * pow(N/volume, 1.0/3.0);
  r_max = sqrt(-log(precision) / eta);
  g_max = 2 * sqrt(-log(precision) * eta);
  
  //cout << "Eta: " << eta << endl;
  //cout << "r_max: " << r_max << endl;
  //cout << "g_max: " << g_max << endl;  
  
  Vector3d box_size_d = res_cell.inverse() * Vector3d(g_max, g_max, g_max);
  for(int i = 0; i < box_size.size(); i++)
    box_size[i] = max(5, int(abs(box_size_d[i])) + 1);
}

double cryst_tools::ewald_sum::get_energy(const Eigen::Vector3d &vd)
{
  double result = 0;
  double real_term = 0;
  double rec_term = 0;
  bool self_iteraction = false;
  for(int i = -box_size[0]; i <= box_size[0]; i++)
  {
    for(int j = -box_size[1]; j <= box_size[1]; j++)
    {
      for(int k = -box_size[2]; k <= box_size[2]; k++)
      {
        Vector3d p(i, j, k);
        double rd = (vd + cell * p).norm();

        // Real-space summation
        if(rd > 1E-3)
          real_term += erfc(sqrt(eta) * rd) / rd;
        else
          self_iteraction = true;

        if( (i != 0) || (j != 0) || (k != 0) )
        {
          // K-space summation
          Vector3d K = res_cell * p;
          double K2 = K.squaredNorm();
          double Kvd = K.transpose() * vd;
          rec_term += cos(Kvd) * exp(-K2 / (4 * eta) )/ K2;
        }  
      }  
    }
  }  
  
  result = 0.5 * real_term + 0.5 * 4 * M_PI / volume * rec_term;
  if(self_iteraction)
    result -= sqrt(eta / M_PI);
  
  return result;
}

Eigen::MatrixXd cryst_tools::ewald_sum::potential_matrix(const std::vector<Eigen::Vector3d> &vd)
{
  Eigen::MatrixXd result;
  
  result.resize(vd.size(), vd.size());
  result.setZero();
  
  for(int i = 0; i < vd.size(); i++)
  {
    for(int j = 0; j < vd.size(); j++)
    {
      result(i, j) = get_energy(vd[i] - vd[j]);
    }  
  }  
  
  return result;
}        
