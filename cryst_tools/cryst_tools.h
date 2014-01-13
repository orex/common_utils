/* 
 * File:   cryst_tools.h
 * Author: kirill
 *
 * Created on December 16, 2013, 4:00 PM
 */

#ifndef CRYST_TOOLS_H
#define	CRYST_TOOLS_H

#include <vector>
#include <Eigen/Core>
#include <Eigen/Geometry>

namespace cryst_tools 
{
  typedef std::vector< std::vector<Eigen::Vector3d> > vc_sets;
  
  
  Eigen::Vector3d min_frac(const Eigen::Vector3d &frac);
  Eigen::Vector3d norm_frac(const Eigen::Vector3d &frac);
  
  std::vector<Eigen::Matrix3d> get_symmetries(const Eigen::Matrix3d &cell, 
                                              const int range = 2);
  
  std::vector<Eigen::Vector3d> get_shifts(const Eigen::Matrix3d &sym_matrix_B,
                                           const Eigen::Matrix3d &cell,  
                                           const std::vector<Eigen::Vector3d> &frac_coords,
                                           double tolerance, bool all_symmetric);

  std::vector<Eigen::Vector3d> shifts_intersect(const Eigen::Matrix3d &cell,  
                                                const std::vector<Eigen::Vector3d> &frac_c1,
                                                const std::vector<Eigen::Vector3d> &frac_c2,
                                                double tolerance);
  
  std::vector<Eigen::Affine3d> get_all_symmetries(const Eigen::Matrix3d &cell,
                                                  const vc_sets &frac_coords,
                                                  const double tol,
                                                  const int range = 2);
};

#endif	/* CRYST_TOOLS_H */

