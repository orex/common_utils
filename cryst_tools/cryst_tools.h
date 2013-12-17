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

namespace cryst_tools 
{
  std::vector<Eigen::Matrix3d> get_symmetries(const Eigen::Matrix3d &cell, 
                                              const int range = 2);
  
  std::vector<Eigen::Vector3d> get_shifts(const Eigen::Matrix3d &sym_matrix_B,
                                           const Eigen::Matrix3d &cell,  
                                           const std::vector<Eigen::Vector3d> &frac_coords,
                                           double tolerance, bool all_symmetric);
  
};

#endif	/* CRYST_TOOLS_H */

