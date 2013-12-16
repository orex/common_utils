/* 
 * File:   obabel_utils.h
 * Author: kirill
 *
 * Created on October 30, 2013, 1:44 PM
 */

#ifndef OBABEL_UTILS_H
#define	OBABEL_UTILS_H

#include <openbabel/mol.h>
#include <vector>

namespace OpenBabel 
{
  OpenBabel::vector3 get_minimal_distance(const OpenBabel::vector3 &dist, 
                                          OpenBabel::OBUnitCell * unitcell);
  
  OpenBabel::vector3 center_mass(const std::vector<OpenBabel::vector3> &atoms_pos,
                                 OpenBabel::OBUnitCell * unitcell,
                                 const double tol);

  std::vector<OpenBabel::vector3> get_image_distances(const OpenBabel::vector3 &dist, 
                                                            OpenBabel::OBUnitCell * unitcell);
  
};  

#endif	/* OBABEL_UTILS_H */

