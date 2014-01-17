/* 
 * File:   obabel_utils.h
 * Author: kirill
 *
 * Created on October 30, 2013, 1:44 PM
 */

#ifndef OBABEL_UTILS_H
#define	OBABEL_UTILS_H

#include <vector>
#include "eigen2babel.h"
#include <cryst_tools/cryst_tools.h>

namespace OpenBabel 
{
  class ob_min_dist : protected cryst_tools::min_dist
  {
    protected:
      using min_dist::operator ();
    public:  
      void set_cell(const OpenBabel::matrix3x3 &cell_v)
      {  cryst_tools::min_dist::set_cell(b2e_matrix<double>(cell_v));  };
                
      OpenBabel::vector3 operator ()(const OpenBabel::vector3 &vd) const
      { return e2b_vector((*this)(b2e_vector<double>(vd))); };
      
      std::vector<OpenBabel::vector3> get_img_dist(const OpenBabel::vector3 &vd, 
                                                   int bx, int by, int bz) const
      {
        std::vector<Eigen::Vector3d> rx = cryst_tools::min_dist::get_img_dist(
                                          b2e_vector<double>(vd), Eigen::Vector3i(bx, by, bz));
        std::vector<OpenBabel::vector3> result;
        for(int i = 0; i < rx.size(); i++)
          result.push_back(e2b_vector(rx[i]));
        
        return result;
      }
    
      OpenBabel::vector3 average_vector(const std::vector<OpenBabel::vector3> &av) const
      { 
        std::vector<Eigen::Vector3d> ve;
        for(int i = 0; i < av.size(); i++)
          ve.push_back(b2e_vector<double>(av[i]));
        return e2b_vector(cryst_tools::min_dist::average_vector(ve));
      };
  };
};  

#endif	/* OBABEL_UTILS_H */

