/* 
 * File:   expc.hpp
 * Author: kirill
 *
 * Created on November 1, 2013, 10:35 AM
 */

#ifndef EXPC_HPP
#define	EXPC_HPP

#include <limits>
#include <complex>

namespace math_functions
{
  template <typename Ta, typename Tx>
  Ta expc(const Ta &a, const Tx &x)
  {
    Ta result;
    //Cubic root of 4! epsilon 
    static Tx const taylor_3_bound = pow(24 * std::numeric_limits<Tx>::epsilon(), 1.0/3.0);
    
    Ta ax = a * x;

    if(abs(ax) > taylor_3_bound )
    {
      result = (exp(ax) - Ta(1)) / x;
    }
    else
    {
      result = a * (Ta(1) + ax / Ta(2) + ax * ax / Ta(6));
    }  
    
    return result;
  }
  
  inline double expcd(const double a, const double x)
  { return expc(a, x); };

  inline std::complex<double> expccd(const std::complex<double> &a, double x)
  { return expc(a, x); };

};

#endif	/* EXPC_HPP */

