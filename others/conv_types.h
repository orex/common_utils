/* 
 * File:   conv_types.h
 * Author: kirill
 *
 * Created on June 5, 2013, 3:13 PM
 */

#ifndef CONV_TYPES_H
#define	CONV_TYPES_H

#include <string>

class conv_types 
{
  public:
  static std::string int2str(const int &value);
  static std::string real2str(const double &value, int digits = 3, bool e_format = false);
};

#endif	/* CONV_TYPES_H */

