/* 
 * File:   conv_types.cpp
 * Author: kirill
 * 
 * Created on June 5, 2013, 3:13 PM
 */

#include "conv_types.h"
#include <sstream>

using namespace std;

std::string conv_types::int2str(const int &value)
{
  stringstream result;
  
  result << value;
  
  return result.str();
}        

std::string conv_types::real2str(const double &value, int digits, bool e_format)
{
  stringstream result;
  
  if(e_format)
    result << std::scientific;
  else
    result << std::fixed;
  
  result.precision(digits);
  
  result << value;
  
  return result.str();  
}
