/* 
 * File:   math_utils.cpp
 * Author: kirill
 * 
 * Created on June 6, 2013, 2:22 PM
 */

#include "math_utils.h"
#include <cmath>
#include <algorithm>

using namespace std;

std::vector<double> math_utils::brd_Lorenzian(const std::vector<double> &data, const double gamma)
{
  std::vector<double> result;
  double gamma2 = pow(gamma, 2);
  
  result.resize(data.size(), 0.0);
  
  for(int i = 0; i < data.size(); i++)  
  {
    if( abs(data[i]) < 1E-8) continue;
    for(int j = 0; j < data.size(); j++)
      result[j] += data[i] / (gamma2 + pow(i - j, 2));
  }

  for(int i = 0; i < data.size(); i++)  
    result[i] *= gamma / M_PI;
    
  return result;
}

std::vector<double> math_utils::brd_Gaussian(const std::vector<double> &data, const double sigma, 
                              const double precision)
{
  std::vector<double> result;
  std::vector<double> precalc_data;
  
  int range = int(sigma * precision) + 1;
  double sigma22 = 2 * pow(sigma, 2);
  
  result.resize(data.size(), 0.0);
  
  precalc_data.resize(range + 5, 0.0);

  for(int i = 0; i < precalc_data.size(); i++)  
  {
    precalc_data[i] = exp(-pow(i,2) / sigma22);
  }
  
  for(int i = 0; i < data.size(); i++)  
  {
    if( abs(data[i]) < 1e-8) continue;
    for(int j = max(0, i - range) ; j < min(int(data.size()), i + range); j++)
      result[j] += data[i] * precalc_data[abs(i-j)];
  }

  for(int i = 0; i < data.size(); i++)  
    result[i] /= sigma * sqrt(2 * M_PI);
    
  return result;
  
}

