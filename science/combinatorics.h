/* 
 * File:   combinatorics.h
 * Author: kirill
 *
 * Created on July 17, 2013, 11:32 AM
 */

#ifndef COMBINATORICS_H
#define	COMBINATORICS_H

#include <vector>
#include <map>
#include <assert.h>
#include <climits>
#include <algorithm>

template <typename T>
int get_sign(const T &v) {
  if( v == 0)
  {
    return 0;
  }
  else if (v > 0)
  {
    return 1;
  }
  else
  {
    return -1;
  }

}

template <typename To, typename Ti>
bool safe_multiplication(const std::vector<Ti> &v, To &result)
{
  result = 0;
  if(v.empty())
    return false;

  int sign = get_sign(v[0]);
  for(int i = 1; i < v.size(); i++) {
    sign *= get_sign(v[i]);
  }

  if(sign == 0) {
    return true;
  }
  else if (sign < 0 && std::numeric_limits<To>::min() == 0)
  {
    return false;
  }

  result = 1;

  bool overflow = false;
  for(int i = 0; i < v.size(); i++)
  {
    To vm = static_cast<To>(std::abs(v[i]));
    if( vm <= std::numeric_limits<To>::max() / result )
    {
      result = result * vm;
    }
    else
    {
      overflow = true;
      result = 0;
      break;
    }
  }
  if(!overflow)
    result = result * static_cast<To>(sign);

  return !overflow;
}

template <typename T>
bool factorial(int k, T &result)
{
  std::vector<int> v;
  v.reserve(k);
  for(int i = 1; i <= k; i++)
    v.push_back(i);
    
  return safe_multiplication(v, result);
}

inline int __delete_vector_by_value(std::vector<int> &nm, int value, int max_operation = -1)
{
  assert(value > 0);
  int result = 0;
  if(max_operation == 0)
    return 0;

  for(int i = 0; i < nm.size(); i++)
  {
    while(nm[i] % value == 0)
    {
      nm[i] /= value;
      result++;
      if(result == max_operation)
        break;
    }
    if(result == max_operation)
      break;
  }
  return result;
}

template <typename T>
bool num_combinations(const std::vector<int> &nm, T &result)
{
  if(nm.empty()) {
    result = 0;
    return true;
  }
  std::vector<int> top_fract;
  std::vector<int> bottom_fract;

  int max_value = 0;
  int counter = 1;
  for(std::vector<int>::const_iterator it = nm.begin(); it != nm.end(); it++)
  {
    for(int i = 1; i <= *it; i++)
    {
      max_value = std::max(max_value, i);
      bottom_fract.push_back(i);
      top_fract.push_back(counter);
      counter++;
    }
  }

  if(max_value == 0) {
    result = 0;
    return true;
  }


  for(int i = max_value; i > 1; i--)
  {
    std::vector<int> bt1, top1;
    top1 = top_fract;
    bt1 = bottom_fract;
    int num_delete_bottom = __delete_vector_by_value(bt1, i);
    int num_delete_top = __delete_vector_by_value(top1, i);
    int num_delete = std::min(num_delete_bottom, num_delete_top);
    num_delete_bottom = __delete_vector_by_value(bottom_fract, i, num_delete);
    num_delete_top = __delete_vector_by_value(top_fract, i, num_delete);

    assert(num_delete_bottom == num_delete_top && num_delete_top == num_delete);
  }
  int bottom_result;
  bool sm_bottom = safe_multiplication(bottom_fract, bottom_result);
  assert(sm_bottom && (bottom_result == 1) );

  return safe_multiplication(top_fract, result);
}

template < typename Tm, typename To>
bool num_combinations(const std::map<Tm, int> &nm, To & result)
{
  typedef std::map<Tm, int> map_type;
  
  std::vector<int> nmv;
  for(typename map_type::const_iterator it = nm.begin(); it != nm.end(); it++)
    nmv.push_back(it->second);
  
  return num_combinations(nmv, result);
}

template < typename T>
std::vector<T> create_start_combination(const std::map<T, int> &nm)
{
  std::vector<T> result;
  for(typename std::map<T, int>::const_iterator it = nm.begin(); it != nm.end(); it++)
  {
    for(int j = 0; j < (*it).second; j++)
      result.push_back((*it).first);
  }
  
  bool prev_perm_exist = std::prev_permutation(result.begin(), result.end());
  bool next_perm_exist = std::next_permutation(result.begin(), result.end());
  
  assert( (!prev_perm_exist) && (!next_perm_exist) );
  
  return result;
}


#endif	/* COMBINATORICS_H */

