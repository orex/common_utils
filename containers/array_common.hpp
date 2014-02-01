/* 
 * File:   array_common.hpp
 * Author: kirill
 *
 * Created on January 13, 2014, 3:13 PM
 */

#ifndef ARRAY_COMMON_HPP
#define	ARRAY_COMMON_HPP

#include <vector>
#include <algorithm>
#include <boost/random.hpp>

namespace array_common
{
  template <typename T>
  void delete_singles(std::vector<T> &data)
  {
    std::sort(data.begin(), data.end());
    int pos_cmb = 0;
    for(int i = 0; i < int(data.size()) - 1; i++)
    {
      if( data[i] == data[i + 1])
      {  
        data[pos_cmb] = data[i];
        pos_cmb++;
      }  
    }
    data.resize(pos_cmb);  
  }
/*
  template<class InputIterator, class RandomBoostGenerator>
  InputIterator random_element(InputIterator first, InputIterator last, RandomBoostGenerator &rnd)
  {
    int size = distance(first, last);
    int elm = 
    for(int i = 0;  )
  }
  
  
  enum thin_type {ttDelFirst, ttDelLast};
  //T - container class
  template <class Container>  
  void thin_to(Container &data, int target_size, thin_type th)
  {
    while (data.size() > target_size) 
    {
      if( th == ttDelFirst )
        data.erase(data.begin());
      else if( th == ttDelLast )
        data.erase(--data.end());
    }        
  }
*/
};


#endif	/* ARRAY_COMMON_HPP */

