/* 
 * File:   rnd_utils.h
 * Author: kokho
 *
 * Created on February 7, 2011, 5:19 PM
 */

#ifndef RND_UTILS_H
#define	RND_UTILS_H

#include <random>
#include <vector>
#include <algorithm>

template <class Container, class RandomGenerator>
void random_thin_to(Container &data, int target_size, RandomGenerator &rnd)
{
  if( target_size >= data.size() )
    return;
  if( target_size <= 0 ) {
    data.clear();
    return;
  }

  while ( data.size() > std::max(target_size, 0) ) {
    std::uniform_int_distribution<int> ds(0, data.size() - 1);
    auto it = data.begin();
    std::advance(it, ds(rnd));
    data.erase(it);
  }
}


#endif	/* RND_UTILS_H */

