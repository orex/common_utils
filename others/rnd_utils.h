/* 
 * File:   rnd_utils.h
 * Author: kokho
 *
 * Created on February 7, 2011, 5:19 PM
 */

#ifndef RND_UTILS_H
#define	RND_UTILS_H

#include <boost/random.hpp>
#include <vector>
#include <algorithm>

boost::mt19937 create_rnd_gen();
double get_rnd_value_in_interval(boost::mt19937 &rng, double min, double max);
int get_rnd_int_value_in_interval(boost::mt19937 &rng, int min, int max);

std::vector<int> get_random_numbers(int count, int max_number);

template <class RandomAccessIterator, class RandomBoostGenerator>
void br_random_shuffle(RandomAccessIterator first, RandomAccessIterator last,
                       RandomBoostGenerator &rnd)
{
  boost::uniform_int<> uni_dist;
  boost::variate_generator<RandomBoostGenerator, boost::uniform_int<> > randomNumber(rnd, uni_dist);
  std::random_shuffle(first, last, randomNumber);   
}

#endif	/* RND_UTILS_H */

