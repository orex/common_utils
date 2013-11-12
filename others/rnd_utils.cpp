/* 
 * File:   rnd_utils.cpp
 * Author: kokho
 * 
 * Created on February 7, 2011, 5:19 PM
 */

#include "rnd_utils.h"
#include <sys/time.h>
#include <cassert>
#include <set>

std::vector<int> get_random_numbers(int count, int max_number)
{
  assert((count >= 0) && (max_number >= 0));
  
  bool hole = count > (max_number / 2);
  
  int count_t = hole ? (max_number + 1 - count) : count;
  
  boost::mt19937 rng = create_rnd_gen();
  
  std::set<int> vc;
  
  while(count_t > 0)
  {
    int num = get_rnd_int_value_in_interval(rng, 0, max_number);
    if(vc.count(num) == 0)
    {
      vc.insert(num);
      count_t--;
    }  
  }
  
  std::vector<int> result;
  
  if(!hole)
  {
    result.resize(vc.size());
    std::copy(vc.begin(), vc.end(), result.begin());
  }  
  else
  {  
    result.clear();
    result.reserve(max_number + 1 - vc.size());
    for(int i = 0; i <= max_number; i++)
    {
      if( vc.count(i) == 0 )
        result.push_back(i);
    }
  }
  
  return result;
}

boost::mt19937 create_rnd_gen()
{
  boost::mt19937 result;
  struct timeval tv;
  gettimeofday(&tv, NULL);
  result.seed(static_cast<unsigned> (tv.tv_usec));

  return result;
}

double get_rnd_value_in_interval(boost::mt19937 &rng, double min, double max)
{
  double result;

  boost::uniform_real<double> uniform(min, max);
  boost::variate_generator< boost::mt19937&, boost::uniform_real<double> > uniform_sampler(rng, uniform);

  result = uniform_sampler();

  return result;
}

int get_rnd_int_value_in_interval(boost::mt19937 &rng, int min, int max)
{
  int result;

  boost::uniform_int<int> uniform(min, max);
  boost::variate_generator< boost::mt19937&, boost::uniform_int<int> > uniform_sampler(rng, uniform);

  result = uniform_sampler();

  return result;
}



