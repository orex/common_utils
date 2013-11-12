#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE random_test

#include <boost/test/unit_test.hpp>

#include <iostream>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <complex>

#include "others/rnd_utils.h"

bool range_and_sort(const std::vector<int> &vc, int max_number)
{
  bool result = true;
  
  for(int i = 0; i < vc.size(); i++)  
  {
    result = result && (vc[i] >= 0) && (vc[i] <= max_number);
    if(i > 0)
      result = result && (vc[i - 1] < vc[i]);
  }  
  
  return result;
}


BOOST_AUTO_TEST_SUITE(RandomTest)

BOOST_AUTO_TEST_CASE(Test_random_array_static)
{
  std::vector<int> vc;
  
  vc = get_random_numbers(10, 100);
  BOOST_CHECK(vc.size() == 10);
  BOOST_CHECK(range_and_sort(vc, 100));

  vc = get_random_numbers(90, 100);
  BOOST_CHECK(vc.size() == 90);
  BOOST_CHECK(range_and_sort(vc, 100));

  vc = get_random_numbers(0, 200);
  BOOST_CHECK(vc.size() == 0);
  BOOST_CHECK(range_and_sort(vc, 200));
  
}

BOOST_AUTO_TEST_CASE(Test_random_array_random)
{
  std::vector<int> vc;
  
  boost::mt19937 rng = create_rnd_gen();
  for(int i = 0; i < 1000; i++)
  {  
    int max_val = get_rnd_int_value_in_interval(rng, 0, 10000);
    int count_val = get_rnd_int_value_in_interval(rng, 0, 10000);
    vc = get_random_numbers(count_val, max_val);
    BOOST_CHECK( (vc.size() == count_val) || (vc.size() == max_val + 1) );
    BOOST_CHECK(range_and_sort(vc, max_val));
  }
}



BOOST_AUTO_TEST_SUITE_END()
