#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_combinatorics

#include <boost/test/unit_test.hpp>

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cstdarg>
#include <ctime>
#include <cmath>
#include <complex>
#include <random>
#include <sys/time.h>

#include <boost/format.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>

#include "science/combinatorics.h"
#include "others/rnd_utils.h"


template <typename T>
void multiplication_check(int v1, int v2, int v3, int v4, bool overflow)
{
  std::vector<int> vm;
  vm.push_back(v1);
  vm.push_back(v2);
  vm.push_back(v3);
  vm.push_back(v4);

  T mult_result;
  bool multiplication_ok = safe_multiplication(vm, mult_result);

  std::string msg = (boost::format("v={%1%, %2%, %3%, %4%}") % v1 % v2 % v3 % v4).str();

  if( !overflow )
  {
    int64_t result_ref =
          static_cast<int64_t>(v1)
        * static_cast<int64_t>(v2)
        * static_cast<int64_t>(v3)
        * static_cast<int64_t>(v4);
    BOOST_CHECK_MESSAGE(multiplication_ok, msg);
    BOOST_CHECK_EQUAL(mult_result, result_ref);
    BOOST_CHECK_MESSAGE(mult_result == result_ref, msg);
  }
  else
  {
    BOOST_CHECK_MESSAGE(!multiplication_ok, "Should be overflow: " + msg);
  }
}

int64_t create_combinations(int v0, ...) {
  va_list ap;
  va_start(ap, v0);
  std::vector<int> vc;
  vc.push_back(v0);
  while(true) {
    int n = va_arg(ap, int);
    if (n < 0)
      break;
    vc.push_back(n);
  }
  va_end(ap);

  int64_t result;
  bool good = num_combinations(vc, result);
  if(good) {
    BOOST_CHECK_GE(result, 0);
  } else
    result = -1;

  return result;
}

BOOST_AUTO_TEST_SUITE(CombinatoricsTest)

BOOST_AUTO_TEST_CASE(Test_safe_mult_simple)
{
  multiplication_check<int>(1, 1, 1, 0, false);

  multiplication_check<int>(1, 1, 1, 1, false);
  multiplication_check<unsigned int>(1, 1, 1, 1, false);

  multiplication_check<int>(1, -1, 1, 1, false);
  multiplication_check<unsigned int>(1, -1, 1, 1, true);
  multiplication_check<unsigned int>(1, -1, -1, 1, false);

  multiplication_check<int64_t>(std::numeric_limits<int>::max(), std::numeric_limits<int>::max(), 1, 1, false);
  multiplication_check<int64_t>(std::numeric_limits<int>::max(), std::numeric_limits<int>::max(), 2, 1, false);
  multiplication_check<int64_t>(std::numeric_limits<int>::max(), std::numeric_limits<int>::max(), 2, 2, true);

  multiplication_check<int>(32767, 32767, 2, 2, true);
  multiplication_check<unsigned int>(32767, 32767, 2, 2, false);
}

BOOST_AUTO_TEST_CASE(Test_safe_mult_random)
{
  std::mt19937 rnd;
  struct timeval tv;
  gettimeofday(&tv, NULL);
  unsigned seed = static_cast<unsigned> (tv.tv_usec);
  rnd.seed(seed);
  BOOST_WARN(seed);

  for(int i = 0; i < 200000; i++) {
    int maxn = std::uniform_int_distribution<int>(0, 10)(rnd);
    std::vector<int> vm;
    vm.reserve(maxn);
    bool all_positive = true;
    for(int j = 0; j < maxn; j++)
    {
      int a = std::uniform_int_distribution<int>(-32767, 32767)(rnd);
      all_positive = all_positive && (a > 0);
      vm.push_back(a);
    }

    int64_t mult_64;
    bool good64 = safe_multiplication(vm, mult_64);

    int mult_32;
    bool good32 = safe_multiplication(vm, mult_32);

    if(!good64 || good32)
      BOOST_CHECK_EQUAL(good32, good64);

    if(good64) {
      int64_t ref = 1;
      for(int j = 0; j < maxn; j++)
        ref *= vm[j];
      if( all_positive)
        BOOST_CHECK_GE(mult_64, 0);
      BOOST_CHECK_EQUAL(mult_64, ref);
      BOOST_CHECK_EQUAL(std::abs(mult_64) <= std::numeric_limits<int>::max(), good32);
      if( good32 && good64 )
        BOOST_CHECK_EQUAL(mult_32, mult_64);
    }
  }
}

BOOST_AUTO_TEST_CASE(Test_combinatorics_simple) {
  BOOST_CHECK_EQUAL(create_combinations(0, -1), 0);
  BOOST_CHECK_EQUAL(create_combinations(1, -1), 1);
  BOOST_CHECK_EQUAL(create_combinations(1, 0, -1), 1);
  BOOST_CHECK_EQUAL(create_combinations(1, 1, -1), 2);
  BOOST_CHECK_EQUAL(create_combinations(100, 1, -1), 101);
  BOOST_CHECK_EQUAL(create_combinations(1000, 2, -1), 1002 * 1001 / 2);
}

BOOST_AUTO_TEST_CASE(Test_combinatorics_ref) {
  std::fstream dat("../tests/data_combinatorics/ctest.dat", std::fstream::in);
  std::string line;
  while(std::getline(dat, line).good()) {
    std::vector<std::string> tks;
    boost::split(tks, line, boost::is_any_of(", "), boost::algorithm::token_compress_on);
    if(tks.size() < 2)
      continue;
    int64_t ref = boost::lexical_cast<int64_t>(tks[0]);
    std::vector<int> vc;
    for(int i = 1; i < tks.size(); i++)
      vc.push_back(boost::lexical_cast<int>(tks[i]));

    int64_t result;
    bool good = num_combinations(vc, result);
    if(!good)
      result = -1;
    BOOST_CHECK_MESSAGE(result == ref, std::string("Data: ") + line);
  }
}

BOOST_AUTO_TEST_CASE(Test_next_permutation_it_rnd) {

  std::mt19937 rnd;
  /*struct timeval tv;
  gettimeofday(&tv, NULL);
  unsigned seed = static_cast<unsigned> (tv.tv_usec);
  rnd.seed(seed);
  BOOST_WARN(seed);*/

  for(int i = 0; i < 12; i++) {
    std::vector<int> v(i, 0);
    for(int j = 0; j < v.size(); j++) {
      v[j] = std::uniform_int_distribution<int>(-1, 4)(rnd);
    }
    std::sort(v.begin(), v.end());
    while(true) {
      auto vn = v;
      auto vl = v;
      auto it = next_permutation_it(vn.begin(), vn.end());
      bool b = std::next_permutation(vl.begin(), vl.end());
      BOOST_CHECK_EQUAL(b, (it != vn.end()) );
      BOOST_CHECK_EQUAL_COLLECTIONS(vn.begin(), vn.end(), vl.begin(), vl.end());
      if( it == vn.end() )
        break;
      int ix = std::distance(vn.begin(), it);
      BOOST_CHECK_NE(v[ix], vn[ix]);
      BOOST_CHECK_EQUAL_COLLECTIONS(vn.begin(), vn.begin() + ix,
                                    vl.begin(), vl.begin() + ix);
      swap(v, vn);
    }
  }
}

BOOST_AUTO_TEST_CASE(Test_permutation_index_rnd) {
  std::mt19937 rnd;
  /*struct timeval tv;
  gettimeofday(&tv, NULL);
  unsigned seed = static_cast<unsigned> (tv.tv_usec);
  rnd.seed(seed);
  BOOST_WARN(seed);*/

  for(int i = 0; i < 12; i++) {
    std::vector<int> v(i, 0);
    for(int j = 0; j < v.size(); j++) {
      v[j] = std::uniform_int_distribution<int>(-1, 4)(rnd);
    }
    std::sort(v.begin(), v.end());
    int64_t count = 0;
    do {
      BOOST_CHECK_EQUAL(permutation_index<int64_t>(v.begin(), v.end()), count);
      count++;
    } while(std::next_permutation(v.begin(), v.end()));
  }
}


BOOST_AUTO_TEST_SUITE_END()

