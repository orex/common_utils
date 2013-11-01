#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE expc_test

#include <boost/test/unit_test.hpp>

#include <iostream>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <complex>

extern "C" {
#include <quadmath.h>
}
  
#include "science/expc.hpp"

double expc_ref_real(double a, double x)
{
  if( fabsq(expq(a * x) - 1) > 0.01)
    return (expq(a * x) - 1) / x;
  
  __float128 sum = 0;
  __float128 p;
  __float128 f = 1;
  int i = 0;

  do
  {  
    p = powq(a * x, i) / f;
    f = f * (i + 2);
    sum += p;            
    i++;
  }  
  while(fabsq(p) > FLT128_EPSILON );  
  
  sum *= a;

  return sum;  
}

double rand_log()
{
  double r1 = -7 + 7 * (double(rand()) / RAND_MAX);
  
  return pow(10, r1);
}

BOOST_AUTO_TEST_SUITE(ExpcTest)

BOOST_AUTO_TEST_CASE(Test_Expc_real)
{
  /*
  BOOST_CHECK(match_wildcard("", ""));
  */

  for(int i = 0; i < 10000; i++)
  {  
    double a = rand_log();
    double x = rand_log();
    
    double dl = abs(expc_ref_real(a, x) - math_functions::expcd(a, x));
    double dp = abs(expc_ref_real(a, x) + math_functions::expcd(a, x));
    
    if( abs(dl / dp) > 10 * std::numeric_limits<double>::epsilon() )
    {  
      std::cout << "a: " << a << "  x: " << x << std::endl;
      std::cout << "ref: " << expc_ref_real(a, x) 
                << "  expc: " << math_functions::expcd(a, x) << std::endl;
      std::cout << dl << std::endl;
    }
  }
}

BOOST_AUTO_TEST_SUITE_END()
