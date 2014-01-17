/* 
 * File:   eigen2babel.h
 * Author: kirill
 *
 * Created on January 9, 2014, 3:02 PM
 */

#ifndef EIGEN2BABEL_H
#define	EIGEN2BABEL_H

#include <openbabel/math/matrix3x3.h>
#include <openbabel/math/transform3d.h>
#include <Eigen/Core>

template <typename T>
OpenBabel::vector3 e2b_vector(const Eigen::Matrix<T, 3, 1> &v)
{
  return OpenBabel::vector3(v[0], v[1], v[2]);
}

template <typename T>
OpenBabel::matrix3x3 e2b_matrix(const Eigen::Matrix<T, 3, 3> &m)
{
  
  return OpenBabel::matrix3x3(e2b_vector(m.row(0)), 
                              e2b_vector(m.row(1)), 
                              e2b_vector(m.row(2)));
}

template <typename T>
Eigen::Matrix<T, 3, 1> b2e_vector(const OpenBabel::vector3 &v)
{
  return Eigen::Matrix<T, 3, 1>(v.GetX(), v.GetY(), v.GetZ());
}

template <typename T>
Eigen::Matrix<T, 3, 3> b2e_matrix(const OpenBabel::matrix3x3 &m)
{
  Eigen::Matrix<T, 3, 3> result;
  
  result.col(0) = b2e_vector<T>(m.GetColumn(0));
  result.col(1) = b2e_vector<T>(m.GetColumn(1));
  result.col(2) = b2e_vector<T>(m.GetColumn(2));  
  
  return result;
}


template <typename T>
Eigen::Transform<T, 3, Eigen::Affine> b2e_affine(const OpenBabel::transform3d &tf)
{
  using namespace OpenBabel;
  
  Eigen::Transform<T, 3, Eigen::Affine> result;
  result.setIdentity();

  vector3 v = tf * vector3(0, 0, 0);
  
  result.linear().col(0) = b2e_vector<T>(tf * vector3(1.0, 0.0, 0.0) - v);
  result.linear().col(1) = b2e_vector<T>(tf * vector3(0.0, 1.0, 0.0) - v);
  result.linear().col(2) = b2e_vector<T>(tf * vector3(0.0, 0.0, 1.0) - v);
  
  result.translation() = b2e_vector<T>(v);
  
  return result;
}

/*
template <typename T>
OpenBabel::transform3d e2b_affine(const Eigen::Transform<T, 3, Eigen::Affine> &tf)
{
  using namespace OpenBabel;
  
  OpenBabel::transform3d result();
  result.setIdentity();

  vector3 v = tf * vector3(0, 0, 0);
  
  result.linear().col(0) = b2e_vector<T>(tf * vector3(1.0, 0.0, 0.0) - v);
  result.linear().col(1) = b2e_vector<T>(tf * vector3(0.0, 1.0, 0.0) - v);
  result.linear().col(2) = b2e_vector<T>(tf * vector3(0.0, 0.0, 1.0) - v);
  
  result.translation() = b2e_vector<T>(v);
  
  return result;
}*/

#endif	/* EIGEN2BABEL_H */

