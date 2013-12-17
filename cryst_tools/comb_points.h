/* 
 * File:   comb_points.h
 * Author: kirill
 *
 * Created on December 16, 2013, 9:27 PM
 */

#ifndef COMB_POINTS_H
#define	COMB_POINTS_H



class cmb_point
{
public:
  virtual double distance(const cmb_point &p) const = 0;
};

class cmb_group
{
public:
  std::vector<cmb_point *> pnts;
};

class comb_points 
{
public:
  comb_points();
  
  comb_points(const comb_points& orig);
  virtual ~comb_points();
private:
};

#endif	/* COMB_POINTS_H */

