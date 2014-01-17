#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE cryst_tool_test

#include <boost/test/unit_test.hpp>

#include <iostream>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <complex>
#include <vector>
#include <boost/filesystem.hpp>

#include "cryst_tools/cryst_tools.h"
#include "openbabel/mol.h"
#include "cryst_tools/comb_points.h"

#include <openbabel/obconversion.h>
#include <obabel/eigen2babel.h>

using namespace std;
using namespace OpenBabel;
using namespace cryst_tools;
using namespace Eigen;

BOOST_AUTO_TEST_SUITE(CrystToolsTest)

vc_sets get_frac_coords(OBMol &mol, OBUnitCell * uc)
{
  vc_sets result;
  result.clear();
  
  typedef map<int, vector<Vector3d> > tm;
  tm mpv;
  
  for(OBAtomIterator it = mol.BeginAtoms(); it != mol.EndAtoms(); ++it)
  {
    int num = (*it)->GetAtomicNum();
    vector3 v3 = uc->CartesianToFractional((*it)->GetVector());
    Vector3d vf(v3.GetX(), v3.GetY(), v3.GetZ());
    mpv[num].push_back(vf);
  }
  
  for(tm::const_iterator it = mpv.begin(); it != mpv.end(); ++it)
  {  
    //cout << etab.GetSymbol(it->first) << ": " << it->second.size() << endl;
    result.push_back(it->second);
  }  
  
  return result;
}

BOOST_AUTO_TEST_CASE(Test_first)
{
  namespace fs = boost::filesystem;
  fs::path test_data_dir(CT_DATA_DIR);
  fs::directory_iterator end_iter;

  
  BOOST_CHECK(fs::exists(test_data_dir));
  BOOST_CHECK(fs::is_directory(test_data_dir));

  for( fs::directory_iterator dir_iter(test_data_dir) ; dir_iter != end_iter ; ++dir_iter)
  {
    if (fs::is_regular_file(dir_iter->status()) )
    {
      if((*dir_iter).path().extension().string() != ".cif")
        continue;
      OBConversion obc;
      OBMol mol_orig, mol_p1;
      obc.SetInAndOutFormats("CIF", "CIF");
      BOOST_CHECK(obc.ReadFile(&mol_orig, (*dir_iter).path().string()));
      BOOST_CHECK(obc.ReadFile(&mol_p1, (*dir_iter).path().string()));      
      
      cout << "Checking file: " << (*dir_iter).path() << endl;
      
      OBUnitCell *orig_unitcell = (OBUnitCell *)mol_orig.GetData(OBGenericDataType::UnitCell);
      OBUnitCell *p1_unitcell   = (OBUnitCell *)mol_p1.GetData(OBGenericDataType::UnitCell);
      
      p1_unitcell->FillUnitCell(&mol_p1);
      
      //BOOST_CHECK(obc.WriteFile(&mol_p1, (*dir_iter).path().string() + "_p1"));
      
      vc_sets frc = get_frac_coords(mol_p1, p1_unitcell);
      
      Matrix3d cell;
      
      cell = b2e_matrix<double>(p1_unitcell->GetCellMatrix().transpose());
      
      vector<Affine3d> syms = get_all_symmetries(cell, frc, 1E-2, 1);
      
      /*
      for(int i = 0; i < syms.size(); i++)
      {
        cout << "Symmetry " << i << ": "<< endl;
        cout << "Matrix: " << endl << syms[i].linear() << endl;
        cout << "Move: " << endl << syms[i].translation().transpose() << endl;
        //cout << "Coord: " << (syms[i] * Vector3d(0, 0, 0)).transpose() << endl;
      }*/
      
      transform3dIterator ti;
      const SpaceGroup* pSG = orig_unitcell->GetSpaceGroup();   
      
      const transform3d *t = pSG->BeginTransform(ti);
      int num_uc_transforms = 0;
      while(t)
      {
        Affine3d et = b2e_affine<double>(*t);
        bool exists = false;
        for(int i = 0; i < syms.size(); i++)
        {
          exists = ((syms[i].linear() - et.linear()).norm() < 1E-4) &&
                   (min_frac(syms[i].translation() - et.translation()).norm() < 1E-4);
          if(exists)        
            break;
        }
        BOOST_CHECK(exists);
        t = pSG->NextTransform(ti);
        num_uc_transforms++;
      }
      BOOST_CHECK_EQUAL(syms.size(), num_uc_transforms);
    }
  }
}

BOOST_AUTO_TEST_SUITE_END()
