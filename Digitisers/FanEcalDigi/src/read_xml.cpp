#include "DD4hep/DetFactoryHelper.h"
#include "XML/Utilities.h"
#include "XML/Layering.h"

using namespace std;
using namespace dd4hep;
using namespace dd4hep::detail;

void read_xml(){
  xml_h e;
  xml_det_t x_det = e;
  xml_dim_t dim = x_det.dimensions();
  string det_name = x_det.nameStr();
  double zhalf = dim.zhalf();
  double rmin  = dim.rmin();
  double rmax  = dim.rmax();
  double alpha = dim.alpha();
  int    nphi  = dim.nphi();
  int    nz    = dim.nz();
  double gap   = dim.gap();
  std::cout << "rmin    = " << rmin/mm << std::endl
            << "rmax    = " << rmax/mm << std::endl
            << "zhalf   = " << zhalf/mm << std::endl
            << "alpha   = " << alpha/degree << std::endl
            << "nphi    = " << nphi << std::endl
            << "nz      = " << nz << std::endl;  
  return;
}
