#ifndef _CALOCROSSBAR_
#define _CALOCROSSBAR_

#include <DD4hep/Objects.h>
#include "TVector3.h"
#include "CaloBar.h"

using namespace dd4hep;
using namespace std;
using namespace dd4hep::detail;

class CaloCrossBar{

public:
  CaloCrossBar(dd4hep::Position _pos, double _E)
  : position(_pos), E(_E) {}; 
  CaloCrossBar() {};

  double getEnergy()   const { return E;     }
  double getScaleEnergy(double f_E)   { return E*f_E; }
  dd4hep::Position getPosition() const { return position; }
  double getPositionR()  const { return position.x(); }
  double getPositionPhi()  const { return position.y(); }
  double getPositionZ()  const { return position.z(); }
 
  void setPosition( dd4hep::Position pos) { position = pos; }
  void setPosition( TVector3 posv3) { position.SetXYZ( posv3.x(), posv3.y(), posv3.z() ); }
  void setEnergy(double _E) { E=_E; }

private:
	dd4hep::Position position;
	double E;      // Q in left readout
};
  
#endif
