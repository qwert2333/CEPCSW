#ifndef _CRD_CALOSTEP_
#define _CRD_CALOSTEP_

namespace CRDEcalDigiEDM {

  class CRDCaloStep{

  public: 
    CRDCaloStep (double _Q, double _T): Q(_Q), T(_T);

    double getQ() const { return Q; }
    double getT() const { return T; }
    inline bool operator < (const CRDCaloStep &x) const { 
      return T <x.T ;
    }

  private: 
    double Q;
    double T;
  };

};
#endif
