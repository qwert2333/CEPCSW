#ifndef _CALOSTEP_H
#define _CALOSTEP_H


class CaloStep{

public: 
  CaloStep (double _Q, double _T): Q(_Q), T(_T) {};
  //CaloStep (double _Q, double _T) { Q=_Q; T=_T; };
  CaloStep() {};

  void setQ(double _Q) { Q =_Q; }
  void setT(double _T) { T =_T; }

  double getQ() const { return Q; }
  double getT() const { return T; }
  inline bool operator < (const CaloStep &x) const { 
    return T <x.T ;
  }

private: 
  double Q;
  double T;
};

#endif
