#ifndef _CRD_CALOHOUGHSPACE_H
#define _CRD_CALOHOUGHSPACE_H
#include "Objects/CRDHoughObject.h"
#include "Objects/CRDCaloHitLongiCluster.h"
#include "TH2.h"
namespace CRDEcalEDM {

  class CRDHoughSpace{
  public: 
    CRDHoughSpace() {};
    ~CRDHoughSpace() { m_sapceMap.Reset(); m_hills.clear(); }
    void Clear() { m_sapceMap.Reset();  m_hills.clear(); }

    //Sub-object classes: Hill. 
	  class HoughHill{
    public:
      HoughHill() {};
      void Clear() {index_alpha.clear(); index_rho.clear(); binSize.clear();}
      std::vector<int> getIndexAlpha() const { return index_alpha; }
      std::vector<int> getIndexRho() const { return index_rho; }

      void AddCell(int _alpha, int _rho, int _size) { index_alpha.push_back(_alpha); index_rho.push_back(_rho); binSize.push_back(_size); }

    private: 
      std::vector<int> index_alpha;
      std::vector<int> index_rho;
      std::vector<int> binSize;

	  };


    //Functions
		TH2F getSpaceMap() const { return m_sapceMap; }
    std::vector<CRDEcalEDM::CRDHoughSpace::HoughHill> getHills() const { return m_hills; }
    double getAlphaBinWidth() const { return m_sapceMap.GetXaxis()->GetBinWidth(0); }
    double getRhoBinWidth() const { return m_sapceMap.GetYaxis()->GetBinWidth(0); }
    double getAlphaLowEdge() const { return m_sapceMap.GetXaxis()->GetBinLowEdge(1); }
    double getAlphaUpEdge() const { return m_sapceMap.GetXaxis()->GetBinUpEdge(m_sapceMap.GetXaxis()->GetNbins()); }
    double getRhoLowEdge() const { return m_sapceMap.GetYaxis()->GetBinLowEdge(1); }
    double getRhoUpEdge() const { return m_sapceMap.GetYaxis()->GetBinUpEdge(m_sapceMap.GetYaxis()->GetNbins()); }

    void SetSpaceMap(TH2F& _map) { m_sapceMap=_map; }
    void SetHoughHills(  std::vector<CRDEcalEDM::CRDHoughSpace::HoughHill> _hillCol ) { m_hills=_hillCol; }

  private: 
    TH2F m_sapceMap; 
		std::vector<CRDEcalEDM::CRDHoughSpace::HoughHill> m_hills; 


};

};
#endif
