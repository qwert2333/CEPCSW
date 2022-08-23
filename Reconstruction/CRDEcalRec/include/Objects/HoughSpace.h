#ifndef CALOHOUGHSPACE_H
#define CALOHOUGHSPACE_H
#include "Objects/HoughObject.h"
#include "Objects/LongiCluster.h"
#include "TH2.h"
namespace PandoraPlus {

  class HoughSpace{
  public: 
    HoughSpace() {};
    ~HoughSpace() { Clear(); }
    void Clear() { m_sapceMap.Reset();  m_hills.clear(); }

    //Sub-object classes: Hill. 
	  class HoughHill{
    public:
      HoughHill() {};
      void Clear() {index_alpha.clear(); index_rho.clear(); binSize.clear();}
      std::vector<int> getIndexAlpha() const { return index_alpha; }
      std::vector<int> getIndexRho() const { return index_rho; }
      std::vector<int> getBinSize() const {return binSize; }

      void AddCell(int _alpha, int _rho, int _size) { index_alpha.push_back(_alpha); index_rho.push_back(_rho); binSize.push_back(_size); }

    private: 
      std::vector<int> index_alpha;
      std::vector<int> index_rho;
      std::vector<int> binSize;

	  };


    //Functions
		TH2F getSpaceMap() const { return m_sapceMap; }
    std::vector<PandoraPlus::HoughSpace::HoughHill> getHills() const { return m_hills; }
    double getAlphaBinWidth() const { return m_sapceMap.GetXaxis()->GetBinWidth(0); }
    double getRhoBinWidth() const { return m_sapceMap.GetYaxis()->GetBinWidth(0); }
    double getAlphaLowEdge() const { return m_sapceMap.GetXaxis()->GetBinLowEdge(1); }
    double getAlphaUpEdge() const { return m_sapceMap.GetXaxis()->GetBinUpEdge(m_sapceMap.GetXaxis()->GetNbins()); }
    double getRhoLowEdge() const { return m_sapceMap.GetYaxis()->GetBinLowEdge(1); }
    double getRhoUpEdge() const { return m_sapceMap.GetYaxis()->GetBinUpEdge(m_sapceMap.GetYaxis()->GetNbins()); }

    void setSpaceMap(TH2F& _map) { m_sapceMap=_map; }
    void setIPBand(TF1& _o, TF1& _i) { IPHoughLine_outer=_o; IPHoughLine_inner=_i; }
    void setHoughHills(  std::vector<PandoraPlus::HoughSpace::HoughHill> _hillCol ) { m_hills=_hillCol; }

    bool isFromIP(PandoraPlus::HoughSpace::HoughHill& _hill) const;


  private: 
    TH2F m_sapceMap; 
		std::vector<PandoraPlus::HoughSpace::HoughHill> m_hills; 
    TF1 IPHoughLine_outer;
    TF1 IPHoughLine_inner;

};

};
#endif
