#ifndef _CRD_CALOHOUGHSPACE_H
#define _CRD_CALOHOUGHSPACE_H
#include "Objects/CRDHoughObject.h"
#include "Objects/CRDCaloHitLongiCluster.h"
#include "TH2.h"
namespace CRDEcalEDM {

  class CRDHoughSpace{
  public: 
    CRDHoughSpace() {};
    ~CRDHoughSpace() { m_sapceMap.Reset();  m_cells.clear(); m_hills.clear(); }
    void Clear() { m_sapceMap.Reset();  m_cells.clear(); m_hills.clear(); }

    //Sub-object classes: Cell and Hill. 
    class HoughCell{
	  public:
      HoughCell() {};

      void Clear() {m_objs.clear(); index_alpha=-1; index_rho=-1; }
      void SetIndex(int _alpha, int _rho) { index_alpha=_alpha; index_rho=_rho; }
      void SetObjects(std::vector<CRDEcalEDM::CRDHoughObject> _objs) { m_objs=_objs; }
      void AddObject(CRDEcalEDM::CRDHoughObject _obj) { m_objs.push_back(_obj); }

    private: 
		  std::vector<CRDEcalEDM::CRDHoughObject> m_objs; 
      int index_alpha;
			int index_rho; 

	  };

	  class HoughHill{
    public:
      HoughHill() {};
      void Clear() {m_cells.clear(); }

      bool isNeighbor( CRDEcalEDM::CRDHoughSpace::HoughCell _cell ) const; 
      void AddCell( CRDEcalEDM::CRDHoughSpace::HoughCell _cell ) { m_cells.push_back(_cell); }
      CRDEcalEDM::CRDCaloHitLongiCluster TransformToCluster(); 


    private: 
		  std::vector<CRDEcalEDM::CRDHoughSpace::HoughCell> m_cells; 

	  };


    //Functions
		TH2F getSpaceMap() const { return m_sapceMap; }
    std::vector<CRDEcalEDM::CRDHoughSpace::HoughHill> getCells() const { return m_Cells; }
    std::vector<CRDEcalEDM::CRDHoughSpace::HoughHill> getHills() const { return m_hills; }

    void SetSpaceMap(TH2F& _map) { m_sapceMap=_map; }
    void SetHoughCells(  std::vector<CRDEcalEDM::CRDHoughSpace::HoughCell> _cellCol ) { m_cells=_cellCol; }
    void SetHoughHills(  std::vector<CRDEcalEDM::CRDHoughSpace::HoughHill> _hillCol ) { m_cells=_hillCol; }

  private: 
    TH2F m_sapceMap; 
    std::vector<CRDEcalEDM::CRDHoughSpace::HoughCell> m_cells; 
		std::vector<CRDEcalEDM::CRDHoughSpace::HoughHill> m_hills; 


};

};
#endif
