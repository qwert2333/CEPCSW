#ifndef PFO_CREATOR_H
#define PFO_CREATOR_H

#include "PandoraPlusDataCol.h"

#include "TVector3.h"
#include "TLorentzVector.h"

class PFOCreator{

public: 

  class Settings{
  public:
    Settings(){};

  };
  
  //initialize a CaloHitCreator
  PFOCreator( Settings& settings ){};
  ~PFOCreator();

  StatusCode ReadSettings( Settings& settings ){ return StatusCode::SUCCESS; };

  StatusCode CreatePFO(PandoraPlusDataCol& dataCol ){ 
    dataCol.ClearPFO();
    //Only have photon now

    std::vector<CRDEcalEDM::CRDCaloHit3DShower> m_clusCol; m_clusCol.clear();
    m_clusCol = dataCol.Clus3DCol;

    std::vector<CRDEcalEDM::PFObject> m_pfoCol; m_pfoCol.clear();
    for(int icl=0;icl<m_clusCol.size();icl++){

      CRDEcalEDM::PFObject m_pfo; m_pfo.Clear(); 
      m_pfo.setPdgID(22);
      m_pfo.setCharge(0); 
      m_pfo.setEcalShower(m_clusCol[icl]);

      TVector3 pfoP3 = m_clusCol[icl].getAxis(); 
      pfoP3.SetMag( m_clusCol[icl].getShowerE() );

      TLorentzVector pfoP4(pfoP3, m_clusCol[icl].getShowerE() );
      m_pfo.setP4(pfoP4);
      m_pfoCol.push_back(m_pfo);
    }
    dataCol.PFOCol = m_pfoCol; 
    return StatusCode::SUCCESS; 
  }


  void Reset(){};

private: 



};
#endif
