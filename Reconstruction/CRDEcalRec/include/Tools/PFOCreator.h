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

//dataCol.Print3DClus();

    std::vector<CRDEcalEDM::CRDCaloHit3DShower> m_clusCol; m_clusCol.clear();
    m_clusCol = dataCol.Clus3DCol;

    std::vector<CRDEcalEDM::PFObject> m_pfoCol; m_pfoCol.clear();
//std::cout<<"Cluster size: "<<m_clusCol.size(); 
    for(int icl=0;icl<m_clusCol.size();icl++){
//std::cout<<"Loop Cluster #"<<icl<<endl;

      CRDEcalEDM::PFObject m_pfo; m_pfo.Clear(); 
      m_pfo.setPdgID(22);
      m_pfo.setCharge(0); 
      m_pfo.setEcalShower(m_clusCol[icl]);

      TVector3 pfoP3 = m_clusCol[icl].getAxis(); 
      pfoP3.SetMag( m_clusCol[icl].getShowerE() );
//printf("PfoP4: (%.2f, %.2f, %.2f, %.2f) \n", pfoP3.x(), pfoP3.y(), pfoP3.z(), m_clusCol[icl].getShowerE());

      TLorentzVector pfoP4(pfoP3, m_clusCol[icl].getShowerE() );
//std::cout<<"here1"<<std::endl;
      m_pfo.setP4(pfoP4);
//std::cout<<"here2"<<std::endl;
      m_pfoCol.push_back(m_pfo);
//std::cout<<"here3"<<std::endl;
    }
//std::cout<<"Loop done. Save in datacol"<<std::endl;
    dataCol.PFOCol = m_pfoCol; 
//std::cout<<"PFOCreator done. Return PandoraPlusPFA"<<std::endl;
    return StatusCode::SUCCESS; 
  }


  void Reset(){};

private: 



};
#endif
