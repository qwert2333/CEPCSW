#ifndef _PFOCREATOR_C
#define _PFOCREATOR_C

#include "Tools/PFOCreator.h"

PFOCreator::PFOCreator(Settings& settings){

};

StatusCode PFOCreator::CreatePFO(PandoraPlusDataCol& dataCol){
  const double Me  = 0.000511;
  const double Mmu = 0.10566;
  const double Mpi = 0.13957;
  const double Mn  = 0.939;
  const double Mgamma = 0.;

  dataCol.ClearPFO();
  //Photon and charged particle now.

  std::vector<CRDEcalEDM::CRDCaloHit3DCluster> m_clusCol = dataCol.Clus3DCol;
  std::vector<CRDEcalEDM::Track>              m_trkCol  = dataCol.TrackCol;


  std::vector<CRDEcalEDM::PFObject> m_pfoCol; m_pfoCol.clear();
  for(int icl=0;icl<m_clusCol.size();icl++){
    m_clusCol[icl].IdentifyCluster();
    CRDEcalEDM::PFObject m_pfo; m_pfo.Clear();
    m_pfo.setEcalShower(m_clusCol[icl]);
/*
    for(int it=0; it<m_trkCol.size(); it++){
      double clus_x = m_clusCol[icl].getShowerCenter().x();
      double var = (m_trkCol[it].getD0()*cos( m_trkCol[it].getPhi0() ) - clus_x)/sin( m_trkCol[it].getPhi0() );
      TVector3 trk_pos = m_trkCol[it].getProspectPos(var);
      bool b_matchtrk = ( ( m_clusCol[icl].getType()!=2 && ((trk_pos - m_clusCol[icl].getShowerCenter()).Mag()<10) ) ||
                          ( m_clusCol[icl].getType()==2 && ((trk_pos - m_clusCol[icl].getShowerCenter()).Mag()<30) )  );
      if( b_matchtrk ){ m_pfo.setTrack( m_trkCol[it] ); break; }
    }

    if(m_pfo.ContainTrack() && m_clusCol[icl].getType()==0  ) { //MIP, Set as Muon.
      if(m_pfo.getTrack().getCharge()>0) m_pfo.setPdgID(-13);
      else m_pfo.setPdgID(13);
    }
    else if(m_pfo.ContainTrack() && m_clusCol[icl].getType()==1){ //Track + EM shower, Set as e+-.
      if(m_pfo.getTrack().getCharge()>0) m_pfo.setPdgID(-11);
      else m_pfo.setPdgID(11);
    }
    else if(m_pfo.ContainTrack() && m_clusCol[icl].getType()==2){ //Track + other shower, set as pi+-
      if(m_pfo.getTrack().getCharge()>0) m_pfo.setPdgID(211);
      else m_pfo.setPdgID(-211);
    }
    else if(!m_pfo.ContainTrack() && m_clusCol[icl].getType()==1) m_pfo.setPdgID(22);
    else m_pfo.setPdgID(2112);
*/

    TLorentzVector pfoP4;
//    if(m_pfo.ContainTrack()){
//      TVector3 pfoP3 = m_pfo.getTrack().getMomentum();
//      if( fabs(m_pfo.getPdgID()==11) )  pfoP4.SetVectM( pfoP3,  Me);
//      if( fabs(m_pfo.getPdgID()==13) )  pfoP4.SetVectM( pfoP3,  Mmu);
//      if( fabs(m_pfo.getPdgID()==211) ) pfoP4.SetVectM( pfoP3,  Mpi);
//    }
//    else{
      TVector3 pfoP3 = m_clusCol[icl].getAxis();
      pfoP3.SetMag( m_clusCol[icl].getShowerE() );
      pfoP4.SetVectM( pfoP3,  Mgamma);
//      if( m_pfo.getPdgID()==22 )   pfoP4.SetVectM( pfoP3,  Mgamma);
//      if( m_pfo.getPdgID()==2112 ) pfoP4.SetVectM( pfoP3,  Mn);
//    }
    m_pfo.setP4(pfoP4);


    m_pfoCol.push_back(m_pfo);
  }
  dataCol.PFOCol = m_pfoCol;
  return StatusCode::SUCCESS;
};


#endif
