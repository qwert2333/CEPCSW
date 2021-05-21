#ifndef _ENERGYTIMEMATCHING_ALG_C
#define _ENERGYTIMEMATCHING_ALG_C

#include "Algorithm/EnergyTimeMatchingAlg.h"

void EnergyTimeMatchingAlg::Settings::SetInitialValue(){
  chi2Wi_E = 1;
  chi2Wi_T = 10;
  th_chi2  = -1;
  sigmaE   = 0.10;
  sigmaPos = 34.89;  //sqrt(10*10/12 + pow((Tres*C/(2*nMat)),2) );;
  nMat = 2.15;

}


EnergyTimeMatchingAlg::EnergyTimeMatchingAlg(){

}

StatusCode Initialize(){
	return StatusCode::SUCCESS;
}

StatusCode EnergyTimeMatchingAlg::RunAlgorithm( EnergyTimeMatchingAlg::Settings& m_settings, PandoraPlusDataCol& m_datacol ){
  settings = m_settings;

  std::vector<CRDEcalEDM::CRDCaloLayer> layers = m_datacol.LayerCol;
  if(layers.size()==0){ std::cout<<"Warning: Empty input in EnergyTimeMatchingAlg! Please check previous algorithm!"<<std::endl; return StatusCode::SUCCESS; }
  std::vector<CRDEcalEDM::CRDCaloHit2DShower> m_2DshowerCol; m_2DshowerCol.clear();

  for(int il=0;il<layers.size();il++){ 

    std::vector<CRDEcalEDM::CRDCaloHit2DShower> m_showersinLayer; m_showersinLayer.clear();

    std::vector<CRDEcalEDM::CRDCaloBarShower> showerXCol = layers[il].barShowerXCol;
    std::vector<CRDEcalEDM::CRDCaloBarShower> showerYCol = layers[il].barShowerYCol;
    if(showerXCol.size()==0 || showerYCol.size()==0) continue;

  
    //For iteration0: Save all combinations for 3D clustering, leave it for next time. 
    if(m_datacol.Flag_Iter==0){
      for(int is=0;is<showerXCol.size();is++){
      for(int js=0;js<showerYCol.size();js++){
        if(showerXCol[is].getBars().size()==0 || showerYCol[js].getBars().size()==0) continue;
        CRDEcalEDM::CRDCaloHit2DShower tmp_shower; tmp_shower.Clear();
        tmp_shower = DigiHitsWithPos(showerXCol[is], showerYCol[js]);
        m_showersinLayer.push_back(tmp_shower);
      }}
    }

    //For next iteration: use several algorithm to solve ghost hit and multiplicity problem.
    else{
      //Case1: 1*N or N*1
      if(showerXCol.size()==1 || showerYCol.size()==1){ 
        std::vector<CRDEcalEDM::CRDCaloBarShower> shower1Col;
        std::vector<CRDEcalEDM::CRDCaloBarShower> showerNCol;
        if(showerXCol.size()==1) { shower1Col = showerXCol; showerNCol = showerYCol; }
        else { shower1Col = showerYCol; showerNCol = showerXCol; }
   
        const int NshY = showerNCol.size();
        double totE_shY = 0;
        double EshY[NshY] = {0};
        for(int is=0;is<NshY;is++){ EshY[is] = showerNCol[is].getE(); totE_shY += EshY[is]; }
        for(int is=0;is<NshY;is++){
          double wi_E = EshY[is]/totE_shY;
          CRDEcalEDM::CRDCaloBarShower m_splitshower1;
   
          CRDEcalEDM::CRDCaloBar m_wiseed = shower1Col[0].getSeed();
          m_wiseed.setQ( wi_E*m_wiseed.getQ1(), wi_E*m_wiseed.getQ2() );
   
          std::vector<CRDEcalEDM::CRDCaloBar> m_wibars; 
          for(int ib=0;ib<shower1Col[0].getBars().size();ib++){
            CRDEcalEDM::CRDCaloBar m_wibar = shower1Col[0].getBars()[ib];
            m_wibar.setQ(wi_E*m_wibar.getQ1(), wi_E*m_wibar.getQ2());
            m_wibars.push_back(m_wibar);
          }
   
          m_splitshower1.setBars( m_wibars );
          m_splitshower1.setSeed( m_wiseed);
   
          CRDEcalEDM::CRDCaloHit2DShower m_shower; m_shower.Clear();
          if(showerXCol.size()==1) m_shower = DigiHitsWithPos( m_splitshower1, showerNCol[is] );
          else m_shower = DigiHitsWithPos( showerNCol[is], m_splitshower1 );
          m_showersinLayer.push_back( m_shower );
        }
      }
      //Case2: N*N
      else if(showerXCol.size()==showerYCol.size()){
        m_showersinLayer = DigiHitsWithMatching( showerXCol, showerYCol );
      }
      //Case3: M*N (M, N!=1)
      else{
        m_showersinLayer = DigiHitsWithMatchingL2( showerXCol, showerYCol );
      }
      if(m_showersinLayer.size()>std::max(showerXCol.size(), showerYCol.size()))  
      printf( "WARNING! In layer #%d: 2Dshower number(%d) is larger than showerX/Y number(%d / %d)! Please check! \n", il, m_showersinLayer.size(), showerXCol.size(), showerYCol.size() );
      }

    m_2DshowerCol.insert( m_2DshowerCol.end(),  m_showersinLayer.begin(), m_showersinLayer.end() );
  }

  m_datacol.Shower2DCol = m_2DshowerCol;
  return StatusCode::SUCCESS;
}



CRDEcalEDM::CRDCaloHit2DShower EnergyTimeMatchingAlg::DigiHitsWithPos( CRDEcalEDM::CRDCaloBarShower& barShowerX,  CRDEcalEDM::CRDCaloBarShower& barShowerY  ){

  CRDEcalEDM::CRDCaloHit2DShower m_2dshower; m_2dshower.Clear();

  std::vector<edm4hep::ConstCalorimeterHit> m_digiCol; m_digiCol.clear();
  int NbarsX = barShowerX.getBars().size();
  int NbarsY = barShowerY.getBars().size();
  if(NbarsX==0 || NbarsY==0){ std::cout<<"WARNING: empty DigiHitsCol returned!"<<std::endl; return m_2dshower;}

  int _module = barShowerX.getBars()[0].getModule(); 
  int _dlayer = barShowerX.getBars()[0].getDlayer();
  int _part   = barShowerX.getBars()[0].getPart();
  int _stave  = barShowerX.getBars()[0].getStave();

  float rotAngle = -_module*PI/4.;

  TVector3 m_vec(0,0,0);
  for(int ibar=0;ibar<NbarsX;ibar++){
    CRDEcalEDM::CRDCaloBar barx = barShowerX.getBars()[ibar];

    m_vec.SetXYZ(barx.getPosition().x(), barx.getPosition().y(), barx.getPosition().z());
    m_vec.RotateZ(rotAngle);
    barx.setPosition( m_vec );

    for(int jbar=0;jbar<NbarsY;jbar++){
      CRDEcalEDM::CRDCaloBar bary = barShowerY.getBars()[jbar];
      m_vec.SetXYZ(bary.getPosition().x(), bary.getPosition().y(), bary.getPosition().z());
      m_vec.RotateZ(rotAngle);
      bary.setPosition( m_vec );

      TVector3 p_hit(bary.getPosition().x(), (barx.getPosition().y()+bary.getPosition().y())/2., barx.getPosition().z() );  
      p_hit.RotateZ(-rotAngle);
      edm4hep::Vector3f m_vec3f(p_hit.x(), p_hit.y(), p_hit.z());
      float m_En = barx.getEnergy()*bary.getEnergy()/barShowerY.getE() + barx.getEnergy()*bary.getEnergy()/barShowerX.getE();
      edm4hep::CalorimeterHit hit;
      hit.setCellID(0);
      hit.setPosition(m_vec3f);
      hit.setEnergy(m_En);
      m_digiCol.push_back(hit);
    }
  }

  m_2dshower.setBarShowers( barShowerX, barShowerY );
  m_2dshower.setCaloHits( m_digiCol );
  m_2dshower.setIDInfo( _module, _stave, _dlayer, _part);
  return m_2dshower;

}


std::vector<CRDEcalEDM::CRDCaloHit2DShower> EnergyTimeMatchingAlg::DigiHitsWithMatching( std::vector<CRDEcalEDM::CRDCaloBarShower>& barShowerXCol, std::vector<CRDEcalEDM::CRDCaloBarShower>& barShowerYCol ){

  std::vector<CRDEcalEDM::CRDCaloHit2DShower> m_showerCol; m_showerCol.clear();

  const int Nshower = barShowerXCol.size();
  double chi2[Nshower][Nshower];
  double chi2_E[Nshower][Nshower];
  double chi2_tx[Nshower][Nshower];
  double chi2_ty[Nshower][Nshower];

  double wi_E = settings.chi2Wi_E/(settings.chi2Wi_E + settings.chi2Wi_T);
  double wi_T = settings.chi2Wi_T/(settings.chi2Wi_E + settings.chi2Wi_T);

  TVector3 m_vec(0,0,0);
  double rotAngle = -(barShowerXCol[0].getBars())[0].getModule()*PI/4.;
  TVector3 Cblock((barShowerXCol[0].getBars())[0].getPosition().x(), (barShowerXCol[0].getBars())[0].getPosition().y(), (barShowerYCol[0].getBars())[0].getPosition().z());
  Cblock.RotateZ(rotAngle);

  for(int ix=0;ix<Nshower;ix++){
  for(int iy=0;iy<Nshower;iy++){
    CRDEcalEDM::CRDCaloBarShower showerX = barShowerXCol[ix];
    CRDEcalEDM::CRDCaloBarShower showerY = barShowerYCol[iy];

    double Ex = showerX.getE();
    double Ey = showerY.getE();
    chi2_E[ix][iy] = pow(fabs(Ex-Ey)/settings.sigmaE, 2);

    double PosTx = C*(showerY.getT1()-showerY.getT2())/(2*settings.nMat) + showerY.getPos().z();
    chi2_tx[ix][iy] = pow( fabs(PosTx-showerX.getPos().z())/settings.sigmaPos, 2 );

    double PosTy = C*(showerX.getT1()-showerX.getT2())/(2*settings.nMat);
    m_vec.SetXYZ(showerY.getPos().x(), showerY.getPos().y(), showerY.getPos().z());
    m_vec.RotateZ(rotAngle);
    chi2_ty[ix][iy] = pow( fabs(PosTy - (m_vec-Cblock).x() )/settings.sigmaPos, 2);

    chi2[ix][iy] = chi2_E[ix][iy]*wi_E + (chi2_tx[ix][iy]+chi2_ty[ix][iy])*wi_T ;

  }}


  int Ncomb=1;
  for(int i=Nshower; i>0; i--) Ncomb = Ncomb*i;

  map<double, vector<pair<int, int>> > matchingMap;
  int num[Nshower];
  int num_init[Nshower];
  for(int i=0;i<Nshower;i++){ num[i]=i; num_init[i]=i;}

  for(int icont=0;icont<Ncomb;icont++){
    vector<pair<int, int>> Index;
    for(int i=0;i<Nshower;i++){
       pair<int, int> p1(num_init[i], num[i]);
       Index.push_back(p1);
    }
    double chi2_tot=0;
    for(int i=0;i<Index.size();i++) chi2_tot += chi2[Index[i].first][Index[i].second];
    matchingMap[chi2_tot] = Index;

    Index.clear();
    if(!next_permutation(num, num+Nshower)) break;
  }

  map<double, vector<pair<int, int>> >::iterator iter = matchingMap.begin();
  vector<pair<int, int>> Index = iter->second;

  for(int i=0;i<Index.size();i++){
    CRDEcalEDM::CRDCaloBarShower showerX = barShowerXCol[Index[i].first];
    CRDEcalEDM::CRDCaloBarShower showerY = barShowerYCol[Index[i].second];

    CRDEcalEDM::CRDCaloHit2DShower tmp_shower; tmp_shower.Clear();
    tmp_shower = DigiHitsWithPos(showerX, showerY);
    m_showerCol.push_back(tmp_shower);
  }

  return m_showerCol;


}


std::vector<CRDEcalEDM::CRDCaloHit2DShower> EnergyTimeMatchingAlg::DigiHitsWithMatchingL2( std::vector<CRDEcalEDM::CRDCaloBarShower>& barShowerXCol, std::vector<CRDEcalEDM::CRDCaloBarShower>& barShowerYCol ){
  std::vector<CRDEcalEDM::CRDCaloHit2DShower> m_showerCol; m_showerCol.clear();

  const int NshowerX = barShowerXCol.size();
  const int NshowerY = barShowerYCol.size();

  double chi2[NshowerX][NshowerY];
  double chi2_E[NshowerX][NshowerY];
  double chi2_tx[NshowerX][NshowerY];
  double chi2_ty[NshowerX][NshowerY];

  double wi_E = settings.chi2Wi_E/(settings.chi2Wi_E + settings.chi2Wi_T);
  double wi_T = settings.chi2Wi_T/(settings.chi2Wi_E + settings.chi2Wi_T);
  double min_chi2 = 999; 

  TVector3 m_vec(0,0,0);
  double rotAngle = -(barShowerXCol[0].getBars())[0].getModule()*PI/4.;
  TVector3 Cblock((barShowerXCol[0].getBars())[0].getPosition().x(), (barShowerXCol[0].getBars())[0].getPosition().y(), (barShowerYCol[0].getBars())[0].getPosition().z());
  Cblock.RotateZ(rotAngle);

  for(int ix=0;ix<NshowerX;ix++){
  for(int iy=0;iy<NshowerY;iy++){
    CRDEcalEDM::CRDCaloBarShower showerX = barShowerXCol[ix];
    CRDEcalEDM::CRDCaloBarShower showerY = barShowerYCol[iy];

    double Ex = showerX.getE();
    double Ey = showerY.getE();
    chi2_E[ix][iy] = pow(fabs(Ex-Ey)/settings.sigmaE, 2);

    double PosTx = C*(showerY.getT1()-showerY.getT2())/(2*settings.nMat) + showerY.getPos().z();
    chi2_tx[ix][iy] = pow( fabs(PosTx-showerX.getPos().z())/settings.sigmaPos, 2);

    double PosTy = C*(showerX.getT1()-showerX.getT2())/(2*settings.nMat);
    m_vec.SetXYZ(showerY.getPos().x(), showerY.getPos().y(), showerY.getPos().z());
    m_vec.RotateZ(rotAngle);
    chi2_ty[ix][iy] = pow( fabs(PosTy - (m_vec-Cblock).x() )/settings.sigmaPos, 2);

    chi2[ix][iy] = chi2_E[ix][iy]*wi_E + (chi2_tx[ix][iy]+chi2_ty[ix][iy])*wi_T ;

    if(chi2[ix][iy]<min_chi2) min_chi2 = chi2[ix][iy];
  }}

  //WIP-- Now: For the hit with minimum chi2, regared it as a 1*1 hit, and remove the corresponding X&Y showers. 
  //           This only work for 2*3 case, and may work for N*(N+1) case after recursion. 
  //      Future: Set a threshold in settings th_chi2, if( hit_chi2<th_chi2 )  regard it as a 1*1 hit and remove X&Y showers. 
  for(int ix=0;ix<NshowerX;ix++){
  for(int iy=0;iy<NshowerY;iy++){

    //if(chi2[ix][iy]<settings.th_chi2){
    if(chi2[ix][iy]==min_chi2){
    CRDEcalEDM::CRDCaloBarShower showerX = barShowerXCol[ix];
    CRDEcalEDM::CRDCaloBarShower showerY = barShowerYCol[iy];
    if(showerX.getBars().size()==0 || showerY.getBars().size()==0) continue;

    CRDEcalEDM::CRDCaloHit2DShower tmp_shower; tmp_shower.Clear();
    tmp_shower = DigiHitsWithPos(showerX, showerY);
    m_showerCol.push_back(tmp_shower);

    CRDEcalEDM::CRDCaloBarShower m_empty; m_empty.Clear();
    barShowerXCol[ix]=m_empty;
    barShowerYCol[iy]=m_empty;
    }
  }}

  for(int is=0;is<barShowerXCol.size();is++){
  for(int js=0;js<barShowerYCol.size();js++){
    if(barShowerXCol[is].getBars().size()==0 || barShowerYCol[js].getBars().size()==0) continue;
    CRDEcalEDM::CRDCaloHit2DShower tmp_shower; tmp_shower.Clear();
    tmp_shower = DigiHitsWithPos(barShowerXCol[is], barShowerYCol[js]);
    m_showerCol.push_back(tmp_shower);
  }}

  return m_showerCol;
}

#endif
