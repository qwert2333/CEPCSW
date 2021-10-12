#ifndef _ENERGYTIMEMATCHING_ALG_C
#define _ENERGYTIMEMATCHING_ALG_C

#include "Algorithm/EnergyTimeMatchingAlg.h"

void EnergyTimeMatchingAlg::Settings::SetInitialValue(){
  chi2Wi_E = 1;
  chi2Wi_T = 10;
  th_chi2  = -1;
  sigmaE   = 0.10;
  fl_useShadowClus = 0;
  fl_UseChi2 = false;
  Debug = 0;
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
  std::vector<CRDEcalEDM::CRDCaloHitTransShower> m_2DshowerCol; m_2DshowerCol.clear();

  //Don't use shadow cluster or chi2 matching, Save all combinations for 3D clustering. 
  if(!settings.fl_useShadowClus && !settings.fl_UseChi2){
    for(int il=0;il<layers.size();il++){
      std::vector<CRDEcalEDM::CRDCaloBarShower> showerXCol = layers[il].barShowerXCol;
      std::vector<CRDEcalEDM::CRDCaloBarShower> showerYCol = layers[il].barShowerYCol;

      std::vector<CRDEcalEDM::CRDCaloHitTransShower> m_showerinlayer; m_showerinlayer.clear(); 
      GetFullMatchedShowers( showerXCol, showerYCol, m_showerinlayer );

      m_2DshowerCol.insert( m_2DshowerCol.end(), m_showerinlayer.begin(), m_showerinlayer.end() );

    }
    m_datacol.Shower2DCol = m_2DshowerCol; 

//m_datacol.PrintLayer();
//m_datacol.PrintShower();
    return StatusCode::SUCCESS;
  }


  //Not use shadow cluster: ET chi2 matching
  if(!settings.fl_useShadowClus && settings.fl_UseChi2){
    for(int il=0;il<layers.size();il++){

      std::vector<CRDEcalEDM::CRDCaloBarShower> showerXCol = layers[il].barShowerXCol;
      std::vector<CRDEcalEDM::CRDCaloBarShower> showerYCol = layers[il].barShowerYCol;
      if(showerXCol.size()==0 || showerYCol.size()==0) continue;

      std::vector<CRDEcalEDM::CRDCaloHitTransShower> m_showerinlayer; m_showerinlayer.clear();
      CRDEcalEDM::CRDCaloHitTransShower tmp_shower; tmp_shower.Clear(); 
      if(showerXCol.size()==1 && showerYCol.size()==1){ XYShowerMatchingL0( showerXCol[0], showerYCol[0], tmp_shower ); m_showerinlayer.push_back(tmp_shower); }
      else if(showerXCol.size()==1) XYShowerMatchingL1( showerXCol[0], showerYCol, m_showerinlayer);
      else if(showerYCol.size()==1) XYShowerMatchingL1( showerYCol[0], showerXCol, m_showerinlayer);
      else if(showerXCol.size()==showerYCol.size()) XYShowerChi2Matching( showerXCol, showerYCol, m_showerinlayer);
      else XYShowerChi2MatchingL1( showerXCol, showerYCol, m_showerinlayer );
      //else GetFullMatchedShowers( showerXCol, showerYCol, m_showerinlayer );

      for(int is=0; is<m_showerinlayer.size(); is++) m_showerinlayer[is].setShadowClusType(0);

      m_2DshowerCol.insert( m_2DshowerCol.end(), m_showerinlayer.begin(), m_showerinlayer.end() );
    }
    m_datacol.Shower2DCol = m_2DshowerCol;

//m_datacol.PrintLayer();
//m_datacol.PrintShower();
    return StatusCode::SUCCESS;
  }


  //Use shadow cluster
//cout<<"EnergyTimeMatchingAlg:: Use shadow cluster to match X/Y."<<endl;
  for(int il=0;il<layers.size();il++){ 
//cout<<"EnergyTimeMatchingAlg: Matching Layer dlayer = "<< layers[il].getDlayer()<<endl;
    std::vector<CRDEcalEDM::CRDCaloBarShower> showerXCol; showerXCol.clear(); 
    std::vector<CRDEcalEDM::CRDCaloBarShower> showerYCol; showerYCol.clear();
    std::vector<CRDEcalEDM::CRDCaloBarShower> m_shX_noCandi; m_shX_noCandi.clear(); 
    std::vector<CRDEcalEDM::CRDCaloBarShower> m_shY_noCandi; m_shX_noCandi.clear(); 

    for(int is=0; is<layers[il].barShowerXCol.size(); is++){
      if( layers[il].barShowerXCol[is].getAllCandiCol().size()==0 ){
        //std::cout<<"WARNING: No shadow cluster in this showerX! "<<std::endl;
        m_shX_noCandi.push_back(layers[il].barShowerXCol[is]);
      }
      else showerXCol.push_back(layers[il].barShowerXCol[is]);
    }
    for(int is=0; is<layers[il].barShowerYCol.size(); is++){
      if( layers[il].barShowerYCol[is].getAllCandiCol().size()==0 ){
        //std::cout<<"WARNING: No shadow cluster in this showerY! "<<std::endl;
        m_shY_noCandi.push_back(layers[il].barShowerYCol[is]);
      }
      else showerYCol.push_back(layers[il].barShowerYCol[is]);
    }

//cout<<"  ShadowClus shower number:  X:"<<showerXCol.size()<<"  Y:"<<showerYCol.size()<<endl;
//cout<<"    ShadowClus size X: "; for(int ic=0; ic<showerXCol.size(); ic++) printf("%d(%d)  ", showerXCol[ic].getTrkCandiCol().size(), showerXCol[ic].getNeuCandiCol().size()); 
//cout<<" | Y: "; for(int ic=0; ic<showerYCol.size(); ic++) printf("%d(%d)  ", showerYCol[ic].getTrkCandiCol().size(), showerYCol[ic].getNeuCandiCol().size() );  cout<<endl;
//cout<<"  No ShadowClus shower number:  X:"<<m_shX_noCandi.size()<<"  Y:"<<m_shY_noCandi.size()<<endl;

    //Pre: matching showers without shadow cluster
    if(m_shX_noCandi.size()!=0 && m_shY_noCandi.size()!=0){
      std::vector<CRDEcalEDM::CRDCaloHitTransShower> m_showerinlayer; m_showerinlayer.clear();
      GetFullMatchedShowers( m_shX_noCandi, m_shY_noCandi, m_showerinlayer );

      m_2DshowerCol.insert( m_2DshowerCol.end(), m_showerinlayer.begin(), m_showerinlayer.end() );
    }
//cout<<"  After Pre-matching. Shower number: "<<showerXCol.size()<<"  "<<showerYCol.size()<<endl;
    if(showerXCol.size()==0 || showerYCol.size()==0) continue;


    //Level 0: pickout 1*1 case.
    std::vector<CRDEcalEDM::CRDCaloHitTransShower> m_lv0showerinlayer; m_lv0showerinlayer.clear(); 
    GetMatchedShowersL0( showerXCol, showerYCol, m_lv0showerinlayer  );

    m_2DshowerCol.insert( m_2DshowerCol.end(), m_lv0showerinlayer.begin(), m_lv0showerinlayer.end() );
//cout<<"  After LV0. Shower number: "<<showerXCol.size()<<"  "<<showerYCol.size()<<endl;
    if(showerXCol.size()==0 || showerYCol.size()==0) continue;


    //Level 1: pickout 1*N case.
    std::vector<CRDEcalEDM::CRDCaloHitTransShower> m_lv1showerinlayer; m_lv1showerinlayer.clear();
    GetMatchedShowersL1( showerXCol, showerYCol, m_lv1showerinlayer );
    GetMatchedShowersL1( showerYCol, showerXCol, m_lv1showerinlayer );

    m_2DshowerCol.insert( m_2DshowerCol.end(), m_lv1showerinlayer.begin(), m_lv1showerinlayer.end() );
//cout<<"  After LV1. Shower number: "<<showerXCol.size()<<"  "<<showerYCol.size()<<endl;
    if(showerXCol.size()==0 || showerYCol.size()==0) continue;


    //Level 2: M*N case. Need to solve a matrix. WIP.
    std::vector<CRDEcalEDM::CRDCaloHitTransShower> m_lv2showerinlayer; m_lv2showerinlayer.clear();
    GetMatchedShowersL2( showerXCol, showerYCol, m_lv2showerinlayer );

//cout<<"  After LV2. Shower number: "<<showerXCol.size()<<"  "<<showerYCol.size()<<endl;
    m_2DshowerCol.insert( m_2DshowerCol.end(), m_lv2showerinlayer.begin(), m_lv2showerinlayer.end() );

  }

  m_datacol.Shower2DCol = m_2DshowerCol;

//m_datacol.PrintLayer();
//m_datacol.PrintShower();

  return StatusCode::SUCCESS;
}


StatusCode EnergyTimeMatchingAlg::GetFullMatchedShowers( 
           std::vector<CRDEcalEDM::CRDCaloBarShower>& barShowerXCol, 
           std::vector<CRDEcalEDM::CRDCaloBarShower>& barShowerYCol, 
           std::vector<CRDEcalEDM::CRDCaloHitTransShower>& outshCol )
{
  outshCol.clear(); 

  for(int is=0;is<barShowerXCol.size();is++){
  for(int js=0;js<barShowerYCol.size();js++){
    if(barShowerXCol[is].getBars().size()==0 || barShowerYCol[js].getBars().size()==0) continue;
    CRDEcalEDM::CRDCaloHitTransShower tmp_shower; tmp_shower.Clear();
    XYShowerMatchingL0(barShowerXCol[is], barShowerYCol[js], tmp_shower);
    tmp_shower.setShadowClusType(0); 
    outshCol.push_back(tmp_shower);
  }}

  return StatusCode::SUCCESS;
}


StatusCode EnergyTimeMatchingAlg::GetMatchedShowersL0(
           std::vector<CRDEcalEDM::CRDCaloBarShower>& barShowerXCol,
           std::vector<CRDEcalEDM::CRDCaloBarShower>& barShowerYCol,
           std::vector<CRDEcalEDM::CRDCaloHitTransShower>& outshCol )
{
//cout<<"  EnergyTimeMatchingAlg::GetMatchedShowersL0"<<endl;

  std::vector< std::pair< CRDEcalEDM::CRDCaloBarShower, CRDEcalEDM::CRDCaloBarShower> > m_showerpairCol; 
  if(barShowerXCol.size()==0 || barShowerYCol.size()==0) return StatusCode::FAILURE; 


  //Save the 1*1 shower into pairVec. 
  for(int is=0; is<barShowerXCol.size(); is++){
    if(barShowerXCol[is].getAllCandiCol().size()!=1) continue;

    for(int js=0; js<barShowerYCol.size(); js++){
      if(barShowerYCol[js].getAllCandiCol().size()!=1) continue;
      if(barShowerXCol[is].getAllCandiCol()[0] == barShowerYCol[js].getAllCandiCol()[0]){
        std::pair<CRDEcalEDM::CRDCaloBarShower, CRDEcalEDM::CRDCaloBarShower> m_shpair( barShowerXCol[is], barShowerYCol[js]);
        m_showerpairCol.push_back(m_shpair);
        barShowerXCol.erase(barShowerXCol.begin()+is); is--;
        barShowerYCol.erase(barShowerYCol.begin()+js); js--;
        break;
      }
    } 
  }
//cout<<"  GetMatchedShowersL0: pair number: "<<m_showerpairCol.size()<<endl;

  //Match in pairVec.
  for(int ip=0; ip<m_showerpairCol.size(); ip++){
    CRDEcalEDM::CRDCaloHitTransShower tmp_shower; tmp_shower.Clear();
    XYShowerMatchingL0(m_showerpairCol[ip].first, m_showerpairCol[ip].second, tmp_shower);   

    tmp_shower.setShadowClus(m_showerpairCol[ip].first.getAllCandiCol()[0]);
    tmp_shower.setShadowClusType( m_showerpairCol[ip].first.getAllCandiCol()[0].Type ); 

    outshCol.push_back(tmp_shower);
  }

  return StatusCode::SUCCESS;
}


StatusCode EnergyTimeMatchingAlg::GetMatchedShowersL1(
           std::vector<CRDEcalEDM::CRDCaloBarShower>& barShowerXCol,
           std::vector<CRDEcalEDM::CRDCaloBarShower>& barShowerYCol,
           std::vector<CRDEcalEDM::CRDCaloHitTransShower>& outshCol )
{
//cout<<"  EnergyTimeMatchingAlg::GetMatchedShowersL1"<<endl;

  std::vector<std::pair< CRDEcalEDM::CRDCaloBarShower, std::vector<CRDEcalEDM::CRDCaloBarShower> >> m_showerpairCol;
  if(barShowerXCol.size()==0 || barShowerYCol.size()==0) return StatusCode::FAILURE;

  for(int is=0; is<barShowerXCol.size(); is++){
    std::vector<CRDEcalEDM::CRDShadowCluster> m_candiCol = barShowerXCol[is].getAllCandiCol();
    std::vector<CRDEcalEDM::CRDCaloBarShower> m_shCol; m_shCol.clear(); 

    for(int js=0; js<barShowerYCol.size(); js++){
      if(barShowerYCol[js].getAllCandiCol().size()!=1) continue;

      CRDEcalEDM::CRDShadowCluster m_Ycandi = barShowerYCol[js].getAllCandiCol()[0];
      bool f_matched = false; 

      for(int ic=0; ic<m_candiCol.size(); ic++)
        if(m_candiCol[ic]==m_Ycandi){ f_matched = true;  break; }

      if(f_matched){
          m_shCol.push_back(barShowerYCol[js]);
          barShowerYCol.erase(barShowerYCol.begin()+js); js--;
          //break;
      }
    } //End loop ShowerYCol. 

    if(m_shCol.size()!=0){
      std::pair< CRDEcalEDM::CRDCaloBarShower, std::vector<CRDEcalEDM::CRDCaloBarShower> > m_showerpair( barShowerXCol[is], m_shCol );
      m_showerpairCol.push_back(m_showerpair);
    }
  }

//cout<<"  GetMatchedShowersL1: pair number: "<<m_showerpairCol.size()<<endl;

  for(int ip=0; ip<m_showerpairCol.size(); ip++){
    std::vector<CRDEcalEDM::CRDCaloHitTransShower> tmp_showerCol; tmp_showerCol.clear();
    XYShowerMatchingL1(m_showerpairCol[ip].first, m_showerpairCol[ip].second, tmp_showerCol);

    for(int is=0; is<tmp_showerCol.size(); is++ ) {
      tmp_showerCol[is].setShadowClus( m_showerpairCol[ip].second[is].getAllCandiCol()[0] );
      tmp_showerCol[is].setShadowClusType( m_showerpairCol[ip].second[is].getAllCandiCol()[0].Type );
    }

    outshCol.insert( outshCol.end(), tmp_showerCol.begin(), tmp_showerCol.end() );
  }

  return StatusCode::SUCCESS;
}


StatusCode EnergyTimeMatchingAlg::GetMatchedShowersL2(
           std::vector<CRDEcalEDM::CRDCaloBarShower>& barShowerXCol,
           std::vector<CRDEcalEDM::CRDCaloBarShower>& barShowerYCol,
           std::vector<CRDEcalEDM::CRDCaloHitTransShower>& outshCol )
{
  outshCol.clear();

  //This algorithm is WIP. Now return all combinations. 
  GetFullMatchedShowers( barShowerXCol, barShowerYCol, outshCol );

  return StatusCode::SUCCESS;
}



StatusCode EnergyTimeMatchingAlg::XYShowerMatchingL0( 
             CRDEcalEDM::CRDCaloBarShower& barShowerX, 
             CRDEcalEDM::CRDCaloBarShower& barShowerY, 
             CRDEcalEDM::CRDCaloHitTransShower& outsh)
{
  CRDEcalEDM::CRDCaloHitTransShower m_2dshower; m_2dshower.Clear();

  std::vector<edm4hep::ConstCalorimeterHit> m_digiCol; m_digiCol.clear();
  int NbarsX = barShowerX.getBars().size();
  int NbarsY = barShowerY.getBars().size();
  if(NbarsX==0 || NbarsY==0){ std::cout<<"WARNING: empty DigiHitsCol returned!"<<std::endl; outsh = m_2dshower;  return StatusCode::SUCCESS; }

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
  outsh = m_2dshower;  

  return StatusCode::SUCCESS;
}


StatusCode EnergyTimeMatchingAlg::XYShowerMatchingL1(
             CRDEcalEDM::CRDCaloBarShower& shower1,
             std::vector<CRDEcalEDM::CRDCaloBarShower>& showerNCol,
             std::vector<CRDEcalEDM::CRDCaloHitTransShower>& outshCol)
{
  outshCol.clear(); 

  int _slayer = shower1.getBars()[0].getSlayer(); 

  const int NshY = showerNCol.size();
  double totE_shY = 0;
  double EshY[NshY] = {0};
  for(int is=0;is<NshY;is++){ EshY[is] = showerNCol[is].getE(); totE_shY += EshY[is]; }
  for(int is=0;is<NshY;is++){
    double wi_E = EshY[is]/totE_shY;
    CRDEcalEDM::CRDCaloBarShower m_splitshower1;

    CRDEcalEDM::CRDCaloBar m_wiseed = shower1.getSeed();
    m_wiseed.setQ( wi_E*m_wiseed.getQ1(), wi_E*m_wiseed.getQ2() );

    std::vector<CRDEcalEDM::CRDCaloBar> m_wibars;
    for(int ib=0;ib<shower1.getBars().size();ib++){
      CRDEcalEDM::CRDCaloBar m_wibar = shower1.getBars()[ib];
      m_wibar.setQ(wi_E*m_wibar.getQ1(), wi_E*m_wibar.getQ2());
      m_wibars.push_back(m_wibar);
    }

    m_splitshower1.setBars( m_wibars );
    m_splitshower1.setSeed( m_wiseed );

    CRDEcalEDM::CRDCaloHitTransShower m_shower; m_shower.Clear();
    if(_slayer==0 ) XYShowerMatchingL0( m_splitshower1, showerNCol[is], m_shower);
    else            XYShowerMatchingL0( showerNCol[is], m_splitshower1, m_shower);
    outshCol.push_back( m_shower );
  }

  return StatusCode::SUCCESS;
}


StatusCode EnergyTimeMatchingAlg::XYShowerChi2Matching(
           std::vector<CRDEcalEDM::CRDCaloBarShower>& barShowerXCol,
           std::vector<CRDEcalEDM::CRDCaloBarShower>& barShowerYCol,
           std::vector<CRDEcalEDM::CRDCaloHitTransShower>& outshCol )
{
  outshCol.clear(); 
  if(barShowerXCol.size() != barShowerYCol.size() ) return StatusCode::FAILURE;

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

    CRDEcalEDM::CRDCaloHitTransShower tmp_shower; tmp_shower.Clear();
    XYShowerMatchingL0(showerX, showerY, tmp_shower);
    outshCol.push_back(tmp_shower);
  }

  return StatusCode::SUCCESS;
}


StatusCode EnergyTimeMatchingAlg::XYShowerChi2MatchingL1(
           std::vector<CRDEcalEDM::CRDCaloBarShower>& barShowerXCol,
           std::vector<CRDEcalEDM::CRDCaloBarShower>& barShowerYCol,
           std::vector<CRDEcalEDM::CRDCaloHitTransShower>& outshCol )
{
  outshCol.clear();

  const int NshowerX = barShowerXCol.size();
  const int NshowerY = barShowerYCol.size();

  double chi2[NshowerX][NshowerY];
  double chi2_E[NshowerX][NshowerY];
  double chi2_tx[NshowerX][NshowerY];
  double chi2_ty[NshowerX][NshowerY];  

  double wi_E = settings.chi2Wi_E/(settings.chi2Wi_E + settings.chi2Wi_T);
  double wi_T = settings.chi2Wi_T/(settings.chi2Wi_E + settings.chi2Wi_T);

  TVector3 m_vec(0,0,0);
  double rotAngle = -(barShowerXCol[0].getBars())[0].getModule()*PI/4.;
  TVector3 Cblock((barShowerXCol[0].getBars())[0].getPosition().x(), (barShowerXCol[0].getBars())[0].getPosition().y(), (barShowerYCol[0].getBars())[0].getPosition().z());
  Cblock.RotateZ(rotAngle);

  map<double, pair<int, int> > m_chi2Map; m_chi2Map.clear(); 

  for(int ix=0;ix<NshowerX;ix++){
  for(int iy=0;iy<NshowerY;iy++){
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

    pair<int, int> p1(ix, iy);
    m_chi2Map[chi2[ix][iy]] = p1;
  }}

  pair<int, int> lastpair; 
  vector<pair<int, int>> indexVec; indexVec.clear(); 
  map<double, pair<int, int> >::iterator iter = m_chi2Map.begin(); 

  for(iter; iter!=m_chi2Map.end(); iter++){
    pair<int, int> indexpair = iter->second; 
    bool inLine = false; 
    bool isLast = false; 
    for(int i=0; i<indexVec.size(); i++){
      if( indexVec.size() == min(NshowerX, NshowerY)-1 ) {lastpair = indexpair; isLast=true;  break; }
      if( indexpair.first == indexVec[i].first || indexpair.second == indexVec[i].second ) { inLine=true; break; }
    }
    if(isLast) break;
    if(inLine) continue; 
    indexVec.push_back(indexpair);
  }
  if( indexVec.size()!=min(NshowerX, NshowerY)-1 ) 
    cout<<"ERROR in XYShowerChi2MatchingL1: found pair size "<<indexVec.size()<<" does not equal to min shower size -1 "<<min(NshowerX, NshowerY)-1<<endl;


  vector<CRDEcalEDM::CRDCaloBarShower> leftShowers; leftShowers.clear(); 
  for(int i=0; i<max(NshowerX, NshowerY); i++){
    bool fl_exist = false; 
    for(int j=0; j<indexVec.size(); j++){
      int m_index = NshowerX>NshowerY ? indexVec[j].first : indexVec[j].second;
      if(i==m_index){ fl_exist = true; break; } 
    }
    if(!fl_exist){
       CRDEcalEDM::CRDCaloBarShower m_shower = NshowerX>NshowerY ? barShowerXCol[i] : barShowerYCol[i] ;
       leftShowers.push_back(m_shower);
    }
  }
  if(leftShowers.size() != fabs( NshowerX-NshowerY )+1 ) 
    cout<<"ERROR in XYShowerChi2MatchingL1: Last pair number "<<leftShowers.size()<<" does not equal to shower difference "<<fabs( NshowerX-NshowerY )+1<<endl;


  for(int ip=0; ip<indexVec.size(); ip++){
    CRDEcalEDM::CRDCaloBarShower showerX = barShowerXCol[indexVec[ip].first];
    CRDEcalEDM::CRDCaloBarShower showerY = barShowerYCol[indexVec[ip].second];

    CRDEcalEDM::CRDCaloHitTransShower tmp_shower; tmp_shower.Clear();
    XYShowerMatchingL0(showerX, showerY, tmp_shower);
    outshCol.push_back(tmp_shower);
  }


  std::vector<CRDEcalEDM::CRDCaloHitTransShower> m_showerinlayer; m_showerinlayer.clear();   
  int ilast = NshowerX<NshowerY ? lastpair.first : lastpair.second; 
  CRDEcalEDM::CRDCaloBarShower m_shower = NshowerX<NshowerY ? barShowerXCol[ilast] : barShowerYCol[ilast] ;
  XYShowerMatchingL1( m_shower, leftShowers, m_showerinlayer);

  outshCol.insert(outshCol.end(), m_showerinlayer.begin(), m_showerinlayer.end());

  return StatusCode::SUCCESS;

}
#endif
