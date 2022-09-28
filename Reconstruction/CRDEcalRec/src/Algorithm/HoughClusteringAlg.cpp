#ifndef HOUGHCLUSTERINGALG_C
#define HOUGHCLUSTERINGALG_C

#include "Algorithm/HoughClusteringAlg.h"
#include <algorithm>
#include "TCanvas.h"
using namespace std;
using namespace TMath;

StatusCode HoughClusteringAlg::ReadSettings(PandoraPlus::Settings& m_settings){
  settings = m_settings;

  if(settings.map_floatPars.find("th_Layers")==settings.map_floatPars.end()) settings.map_floatPars["th_Layers"] = 10;  
  if(settings.map_floatPars.find("th_peak")==settings.map_floatPars.end()) settings.map_floatPars["th_peak"] = 3;  
  if(settings.map_floatPars.find("Nbins_alpha")==settings.map_floatPars.end()) settings.map_floatPars["Nbins_alpha"] = 800;  
  if(settings.map_floatPars.find("Nbins_rho")==settings.map_floatPars.end()) settings.map_floatPars["Nbins_rho"] = 800;  
  if(settings.map_floatPars.find("Range_rho")==settings.map_floatPars.end()) settings.map_floatPars["Range_rho"] = 130.; //[-range, range] 
  if(settings.map_floatPars.find("HoughBinDivide")==settings.map_floatPars.end()) settings.map_floatPars["HoughBinDivide"] = 4;  
  if(settings.map_floatPars.find("th_continuetrkN")==settings.map_floatPars.end()) settings.map_floatPars["th_continuetrkN"] = 3;  
  if(settings.map_boolPars.find("passIP")==settings.map_boolPars.end()) settings.map_boolPars["passIP"] = true;
  if(settings.map_floatPars.find("IPBandSmear")==settings.map_floatPars.end()) settings.map_floatPars["IPBandSmear"] = 200.;
  

  if(settings.map_floatPars.find("th_breakLaer")==settings.map_floatPars.end()) settings.map_floatPars["th_breakLaer"] = 2;  
  if(settings.map_floatPars.find("th_AxisE")==settings.map_floatPars.end()) settings.map_floatPars["th_AxisE"] = 1.;  
  if(settings.map_floatPars.find("th_intercept")==settings.map_floatPars.end()) settings.map_floatPars["th_intercept"] = 10e7;  
  if(settings.map_floatPars.find("th_dAlpha1")==settings.map_floatPars.end()) settings.map_floatPars["th_dAlpha1"] = 0.3;  
  if(settings.map_floatPars.find("th_dAlpha2")==settings.map_floatPars.end()) settings.map_floatPars["th_dAlpha2"] = 0.8;  
  if(settings.map_floatPars.find("th_dRho1")==settings.map_floatPars.end()) settings.map_floatPars["th_dRho1"] = 20;  
  if(settings.map_floatPars.find("th_dRho2")==settings.map_floatPars.end()) settings.map_floatPars["th_dRho2"] = 50;  
  if(settings.map_floatPars.find("th_overlapE1")==settings.map_floatPars.end()) settings.map_floatPars["th_overlapE1"] = 0.7;  
  if(settings.map_floatPars.find("th_overlapE2")==settings.map_floatPars.end()) settings.map_floatPars["th_overlapE2"] = 0.9;  

  if(settings.map_stringPars.find("ReadinLocalMaxName")==settings.map_stringPars.end())  settings.map_stringPars["ReadinLocalMaxName"] = "AllLocalMax";
  if(settings.map_stringPars.find("LeftLocalMaxName")==settings.map_stringPars.end())    settings.map_stringPars["LeftLocalMaxName"] = "LeftLocalMax";
  if(settings.map_stringPars.find("OutputLongiClusName")==settings.map_stringPars.end()) settings.map_stringPars["OutputLongiClusName"] = "EMLongiCluster"; 

  return StatusCode::SUCCESS;
};

StatusCode HoughClusteringAlg::Initialize(){

  return StatusCode::SUCCESS;
};

StatusCode HoughClusteringAlg::RunAlgorithm( PandoraPlusDataCol& m_datacol ){

  std::vector<PandoraPlus::Calo3DCluster*>* p_3DClusters = &(m_datacol.Cluster3DCol); 
  if(!p_3DClusters){ std::cout<<"ERROR: No 3DCluster in present data collection! "<<std::endl; return StatusCode::FAILURE; }
cout<<"3DCluster size: "<<p_3DClusters->size()<<endl;

  for(int it=0; it<p_3DClusters->size(); it++){
printf("  In 3DCluster #%d: block size = %d \n", it, p_3DClusters->at(it)->getCluster().size());

    std::vector<const PandoraPlus::Calo1DCluster*> m_localMaxUCol; m_localMaxUCol.clear();  
    std::vector<const PandoraPlus::Calo1DCluster*> m_localMaxVCol; m_localMaxVCol.clear();
    std::vector<const PandoraPlus::Calo1DCluster*> tmp_localMaxUCol = p_3DClusters->at(it)->getLocalMaxUCol(settings.map_stringPars["ReadinLocalMaxName"]); 
    std::vector<const PandoraPlus::Calo1DCluster*> tmp_localMaxVCol = p_3DClusters->at(it)->getLocalMaxVCol(settings.map_stringPars["ReadinLocalMaxName"]);

    for(int il=0; il<tmp_localMaxUCol.size(); il++)
      if(tmp_localMaxUCol[il]->getDlayer()<=settings.map_floatPars["th_Layers"]) m_localMaxUCol.push_back(tmp_localMaxUCol[il]);
    for(int il=0; il<tmp_localMaxVCol.size(); il++)
      if(tmp_localMaxVCol[il]->getDlayer()<=settings.map_floatPars["th_Layers"]) m_localMaxVCol.push_back(tmp_localMaxVCol[il]);


    if(m_localMaxUCol.size()==0 && m_localMaxVCol.size()==0) continue; 
cout<<"  Local maximum size: barX "<<m_localMaxUCol.size()<<"  barY "<<m_localMaxVCol.size()<<endl;

    std::vector<PandoraPlus::HoughObject> m_HoughObjectsU; m_HoughObjectsU.clear(); 
    std::vector<PandoraPlus::HoughObject> m_HoughObjectsV; m_HoughObjectsV.clear(); 
    for(int il=0; il<m_localMaxUCol.size(); il++){
      PandoraPlus::HoughObject m_obj; m_obj.Clear();
      m_obj.addLocalMax( m_localMaxUCol[il] );
      m_HoughObjectsU.push_back(m_obj);
    }
    for(int il=0; il<m_localMaxVCol.size(); il++){
      PandoraPlus::HoughObject m_obj; m_obj.Clear();
      m_obj.addLocalMax( m_localMaxVCol[il] );
      m_HoughObjectsV.push_back(m_obj);
    }

cout<<"  HoughClusteringAlg: Conformal transformation"<<endl;
    //Conformal transformation  
    ConformalTransformation(m_HoughObjectsU);
    ConformalTransformation(m_HoughObjectsV);

/*
cout<<"    Local max (HoughObjectX): "<<endl;
for(int i=0; i<m_HoughObjectsU.size(); i++) printf("(%.2f, %.2f, %.2f) \t", m_HoughObjectsU[i].getLocalMax()[0]->getPos().x(), m_HoughObjectsU[i].getLocalMax()[0]->getPos().y(), m_HoughObjectsU[i].getLocalMax()[0]->getPos().z() );
cout<<endl;
cout<<"    Conformal Point (HoughObjectX): "<<endl;
for(int i=0; i<m_HoughObjectsU.size(); i++) printf("(%.2f, %.2f) \t", m_HoughObjectsU[i].getConformPointUR().X(), m_HoughObjectsU[i].getConformPointUR().Y());
cout<<endl;

cout<<"    Local max (HoughObjectY): "<<endl;
for(int i=0; i<m_HoughObjectsV.size(); i++) printf("(%.2f, %.2f, %.2f) \t", m_HoughObjectsV[i].getLocalMax()[0]->getPos().x(), m_HoughObjectsV[i].getLocalMax()[0]->getPos().y(), m_HoughObjectsV[i].getLocalMax()[0]->getPos().z() );
cout<<endl;
cout<<"    Conformal Point (HoughObjectY): "<<endl;
for(int i=0; i<m_HoughObjectsV.size(); i++) printf("(%.2f, %.2f) \t", m_HoughObjectsV[i].getConformPointUR().X(), m_HoughObjectsV[i].getConformPointUR().Y());
cout<<endl;
*/


//cout<<"  HoughClusteringAlg: Hough transformation"<<endl;

    //Hough transformation
    HoughTransformation(m_HoughObjectsU);
    HoughTransformation(m_HoughObjectsV);

/*
cout<<"  Check all Hough Objects: "<<endl;

cout<<"  Transformed Hough lines X: "<<endl;
for(int i=0; i<m_HoughObjectsU.size(); i++){
  cout<<"    HoughObj #"<<i<<": Origin localMax position ";
  printf(" (%.3f, %.3f, %.3f) \n", m_HoughObjectsU[i].getLocalMax()->getPos().X(), m_HoughObjectsU[i].getLocalMax()->getPos().Z(), m_HoughObjectsU[i].getLocalMax()->getEnergy());
  printf("    Conformal Point: (%.3f, %.3f) \n", m_HoughObjectsU[i].getConformPointUR().X()-5, m_HoughObjectsU[i].getConformPointUR().Y()-5 );
  cout<<"    LineUR pars: "<<m_HoughObjectsU[i].getHoughLineUR().GetParameter(0)<<"  "<<m_HoughObjectsU[i].getHoughLineUR().GetParameter(1)<<", range: ["<<m_HoughObjectsU[i].getHoughLineUR().GetMaximum()<<", "<<m_HoughObjectsU[i].getHoughLineUR().GetMinimum()<<"]"<<endl;
  cout<<"    LineUL pars: "<<m_HoughObjectsU[i].getHoughLineUL().GetParameter(0)<<"  "<<m_HoughObjectsU[i].getHoughLineUL().GetParameter(1)<<", range: ["<<m_HoughObjectsU[i].getHoughLineUL().GetMaximum()<<", "<<m_HoughObjectsU[i].getHoughLineUL().GetMinimum()<<"]"<<endl;
  cout<<"    LineDR pars: "<<m_HoughObjectsU[i].getHoughLineDR().GetParameter(0)<<"  "<<m_HoughObjectsU[i].getHoughLineDR().GetParameter(1)<<", range: ["<<m_HoughObjectsU[i].getHoughLineDR().GetMaximum()<<", "<<m_HoughObjectsU[i].getHoughLineDR().GetMinimum()<<"]"<<endl;
  cout<<"    LineDL pars: "<<m_HoughObjectsU[i].getHoughLineDL().GetParameter(0)<<"  "<<m_HoughObjectsU[i].getHoughLineDL().GetParameter(1)<<", range: ["<<m_HoughObjectsU[i].getHoughLineDL().GetMaximum()<<", "<<m_HoughObjectsU[i].getHoughLineDL().GetMinimum()<<"]"<<endl;
  cout<<endl;
}


cout<<"  Transformed Hough lines Y: "<<endl;
for(int i=0; i<m_HoughObjectsV.size(); i++){
  cout<<"    HoughObj #"<<i<<": Origin localMax position ";
  printf(" (%.3f, %.3f, %.3f) \n", m_HoughObjectsV[i].getLocalMax()->getPos().X(), m_HoughObjectsV[i].getLocalMax()->getPos().Z(), m_HoughObjectsV[i].getLocalMax()->getEnergy());
  printf("    Conformal Point: (%.3f, %.3f) \n", m_HoughObjectsV[i].getConformPointUR().X()-5, m_HoughObjectsV[i].getConformPointUR().Y()-5 );
  cout<<"    LineUR pars: "<<m_HoughObjectsV[i].getHoughLineUR().GetParameter(0)<<"  "<<m_HoughObjectsV[i].getHoughLineUR().GetParameter(1)<<", range: ["<<m_HoughObjectsV[i].getHoughLineUR().GetMaximum()<<", "<<m_HoughObjectsV[i].getHoughLineUR().GetMinimum()<<"]"<<endl;
  cout<<"    LineUL pars: "<<m_HoughObjectsV[i].getHoughLineUL().GetParameter(0)<<"  "<<m_HoughObjectsV[i].getHoughLineUL().GetParameter(1)<<", range: ["<<m_HoughObjectsV[i].getHoughLineUL().GetMaximum()<<", "<<m_HoughObjectsV[i].getHoughLineUL().GetMinimum()<<"]"<<endl;
  cout<<"    LineDR pars: "<<m_HoughObjectsV[i].getHoughLineDR().GetParameter(0)<<"  "<<m_HoughObjectsV[i].getHoughLineDR().GetParameter(1)<<", range: ["<<m_HoughObjectsV[i].getHoughLineDR().GetMaximum()<<", "<<m_HoughObjectsV[i].getHoughLineDR().GetMinimum()<<"]"<<endl;
  cout<<"    LineDL pars: "<<m_HoughObjectsV[i].getHoughLineDL().GetParameter(0)<<"  "<<m_HoughObjectsV[i].getHoughLineDL().GetParameter(1)<<", range: ["<<m_HoughObjectsV[i].getHoughLineDL().GetMaximum()<<", "<<m_HoughObjectsV[i].getHoughLineDL().GetMinimum()<<"]"<<endl;
  cout<<endl;
}
*/

   
cout<<"  HoughClusteringAlg: Fill bins to get the hills"<<endl;
    //Fill bins to get the hills
    PandoraPlus::HoughSpace m_HoughSpaceU; 
    PandoraPlus::HoughSpace m_HoughSpaceV; 
    FillHoughSpace(m_HoughObjectsU, m_HoughSpaceU);
    FillHoughSpace(m_HoughObjectsV, m_HoughSpaceV);

/*
cout<<"  Save Hough space map"<<endl;
TCanvas *c1 = new TCanvas("c1", "c1", 700,500);
c1->cd();
TH2F *h2_mapXY = new TH2F();
*h2_mapXY = m_HoughSpaceV.getSpaceMap();
cout<<" If space map is empty: "<<h2_mapXY->Integral()<<endl;
h2_mapXY->Draw("colz");
//c1->Draw();
c1->SaveAs("/cefs/higgs/guofy/CEPCSW_v203/run/HoughMapXY.png");
c1->SaveAs("/cefs/higgs/guofy/CEPCSW_v203/run/HoughMapXY.C");
delete c1, h2_mapXY;
*/

cout<<"  HoughClusteringAlg: Find hills"<<endl;
    //Find hills
    FindingHills(m_HoughSpaceU);
    FindingHills(m_HoughSpaceV);

/*
cout<<"  Print Hough hills in HoughSpaceU: "<<endl;
cout<<"  Hill size: "<<m_HoughSpaceU.getHills().size()<<endl;
for(int ih=0; ih<m_HoughSpaceU.getHills().size(); ih++){
  int m_size = m_HoughSpaceU.getHills()[ih].getIndexAlpha().size(); 
  std::vector<int> vec_alpha = m_HoughSpaceU.getHills()[ih].getIndexAlpha();
  std::vector<int> vec_rho = m_HoughSpaceU.getHills()[ih].getIndexRho();

  printf("    In Hill #%d: cell size = %d, alpha range [%d, %d], rho range [%d, %d]. \n", 
          ih, m_size, *min_element(vec_alpha.begin(), vec_alpha.end()), *max_element(vec_alpha.begin(), vec_alpha.end()), 
          *min_element(vec_rho.begin(), vec_rho.end()), *max_element(vec_rho.begin(), vec_rho.end())  );
}
*/

cout<<"  HoughClusteringAlg: Create output HoughClusters. "<<endl;
    //Create output HoughClusters 
    std::vector<const PandoraPlus::LongiCluster*> m_longiClusUCol; m_longiClusUCol.clear();
    std::vector<const PandoraPlus::LongiCluster*> m_longiClusVCol; m_longiClusVCol.clear(); 
    Transform2Clusters(m_HoughSpaceU, m_HoughObjectsU, m_longiClusUCol);
    Transform2Clusters(m_HoughSpaceV, m_HoughObjectsV, m_longiClusVCol);

    for(int ic=0; ic<m_longiClusUCol.size(); ic++) m_datacol.bk_LongiClusCol.push_back( const_cast<PandoraPlus::LongiCluster*>(m_longiClusUCol[ic]) );
    for(int ic=0; ic<m_longiClusVCol.size(); ic++) m_datacol.bk_LongiClusCol.push_back( const_cast<PandoraPlus::LongiCluster*>(m_longiClusVCol[ic]) );

    std::vector<const PandoraPlus::Calo1DCluster*> left_localMaxUCol; left_localMaxUCol.clear(); 
    std::vector<const PandoraPlus::Calo1DCluster*> left_localMaxVCol; left_localMaxVCol.clear(); 


//cout<<"Total localMax U address: "<<endl;
//for(int i=0; i<tmp_localMaxUCol.size(); i++) 
//  printf("  #%d: (%.3f, %.3f, %.3f), %p \n", i, tmp_localMaxUCol[i]->getPos().x(), tmp_localMaxUCol[i]->getPos().y(), tmp_localMaxUCol[i]->getPos().z(), tmp_localMaxUCol[i]);
//cout<<endl;

    for(int is=0; is<tmp_localMaxUCol.size(); is++){
      bool fl_incluster = false; 
      for(int ic=0; ic<m_longiClusUCol.size(); ic++){
        std::vector<const PandoraPlus::Calo1DCluster*> p_showers = m_longiClusUCol[ic]->getBarShowers();
        if( find(p_showers.begin(), p_showers.end(), tmp_localMaxUCol[is])!=p_showers.end() ) { fl_incluster = true; break; }
      }
      if(!fl_incluster && find(left_localMaxUCol.begin(), left_localMaxUCol.end(), tmp_localMaxUCol[is])==left_localMaxUCol.end() ) left_localMaxUCol.push_back(tmp_localMaxUCol[is]);
    }
    for(int is=0; is<tmp_localMaxVCol.size(); is++){
      bool fl_incluster = false;
      for(int ic=0; ic<m_longiClusVCol.size(); ic++){
        std::vector<const PandoraPlus::Calo1DCluster*> p_showers = m_longiClusVCol[ic]->getBarShowers();
        if( find(p_showers.begin(), p_showers.end(), tmp_localMaxVCol[is])!=p_showers.end() ) { fl_incluster = true; break; }
      }
      if(!fl_incluster && find(left_localMaxVCol.begin(), left_localMaxVCol.end(), tmp_localMaxVCol[is])==left_localMaxVCol.end() ) left_localMaxVCol.push_back(tmp_localMaxVCol[is]);
    }


/*
    for(int ic=0; ic<m_longiClusUCol.size(); ic++){
      std::vector<const PandoraPlus::Calo1DCluster*> p_showers = m_longiClusUCol[ic]->getBarShowers();
cout<<"LocalMax in LongiCluster #"<<ic<<endl;
for(int i=0; i<p_showers.size(); i++)
  printf("  #%d: (%.3f, %.3f, %.3f), %p \n", i, p_showers[i]->getPos().x(), p_showers[i]->getPos().y(), p_showers[i]->getPos().z(), p_showers[i]);
cout<<endl;
      for(int is=0; is<tmp_localMaxUCol.size(); is++){
printf("  Looping for localMax #%d: (%.3f, %.3f, %.3f), %p \n", is, tmp_localMaxUCol[is]->getPos().x(), tmp_localMaxUCol[is]->getPos().y(), tmp_localMaxUCol[is]->getPos().z(), tmp_localMaxUCol[is]);
cout<<"  Current left collection size: "<<left_localMaxUCol.size()<<endl;

        if( find(p_showers.begin(), p_showers.end(), tmp_localMaxUCol[is])==p_showers.end() && 
            find(left_localMaxUCol.begin(), left_localMaxUCol.end(), tmp_localMaxUCol[is])==left_localMaxUCol.end() ){ 
cout<<"    Not found in LongiCluster or current left collection"<<endl;
          left_localMaxUCol.push_back(tmp_localMaxUCol[is]);
        }
      }
    }

    for(int ic=0; ic<m_longiClusVCol.size(); ic++){
      std::vector<const PandoraPlus::Calo1DCluster*> p_showers = m_longiClusVCol[ic]->getBarShowers();
      for(int is=0; is<tmp_localMaxVCol.size(); is++)
        if( find(p_showers.begin(), p_showers.end(), tmp_localMaxVCol[is])==p_showers.end() && 
            find(left_localMaxVCol.begin(), left_localMaxVCol.end(), tmp_localMaxVCol[is])==left_localMaxVCol.end() ) left_localMaxVCol.push_back(tmp_localMaxVCol[is]);
    }
*/

    p_3DClusters->at(it)->setLocalMax( settings.map_stringPars["LeftLocalMaxName"], left_localMaxUCol, 
                                       settings.map_stringPars["LeftLocalMaxName"], left_localMaxVCol  );

//cout<<"Left LocalMax U: "<<endl;
//for(int i=0; i<left_localMaxUCol.size(); i++)
//  printf("  #%d: (%.3f, %.3f, %.3f) \n ", i, left_localMaxUCol[i]->getPos().x(), left_localMaxUCol[i]->getPos().y(), left_localMaxUCol[i]->getPos().z());
//cout<<"Left LocalMax V: "<<endl;
//for(int i=0; i<left_localMaxVCol.size(); i++)
//  printf("  #%d: (%.3f, %.3f, %.3f) \n ", i, left_localMaxVCol[i]->getPos().x(), left_localMaxVCol[i]->getPos().y(), left_localMaxVCol[i]->getPos().z());


/*
cout<<"  HoughCluster size: "<<m_longiClusUCol.size()<<" / "<<m_longiClusVCol.size()<<endl;

cout<<"  Print LongiClusterX: size = "<<m_longiClusUCol.size()<<endl;
for(int il=0; il<m_longiClusUCol.size(); il++){
  printf("    Clus#%d: Hough Par = (%.3f, %.3f, %.3f), shower layers: \n", il, m_longiClusUCol[il]->getHoughAlpha(), m_longiClusUCol[il]->getHoughRho(), m_longiClusUCol[il]->getHoughIntercept());
  for(int is=0; is<m_longiClusUCol[il]->getBarShowers().size(); is++) 
    printf("      Dlayer = %d, Pos/E = (%.2f, %.2f, %.2f, %.3f) \n", m_longiClusUCol[il]->getBarShowers()[is]->getDlayer(), 
                                                                     m_longiClusUCol[il]->getBarShowers()[is]->getPos().x(),
                                                                     m_longiClusUCol[il]->getBarShowers()[is]->getPos().y(),
                                                                     m_longiClusUCol[il]->getBarShowers()[is]->getPos().z(),
                                                                     m_longiClusUCol[il]->getBarShowers()[is]->getEnergy() );
  cout<<endl;
}

cout<<"  Print LongiClusterY: size = "<<m_longiClusVCol.size()<<endl;
for(int il=0; il<m_longiClusVCol.size(); il++){
  printf("    Clus#%d: Hough Par = (%.3f, %.3f, %.3f), shower layers: \n", il, m_longiClusVCol[il]->getHoughAlpha(), m_longiClusVCol[il]->getHoughRho(), m_longiClusVCol[il]->getHoughIntercept());
  for(int is=0; is<m_longiClusVCol[il]->getBarShowers().size(); is++)
    printf("      Dlayer = %d, Pos/E = (%.2f, %.2f, %.2f, %.3f) \n", m_longiClusVCol[il]->getBarShowers()[is]->getDlayer(),
                                                                     m_longiClusVCol[il]->getBarShowers()[is]->getPos().x(),
                                                                     m_longiClusVCol[il]->getBarShowers()[is]->getPos().y(),
                                                                     m_longiClusVCol[il]->getBarShowers()[is]->getPos().z(),
                                                                     m_longiClusVCol[il]->getBarShowers()[is]->getEnergy() );
  cout<<endl;
}
*/
    p_3DClusters->at(it)->setLongiClusters( settings.map_stringPars["OutputLongiClusName"], m_longiClusUCol, 
                                            settings.map_stringPars["OutputLongiClusName"], m_longiClusVCol); 

  }//End loop 3D cluster

  //m_datacol.TowerCol = p_3DClusters;
//cout<<"End in HoughClusteringAlg"<<endl;
  p_3DClusters = nullptr;
  return StatusCode::SUCCESS;
};



StatusCode HoughClusteringAlg::ClearAlgorithm(){

  return StatusCode::SUCCESS;
};


StatusCode HoughClusteringAlg::ConformalTransformation(std::vector<PandoraPlus::HoughObject>& m_Hobjects){
  if(m_Hobjects.size()==0) return StatusCode::SUCCESS;

  int module = -1;
  for(int io=0; io<m_Hobjects.size(); io++) 
    if(m_Hobjects[io].getLocalMax().size()!=0) { module = m_Hobjects[io].getLocalMax()[0]->getTowerID()[0][0]; break; }
  double rotAngle = -module*PI/4.;

  for(int io=0; io<m_Hobjects.size(); io++){
    if(m_Hobjects[io].getLocalMax().size()==0) continue;
    PandoraPlus::Calo1DCluster m_localMax = *(m_Hobjects[io].getLocalMax()[0]);

//cout<<"  Object #"<<io<<": local max position: ";
//printf("(%.2f, %.2f, %.2f) \n", m_localMax.getPos().x(), m_localMax.getPos().y(), m_localMax.getPos().z());

    int slayer = m_localMax.getSlayer();
    bool f_isXclus = (slayer==0 ? true : false);
//cout<<"  Objcet #"<<io<<": slayer="<<slayer<<", module="<<module<<", isX="<<f_isXclus<<", rotAngle="<<rotAngle<<endl;

    TVector3 p_localMax(0., 0., 0.);
    p_localMax.SetXYZ(m_localMax.getPos().x(), m_localMax.getPos().y(), m_localMax.getPos().z());
    p_localMax.RotateZ(rotAngle);

//cout<<"  Objcet #"<<io<<": position after rotation: ";
//printf("(%.2f, %.2f, %.2f) \n", p_localMax.X(), p_localMax.Y(), p_localMax.Z());
//cout<<endl;

    TVector2 pos;
    pos.SetX(p_localMax.Y());
    if(f_isXclus){
      pos.SetY(p_localMax.Z());
    }
    else{
      pos.SetY(p_localMax.X());
    }

    //double mod2_u = pos_u.Mod2();
    //double mod2_d = pos_d.Mod2();
    //pos_u.SetX(2*pos_u.X()/mod2_u); pos_u.SetY(2*pos_u.Y()/mod2_u);
    //pos_d.SetX(2*pos_d.X()/mod2_d); pos_d.SetY(2*pos_d.Y()/mod2_d);

    m_Hobjects[io].setCellSize(10.);
    m_Hobjects[io].setSlayer(slayer);
    m_Hobjects[io].setConformalPoint(pos);
  }

//cout<<"  ConformalTransformation: Hough Object size "<<m_Hobjects.size()<<endl;

  for(int io=0; io<m_Hobjects.size() && m_Hobjects.size()>1; io++){
    for(int jo=io+1; jo<m_Hobjects.size(); jo++){
//printf("    Compare (%d, %d): point (%f, %f) vs (%f, %f) \n", io, jo, 
//m_Hobjects[io].getConformPointUR().X(), m_Hobjects[io].getConformPointUR().Y(), m_Hobjects[jo].getConformPointUR().X(), m_Hobjects[jo].getConformPointUR().Y());

      if( fabs(m_Hobjects[io].getConformPointUR().X()-m_Hobjects[jo].getConformPointUR().X())<10e-3 &&
          fabs(m_Hobjects[io].getConformPointUR().Y()-m_Hobjects[jo].getConformPointUR().Y())<10e-3 ){
//cout<<"    --These 2 points are the same! "<<endl;
        for(int il=0; il<m_Hobjects[jo].getLocalMax().size(); il++) m_Hobjects[io].addLocalMax( m_Hobjects[jo].getLocalMax()[il] );
        m_Hobjects.erase(m_Hobjects.begin()+jo);
        jo--;
        if(jo<io) jo=io; 
      }
  }}

  return StatusCode::SUCCESS;
}


StatusCode HoughClusteringAlg::HoughTransformation(std::vector<PandoraPlus::HoughObject>& m_Hobjects){

  if(m_Hobjects.size()==0) return StatusCode::SUCCESS;
  
	double m_minX = 9999;
	double m_maxX = -9999;
	double m_minY = 9999;
	double m_maxY = -9999;
	for(int iobj=0; iobj<m_Hobjects.size(); iobj++){
    if(m_Hobjects[iobj].getConformPointUR().X()<m_minX) m_minX = m_Hobjects[iobj].getConformPointUR().X();
    if(m_Hobjects[iobj].getConformPointUR().X()>m_maxX) m_maxX = m_Hobjects[iobj].getConformPointUR().X();
    if(m_Hobjects[iobj].getConformPointUR().Y()<m_minY) m_minY = m_Hobjects[iobj].getConformPointUR().Y();
    if(m_Hobjects[iobj].getConformPointUR().Y()>m_maxY) m_maxY = m_Hobjects[iobj].getConformPointUR().Y();
	}
  double centerX = (m_minX+m_maxX)/2.; 
	double centerY = (m_minY+m_maxY)/2.; 

cout<<"Origin point in conformal space: (X, Y) = ("<<centerX<<", "<<centerY<<")"<<endl;

  for(int iobj=0; iobj<m_Hobjects.size(); iobj++){
    //UR
    TF1 line1("line_ur", "[0]*cos(x)+[1]*sin(x)", -0.1, PI/2);  
    line1.SetParameters(m_Hobjects[iobj].getConformPointUR().X()-centerX, m_Hobjects[iobj].getConformPointUR().Y()-centerY);

    //UL
    TF1 line2("line_ul", "[0]*cos(x)+[1]*sin(x)", PI/2, PI);
    line2.SetParameters(m_Hobjects[iobj].getConformPointUL().X()-centerX, m_Hobjects[iobj].getConformPointUL().Y()-centerY);

    //DR
    TF1 line3("line_dr", "[0]*cos(x)+[1]*sin(x)", PI/2, PI);
    line3.SetParameters(m_Hobjects[iobj].getConformPointDR().X()-centerX, m_Hobjects[iobj].getConformPointDR().Y()-centerY);

    //DL   
    TF1 line4("line_dl", "[0]*cos(x)+[1]*sin(x)", -0.1, PI/2);  
    line4.SetParameters(m_Hobjects[iobj].getConformPointDL().X()-centerX, m_Hobjects[iobj].getConformPointDL().Y()-centerY);

		m_Hobjects[iobj].setHoughLine(line1, line2, line3, line4);
  }

  return StatusCode::SUCCESS;
}


StatusCode HoughClusteringAlg::FillHoughSpace(std::vector<PandoraPlus::HoughObject>& m_Hobjects, PandoraPlus::HoughSpace& m_Hspace){
  if(m_Hobjects.size()==0) return StatusCode::SUCCESS;

  m_Hspace.Clear(); 

  //Get central conform point. 
  double m_minX = 9999;
  double m_maxX = -9999;
  double m_minY = 9999;
  double m_maxY = -9999;
  for(int iobj=0; iobj<m_Hobjects.size(); iobj++){
    if(m_Hobjects[iobj].getConformPointUR().X()<m_minX) m_minX = m_Hobjects[iobj].getConformPointUR().X();
    if(m_Hobjects[iobj].getConformPointUR().X()>m_maxX) m_maxX = m_Hobjects[iobj].getConformPointUR().X();
    if(m_Hobjects[iobj].getConformPointUR().Y()<m_minY) m_minY = m_Hobjects[iobj].getConformPointUR().Y();
    if(m_Hobjects[iobj].getConformPointUR().Y()>m_maxY) m_maxY = m_Hobjects[iobj].getConformPointUR().Y();
  }
  double centerX = (m_minX+m_maxX)/2.;
  double centerY = (m_minY+m_maxY)/2.;


  //Create Hough space
  TH2F m_houghMap("","",(int)settings.map_floatPars["Nbins_alpha"], -0.1, PI, (int)settings.map_floatPars["Nbins_rho"], -settings.map_floatPars["Range_rho"], settings.map_floatPars["Range_rho"]);
  double width_alpha = (PI+0.1)/(double)settings.map_floatPars["Nbins_alpha"];
  double width_rho = 2.*settings.map_floatPars["Range_rho"]/(double)settings.map_floatPars["Nbins_rho"]; 

  //Fill Hough space
  for(int io=0; io<m_Hobjects.size(); io++){
//cout<<"  In HoughObj #"<<io<<endl;

    TF1 line_ur = m_Hobjects[io].getHoughLineUR();
    TF1 line_ul = m_Hobjects[io].getHoughLineUL();
    TF1 line_dr = m_Hobjects[io].getHoughLineDR();
    TF1 line_dl = m_Hobjects[io].getHoughLineDL();


    //Loop for alpha bins
		for(int ibin=0; ibin<m_houghMap.GetNbinsX(); ibin++){
      double m_alphaL = ibin*width_alpha-0.1; 
      if(line_ur.Eval(m_alphaL)>settings.map_floatPars["Range_rho"]-2*width_rho
         || line_ur.Eval(m_alphaL)<-settings.map_floatPars["Range_rho"]+2*width_rho
         || line_dl.Eval(m_alphaL)>settings.map_floatPars["Range_rho"]-2*width_rho
         || line_dl.Eval(m_alphaL)<-settings.map_floatPars["Range_rho"]+2*width_rho) continue;


      int divide = settings.map_floatPars["HoughBinDivide"];
      int nbin_rho_u, nbin_rho_d; 
      double rho_min = 9999.;
      double rho_max = -9999.; 

      for(int id=0; id<=divide; id++){

        double alpha = ibin*width_alpha-0.1 + (double)id*width_alpha/divide; 

        double rho_u = alpha<((PI+0.1)/2.-0.1) ? line_ur.Eval(alpha) : line_ul.Eval(alpha);
        double rho_d = alpha<((PI+0.1)/2.-0.1) ? line_dl.Eval(alpha) : line_dr.Eval(alpha);
        //double rho_u_p = (alpha+width_alpha)<((PI+0.1)/2.-0.1) ? line_ur.Eval(alpha+width_alpha) : line_ul.Eval(alpha+width_alpha);
        //double rho_d_p = (alpha+width_alpha)<((PI+0.1)/2.-0.1) ? line_dl.Eval(alpha+width_alpha) : line_dr.Eval(alpha+width_alpha);

        //if(rho_min>min(rho_u, rho_d)) rho_min = min(rho_u, rho_d); 
        //if(rho_max<max(rho_u, rho_d)) rho_max = max(rho_u, rho_d);
        if(rho_min>rho_u) rho_min = rho_u;
        if(rho_max<rho_u) rho_max = rho_u;
        if(rho_min>rho_d) rho_min = rho_d;
        if(rho_max<rho_d) rho_max = rho_d;
//if(io==0) printf("      In divied %d: rho = %.2f(%.2f), min_rho = %.2f, max_rho = %.2f \n", id, rho_u,rho_d,rho_min,rho_max);
      }
      nbin_rho_d = floor( (rho_min+settings.map_floatPars["Range_rho"])/width_rho );
      nbin_rho_u = floor( (rho_max+settings.map_floatPars["Range_rho"])/width_rho );

//if(io==0) printf("    Loop in bin #%d: alpha range [%.3f, %.3f], rho bin range [%d, %d] \n",
//                 ibin, ibin*width_alpha-0.1, (ibin+1)*width_alpha-0.1, nbin_rho_d, nbin_rho_u); 

      for(int irho=nbin_rho_d; irho<=nbin_rho_u; irho++){
        //double tmp_alpha = ibin*width_alpha-0.1+0.001; 
        //double tmp_rho = irho*width_rho-130.; 
        //m_houghMap.Fill(tmp_alpha, tmp_rho);
        m_houghMap.SetBinContent(ibin+1, irho+1, m_houghMap.GetBinContent(ibin+1, irho+1)+1 );
//if(io==0) printf("      Fill HoughMap bin: (%d, %d, %.1f) \n", ibin, irho, m_houghMap.GetBinContent(ibin+1, irho+1) );
      }

    }
	}//End loop Hough Objects

  //Create IP Hough band

  TF1 m_IPline_outer("line_IPo", "[0]*cos(x)+[1]*sin(x)", -0.1, PI);
  TF1 m_IPline_inner("line_IPi", "[0]*cos(x)+[1]*sin(x)", -0.1, PI);
  m_IPline_outer.SetParameters(-centerX, -centerY+settings.map_floatPars["IPBandSmear"]);
  m_IPline_inner.SetParameters(-centerX, -centerY-settings.map_floatPars["IPBandSmear"]);


  m_Hspace.setSpaceMap(m_houghMap);
  m_Hspace.setIPBand( m_IPline_outer, m_IPline_inner ); //central, outer, inner. 

  return StatusCode::SUCCESS;
}


StatusCode HoughClusteringAlg::FindingHills(PandoraPlus::HoughSpace& m_Hspace){

  std::vector<PandoraPlus::HoughSpace::HoughHill> m_hillCol; m_hillCol.clear(); 
  TH2F m_houghMap = m_Hspace.getSpaceMap();
  
  for(int ia=0; ia<settings.map_floatPars["Nbins_alpha"]; ia++){
  for(int ir=0; ir<settings.map_floatPars["Nbins_alpha"]; ir++){
/*
cout<<"  Looping for cell: "<<ia<<", "<<ir<<endl;
cout<<"  cell and neighbor: "<<endl;
if (!(ia==0 || ia==settings.Nbins_alpha-1 || ir==0 || ir==settings.Nbins_rho-1)){
cout<<m_houghMap.GetBinContent(ia, ir)<<'\t'<<m_houghMap.GetBinContent(ia+1, ir)<<'\t'<<m_houghMap.GetBinContent(ia+2, ir)<<endl;
cout<<m_houghMap.GetBinContent(ia, ir+1)<<'\t'<<m_houghMap.GetBinContent(ia+1, ir+1)<<'\t'<<m_houghMap.GetBinContent(ia+2, ir+1)<<endl;
cout<<m_houghMap.GetBinContent(ia, ir+2)<<'\t'<<m_houghMap.GetBinContent(ia+1, ir+2)<<'\t'<<m_houghMap.GetBinContent(ia+2, ir+2)<<endl;
}
*/
    if(m_houghMap.GetBinContent(ia+1, ir+1)>=settings.map_floatPars["th_peak"]){
      PandoraPlus::HoughSpace::HoughHill m_hill; m_hill.Clear();
      m_hill.AddCell(ia, ir, m_houghMap.GetBinContent(ia+1, ir+1));
      int flag_expend = ExpandingPeak( m_houghMap, ia, ir, m_hill );
//cout<<"  expand flag: "<<flag_expend<<endl;
      if(flag_expend==1) continue;
      //else m_hillCol.push_back(m_hill);
      else if( m_Hspace.isFromIP(m_hill) ) { m_hillCol.push_back(m_hill); }
    }
  }}

  m_Hspace.setHoughHills( m_hillCol );
  return StatusCode::SUCCESS;
}

int HoughClusteringAlg::ExpandingPeak(TH2 &houghMap, int index_a, int index_b, PandoraPlus::HoughSpace::HoughHill& hill){

  if(index_a<0 || index_b<0 || index_a>=settings.map_floatPars["Nbins_alpha"] || index_b>=settings.map_floatPars["Nbins_alpha"]) return 0;
  houghMap.SetBinContent(index_a+1, index_b+1, -1.*houghMap.GetBinContent(index_a+1, index_b+1) );
  int count = 0;
  for(int fl1=-1; fl1<=1; fl1++){
  for(int fl2=-1; fl2<=1; fl2++){
    if(fl1!=0 || fl2!=0){
      if( Abs(houghMap.GetBinContent(index_a+fl1+1, index_b+fl2+1)) > Abs( houghMap.GetBinContent(index_a+1, index_b+1))){
        houghMap.SetBinContent(index_a+1, index_b+1, -1.*houghMap.GetBinContent(index_a+1, index_b+1) );
        return 1;
      }
      if(houghMap.GetBinContent(index_a+fl1+1, index_b+fl2+1) == Abs(houghMap.GetBinContent(index_a+1, index_b+1)) ){
        hill.AddCell(index_a+fl1, index_b+fl2, houghMap.GetBinContent(index_a+fl1+1, index_b+fl2+1) );
        count += ExpandingPeak(houghMap, index_a+fl1, index_b+fl2, hill );
        if(count>0){ 
          houghMap.SetBinContent(index_a+1, index_b+1, -1.*houghMap.GetBinContent(index_a+1, index_b+1) );
          return 1;
        }
      }
    }
  }}
  return 0;
}

StatusCode HoughClusteringAlg::Transform2Clusters( PandoraPlus::HoughSpace& m_Hspace, 
                                                   std::vector<PandoraPlus::HoughObject>& m_Hobjects,
                                                   std::vector<const PandoraPlus::LongiCluster*>& m_longiClusCol )
{
  if(m_Hobjects.size()==0) return StatusCode::SUCCESS;

/*
cout<<"    Print Local max (HoughObject): "<<endl;
for(int i=0; i<m_Hobjects.size(); i++){
  for(int il=0; il<m_Hobjects[i].getLocalMax().size(); il++){
    printf(" (%.2f, %.2f, %.2f), %p \t", m_Hobjects[i].getLocalMax()[il]->getPos().x(), m_Hobjects[i].getLocalMax()[il]->getPos().y(), m_Hobjects[i].getLocalMax()[il]->getPos().z(), m_Hobjects[i].getLocalMax()[il] );
  }
  cout<<endl;
}
cout<<endl;
*/

  //Get center of hits, for calculating intercept. 
  double m_minX = 9999;
  double m_maxX = -9999;
  double m_minY = 9999;
  double m_maxY = -9999;
  for(int iobj=0; iobj<m_Hobjects.size(); iobj++){
    if(m_Hobjects[iobj].getConformPointUR().X()<m_minX) m_minX = m_Hobjects[iobj].getConformPointUR().X();
    if(m_Hobjects[iobj].getConformPointUR().X()>m_maxX) m_maxX = m_Hobjects[iobj].getConformPointUR().X();
    if(m_Hobjects[iobj].getConformPointUR().Y()<m_minY) m_minY = m_Hobjects[iobj].getConformPointUR().Y();
    if(m_Hobjects[iobj].getConformPointUR().Y()>m_maxY) m_maxY = m_Hobjects[iobj].getConformPointUR().Y();
  }
  double centerX = (m_minX+m_maxX)/2.;
  double centerY = (m_minY+m_maxY)/2.;


  std::vector<PandoraPlus::LongiCluster*> m_clusCol; m_clusCol.clear(); 
  std::vector<PandoraPlus::HoughSpace::HoughHill> m_hills = m_Hspace.getHills(); 
  double alpha_min = m_Hspace.getAlphaLowEdge();
  double alpha_max = m_Hspace.getAlphaUpEdge();
  double rho_min = m_Hspace.getRhoLowEdge();
  double rho_max = m_Hspace.getRhoUpEdge();
  double alpha_width = m_Hspace.getAlphaBinWidth();
  double rho_width = m_Hspace.getRhoBinWidth();


std::vector<const PandoraPlus::LongiCluster*> tmp_clusCol; tmp_clusCol.clear(); 
cout<<"  Transform2Clusters: input hill size: "<<m_hills.size()<<endl;

  for(int ih=0; ih<m_hills.size(); ih++){
    PandoraPlus::LongiCluster* m_clus = new PandoraPlus::LongiCluster();

    std::vector<int> index_alpha = m_hills[ih].getIndexAlpha();
    std::vector<int> index_rho = m_hills[ih].getIndexRho();
    double sum_alpha=0; double sum_rho=0; 
    for(int i=0; i<index_alpha.size(); i++) sum_alpha += index_alpha[i];
    for(int i=0; i<index_rho.size(); i++) sum_rho += index_rho[i];

    double ave_alpha = alpha_min + (sum_alpha/index_alpha.size()+0.5)*alpha_width;
    double ave_rho = rho_min + (sum_rho/index_rho.size()+0.5)*rho_width;

//cout<<"    Hill #"<<ih<<": alpha = "<<ave_alpha<<", rho = "<<ave_rho<<endl;

    if(ave_alpha>alpha_min-alpha_width && ave_alpha <= alpha_min+(alpha_max-alpha_min)/2.) {
    for(int io=0; io<m_Hobjects.size(); io++){
        if( ( m_Hobjects[io].getHoughLineDL().Eval(ave_alpha)<ave_rho+rho_width &&  
              m_Hobjects[io].getHoughLineUR().Eval(ave_alpha)>ave_rho-rho_width) ||
            ( m_Hobjects[io].getHoughLineDL().Eval(ave_alpha)>ave_rho-rho_width &&
              m_Hobjects[io].getHoughLineUR().Eval(ave_alpha)<ave_rho+rho_width) )
          { 
          for(int il=0; il<m_Hobjects[io].getLocalMax().size(); il++) m_clus->addBarShower( m_Hobjects[io].getLocalMax()[il], 0 ); 
          m_clus->setHoughPars(ave_alpha, ave_rho ); 
          }
    }}

    else if(ave_alpha>alpha_min+(alpha_max-alpha_min)/2. && ave_alpha<alpha_max+alpha_width){
    for(int io=0; io<m_Hobjects.size(); io++){
        if( ( m_Hobjects[io].getHoughLineDR().Eval(ave_alpha)<ave_rho+rho_width &&
              m_Hobjects[io].getHoughLineUL().Eval(ave_alpha)>ave_rho-rho_width) ||
            ( m_Hobjects[io].getHoughLineDR().Eval(ave_alpha)>ave_rho-rho_width &&
              m_Hobjects[io].getHoughLineUL().Eval(ave_alpha)<ave_rho+rho_width) )

          {   
          for(int il=0; il<m_Hobjects[io].getLocalMax().size(); il++) m_clus->addBarShower( m_Hobjects[io].getLocalMax()[il], 0 ); 
          m_clus->setHoughPars(ave_alpha, ave_rho );
          }
    }}
tmp_clusCol.push_back(m_clus);

    //if(m_clus.getBarShowers().size()<settings.th_peak) continue;
    //if(settings.fl_continuetrk && !m_clus.isContinue()) continue;
    if(!m_clus->isContinueN(settings.map_floatPars["th_continuetrkN"])) continue; 

    double sin_alpha = sin(m_clus->getHoughAlpha()); 
    double cos_alpha = cos(m_clus->getHoughAlpha());
    double intercept = centerY + cos_alpha*centerX/sin_alpha + m_clus->getHoughRho()/sin_alpha; 
    m_clus->setIntercept(intercept); 
    //m_clus.FitAxis();
    m_clusCol.push_back(m_clus);

  }

cout<<"  Transform2Cluster: track removal cutflow: "<<endl;
cout<<"    Initial: "<<tmp_clusCol.size()<<endl;
cout<<"    After >= 3 hits "<<m_clusCol.size()<<endl;
/*
cout<<"  Print clusters: "<<endl;
for(int ic=0; ic<tmp_clusCol.size(); ic++){
  printf("    Nhit %d, Energy %.2f, HoughPar (%.2f, %.2f).  ", 
          tmp_clusCol[ic].getBarShowers().size(), tmp_clusCol[ic].getEnergy(),
          tmp_clusCol[ic].getHoughAlpha(), tmp_clusCol[ic].getHoughRho() );
  for(int ihit=0; ihit<tmp_clusCol[ic].getBarShowers().size(); ihit++) cout<<tmp_clusCol[ic].getBarShowers()[ihit].getDlayer()<<"  ";
  cout<<endl;
}
cout<<endl;
*/
cout<<"  Print clusters: "<<endl;
for(int ic=0; ic<m_clusCol.size(); ic++)
  printf("    Nhit %d, Energy %.2f, HoughPar (%.2f, %.2f) \n",
          m_clusCol[ic]->getBarShowers().size(), m_clusCol[ic]->getEnergy(), 
          m_clusCol[ic]->getHoughAlpha(), m_clusCol[ic]->getHoughRho() );
cout<<endl;


  CleanClusters(m_clusCol);

  m_longiClusCol.insert(m_longiClusCol.end(), m_clusCol.begin(), m_clusCol.end());

  return StatusCode::SUCCESS;
}


StatusCode HoughClusteringAlg::CleanClusters( std::vector<PandoraPlus::LongiCluster*>& m_longiClusCol ){
  if(m_longiClusCol.size()==0) return StatusCode::SUCCESS;

  //Depart the HoughCluster to 2 sub-clusters if it has blank in middle.
  for(int ic=0; ic<m_longiClusCol.size(); ic++){
    int m_nhit = m_longiClusCol[ic]->getBarShowers().size(); 
    m_longiClusCol[ic]->SortBarShowersByLayer();

//cout<<"In clus #"<<ic<<": hit size = "<<m_nhit<<endl;
//cout<<"Print the clus layer: ";
//for(int ah = 0; ah<m_nhit; ah++) cout<<m_longiClusCol[ic]->getBarShowers()[ah]->getDlayer()<<'\t';
//cout<<endl;

    for(int ih=0; ih<m_nhit-1; ih++){
    if(m_longiClusCol[ic]->getBarShowers()[ih+1]->getDlayer() - m_longiClusCol[ic]->getBarShowers()[ih]->getDlayer() > settings.map_floatPars["th_breakLaer"]){
//cout<<"  Need to depart in hit #"<<ih<<endl;

      PandoraPlus::LongiCluster* clus_head = new PandoraPlus::LongiCluster();
      PandoraPlus::LongiCluster* clus_tail = new PandoraPlus::LongiCluster();

      clus_head->setHoughPars( m_longiClusCol[ic]->getHoughAlpha(), m_longiClusCol[ic]->getHoughRho() );
      clus_head->setIntercept( m_longiClusCol[ic]->getHoughIntercept() );
      clus_tail->setHoughPars( m_longiClusCol[ic]->getHoughAlpha(), m_longiClusCol[ic]->getHoughRho() );
      clus_tail->setIntercept( m_longiClusCol[ic]->getHoughIntercept() );

      for(int jh=0; jh<=ih; jh++) clus_head->addBarShower( m_longiClusCol[ic]->getBarShowers()[jh], 0 );
      for(int jh=ih+1; jh<m_nhit; jh++) clus_tail->addBarShower( m_longiClusCol[ic]->getBarShowers()[jh], 0 );

//cout<<"  Head cluster size: "<<clus_head->getBarShowers().size()<<", isContinueN: "<<clus_head->isContinueN(settings.map_floatPars["th_continuetrkN"])<<endl;
//cout<<"  Tail cluster size: "<<clus_tail->getBarShowers().size()<<", isContinueN: "<<clus_tail->isContinueN(settings.map_floatPars["th_continuetrkN"])<<endl;

      if( clus_head->isContinueN(settings.map_floatPars["th_continuetrkN"]) ) m_longiClusCol.push_back(clus_head);
      if( clus_tail->isContinueN(settings.map_floatPars["th_continuetrkN"]) ) m_longiClusCol.push_back(clus_tail);

      delete m_longiClusCol[ic]; m_longiClusCol[ic] = NULL;
      m_longiClusCol.erase(m_longiClusCol.begin()+ic);
      ic--;
      break;
    }}

  }


cout<<"    After isolated hits removal: "<<m_longiClusCol.size()<<endl;
for(int il=0; il<m_longiClusCol.size(); il++){
  printf("    Clus#%d: Hough Par = (%.3f, %.3f, %.3f), shower layers: \n", il, m_longiClusCol[il]->getHoughAlpha(), m_longiClusCol[il]->getHoughRho(), m_longiClusCol[il]->getHoughIntercept());
  for(int is=0; is<m_longiClusCol[il]->getBarShowers().size(); is++)
    printf("      Dlayer = %d, Pos/E = (%.2f, %.2f, %.2f, %.3f), address %p \n", m_longiClusCol[il]->getBarShowers()[is]->getDlayer(),
                                                                     m_longiClusCol[il]->getBarShowers()[is]->getPos().x(),
                                                                     m_longiClusCol[il]->getBarShowers()[is]->getPos().y(),
                                                                     m_longiClusCol[il]->getBarShowers()[is]->getPos().z(),
                                                                     m_longiClusCol[il]->getBarShowers()[is]->getEnergy(),
                                                                     m_longiClusCol[il]->getBarShowers()[is] );
  cout<<endl;
}


  //Remove the repeated tracks
  for(int ic=0; ic<m_longiClusCol.size(); ic++){
  for(int jc=0; jc<m_longiClusCol.size(); jc++){

    if(ic>=m_longiClusCol.size()) ic--;
    if(ic==jc) continue;
//printf("   Check clus pair: (%d, %d): isSubset: %d \n", ic, jc, m_longiClusCol[ic]->isSubset(m_longiClusCol[jc]));

    if( m_longiClusCol[ic]->isSubset(m_longiClusCol[jc]) ){  //jc is the subset of ic. remove jc. 
      delete m_longiClusCol[jc]; m_longiClusCol[jc] = NULL;
      m_longiClusCol.erase(m_longiClusCol.begin()+jc );
      jc--;
      if(ic>jc+1) ic--;
    }
  }}


cout<<"    After subset removal: "<<m_longiClusCol.size()<<endl;
for(int il=0; il<m_longiClusCol.size(); il++){
  printf("    Clus#%d: Hough Par = (%.3f, %.3f, %.3f), E = %.3f, shower layers: \n", il, m_longiClusCol[il]->getHoughAlpha(), m_longiClusCol[il]->getHoughRho(), m_longiClusCol[il]->getHoughIntercept(), m_longiClusCol[il]->getEnergy());
  for(int is=0; is<m_longiClusCol[il]->getBarShowers().size(); is++)
    printf("      Dlayer = %d, Pos/E = (%.2f, %.2f, %.2f, %.3f) \n", m_longiClusCol[il]->getBarShowers()[is]->getDlayer(),
                                                                     m_longiClusCol[il]->getBarShowers()[is]->getPos().x(),
                                                                     m_longiClusCol[il]->getBarShowers()[is]->getPos().y(),
                                                                     m_longiClusCol[il]->getBarShowers()[is]->getPos().z(),
                                                                     m_longiClusCol[il]->getBarShowers()[is]->getEnergy() );
  cout<<endl;
}


  //Cut on energy and intercept
  for(int ic=0; ic<m_longiClusCol.size(); ic++){
    if( fabs(m_longiClusCol[ic]->getHoughIntercept())>=settings.map_floatPars["th_intercept"] || m_longiClusCol[ic]->getEnergy()<settings.map_floatPars["th_AxisE"]){
      delete m_longiClusCol[ic]; m_longiClusCol[ic]=NULL;
      m_longiClusCol.erase(m_longiClusCol.begin()+ic );
      ic--;
    }
  }
cout<<"    After energy and intercetp cut: "<<m_longiClusCol.size()<<endl;


  //Overlap with other clusters: Iter 1. 
  for(int ic=0; m_longiClusCol.size()>1 && ic<m_longiClusCol.size()-1; ic++){
  for(int jc=ic+1; m_longiClusCol.size()>1 && jc<m_longiClusCol.size(); jc++){
    if(ic>=m_longiClusCol.size()) ic--;

    if( fabs(m_longiClusCol[ic]->getHoughAlpha() - m_longiClusCol[jc]->getHoughAlpha()) > settings.map_floatPars["th_dAlpha1"] || 
        fabs(m_longiClusCol[ic]->getHoughRho() - m_longiClusCol[jc]->getHoughRho()) > settings.map_floatPars["th_dRho1"]  ) 
      continue;

    double m_ratio1 = m_longiClusCol[ic]->OverlapRatioE(m_longiClusCol[jc]);
    double m_ratio2 = m_longiClusCol[jc]->OverlapRatioE(m_longiClusCol[ic]);

//printf("      Tag branch: Clus #%d: (R, E) = (%.2f, %.2f). Clus #%d: (%.2f, %.2f). \n", ic, m_ratio1, m_longiClusCol[ic]->getEnergy(), jc, m_ratio2, m_longiClusCol[jc]->getEnergy() );
    
    if(m_ratio1>settings.map_floatPars["th_overlapE1"] && m_longiClusCol[ic]->getEnergy()<m_longiClusCol[jc]->getEnergy()){
      delete m_longiClusCol[ic]; m_longiClusCol[ic] = NULL;
      m_longiClusCol.erase( m_longiClusCol.begin()+ic );
      ic--;
      break;
    }

    if(m_ratio2>settings.map_floatPars["th_overlapE1"] && m_longiClusCol[jc]->getEnergy()<m_longiClusCol[ic]->getEnergy()){
      delete m_longiClusCol[jc]; m_longiClusCol[jc] = NULL;
      m_longiClusCol.erase( m_longiClusCol.begin()+jc );
      jc--;
    }

  }}
cout<<"    After 1st iter tagging: "<<m_longiClusCol.size()<<endl;


  //Overlap with other clusters: Iter 2. 
  for(int ic=0; ic<m_longiClusCol.size()-1; ic++){
  for(int jc=ic+1; jc<m_longiClusCol.size(); jc++){
    if(ic>=m_longiClusCol.size()) ic--;

    if( fabs(m_longiClusCol[ic]->getHoughAlpha() - m_longiClusCol[jc]->getHoughAlpha()) > settings.map_floatPars["th_dAlpha2"] ||
        fabs(m_longiClusCol[ic]->getHoughRho() - m_longiClusCol[jc]->getHoughRho()) > settings.map_floatPars["th_dRho2"]  )
      continue;

    double m_ratio1 = m_longiClusCol[ic]->OverlapRatioE(m_longiClusCol[jc]);
    double m_ratio2 = m_longiClusCol[jc]->OverlapRatioE(m_longiClusCol[ic]);

//printf("      Tag branch: Clus #%d: (R, E) = (%.2f, %.2f). Clus #%d: (%.2f, %.2f). \n", ic, m_ratio1, m_longiClusCol[ic]->getEnergy(), jc, m_ratio2, m_longiClusCol[jc]->getEnergy() );

    if(m_ratio1>settings.map_floatPars["th_overlapE2"] && m_longiClusCol[ic]->getEnergy()<m_longiClusCol[jc]->getEnergy()){
      delete m_longiClusCol[ic]; m_longiClusCol[ic]=NULL;
      m_longiClusCol.erase( m_longiClusCol.begin()+ic );
      ic--;
      break;
    }

    if(m_ratio2>settings.map_floatPars["th_overlapE2"] && m_longiClusCol[jc]->getEnergy()<m_longiClusCol[ic]->getEnergy()){
      delete m_longiClusCol[jc]; m_longiClusCol[jc]=NULL;
      m_longiClusCol.erase( m_longiClusCol.begin()+jc );
      jc--;
    }

  }}
cout<<"    After 2nd iter tagging: "<<m_longiClusCol.size()<<endl;

/*
cout<<"  Print clusters: "<<endl;
for(int ic=0; ic<m_longiClusCol.size(); ic++)
  printf("    Nhit %d, Energy %.2f, HoughPar (%.2f, %.2f, %.2f) \n",
          m_longiClusCol[ic].getBarShowers().size(), m_longiClusCol[ic].getEnergy(), 
          m_longiClusCol[ic].getHoughAlpha(), m_longiClusCol[ic].getHoughRho(), m_longiClusCol[ic].getHoughIntercept() );
cout<<endl;
*/
  return StatusCode::SUCCESS;
};

#endif
