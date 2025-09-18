#ifndef HOUGHCLUSTERINGALG_C
#define HOUGHCLUSTERINGALG_C

#include "Algorithm/HoughClusteringAlg.h"
#include <algorithm>
#include <cmath>
#include <set>

using namespace std;

StatusCode HoughClusteringAlg::ReadSettings(Cyber::Settings& m_settings){
  settings = m_settings;

  // ECAL geometry settings
  //if(settings.map_floatPars.find("cell_size")==settings.map_floatPars.end())
  //  settings.map_floatPars["cell_size"] = 10; // unit: mm
  //if(settings.map_floatPars.find("ecal_inner_radius")==settings.map_floatPars.end())
  //  settings.map_floatPars["ecal_inner_radius"] = 1860; // unit: mm

  // Hough space settings
  // alpha in V plane (bars parallel to z axis)
  if(settings.map_floatPars.find("alpha_lowV")==settings.map_floatPars.end())
    settings.map_floatPars["alpha_lowV"] = -0.1;  
  if(settings.map_floatPars.find("alpha_highV")==settings.map_floatPars.end())
    settings.map_floatPars["alpha_highV"] = 2.*TMath::Pi(); 
  if(settings.map_intPars.find("Nbins_alphaV")==settings.map_intPars.end())
    settings.map_intPars["Nbins_alphaV"] = 3000; 
  if(settings.map_floatPars.find("bin_width_alphaV")==settings.map_floatPars.end())
    settings.map_floatPars["bin_width_alphaV"] = ( settings.map_floatPars["alpha_highV"] - settings.map_floatPars["alpha_lowV"] ) / (double)settings.map_intPars["Nbins_alphaV"];
  // double bin_width_alphaV = (alpha_highV - alpha_lowV) / (double)Nbins_alphaV;
  
  // alpha in U plane (bars perpendicular to z axis)
  if(settings.map_floatPars.find("alpha_lowU")==settings.map_floatPars.end())
    settings.map_floatPars["alpha_lowU"] = 0.;
  if(settings.map_floatPars.find("alpha_highU")==settings.map_floatPars.end())
    settings.map_floatPars["alpha_highU"] = TMath::Pi();  
  if(settings.map_intPars.find("Nbins_alphaU")==settings.map_intPars.end())
    settings.map_intPars["Nbins_alphaU"] = 5000;  
  if(settings.map_floatPars.find("bin_width_alphaU")==settings.map_floatPars.end())
    settings.map_floatPars["bin_width_alphaU"] = ( settings.map_floatPars["alpha_highU"] - settings.map_floatPars["alpha_lowU"] ) / (double)settings.map_intPars["Nbins_alphaU"];
  // double bin_width_alphaU = (alpha_highU - alpha_lowU) / (double)Nbins_alphaU;

  // alpha in endcap. Works for both V and U planes
  if(settings.map_floatPars.find("alpha_low_endcap")==settings.map_floatPars.end())
    settings.map_floatPars["alpha_low_endcap"] = 0.;
  if(settings.map_floatPars.find("alpha_high_endcap")==settings.map_floatPars.end())
    settings.map_floatPars["alpha_high_endcap"] = TMath::Pi();
  if(settings.map_intPars.find("Nbins_alpha_endcap")==settings.map_intPars.end())
    settings.map_intPars["Nbins_alpha_endcap"] = 3000;
  if(settings.map_floatPars.find("bin_width_alpha_endcap")==settings.map_floatPars.end())
    settings.map_floatPars["bin_width_alpha_endcap"] = ( settings.map_floatPars["alpha_high_endcap"] - settings.map_floatPars["alpha_low_endcap"] ) / (double)settings.map_intPars["Nbins_alpha_endcap"]; 

  // rho
  if(settings.map_floatPars.find("rho_low")==settings.map_floatPars.end())
    settings.map_floatPars["rho_low"] = -50.;
  if(settings.map_floatPars.find("rho_high")==settings.map_floatPars.end())
    settings.map_floatPars["rho_high"] = 50.;
  if(settings.map_intPars.find("Nbins_rho")==settings.map_intPars.end())
    settings.map_intPars["Nbins_rho"] = 20;   // (rho_high - rho_low)/5
  if(settings.map_floatPars.find("bin_width_rho")==settings.map_floatPars.end())
    settings.map_floatPars["bin_width_rho"] = ( settings.map_floatPars["rho_high"] - settings.map_floatPars["rho_low"] ) / (double)settings.map_intPars["Nbins_rho"];

  // Algorithm parameter settings
  if(settings.map_intPars.find("th_Layers")==settings.map_intPars.end())
    settings.map_intPars["th_Layers"] = 7;
  if(settings.map_intPars.find("th_peak")==settings.map_intPars.end())
    settings.map_intPars["th_peak"] = 3;
  if(settings.map_intPars.find("th_continueN")==settings.map_intPars.end())
    settings.map_intPars["th_continueN"] = 3;
  if(settings.map_floatPars.find("th_AxisE")==settings.map_floatPars.end())
    settings.map_floatPars["th_AxisE"] = 0.15; // unit: GeV
  if(settings.map_floatPars.find("th_overlapE")==settings.map_floatPars.end())
    settings.map_floatPars["th_overlapE"] = 0.5;
  if(settings.map_floatPars.find("th_dAlpha1")==settings.map_floatPars.end())
    settings.map_floatPars["th_dAlpha1"] = 0.1;
  if(settings.map_floatPars.find("th_dAlpha2")==settings.map_floatPars.end())
    settings.map_floatPars["th_dAlpha2"] = 0.05;
  if(settings.map_floatPars.find("th_ERatio")==settings.map_floatPars.end())
    settings.map_floatPars["th_ERatio"] = 0.04;

  if(settings.map_stringPars.find("ReadinLocalMaxName")==settings.map_stringPars.end())  
    settings.map_stringPars["ReadinLocalMaxName"] = "AllLocalMax";
  if(settings.map_stringPars.find("LeftLocalMaxName")==settings.map_stringPars.end())    
    settings.map_stringPars["LeftLocalMaxName"] = "LeftLocalMax";
  if(settings.map_stringPars.find("OutputLongiClusName")==settings.map_stringPars.end()) 
    settings.map_stringPars["OutputLongiClusName"] = "HoughAxis"; 

  return StatusCode::SUCCESS;  
}

StatusCode HoughClusteringAlg::Initialize( CyberDataCol& m_datacol ){
  barrel_HalfClusterU.clear(); 
  barrel_HalfClusterV.clear(); 
  endcap0_HalfClusterV.clear();
  endcap0_HalfClusterU.clear();
  endcap1_HalfClusterV.clear();
  endcap1_HalfClusterU.clear();

  for(int ih=0; ih<m_datacol.map_HalfCluster["HalfClusterColU"].size(); ih++){
    if(m_datacol.map_HalfCluster["HalfClusterColU"][ih].get()->getBars()[0]->getSystem()==Cyber::CaloUnit::System_Barrel){
      barrel_HalfClusterU.push_back( m_datacol.map_HalfCluster["HalfClusterColU"][ih].get() );
    }
    else if ((m_datacol.map_HalfCluster["HalfClusterColU"][ih].get()->getBars()[0]->getSystem()==Cyber::CaloUnit::System_Endcap) 
              && (m_datacol.map_HalfCluster["HalfClusterColU"][ih].get()->getBars()[0]->getModule()==0)){
      endcap0_HalfClusterU.push_back( m_datacol.map_HalfCluster["HalfClusterColU"][ih].get() );
    }
    else if ((m_datacol.map_HalfCluster["HalfClusterColU"][ih].get()->getBars()[0]->getSystem()==Cyber::CaloUnit::System_Endcap) 
              && (m_datacol.map_HalfCluster["HalfClusterColU"][ih].get()->getBars()[0]->getModule()==1)){
      endcap1_HalfClusterU.push_back( m_datacol.map_HalfCluster["HalfClusterColU"][ih].get() );
    }
  }
  for(int ih=0; ih<m_datacol.map_HalfCluster["HalfClusterColV"].size(); ih++){
    if(m_datacol.map_HalfCluster["HalfClusterColV"][ih].get()->getBars()[0]->getSystem()==Cyber::CaloUnit::System_Barrel){
      barrel_HalfClusterV.push_back( m_datacol.map_HalfCluster["HalfClusterColV"][ih].get() );
    }
    else if ((m_datacol.map_HalfCluster["HalfClusterColV"][ih].get()->getBars()[0]->getSystem()==Cyber::CaloUnit::System_Endcap) 
              && (m_datacol.map_HalfCluster["HalfClusterColV"][ih].get()->getBars()[0]->getModule()==0)){
      endcap0_HalfClusterV.push_back( m_datacol.map_HalfCluster["HalfClusterColV"][ih].get() );
    }
    else if ((m_datacol.map_HalfCluster["HalfClusterColV"][ih].get()->getBars()[0]->getSystem()==Cyber::CaloUnit::System_Endcap) 
              && (m_datacol.map_HalfCluster["HalfClusterColV"][ih].get()->getBars()[0]->getModule()==1)){
      endcap1_HalfClusterV.push_back( m_datacol.map_HalfCluster["HalfClusterColV"][ih].get() );
    }
  }

//cout << "yyy: Number of barrel_HalfClusterU: " << barrel_HalfClusterU.size() << endl;
//cout << "yyy: Number of barrel_HalfClusterV: " << barrel_HalfClusterV.size() << endl;


  return StatusCode::SUCCESS;
}

StatusCode HoughClusteringAlg::RunAlgorithm( CyberDataCol& m_datacol ){
  if(barrel_HalfClusterV.size()==0 && endcap0_HalfClusterV.size()==0 && endcap1_HalfClusterV.size()==0){ 
    std::cout<<"  HoughClusteringAlg: No HalfClusterV in present data collection! "<<std::endl; 
  }
  if(barrel_HalfClusterU.size()==0 && endcap0_HalfClusterU.size()==0 && endcap1_HalfClusterU.size()==0){ 
    std::cout<<"  HoughClusteringAlg: No HalfClusterU in present data collection! "<<std::endl; 
  }

  // Processing V bar in barrel (V bar parallel to z axis)
  for(int it=0; it<barrel_HalfClusterV.size(); it++){ // process each HalfCluster respectively
    m_localMaxVCol.clear();
    std::vector<const Cyber::Calo1DCluster*> tmp_localMaxVCol = barrel_HalfClusterV[it]->getLocalMaxCol(settings.map_stringPars["ReadinLocalMaxName"]);

    for(int il=0; il<tmp_localMaxVCol.size(); il++){
      if(tmp_localMaxVCol[il]->getDlayer()<=settings.map_intPars["th_Layers"]) 
        m_localMaxVCol.push_back(tmp_localMaxVCol[il]);
    }

    if(m_localMaxVCol.size()<settings.map_intPars["th_peak"]){
      continue; 
    } 

    std::vector<Cyber::HoughObject> m_HoughObjectsV; m_HoughObjectsV.clear(); 
    for(int il=0; il<m_localMaxVCol.size(); il++){
      Cyber::HoughObject m_obj(m_localMaxVCol[il], Cyber::CaloUnit::barsize, Cyber::CaloUnit::ecal_innerR);
      m_HoughObjectsV.push_back(m_obj);
    }

    HoughTransformation(m_HoughObjectsV);

    // cout<<"  HoughClusteringAlg: Creating hough_spaceV"<<endl;
    Cyber::HoughSpace hough_spaceV( settings.map_floatPars["alpha_lowV"], settings.map_floatPars["alpha_highV"], 
                                    settings.map_floatPars["bin_width_alphaV"], settings.map_intPars["Nbins_alphaV"], 
                                    settings.map_floatPars["rho_low"], settings.map_floatPars["rho_high"], 
                                    settings.map_floatPars["bin_width_rho"], settings.map_intPars["Nbins_rho"]);

    FillHoughSpace(m_HoughObjectsV, hough_spaceV);

    //Create output HoughClusters
    m_longiClusVCol.clear(); 
    ClusterFinding(m_HoughObjectsV, hough_spaceV, m_longiClusVCol);

    CleanClusters(m_longiClusVCol);
    m_datacol.map_HalfCluster["bkHalfCluster"].insert( m_datacol.map_HalfCluster["bkHalfCluster"].end(), m_longiClusVCol.begin(), m_longiClusVCol.end() );

    std::vector<const Cyber::CaloHalfCluster*> m_constHoughCluster; m_constHoughCluster.clear();
    std::vector<const Cyber::Calo1DCluster*> left_localMaxVCol; left_localMaxVCol.clear();
    std::vector<const Cyber::Calo1DCluster*> m_houghMax; m_houghMax.clear(); 
    for(int is=0; is<tmp_localMaxVCol.size(); is++){
      bool fl_incluster = false; 
      for(int ic=0; ic<m_longiClusVCol.size(); ic++){
        std::vector<const Cyber::Calo1DCluster*> p_showers = m_longiClusVCol[ic]->getCluster();
        if( find(p_showers.begin(), p_showers.end(), tmp_localMaxVCol[is])!=p_showers.end() ) { fl_incluster = true; break; }
      }
      if(!fl_incluster && find(left_localMaxVCol.begin(), left_localMaxVCol.end(), tmp_localMaxVCol[is])==left_localMaxVCol.end() ) left_localMaxVCol.push_back(tmp_localMaxVCol[is]);
      m_houghMax.push_back( tmp_localMaxVCol[is] );
    }
    for(int ic=0; ic<m_longiClusVCol.size(); ic++)
      m_constHoughCluster.push_back(m_longiClusVCol[ic].get());


    barrel_HalfClusterV[it]->setLocalMax("HoughLocalMax", m_houghMax);
    barrel_HalfClusterV[it]->setLocalMax(settings.map_stringPars["LeftLocalMaxName"], left_localMaxVCol);
    barrel_HalfClusterV[it]->setHalfClusters(settings.map_stringPars["OutputLongiClusName"], m_constHoughCluster);
    m_houghMax.clear();
    left_localMaxVCol.clear();

  }

  // Processing U bar in barrel (U bar perpendicular to z axis)
  for(int it=0; it<barrel_HalfClusterU.size(); it++){ // process each HalfCluster respectively
    m_localMaxUCol.clear();
    std::vector<const Cyber::Calo1DCluster*> tmp_localMaxUCol = barrel_HalfClusterU[it]->getLocalMaxCol(settings.map_stringPars["ReadinLocalMaxName"]);

    for(int il=0; il<tmp_localMaxUCol.size(); il++){
      if(tmp_localMaxUCol[il]->getDlayer()<=settings.map_intPars["th_Layers"])
        m_localMaxUCol.push_back(tmp_localMaxUCol[il]);
    }

    if(m_localMaxUCol.size()<settings.map_intPars["th_peak"]){
      continue;
    }


    std::vector<Cyber::HoughObject> m_HoughObjectsU; m_HoughObjectsU.clear();
    for(int il=0; il<m_localMaxUCol.size(); il++){
      Cyber::HoughObject m_obj(m_localMaxUCol[il], Cyber::CaloUnit::barsize, Cyber::CaloUnit::ecal_innerR);
      m_HoughObjectsU.push_back(m_obj);
    }

    HoughTransformation(m_HoughObjectsU);

    // cout<<"  HoughClusteringAlg: Creating hough_spaceU"<<endl;
    Cyber::HoughSpace hough_spaceU(settings.map_floatPars["alpha_lowU"], settings.map_floatPars["alpha_highU"],
                                         settings.map_floatPars["bin_width_alphaU"], settings.map_intPars["Nbins_alphaU"],
                                         settings.map_floatPars["rho_low"], settings.map_floatPars["rho_high"],
                                         settings.map_floatPars["bin_width_rho"], settings.map_intPars["Nbins_rho"]);
    FillHoughSpace(m_HoughObjectsU, hough_spaceU);

    //Create output HoughClusters
    m_longiClusUCol.clear();
    ClusterFinding(m_HoughObjectsU, hough_spaceU, m_longiClusUCol  );
    CleanClusters(m_longiClusUCol);

    m_datacol.map_HalfCluster["bkHalfCluster"].insert( m_datacol.map_HalfCluster["bkHalfCluster"].end(), m_longiClusUCol.begin(), m_longiClusUCol.end() );

    std::vector<const Cyber::CaloHalfCluster*> m_constHoughCluster; m_constHoughCluster.clear();
    std::vector<const Cyber::Calo1DCluster*> left_localMaxUCol; left_localMaxUCol.clear();
    std::vector<const Cyber::Calo1DCluster*> m_houghMax; m_houghMax.clear();
    for(int is=0; is<tmp_localMaxUCol.size(); is++){
      bool fl_incluster = false;
      for(int ic=0; ic<m_longiClusUCol.size(); ic++){
        std::vector<const Cyber::Calo1DCluster*> p_showers = m_longiClusUCol[ic]->getCluster();
        if( find(p_showers.begin(), p_showers.end(), tmp_localMaxUCol[is])!=p_showers.end() ) { fl_incluster = true; break; }
      }
      if(!fl_incluster && find(left_localMaxUCol.begin(), left_localMaxUCol.end(), tmp_localMaxUCol[is])==left_localMaxUCol.end() ) left_localMaxUCol.push_back(tmp_localMaxUCol[is]);
      m_houghMax.push_back( tmp_localMaxUCol[is] );
    }
    for(int ic=0; ic<m_longiClusUCol.size(); ic++)
      m_constHoughCluster.push_back(m_longiClusUCol[ic].get());

    barrel_HalfClusterU[it]->setLocalMax("HoughLocalMax", m_houghMax);
    barrel_HalfClusterU[it]->setLocalMax(settings.map_stringPars["LeftLocalMaxName"], left_localMaxUCol);
    barrel_HalfClusterU[it]->setHalfClusters(settings.map_stringPars["OutputLongiClusName"], m_constHoughCluster);
    m_houghMax.clear();
    left_localMaxUCol.clear();

  }  // end of U plane

  // Processing V bar in endcap 0 (V bar parallel to x axis, endcap 0 at z~-2900mm)
  for(int it=0; it<endcap0_HalfClusterV.size(); it++){ // process each HalfCluster respectively
    m_localMaxVCol.clear();
    std::vector<const Cyber::Calo1DCluster*> tmp_localMaxVCol = endcap0_HalfClusterV[it]->getLocalMaxCol(settings.map_stringPars["ReadinLocalMaxName"]);

    for(int il=0; il<tmp_localMaxVCol.size(); il++){
      if(tmp_localMaxVCol[il]->getDlayer()<=settings.map_intPars["th_Layers"])
        m_localMaxVCol.push_back(tmp_localMaxVCol[il]);
    }

    if(m_localMaxVCol.size()<settings.map_intPars["th_peak"]){
      continue;
    }

    std::vector<Cyber::HoughObject> m_HoughObjectsV; m_HoughObjectsV.clear();
    for(int il=0; il<m_localMaxVCol.size(); il++){
      Cyber::HoughObject m_obj(m_localMaxVCol[il], Cyber::CaloUnit::barsize, Cyber::CaloUnit::ecal_innerR);
      m_HoughObjectsV.push_back(m_obj);
    }

    HoughTransformation(m_HoughObjectsV);

    // cout<<"  HoughClusteringAlg: Creating hough_spaceV"<<endl;
    Cyber::HoughSpace hough_spaceV(settings.map_floatPars["alpha_low_endcap"], settings.map_floatPars["alpha_high_endcap"],
                                    settings.map_floatPars["bin_width_alpha_endcap"], settings.map_intPars["Nbins_alpha_endcap"],
                                    settings.map_floatPars["rho_low"], settings.map_floatPars["rho_high"],
                                    settings.map_floatPars["bin_width_rho"], settings.map_intPars["Nbins_rho"]);
    
    FillHoughSpace(m_HoughObjectsV, hough_spaceV);

    //Create output HoughClusters
    m_longiClusVCol.clear();
    ClusterFinding(m_HoughObjectsV, hough_spaceV, m_longiClusVCol);

    CleanClusters(m_longiClusVCol);
    m_datacol.map_HalfCluster["bkHalfCluster"].insert( m_datacol.map_HalfCluster["bkHalfCluster"].end(), m_longiClusVCol.begin(), m_longiClusVCol.end() );

    std::vector<const Cyber::CaloHalfCluster*> m_constHoughCluster; m_constHoughCluster.clear();
    std::vector<const Cyber::Calo1DCluster*> left_localMaxVCol; left_localMaxVCol.clear();
    std::vector<const Cyber::Calo1DCluster*> m_houghMax; m_houghMax.clear();
    for(int is=0; is<tmp_localMaxVCol.size(); is++){
      bool fl_incluster = false;
      for(int ic=0; ic<m_longiClusVCol.size(); ic++){
        std::vector<const Cyber::Calo1DCluster*> p_showers = m_longiClusVCol[ic]->getCluster();
        if( find(p_showers.begin(), p_showers.end(), tmp_localMaxVCol[is])!=p_showers.end() ) { fl_incluster = true; break; }
      }
      if(!fl_incluster && find(left_localMaxVCol.begin(), left_localMaxVCol.end(), tmp_localMaxVCol[is])==left_localMaxVCol.end() ) left_localMaxVCol.push_back(tmp_localMaxVCol[is]);
      m_houghMax.push_back( tmp_localMaxVCol[is] );
    }
    for(int ic=0; ic<m_longiClusVCol.size(); ic++)
      m_constHoughCluster.push_back(m_longiClusVCol[ic].get());

    endcap0_HalfClusterV[it]->setLocalMax("HoughLocalMax", m_houghMax);
    endcap0_HalfClusterV[it]->setLocalMax(settings.map_stringPars["LeftLocalMaxName"], left_localMaxVCol);
    endcap0_HalfClusterV[it]->setHalfClusters(settings.map_stringPars["OutputLongiClusName"], m_constHoughCluster);
    m_houghMax.clear();
    left_localMaxVCol.clear();
  }


  // Processing U bar in endcap 0 (U bar parallel to y axis, endcap 0 at z~-2900mm)
  for(int it=0; it<endcap0_HalfClusterU.size(); it++){ // process each HalfCluster respectively
    m_localMaxUCol.clear();
    std::vector<const Cyber::Calo1DCluster*> tmp_localMaxUCol = endcap0_HalfClusterU[it]->getLocalMaxCol(settings.map_stringPars["ReadinLocalMaxName"]);

    for(int il=0; il<tmp_localMaxUCol.size(); il++){
      if(tmp_localMaxUCol[il]->getDlayer()<=settings.map_intPars["th_Layers"])
        m_localMaxUCol.push_back(tmp_localMaxUCol[il]);
    }

    if(m_localMaxUCol.size()<settings.map_intPars["th_peak"]){
      continue;
    }

    std::vector<Cyber::HoughObject> m_HoughObjectsU; m_HoughObjectsU.clear();
    for(int il=0; il<m_localMaxUCol.size(); il++){
      Cyber::HoughObject m_obj(m_localMaxUCol[il], Cyber::CaloUnit::barsize, Cyber::CaloUnit::ecal_innerR);
      m_HoughObjectsU.push_back(m_obj);
    }

    HoughTransformation(m_HoughObjectsU);

    // cout<<"  HoughClusteringAlg: Creating hough_spaceU"<<endl;
    Cyber::HoughSpace hough_spaceU(settings.map_floatPars["alpha_low_endcap"], settings.map_floatPars["alpha_high_endcap"],
                                    settings.map_floatPars["bin_width_alpha_endcap"], settings.map_intPars["Nbins_alpha_endcap"],
                                    settings.map_floatPars["rho_low"], settings.map_floatPars["rho_high"],
                                    settings.map_floatPars["bin_width_rho"], settings.map_intPars["Nbins_rho"]);

    FillHoughSpace(m_HoughObjectsU, hough_spaceU);

    //Create output HoughClusters
    m_longiClusUCol.clear();
    ClusterFinding(m_HoughObjectsU, hough_spaceU, m_longiClusUCol);

    CleanClusters(m_longiClusUCol);
    m_datacol.map_HalfCluster["bkHalfCluster"].insert( m_datacol.map_HalfCluster["bkHalfCluster"].end(), m_longiClusUCol.begin(), m_longiClusUCol.end() );

    std::vector<const Cyber::CaloHalfCluster*> m_constHoughCluster; m_constHoughCluster.clear();
    std::vector<const Cyber::Calo1DCluster*> left_localMaxUCol; left_localMaxUCol.clear();
    std::vector<const Cyber::Calo1DCluster*> m_houghMax; m_houghMax.clear();
    for(int is=0; is<tmp_localMaxUCol.size(); is++){
      bool fl_incluster = false;
      for(int ic=0; ic<m_longiClusUCol.size(); ic++){
        std::vector<const Cyber::Calo1DCluster*> p_showers = m_longiClusUCol[ic]->getCluster();
        if( find(p_showers.begin(), p_showers.end(), tmp_localMaxUCol[is])!=p_showers.end() ) { fl_incluster = true; break; }
      }
      if(!fl_incluster && find(left_localMaxUCol.begin(), left_localMaxUCol.end(), tmp_localMaxUCol[is])==left_localMaxUCol.end() ) left_localMaxUCol.push_back(tmp_localMaxUCol[is]);
      m_houghMax.push_back( tmp_localMaxUCol[is] );
    }
    for(int ic=0; ic<m_longiClusUCol.size(); ic++)
      m_constHoughCluster.push_back(m_longiClusUCol[ic].get());

    endcap0_HalfClusterU[it]->setLocalMax("HoughLocalMax", m_houghMax);
    endcap0_HalfClusterU[it]->setLocalMax(settings.map_stringPars["LeftLocalMaxName"], left_localMaxUCol);
    endcap0_HalfClusterU[it]->setHalfClusters(settings.map_stringPars["OutputLongiClusName"], m_constHoughCluster);
    m_houghMax.clear();
    left_localMaxUCol.clear();
  }


  // Processing V bar in endcap 1 (V bar parallel to x axis, endcap 1 at z~2900mm)
  for(int it=0; it<endcap1_HalfClusterV.size(); it++){ // process each HalfCluster respectively
    m_localMaxVCol.clear();
    std::vector<const Cyber::Calo1DCluster*> tmp_localMaxVCol = endcap1_HalfClusterV[it]->getLocalMaxCol(settings.map_stringPars["ReadinLocalMaxName"]);

    for(int il=0; il<tmp_localMaxVCol.size(); il++){
      if(tmp_localMaxVCol[il]->getDlayer()<=settings.map_intPars["th_Layers"])
        m_localMaxVCol.push_back(tmp_localMaxVCol[il]);
    }
    //cout << "yyy: Number of localMaxVCol: " << m_localMaxVCol.size() << endl;

    if(m_localMaxVCol.size()<settings.map_intPars["th_peak"]){
      continue;
    }

    std::vector<Cyber::HoughObject> m_HoughObjectsV; m_HoughObjectsV.clear();
    for(int il=0; il<m_localMaxVCol.size(); il++){
      Cyber::HoughObject m_obj(m_localMaxVCol[il], Cyber::CaloUnit::barsize, Cyber::CaloUnit::ecal_innerR);
      m_HoughObjectsV.push_back(m_obj);
    }
    //cout << "yyy: Number of HoughObjectsV: " << m_HoughObjectsV.size() << endl;

    HoughTransformation(m_HoughObjectsV);

    // cout<<"  HoughClusteringAlg: Creating hough_spaceV"<<endl;
    Cyber::HoughSpace hough_spaceV(settings.map_floatPars["alpha_low_endcap"], settings.map_floatPars["alpha_high_endcap"],
                                    settings.map_floatPars["bin_width_alpha_endcap"], settings.map_intPars["Nbins_alpha_endcap"],
                                    settings.map_floatPars["rho_low"], settings.map_floatPars["rho_high"],
                                    settings.map_floatPars["bin_width_rho"], settings.map_intPars["Nbins_rho"]);

    FillHoughSpace(m_HoughObjectsV, hough_spaceV);

    //Create output HoughClusters
    m_longiClusVCol.clear();
    ClusterFinding(m_HoughObjectsV, hough_spaceV, m_longiClusVCol);
    //cout << "yyy: Number of longiClusVCol: " << m_longiClusVCol.size() << endl;

    CleanClusters(m_longiClusVCol);
    //cout << "yyy: Number of longiClusVCol after cleaning: " << m_longiClusVCol.size() << endl;
    m_datacol.map_HalfCluster["bkHalfCluster"].insert( m_datacol.map_HalfCluster["bkHalfCluster"].end(), m_longiClusVCol.begin(), m_longiClusVCol.end() );

    std::vector<const Cyber::CaloHalfCluster*> m_constHoughCluster; m_constHoughCluster.clear();
    std::vector<const Cyber::Calo1DCluster*> left_localMaxVCol; left_localMaxVCol.clear();
    std::vector<const Cyber::Calo1DCluster*> m_houghMax; m_houghMax.clear();
    for(int is=0; is<tmp_localMaxVCol.size(); is++){
      bool fl_incluster = false;
      for(int ic=0; ic<m_longiClusVCol.size(); ic++){
        std::vector<const Cyber::Calo1DCluster*> p_showers = m_longiClusVCol[ic]->getCluster();
        if( find(p_showers.begin(), p_showers.end(), tmp_localMaxVCol[is])!=p_showers.end() ) { fl_incluster = true; break; }
      }
      if(!fl_incluster && find(left_localMaxVCol.begin(), left_localMaxVCol.end(), tmp_localMaxVCol[is])==left_localMaxVCol.end() ) left_localMaxVCol.push_back(tmp_localMaxVCol[is]);
      m_houghMax.push_back( tmp_localMaxVCol[is] );
    }
    for(int ic=0; ic<m_longiClusVCol.size(); ic++)
      m_constHoughCluster.push_back(m_longiClusVCol[ic].get());

    endcap1_HalfClusterV[it]->setLocalMax("HoughLocalMax", m_houghMax);
    endcap1_HalfClusterV[it]->setLocalMax(settings.map_stringPars["LeftLocalMaxName"], left_localMaxVCol);
    endcap1_HalfClusterV[it]->setHalfClusters(settings.map_stringPars["OutputLongiClusName"], m_constHoughCluster);
    m_houghMax.clear();
    left_localMaxVCol.clear();
  }


  // Processing U bar in endcap 1 (U bar parallel to y axis, endcap 1 at z~2900mm)
  for(int it=0; it<endcap1_HalfClusterU.size(); it++){ // process each HalfCluster respectively
    m_localMaxUCol.clear();
    std::vector<const Cyber::Calo1DCluster*> tmp_localMaxUCol = endcap1_HalfClusterU[it]->getLocalMaxCol(settings.map_stringPars["ReadinLocalMaxName"]);

    for(int il=0; il<tmp_localMaxUCol.size(); il++){
      if(tmp_localMaxUCol[il]->getDlayer()<=settings.map_intPars["th_Layers"])
        m_localMaxUCol.push_back(tmp_localMaxUCol[il]);
    }

    if(m_localMaxUCol.size()<settings.map_intPars["th_peak"]){
      continue;
    }

    std::vector<Cyber::HoughObject> m_HoughObjectsU; m_HoughObjectsU.clear();
    for(int il=0; il<m_localMaxUCol.size(); il++){
      Cyber::HoughObject m_obj(m_localMaxUCol[il], Cyber::CaloUnit::barsize, Cyber::CaloUnit::ecal_innerR);
      m_HoughObjectsU.push_back(m_obj);
    }

    HoughTransformation(m_HoughObjectsU);

    // cout<<"  HoughClusteringAlg: Creating hough_spaceU"<<endl;
    Cyber::HoughSpace hough_spaceU(settings.map_floatPars["alpha_low_endcap"], settings.map_floatPars["alpha_high_endcap"],
                                    settings.map_floatPars["bin_width_alpha_endcap"], settings.map_intPars["Nbins_alpha_endcap"],
                                    settings.map_floatPars["rho_low"], settings.map_floatPars["rho_high"],
                                    settings.map_floatPars["bin_width_rho"], settings.map_intPars["Nbins_rho"]);

    FillHoughSpace(m_HoughObjectsU, hough_spaceU);

    //Create output HoughClusters
    m_longiClusUCol.clear();
    ClusterFinding(m_HoughObjectsU, hough_spaceU, m_longiClusUCol);

    CleanClusters(m_longiClusUCol);
    m_datacol.map_HalfCluster["bkHalfCluster"].insert( m_datacol.map_HalfCluster["bkHalfCluster"].end(), m_longiClusUCol.begin(), m_longiClusUCol.end() );

    std::vector<const Cyber::CaloHalfCluster*> m_constHoughCluster; m_constHoughCluster.clear();
    std::vector<const Cyber::Calo1DCluster*> left_localMaxUCol; left_localMaxUCol.clear();
    std::vector<const Cyber::Calo1DCluster*> m_houghMax; m_houghMax.clear();
    for(int is=0; is<tmp_localMaxUCol.size(); is++){
      bool fl_incluster = false;
      for(int ic=0; ic<m_longiClusUCol.size(); ic++){
        std::vector<const Cyber::Calo1DCluster*> p_showers = m_longiClusUCol[ic]->getCluster();
        if( find(p_showers.begin(), p_showers.end(), tmp_localMaxUCol[is])!=p_showers.end() ) { fl_incluster = true; break; }
      }
      if(!fl_incluster && find(left_localMaxUCol.begin(), left_localMaxUCol.end(), tmp_localMaxUCol[is])==left_localMaxUCol.end() ) left_localMaxUCol.push_back(tmp_localMaxUCol[is]);
      m_houghMax.push_back( tmp_localMaxUCol[is] );
    }
    for(int ic=0; ic<m_longiClusUCol.size(); ic++)
      m_constHoughCluster.push_back(m_longiClusUCol[ic].get());

    endcap1_HalfClusterU[it]->setLocalMax("HoughLocalMax", m_houghMax);
    endcap1_HalfClusterU[it]->setLocalMax(settings.map_stringPars["LeftLocalMaxName"], left_localMaxUCol);
    endcap1_HalfClusterU[it]->setHalfClusters(settings.map_stringPars["OutputLongiClusName"], m_constHoughCluster);
    m_houghMax.clear();
    left_localMaxUCol.clear();
  }
/*
int Ncl_hough = 0;
int Ncl_trk = 0;
int Naxis_hough = 0;
int Naxis_trk = 0;
double Etot_hough = 0.;
double Etot_trk = 0.;
for(int i=0; i<barrel_HalfClusterU.size(); i++){
  //printf("  HalfClusterU #%d: energy %.4f, Hough axis size %d \n", i, barrel_HalfClusterU[i]->getEnergy(), barrel_HalfClusterU[i]->getHalfClusterCol(settings.map_stringPars["OutputLongiClusName"]).size()  );
  if(barrel_HalfClusterU[i]->getHalfClusterCol(settings.map_stringPars["OutputLongiClusName"]).size()!=0){
    Ncl_hough++;
    Naxis_hough += barrel_HalfClusterU[i]->getHalfClusterCol(settings.map_stringPars["OutputLongiClusName"]).size();
    Etot_hough += barrel_HalfClusterU[i]->getEnergy(); 
  }

  if(barrel_HalfClusterU[i]->getHalfClusterCol("TrackAxis").size()!=0){
    Ncl_trk++;
    Naxis_trk += barrel_HalfClusterU[i]->getHalfClusterCol("TrackAxis").size();
    Etot_trk += barrel_HalfClusterU[i]->getEnergy();
  }
}
cout<<"ClusterU number with Hough axis: "<<Ncl_hough<<", total Hough axis size: "<<Naxis_hough<<", total energy of cluster with Hough axis: "<<Etot_hough<<endl;
cout<<"ClusterU number with track axis: "<<Ncl_trk<<", total Track axis size: "<<Naxis_trk<<", total energy of cluster with Track axis: "<<Etot_trk<<endl;
Ncl_hough = 0;
Ncl_trk = 0;
Naxis_hough = 0;
Naxis_trk = 0;
Etot_hough = 0.;
Etot_trk = 0.;
for(int i=0; i<barrel_HalfClusterV.size(); i++){
  //printf("HalfClusterV #%d: energy %.4f, Hough axis size %d \n", i, barrel_HalfClusterV[i]->getEnergy(), barrel_HalfClusterV[i]->getHalfClusterCol(settings.map_stringPars["OutputLongiClusName"]).size()  );
  if(barrel_HalfClusterV[i]->getHalfClusterCol(settings.map_stringPars["OutputLongiClusName"]).size()!=0){
    Ncl_hough++;
    Naxis_hough += barrel_HalfClusterV[i]->getHalfClusterCol(settings.map_stringPars["OutputLongiClusName"]).size();
    Etot_hough += barrel_HalfClusterV[i]->getEnergy();
  }

  if(barrel_HalfClusterV[i]->getHalfClusterCol("TrackAxis").size()!=0){
    Ncl_trk++;
    Naxis_trk += barrel_HalfClusterV[i]->getHalfClusterCol("TrackAxis").size();
    Etot_trk += barrel_HalfClusterV[i]->getEnergy();
  }
}
cout<<"ClusterV number with Hough axis: "<<Ncl_hough<<", total Hough axis size: "<<Naxis_hough<<", total energy of cluster with Hough axis: "<<Etot_hough<<endl;
cout<<"ClusterV number with track axis: "<<Ncl_trk<<", total Track axis size: "<<Naxis_trk<<", total energy of cluster with Track axis: "<<Etot_trk<<endl;
*/
  return StatusCode::SUCCESS;
}

StatusCode HoughClusteringAlg::ClearAlgorithm(){
  barrel_HalfClusterV.clear();
  barrel_HalfClusterU.clear(); 
  endcap0_HalfClusterV.clear();
  endcap0_HalfClusterU.clear();
  endcap1_HalfClusterV.clear();
  endcap1_HalfClusterU.clear();
  m_localMaxVCol.clear();
  m_localMaxUCol.clear(); 
  m_longiClusVCol.clear();
  m_longiClusUCol.clear();

  return StatusCode::SUCCESS;
};


StatusCode HoughClusteringAlg::HoughTransformation(std::vector<Cyber::HoughObject>& Hobjects){
  if(Hobjects.size()<settings.map_intPars["th_peak"]) return StatusCode::SUCCESS;

  for(int iobj=0; iobj<Hobjects.size(); iobj++){
    int t_slayer = Hobjects[iobj].getSlayer();
    int t_system = Hobjects[iobj].getSystem();
    double point_Phi = Hobjects[iobj].getCenterPoint().Phi();
    double alpha_min, alpha_max;

    // Set range of alpha for Hough band. The range is different for different systems(Barrel/Endcap) and slayers(U/V)
    // Barrel ECAL
    if(t_system==Cyber::CaloUnit::System_Barrel){
      // U plane
      if(t_slayer==0){
        alpha_min = 0;
        alpha_max = TMath::Pi();
      }
      // V plane
      else{
        if( point_Phi < 5*TMath::Pi()/8. && point_Phi >= TMath::PiOver2() ){
          alpha_min = -0.1;
          alpha_max = TMath::PiOver4();
        }
        else if(point_Phi>3*TMath::Pi()/8. && point_Phi<TMath::PiOver2()){
          alpha_min = 7.*TMath::Pi()/4.;
          alpha_max = 2*TMath::Pi();
        }
        else{
          alpha_min = floor(4*point_Phi/TMath::Pi() - 1.5)*TMath::PiOver4() - TMath::PiOver4();
          alpha_max = floor(4*point_Phi/TMath::Pi() - 1.5)*TMath::PiOver4() + TMath::PiOver4();
        }

        if(alpha_min<=0 && alpha_max<=0){
          alpha_min += 2*TMath::Pi();
          alpha_max += 2*TMath::Pi();
        }
      }
    }
    // Endcap ECAL
    else if(t_system==Cyber::CaloUnit::System_Endcap){
      alpha_min = 0;
      alpha_max = TMath::Pi();
    }
    else{
      cout << "  HoughClusteringAlg: Unknown system ID: " << t_system << endl;
      return StatusCode::FAILURE;
    }

    // Create Hough lines. The range between
    TF1 line1("line1", "[0]*cos(x)+[1]*sin(x)", alpha_min, alpha_max);
    TF1 line2("line2", "[0]*cos(x)+[1]*sin(x)", alpha_min, alpha_max);

    // Set the parameters of the lines.
    line1.SetParameters( Hobjects[iobj].getUpperPoint().X(), Hobjects[iobj].getUpperPoint().Y() );
    line2.SetParameters( Hobjects[iobj].getLowerPoint().X(), Hobjects[iobj].getLowerPoint().Y() );
    
    Hobjects[iobj].setHoughLine(line1, line2);  

    //cout << "  yyy: HoughTransformation: HoughObject " << iobj << " has been transformed." << endl;
    //cout << "       alpha_min: " << alpha_min << " alpha_max: " << alpha_max << endl;
    //cout << "       line1: " << line1.GetParameter(0) << " " << line1.GetParameter(1) << endl;
    //cout << "       line2: " << line2.GetParameter(0) << " " << line2.GetParameter(1) << endl;
  } 

  return StatusCode::SUCCESS;
} 

StatusCode HoughClusteringAlg::FillHoughSpace(vector<Cyber::HoughObject>& Hobjects, Cyber::HoughSpace& Hspace){  
  for(int ih=0; ih<Hobjects.size(); ih++){
    TF1 line1 = Hobjects[ih].getHoughLine1();
    TF1 line2 = Hobjects[ih].getHoughLine2();

    double range_min, range_max;
    line1.GetRange(range_min, range_max);
    
    // Get bin num in alpha axis
    int bin_min = Hspace.getAlphaBin(range_min);
    int bin_max = Hspace.getAlphaBin(range_max);

    // Loop for alpha bins, line1 and line2
    for(int ialpha=bin_min; ialpha<=bin_max; ialpha++) {
      // The lines should be monotone at this range
      double line1_rho1 = line1.Eval( Hspace.getAlphaBinLowEdge(ialpha) );
      double line1_rho2 = line1.Eval( Hspace.getAlphaBinUpEdge(ialpha)  );
      double line2_rho1 = line2.Eval( Hspace.getAlphaBinLowEdge(ialpha) );
      double line2_rho2 = line2.Eval( Hspace.getAlphaBinUpEdge(ialpha)  );

      double line1_rho_min = TMath::Min(line1_rho1, line1_rho2);
      double line1_rho_max = TMath::Max(line1_rho1, line1_rho2);
      double line2_rho_min = TMath::Min(line2_rho1, line2_rho2);
      double line2_rho_max = TMath::Max(line2_rho1, line2_rho2);

      if(line1_rho_min>line1_rho_max || line2_rho_min>line2_rho_max){
        cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
      }

      double rho_min = TMath::Min( line1_rho_min, line2_rho_min);
      double rho_max = TMath::Max( line1_rho_max, line2_rho_max);

      if(rho_max<settings.map_floatPars["rho_low"] || rho_min>settings.map_floatPars["rho_high"]) continue;

      int nbin_rho_min = TMath::Max( int(ceil( (rho_min-settings.map_floatPars["rho_low"]) / settings.map_floatPars["bin_width_rho"] )), 1 );
      int nbin_rho_max = TMath::Min( int(ceil( (rho_max-settings.map_floatPars["rho_low"]) / settings.map_floatPars["bin_width_rho"] )), settings.map_intPars["Nbins_rho"] );


      for(int irho=nbin_rho_min; irho<=nbin_rho_max; irho++){
        Hspace.AddBinHobj(ialpha, irho, ih);
      }
    }  

  } 

  return StatusCode::SUCCESS;
} 


StatusCode HoughClusteringAlg::ClusterFinding(vector<Cyber::HoughObject>& Hobjects, Cyber::HoughSpace& Hspace, 
                                              std::vector<std::shared_ptr<Cyber::CaloHalfCluster>>& m_longiClusCol ){
  if(Hobjects.size()==0) return StatusCode::SUCCESS;

  map< pair<int, int>, set<int> > Hough_bins = Hspace.getHoughBins();
  // transform candidate to longicluster
  vector<std::shared_ptr<Cyber::CaloHalfCluster>> m_clusCol; m_clusCol.clear();
  for(auto ihb : Hough_bins){
    if(ihb.second.size()<settings.map_intPars["th_peak"]) continue;
    
    std::shared_ptr<Cyber::CaloHalfCluster> m_clus = std::make_shared<Cyber::CaloHalfCluster>();
    for(auto it = (ihb.second).begin(); it!=(ihb.second).end(); it++){
      m_clus->addUnit(Hobjects[*it].getLocalMax());
    }

    double t_alpha = Hspace.getAlphaBinCenter((ihb.first).first);
    double t_rho = Hspace.getRhoBinCenter((ihb.first).second);
    m_clus->setHoughPars(t_alpha, t_rho);

    if( !m_clus->isContinueN(settings.map_intPars["th_continueN"]) ){
      continue;
    } 
    m_clus->setType(100); //EM-type axis.     
    m_clus->getLinkedMCPfromUnit();
    m_clusCol.push_back(m_clus);
  }
  // Clean cluster
  //CleanClusters(m_clusCol);

  //bk_HFclus.insert( bk_HFclus.end(), m_clusCol.begin(), m_clusCol.end() );
  m_longiClusCol.insert( m_longiClusCol.end(), m_clusCol.begin(), m_clusCol.end() );
  m_clusCol.clear();

  return StatusCode::SUCCESS;
}  // ClusterFinding() end


StatusCode HoughClusteringAlg::CleanClusters( std::vector<std::shared_ptr<Cyber::CaloHalfCluster>>& m_longiClusCol){

  if(m_longiClusCol.size()==0)  return StatusCode::SUCCESS;

  // Remove repeated tracks
  for(int ic=0; ic<m_longiClusCol.size(); ic++){
  for(int jc=0; jc<m_longiClusCol.size(); jc++){
    if(ic>=m_longiClusCol.size()) ic--;
    if(ic==jc) continue;

    if( m_longiClusCol[ic].get()->isSubset(m_longiClusCol[jc].get()) ){  //jc is the subset of ic. remove jc. 
      //delete m_longiClusCol[jc]; m_longiClusCol[jc] = NULL;
      m_longiClusCol.erase(m_longiClusCol.begin()+jc );
      jc--;
      if(ic>jc+1) ic--;
    }
  }}

  //Depart the HoughCluster to 2 sub-clusters if it has blank in middle.
  for(int ic=0; ic<m_longiClusCol.size(); ic++){
    int m_nhit = m_longiClusCol[ic].get()->getCluster().size(); 
    m_longiClusCol[ic].get()->sortBarShowersByLayer();

    for(int ih=0; ih<m_nhit-1; ih++){
      if(m_longiClusCol[ic].get()->getCluster()[ih+1]->getDlayer() - m_longiClusCol[ic].get()->getCluster()[ih]->getDlayer() > 2){
        //Cyber::CaloHalfCluster* clus_head = new Cyber::CaloHalfCluster();
        //Cyber::CaloHalfCluster* clus_tail = new Cyber::CaloHalfCluster();
        std::shared_ptr<Cyber::CaloHalfCluster> clus_head = std::make_shared<Cyber::CaloHalfCluster>();
        std::shared_ptr<Cyber::CaloHalfCluster> clus_tail = std::make_shared<Cyber::CaloHalfCluster>();

        for(int jh=0; jh<=ih; jh++) 
          clus_head->addUnit( m_longiClusCol[ic].get()->getCluster()[jh]);
        for(int jh=ih+1; jh<m_nhit; jh++) 
          clus_tail->addUnit( m_longiClusCol[ic].get()->getCluster()[jh]);
        
        if( clus_head->isContinueN(settings.map_intPars["th_continueN"]) ) {
            clus_head->setType(100); 
            clus_head->setHoughPars(m_longiClusCol[ic].get()->getHoughAlpha(), m_longiClusCol[ic].get()->getHoughRho());
            clus_head->getLinkedMCPfromUnit();
            m_longiClusCol.push_back(clus_head);
        }
        //else{
        //  delete clus_head;
        //}
        if( clus_tail->isContinueN(settings.map_intPars["th_continueN"]) ) {
          clus_tail->setHoughPars(m_longiClusCol[ic].get()->getHoughAlpha(), m_longiClusCol[ic].get()->getHoughRho());
          clus_tail->setType(100);
          clus_tail->getLinkedMCPfromUnit();
          m_longiClusCol.push_back(clus_tail);
        }
        //else{
        //  delete clus_tail;
        //}
            

        //delete m_longiClusCol[ic]; m_longiClusCol[ic] = NULL;
        m_longiClusCol.erase(m_longiClusCol.begin()+ic);
        ic--;
        break;
      }
    }
  }

  // Remove repeated axes
  for(int ic=0; ic<m_longiClusCol.size(); ic++){
  for(int jc=0; jc<m_longiClusCol.size(); jc++){
    if(ic>=m_longiClusCol.size()) ic--;
    if(ic==jc) continue;

    if( m_longiClusCol[ic].get()->isSubset(m_longiClusCol[jc].get()) ){  //jc is the subset of ic. remove jc. 
      //delete m_longiClusCol[jc]; m_longiClusCol[jc] = NULL;
      m_longiClusCol.erase(m_longiClusCol.begin()+jc );
      jc--;
      if(ic>jc+1) ic--;
    }
  }}

  // Cut energy
  for(int ic=0; ic<m_longiClusCol.size(); ic++){
    if(m_longiClusCol[ic].get()->getEnergy()<settings.map_floatPars["th_AxisE"]){
      //delete m_longiClusCol[ic]; m_longiClusCol[ic]=NULL;
      m_longiClusCol.erase(m_longiClusCol.begin()+ic );
      ic--;
    }
  }

  // Overlap with other clusters: 
  if(m_longiClusCol.size()>=2){
    for(int ic=0; ic<m_longiClusCol.size()-1; ic++){
    for(int jc=ic+1; jc<m_longiClusCol.size(); jc++){
      if(ic>=m_longiClusCol.size()) ic--;

      double delta_alpha = TMath::Abs(m_longiClusCol[ic].get()->getHoughAlpha() -  m_longiClusCol[jc].get()->getHoughAlpha());
      if( (delta_alpha > settings.map_floatPars["th_dAlpha1"])
          && (delta_alpha < 2*TMath::Pi()-settings.map_floatPars["th_dAlpha1"]) ) continue;

      double m_ratio1 = m_longiClusCol[ic].get()->OverlapRatioE(m_longiClusCol[jc].get());
      double m_ratio2 = m_longiClusCol[jc].get()->OverlapRatioE(m_longiClusCol[ic].get());

      if(m_ratio1>settings.map_floatPars["th_overlapE"] && m_longiClusCol[ic].get()->getEnergy()<m_longiClusCol[jc].get()->getEnergy()){
        //delete m_longiClusCol[ic]; m_longiClusCol[ic] = NULL;
        m_longiClusCol.erase( m_longiClusCol.begin()+ic );
        ic--;
        break;
      }

      if(m_ratio2>settings.map_floatPars["th_overlapE"] && m_longiClusCol[jc].get()->getEnergy()<m_longiClusCol[ic].get()->getEnergy()){
        //delete m_longiClusCol[jc]; m_longiClusCol[jc] = NULL;
        m_longiClusCol.erase( m_longiClusCol.begin()+jc );
        jc--;
      }
    }}
  }

  // If two cluster are close to each other, and E_small/E_large < threshold, delete the small ones
  if(m_longiClusCol.size()>=2){
    for(int ic=0; ic<m_longiClusCol.size()-1; ic++){
    for(int jc=ic+1; jc<m_longiClusCol.size(); jc++){
      if(ic>=m_longiClusCol.size()) ic--;

      double delta_alpha = TMath::Abs(m_longiClusCol[ic].get()->getHoughAlpha() -  m_longiClusCol[jc].get()->getHoughAlpha());
      if( delta_alpha > settings.map_floatPars["th_dAlpha2"] ) continue;

      double E_ratio1 = m_longiClusCol[ic].get()->getEnergy() / m_longiClusCol[jc].get()->getEnergy();
      double E_ratio2 = m_longiClusCol[jc].get()->getEnergy() / m_longiClusCol[ic].get()->getEnergy();

      if( E_ratio1 < settings.map_floatPars["th_ERatio"] ){
        //delete m_longiClusCol[ic]; m_longiClusCol[ic] = NULL;
        m_longiClusCol.erase( m_longiClusCol.begin()+ic );
        ic--;
        break;
      }
      else if( E_ratio2 < settings.map_floatPars["th_ERatio"] ){
        //delete m_longiClusCol[jc]; m_longiClusCol[jc] = NULL;
        m_longiClusCol.erase( m_longiClusCol.begin()+jc );
        jc--;
      }
    }}
  }

  return StatusCode::SUCCESS;
}

#endif
