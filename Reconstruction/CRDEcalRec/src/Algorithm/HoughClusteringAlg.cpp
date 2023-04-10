#ifndef HOUGHCLUSTERINGALG_C
#define HOUGHCLUSTERINGALG_C

#include "Algorithm/HoughClusteringAlg.h"
#include <algorithm>
#include <set>
#include "TCanvas.h"
using namespace std;

StatusCode HoughClusteringAlg::ReadSettings(PandoraPlus::Settings& m_settings){
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
    settings.map_intPars["th_Layers"] = 10;
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

StatusCode HoughClusteringAlg::Initialize( PandoraPlusDataCol& m_datacol ){
  p_HalfClusterU.clear(); 
  p_HalfClusterV.clear(); 


  p_HalfClusterU = m_datacol.map_HalfCluster["HalfClusterColU"];
  p_HalfClusterV = m_datacol.map_HalfCluster["HalfClusterColV"];

  return StatusCode::SUCCESS;
}

StatusCode HoughClusteringAlg::RunAlgorithm( PandoraPlusDataCol& m_datacol ){

  //if( (p_HalfClusterU.size()+p_HalfClusterV.size())<1 ){
  //  std::cout << "HoughClusteringAlg: No HalfCluster input"<<std::endl;
  //  return StatusCode::SUCCESS;
  //}

  
  if(p_HalfClusterV.size()==0){ std::cout<<"  HoughClusteringAlg: No HalfClusterV in present data collection! "<<std::endl; }
  if(p_HalfClusterU.size()==0){ std::cout<<"  HoughClusteringAlg: No HalfClusterU in present data collection! "<<std::endl; }

  // Processing V(xy) plane
  for(int it=0; it<p_HalfClusterV.size(); it++){ // process each HalfCluster respectively
    m_localMaxVCol.clear();
    std::vector<const PandoraPlus::Calo1DCluster*> tmp_localMaxVCol = p_HalfClusterV[it]->getLocalMaxCol(settings.map_stringPars["ReadinLocalMaxName"]);

    for(int il=0; il<tmp_localMaxVCol.size(); il++){
      if(tmp_localMaxVCol[il]->getDlayer()<=settings.map_intPars["th_Layers"]) 
        m_localMaxVCol.push_back(tmp_localMaxVCol[il]);
    }

    if(m_localMaxVCol.size()<settings.map_intPars["th_peak"]){
      //std::cout << "    yyy: m_localMaxVCol.size()<th_peak, continue" << std::endl;
      continue; 
    } 

    // cout<<"  HoughClusteringAlg: Local maximum size V = "<<m_localMaxVCol.size()<<endl;
    
    // cout<<"  HoughClusteringAlg: Creating m_HoughObjectsV"<<endl;
    std::vector<PandoraPlus::HoughObject> m_HoughObjectsV; m_HoughObjectsV.clear(); 
    for(int il=0; il<m_localMaxVCol.size(); il++){
      PandoraPlus::HoughObject m_obj(m_localMaxVCol[il], PandoraPlus::CaloUnit::barsize, PandoraPlus::CaloUnit::ecal_innerR);
      m_HoughObjectsV.push_back(m_obj);
    }

    // cout<<"  HoughClusteringAlg: Hough transformation"<<endl;
    HoughTransformation(m_HoughObjectsV);

    // cout<<"  HoughClusteringAlg: Creating hough_spaceV"<<endl;
    PandoraPlus::HoughSpace hough_spaceV(settings.map_floatPars["alpha_lowV"], settings.map_floatPars["alpha_highV"], 
                                         settings.map_floatPars["bin_width_alphaV"], settings.map_intPars["Nbins_alphaV"], 
                                         settings.map_floatPars["rho_low"], settings.map_floatPars["rho_high"], 
                                         settings.map_floatPars["bin_width_rho"], settings.map_intPars["Nbins_rho"]);

    // cout<<"  HoughClusteringAlg: Filling hough_spaceV"<<endl;
    FillHoughSpace(m_HoughObjectsV, hough_spaceV);

    // cout<<"  HoughClusteringAlg: Finding clusters from Hough space"<<endl;

    //Create output HoughClusters
    m_longiClusVCol.clear(); 
    ClusterFinding(m_HoughObjectsV, hough_spaceV, m_longiClusVCol);

    cout << "  HoughClusteringAlg: final output m_longiClusVCol.size() = " << m_longiClusVCol.size() << endl;

    for(int ic=0; ic<m_longiClusVCol.size(); ic++)
      m_datacol.bk_ClusterHalfCol.push_back( const_cast<PandoraPlus::CaloHalfCluster*>(m_longiClusVCol[ic]) );

    std::vector<const PandoraPlus::Calo1DCluster*> left_localMaxVCol; left_localMaxVCol.clear();
    std::vector<const PandoraPlus::Calo1DCluster*> m_houghMax; m_houghMax.clear(); 
    for(int is=0; is<tmp_localMaxVCol.size(); is++){
      bool fl_incluster = false; 
      for(int ic=0; ic<m_longiClusVCol.size(); ic++){
        std::vector<const PandoraPlus::Calo1DCluster*> p_showers = m_longiClusVCol[ic]->getCluster();
        if( find(p_showers.begin(), p_showers.end(), tmp_localMaxVCol[is])!=p_showers.end() ) { fl_incluster = true; break; }
      }
      if(!fl_incluster && find(left_localMaxVCol.begin(), left_localMaxVCol.end(), tmp_localMaxVCol[is])==left_localMaxVCol.end() ) left_localMaxVCol.push_back(tmp_localMaxVCol[is]);
      m_houghMax.push_back( tmp_localMaxVCol[is] );
    }

    p_HalfClusterV[it]->setLocalMax("HoughLocalMax", m_houghMax);
    p_HalfClusterV[it]->setLocalMax(settings.map_stringPars["LeftLocalMaxName"], left_localMaxVCol);
    p_HalfClusterV[it]->setHalfClusters(settings.map_stringPars["OutputLongiClusName"], m_longiClusVCol);
    m_houghMax.clear();
    left_localMaxVCol.clear();

  }  // end of V plane

  // Processing U(r-phi) plane
  for(int it=0; it<p_HalfClusterU.size(); it++){
    m_localMaxUCol.clear();
    std::vector<const PandoraPlus::Calo1DCluster*> tmp_localMaxUCol = p_HalfClusterU[it]->getLocalMaxCol(settings.map_stringPars["ReadinLocalMaxName"]);

    for(int il=0; il<tmp_localMaxUCol.size(); il++){
      if(tmp_localMaxUCol[il]->getDlayer()<=settings.map_intPars["th_Layers"]) 
        m_localMaxUCol.push_back(tmp_localMaxUCol[il]);
    }

    if(m_localMaxUCol.size()<settings.map_intPars["th_peak"]){
//std::cout << "    yyy: m_localMaxUCol.size()<th_peak, continue" << std::endl;
      continue; 
    } 
//cout << "  HoughClusteringAlg: Local maximum size U = " << m_localMaxUCol.size() << endl;

//cout << "  HoughClusteringAlg: Creating m_HoughObjectsU" << endl;
    std::vector<PandoraPlus::HoughObject> m_HoughObjectsU; m_HoughObjectsU.clear(); 
    for(int il=0; il<m_localMaxUCol.size(); il++){
      PandoraPlus::HoughObject m_obj(m_localMaxUCol[il], PandoraPlus::CaloUnit::barsize, PandoraPlus::CaloUnit::ecal_innerR);
      m_HoughObjectsU.push_back(m_obj);
    }

    // cout<<"  HoughClusteringAlg: Hough transformation"<<endl;
    HoughTransformation(m_HoughObjectsU);

    // cout<<"  HoughClusteringAlg: Divide m_HoughObjectsU module by module"<<endl;
    set<int> modules;
    for(int il=0; il<m_localMaxUCol.size(); il++){
      modules.insert((m_localMaxUCol[il]->getTowerID())[0][0]);
    }
// cout << "  yyy: modules.size() = " << modules.size() << endl;
    std::vector< std::vector<PandoraPlus::HoughObject> > m_HoughObjectsU_modules; m_HoughObjectsU_modules.clear();
    for(auto im=modules.begin(); im!=modules.end(); im++){
      std::vector<PandoraPlus::HoughObject> hobj;
      for(int ih=0; ih<m_HoughObjectsU.size(); ih++){
        if(m_HoughObjectsU[ih].getModule()==*im)
          hobj.push_back(m_HoughObjectsU[ih]);
      }
      if (hobj.size() < settings.map_intPars["th_peak"]) continue;
      m_HoughObjectsU_modules.push_back(hobj);
    }

// cout << "  HoughClusteringAlg: Creating hough_spacesU module by module" << endl;
    std::vector<PandoraPlus::HoughSpace> hough_spacesU;
    for(int im=0; im<m_HoughObjectsU_modules.size(); im++){
      PandoraPlus::HoughSpace hspaceU(settings.map_floatPars["alpha_lowU"], settings.map_floatPars["alpha_highU"], 
                                         settings.map_floatPars["bin_width_alphaU"], settings.map_intPars["Nbins_alphaU"], 
                                         settings.map_floatPars["rho_low"], settings.map_floatPars["rho_high"], 
                                         settings.map_floatPars["bin_width_rho"], settings.map_intPars["Nbins_rho"]);
      FillHoughSpace(m_HoughObjectsU_modules[im], hspaceU);
      hough_spacesU.push_back(hspaceU);
    }

// cout<<"  HoughClusteringAlg: Finding clusters from Hough spaces"<<endl;
    m_longiClusUCol.clear(); 
    for(int ih=0; ih<hough_spacesU.size(); ih++){
      ClusterFinding(m_HoughObjectsU_modules[ih], hough_spacesU[ih], m_longiClusUCol);
    }

 cout << "  HoughClusteringAlg: final output m_longiClusUCol.size() = " << m_longiClusUCol.size() << endl;
    for(int ic=0; ic<m_longiClusUCol.size(); ic++)
      m_datacol.bk_ClusterHalfCol.push_back( const_cast<PandoraPlus::CaloHalfCluster*>(m_longiClusUCol[ic]) );
    
    std::vector<const PandoraPlus::Calo1DCluster*> left_localMaxUCol; left_localMaxUCol.clear();
    std::vector<const PandoraPlus::Calo1DCluster*> m_houghMax; m_houghMax.clear();
    for(int is=0; is<tmp_localMaxUCol.size(); is++){
      bool fl_incluster = false; 
      for(int ic=0; ic<m_longiClusUCol.size(); ic++){
        std::vector<const PandoraPlus::Calo1DCluster*> p_showers = m_longiClusUCol[ic]->getCluster();
        if( find(p_showers.begin(), p_showers.end(), tmp_localMaxUCol[is])!=p_showers.end() ) { fl_incluster = true; break; }
      }
      if(!fl_incluster && find(left_localMaxUCol.begin(), left_localMaxUCol.end(), tmp_localMaxUCol[is])==left_localMaxUCol.end() ) left_localMaxUCol.push_back(tmp_localMaxUCol[is]);
      else m_houghMax.push_back( tmp_localMaxUCol[is] );
    }

    p_HalfClusterU[it]->setLocalMax("HoughLocalMax", m_houghMax);
    p_HalfClusterU[it]->setLocalMax(settings.map_stringPars["LeftLocalMaxName"], left_localMaxUCol);
    p_HalfClusterU[it]->setHalfClusters(settings.map_stringPars["OutputLongiClusName"], m_longiClusUCol);
    m_houghMax.clear();
    left_localMaxUCol.clear(); 

  }  // end of U plane


  return StatusCode::SUCCESS;
}

StatusCode HoughClusteringAlg::ClearAlgorithm(){
  p_HalfClusterV.clear();
  p_HalfClusterU.clear(); 
  m_localMaxVCol.clear();
  m_localMaxUCol.clear(); 
  m_longiClusVCol.clear();
  m_longiClusUCol.clear();

  return StatusCode::SUCCESS;
};


StatusCode HoughClusteringAlg::HoughTransformation(std::vector<PandoraPlus::HoughObject>& Hobjects){
  if(Hobjects.size()<settings.map_intPars["th_peak"]) return StatusCode::SUCCESS;

  // range of alpha of different lines
  double range12[2] = {0, 0};
  double range34[2] = {0, 0};

  for(int iobj=0; iobj<Hobjects.size(); iobj++){
    int t_module = Hobjects[iobj].getModule();
    int t_slayer = Hobjects[iobj].getSlayer();
    SetLineRange(t_module, t_slayer, range12, range34);

    TF1 line1("line1", "[0]*cos(x)+[1]*sin(x)", range12[0], range12[1]);
    TF1 line2("line2", "[0]*cos(x)+[1]*sin(x)", range12[0], range12[1]);
    TF1 line3("line3", "[0]*cos(x)+[1]*sin(x)", range34[0], range34[1]);
    TF1 line4("line4", "[0]*cos(x)+[1]*sin(x)", range34[0], range34[1]);

    if(t_slayer==0){
      line1.SetParameters( Hobjects[iobj].getPointUR().X(), Hobjects[iobj].getPointUR().Y() );
      line2.SetParameters( Hobjects[iobj].getPointDL().X(), Hobjects[iobj].getPointDL().Y() );
      line3.SetParameters( Hobjects[iobj].getPointUL().X(), Hobjects[iobj].getPointUL().Y() );
      line4.SetParameters( Hobjects[iobj].getPointDR().X(), Hobjects[iobj].getPointDR().Y() );
    }
    else if(t_slayer==1){
      if(t_module % 2 == 0){
        line1.SetParameters( Hobjects[iobj].getPointUR().X(), Hobjects[iobj].getPointUR().Y() );
        line2.SetParameters( Hobjects[iobj].getPointDL().X(), Hobjects[iobj].getPointDL().Y() );
        line3.SetParameters( Hobjects[iobj].getPointUL().X(), Hobjects[iobj].getPointUL().Y() );
        line4.SetParameters( Hobjects[iobj].getPointDR().X(), Hobjects[iobj].getPointDR().Y() );
      }else{
        line1.SetParameters( Hobjects[iobj].getPointU().X(), Hobjects[iobj].getPointU().Y() );
        line2.SetParameters( Hobjects[iobj].getPointD().X(), Hobjects[iobj].getPointD().Y() );
        line3.SetParameters( Hobjects[iobj].getPointL().X(), Hobjects[iobj].getPointL().Y() );
        line4.SetParameters( Hobjects[iobj].getPointR().X(), Hobjects[iobj].getPointR().Y() );
      }
    }
    
    Hobjects[iobj].setHoughLine(t_module, line1, line2, line3, line4);  

  }

  return StatusCode::SUCCESS;
}  // HoughTransformation() end


StatusCode HoughClusteringAlg::SetLineRange(int module, int slayer, double *range12, double* range34){
  // range12: ur, dl, u, d
  // range34: ul, dr, l, r
  if(slayer == 0){
    range12[0] = 0.;
    range12[1] = TMath::Pi()/2.;
    range34[0] = TMath::Pi()/2.;
    range34[1] = TMath::Pi();
  }
  else if(slayer == 1){
    switch(module){
      case 0:{
        range12[0] = -0.1;
        range12[1] = TMath::Pi()/4.;
        range34[0] = 7.*TMath::Pi()/4.;
        range34[1] = range34[0] + TMath::Pi()/4.;
        break;
      }
      case 1:{
        range12[0] = TMath::Pi()/4.;
        range12[1] = range12[0] + TMath::Pi()/4.;
        range34[0] = 0;
        range34[1] = range34[0] + TMath::Pi()/4.;
        break;
      }
      case 2:{
        range12[0] = TMath::Pi()/4.;
        range12[1] = range12[0] + TMath::Pi()/4.;
        range34[0] = TMath::Pi()/2;
        range34[1] = range34[0] + TMath::Pi()/4.;
        break;
      }
      case 3:{
        range12[0] = TMath::Pi()/2.;
        range12[1] = range12[0] + TMath::Pi()/4.;
        range34[0] = 3.*TMath::Pi()/4.;
        range34[1] = range34[0] + TMath::Pi()/4.;
        break;
      }
      case 4:{
        range12[0] = TMath::Pi();
        range12[1] = range12[0] + TMath::Pi()/4.;
        range34[0] = 3.*TMath::Pi()/4.;
        range34[1] = range34[0] + TMath::Pi()/4.;
        break;
      }
      case 5:{
        range12[0] = 5.*TMath::Pi()/4.;
        range12[1] = range12[0] + TMath::Pi()/4.;
        range34[0] = TMath::Pi();
        range34[1] = range34[0] + TMath::Pi()/4.;
        break;
      }
      case 6:{
        range12[0] = 5.*TMath::Pi()/4.;
        range12[1] = range12[0] + TMath::Pi()/4.;
        range34[0] = 3.*TMath::Pi()/2.;
        range34[1] = range34[0] + TMath::Pi()/4.;
        break;
      }
      case 7:{
        range12[0] = 3.*TMath::Pi()/2.;
        range12[1] = range12[0] + TMath::Pi()/4.;
        range34[0] = 7.*TMath::Pi()/4.;
        range34[1] = range34[0] + TMath::Pi()/4.;
        break;
      }
      default:{
        cout << "Wrong module: module = " << module << endl;
      }

    }
  }

  return StatusCode::SUCCESS;
}  // SetLineRange() end


StatusCode HoughClusteringAlg::FillHoughSpace(vector<PandoraPlus::HoughObject>& Hobjects, PandoraPlus::HoughSpace& Hspace){  
  // Fill Hough space
  // Loop Hough objects
  for(int ih=0; ih<Hobjects.size(); ih++){
    TF1 line1 = Hobjects[ih].getHoughLine1();
    TF1 line2 = Hobjects[ih].getHoughLine2();
    TF1 line3 = Hobjects[ih].getHoughLine3();
    TF1 line4 = Hobjects[ih].getHoughLine4();

    // line1 and line2 share the same range in alpha, so does line3 and line4
    double range_12_min, range_12_max, range_34_min, range_34_max;
    line1.GetRange(range_12_min, range_12_max);
    line3.GetRange(range_34_min, range_34_max);
    
    // Get bin num in alpha axis
    int bin_12_min = Hspace.getAlphaBin(range_12_min);
    int bin_12_max = Hspace.getAlphaBin(range_12_max);
    int bin_34_min = Hspace.getAlphaBin(range_34_min);
    int bin_34_max = Hspace.getAlphaBin(range_34_max);
    if (bin_12_max == bin_34_min) bin_34_min ++;
    if (bin_34_max == bin_12_min) bin_12_min ++;


    // Loop for alpha bins, line1 and line2
    for(int ialpha=bin_12_min; ialpha<=bin_12_max; ialpha++) {
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
    }  // end loop alpha bin, line1 and line2
    // Loop for alpha bins, line3 and line4
    for(int ialpha=bin_34_min; ialpha<=bin_34_max; ialpha++) {
      // The lines should be monotone at this range
      double line3_rho1 = line3.Eval( Hspace.getAlphaBinLowEdge(ialpha) );
      double line3_rho2 = line3.Eval( Hspace.getAlphaBinUpEdge(ialpha)  );
      double line4_rho1 = line4.Eval( Hspace.getAlphaBinLowEdge(ialpha) );
      double line4_rho2 = line4.Eval( Hspace.getAlphaBinUpEdge(ialpha)  );

      double line3_rho_min = TMath::Min(line3_rho1, line3_rho2);
      double line3_rho_max = TMath::Max(line3_rho1, line3_rho2);;
      double line4_rho_min = TMath::Min(line4_rho1, line4_rho2);
      double line4_rho_max = TMath::Max(line4_rho1, line4_rho2);

      if(line3_rho_min>line3_rho_max || line4_rho_min>line4_rho_max){
        cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
      }

      double rho_min = TMath::Min( line3_rho_min, line4_rho_min);
      double rho_max = TMath::Max( line3_rho_max, line4_rho_max);

      if(rho_max<settings.map_floatPars["rho_low"] || rho_min>settings.map_floatPars["rho_high"]) continue;

      int nbin_rho_min = TMath::Max( int(ceil( (rho_min-settings.map_floatPars["rho_low"]) / settings.map_floatPars["bin_width_rho"] )), 1 );
      int nbin_rho_max = TMath::Min( int(ceil( (rho_max-settings.map_floatPars["rho_low"]) / settings.map_floatPars["bin_width_rho"] )), settings.map_intPars["Nbins_rho"] );

      for(int irho=nbin_rho_min; irho<=nbin_rho_max; irho++){
        Hspace.AddBinHobj(ialpha, irho, ih);
      }
    }  // end loop alpha bin, line3 and line4

  }  // End loop Hough objects

  return StatusCode::SUCCESS;
}  // FillHoughSpace() end


StatusCode HoughClusteringAlg::ClusterFinding(vector<PandoraPlus::HoughObject>& Hobjects, PandoraPlus::HoughSpace& Hspace, 
                                              vector<const PandoraPlus::CaloHalfCluster*>& m_longiClusCol){
  if(Hobjects.size()==0) return StatusCode::SUCCESS;

  map< pair<int, int>, set<int> > Hough_bins = Hspace.getHoughBins();

  // transform candidate to longicluster
  vector<PandoraPlus::CaloHalfCluster*> m_clusCol; m_clusCol.clear();
  // for(auto ihb=Hough_bins.begin(); ihb!=Hough_bins.end(); ihb++){
    for(auto ihb : Hough_bins){
    if(ihb.second.size()<settings.map_intPars["th_peak"]) continue;
    
    PandoraPlus::CaloHalfCluster* m_clus = new PandoraPlus::CaloHalfCluster();
    for(auto it = (ihb.second).begin(); it!=(ihb.second).end(); it++){
      m_clus->addUnit(Hobjects[*it].getLocalMax());
    }

    double t_alpha = Hspace.getAlphaBinCenter((ihb.first).first);
    double t_rho = Hspace.getRhoBinCenter((ihb.first).second);
    m_clus->setHoughPars(t_alpha, t_rho);

    if( !m_clus->isContinueN(settings.map_intPars["th_continueN"]) ){
      delete m_clus;
      continue;
    } 
    m_clus->setType(1); //EM-type axis.     
    m_clusCol.push_back(m_clus);
  }

  // Clean cluster
  CleanClusters(m_clusCol);

  m_longiClusCol.insert( m_longiClusCol.end(), m_clusCol.begin(), m_clusCol.end() );

  return StatusCode::SUCCESS;
}  // ClusterFinding() end


StatusCode HoughClusteringAlg::CleanClusters(vector<PandoraPlus::CaloHalfCluster*>& m_longiClusCol){
  if(m_longiClusCol.size()==0)  return StatusCode::SUCCESS;

  // Remove repeated tracks
  for(int ic=0; ic<m_longiClusCol.size(); ic++){
  for(int jc=0; jc<m_longiClusCol.size(); jc++){
    if(ic>=m_longiClusCol.size()) ic--;
    if(ic==jc) continue;

    if( m_longiClusCol[ic]->isSubset(m_longiClusCol[jc]) ){  //jc is the subset of ic. remove jc. 
      delete m_longiClusCol[jc]; m_longiClusCol[jc] = NULL;
      m_longiClusCol.erase(m_longiClusCol.begin()+jc );
      jc--;
      if(ic>jc+1) ic--;
    }
  }}

  //Depart the HoughCluster to 2 sub-clusters if it has blank in middle.
  for(int ic=0; ic<m_longiClusCol.size(); ic++){
    int m_nhit = m_longiClusCol[ic]->getCluster().size(); 
    m_longiClusCol[ic]->sortBarShowersByLayer();

    for(int ih=0; ih<m_nhit-1; ih++){
      if(m_longiClusCol[ic]->getCluster()[ih+1]->getDlayer() - m_longiClusCol[ic]->getCluster()[ih]->getDlayer() > 2){
        PandoraPlus::CaloHalfCluster* clus_head = new PandoraPlus::CaloHalfCluster();
        PandoraPlus::CaloHalfCluster* clus_tail = new PandoraPlus::CaloHalfCluster();
        for(int jh=0; jh<=ih; jh++) 
          clus_head->addUnit( m_longiClusCol[ic]->getCluster()[jh]);
        for(int jh=ih+1; jh<m_nhit; jh++) 
          clus_tail->addUnit( m_longiClusCol[ic]->getCluster()[jh]);
        
        if( clus_head->isContinueN(settings.map_intPars["th_continueN"]) ) {
            clus_head->setType(1); 
            clus_head->setHoughPars(m_longiClusCol[ic]->getHoughAlpha(), m_longiClusCol[ic]->getHoughRho());
            m_longiClusCol.push_back(clus_head);
        }
        else{
          delete clus_head;
        }
        if( clus_tail->isContinueN(settings.map_intPars["th_continueN"]) ) {
          clus_tail->setHoughPars(m_longiClusCol[ic]->getHoughAlpha(), m_longiClusCol[ic]->getHoughRho());
          clus_tail->setType(1);
          m_longiClusCol.push_back(clus_tail);
        }
        else{
          delete clus_tail;
        }
            

        delete m_longiClusCol[ic]; m_longiClusCol[ic] = NULL;
        m_longiClusCol.erase(m_longiClusCol.begin()+ic);
        ic--;
        break;
      }
    }
  }

  // Remove repeated tracks
  for(int ic=0; ic<m_longiClusCol.size(); ic++){
  for(int jc=0; jc<m_longiClusCol.size(); jc++){
    if(ic>=m_longiClusCol.size()) ic--;
    if(ic==jc) continue;

    if( m_longiClusCol[ic]->isSubset(m_longiClusCol[jc]) ){  //jc is the subset of ic. remove jc. 
      delete m_longiClusCol[jc]; m_longiClusCol[jc] = NULL;
      m_longiClusCol.erase(m_longiClusCol.begin()+jc );
      jc--;
      if(ic>jc+1) ic--;
    }
  }}

  // Cut energy
  for(int ic=0; ic<m_longiClusCol.size(); ic++){
    if(m_longiClusCol[ic]->getEnergy()<settings.map_floatPars["th_AxisE"]){
      delete m_longiClusCol[ic]; m_longiClusCol[ic]=NULL;
      m_longiClusCol.erase(m_longiClusCol.begin()+ic );
      ic--;
    }
  }

  // Overlap with other clusters: 
    if(m_longiClusCol.size()>=2){
      for(int ic=0; ic<m_longiClusCol.size()-1; ic++){
      for(int jc=ic+1; jc<m_longiClusCol.size(); jc++){
        if(ic>=m_longiClusCol.size()) ic--;

        double delta_alpha = TMath::Abs(m_longiClusCol[ic]->getHoughAlpha() -  m_longiClusCol[jc]->getHoughAlpha());
        if( delta_alpha > settings.map_floatPars["th_dAlpha1"] ) continue;

        double m_ratio1 = m_longiClusCol[ic]->OverlapRatioE(m_longiClusCol[jc]);
        double m_ratio2 = m_longiClusCol[jc]->OverlapRatioE(m_longiClusCol[ic]);

        if(m_ratio1>settings.map_floatPars["th_overlapE"] && m_longiClusCol[ic]->getEnergy()<m_longiClusCol[jc]->getEnergy()){
          delete m_longiClusCol[ic]; m_longiClusCol[ic] = NULL;
          m_longiClusCol.erase( m_longiClusCol.begin()+ic );
          ic--;
          break;
        }

        if(m_ratio2>settings.map_floatPars["th_overlapE"] && m_longiClusCol[jc]->getEnergy()<m_longiClusCol[ic]->getEnergy()){
          delete m_longiClusCol[jc]; m_longiClusCol[jc] = NULL;
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

      double delta_alpha = TMath::Abs(m_longiClusCol[ic]->getHoughAlpha() -  m_longiClusCol[jc]->getHoughAlpha());
      if( delta_alpha > settings.map_floatPars["th_dAlpha2"] ) continue;

      double E_ratio1 = m_longiClusCol[ic]->getEnergy() / m_longiClusCol[jc]->getEnergy();
      double E_ratio2 = m_longiClusCol[jc]->getEnergy() / m_longiClusCol[ic]->getEnergy();

      if( E_ratio1 < settings.map_floatPars["th_ERatio"] ){
        delete m_longiClusCol[ic]; m_longiClusCol[ic] = NULL;
        m_longiClusCol.erase( m_longiClusCol.begin()+ic );
        ic--;
        break;
      }
      else if( E_ratio2 < settings.map_floatPars["th_ERatio"] ){
        delete m_longiClusCol[jc]; m_longiClusCol[jc] = NULL;
        m_longiClusCol.erase( m_longiClusCol.begin()+jc );
        jc--;
      }
    }}
  }

  return StatusCode::SUCCESS;
}  // CleanClusters() end

#endif
