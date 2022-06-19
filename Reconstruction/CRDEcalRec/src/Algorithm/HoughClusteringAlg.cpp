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
  if(settings.map_floatPars.find("HoughBinDivide")==settings.map_floatPars.end()) settings.map_floatPars["HoughBinDivide"] = 4;  
  if(settings.map_floatPars.find("th_continuetrkN")==settings.map_floatPars.end()) settings.map_floatPars["th_continuetrkN"] = 3;  
  if(settings.map_floatPars.find("th_AxisE")==settings.map_floatPars.end()) settings.map_floatPars["th_AxisE"] = 1.;  
  if(settings.map_floatPars.find("th_intercept")==settings.map_floatPars.end()) settings.map_floatPars["th_intercept"] = 300;  
  if(settings.map_floatPars.find("th_dAlpha1")==settings.map_floatPars.end()) settings.map_floatPars["th_dAlpha1"] = 0.3;  
  if(settings.map_floatPars.find("th_dAlpha2")==settings.map_floatPars.end()) settings.map_floatPars["th_dAlpha2"] = 0.8;  
  if(settings.map_floatPars.find("th_dRho1")==settings.map_floatPars.end()) settings.map_floatPars["th_dRho1"] = 20;  
  if(settings.map_floatPars.find("th_dRho2")==settings.map_floatPars.end()) settings.map_floatPars["th_dRho2"] = 50;  
  if(settings.map_floatPars.find("th_overlapE1")==settings.map_floatPars.end()) settings.map_floatPars["th_overlapE1"] = 0.7;  
  if(settings.map_floatPars.find("th_overlapE2")==settings.map_floatPars.end()) settings.map_floatPars["th_overlapE2"] = 0.9;  
  return StatusCode::SUCCESS;
};

StatusCode HoughClusteringAlg::Initialize(){

  return StatusCode::SUCCESS;
};

StatusCode HoughClusteringAlg::RunAlgorithm( PandoraPlusDataCol& m_datacol ){

  std::vector<PandoraPlus::CaloTower*>* p_towers = &(m_datacol.TowerCol); 
cout<<"Tower size: "<<p_towers->size()<<endl;
  for(int it=0; it<p_towers->size(); it++){
printf("  In tower #%d: block size = %d \n", it, p_towers->at(it)->getBlocks().size());

    std::vector<const PandoraPlus::CaloBarShower*> m_localMaxXCol; m_localMaxXCol.clear();
    std::vector<const PandoraPlus::CaloBarShower*> m_localMaxYCol; m_localMaxYCol.clear();
    for(int ib=0; ib<p_towers->at(it)->getBlocks().size(); ib++){
      if(settings.map_floatPars["th_Layers"]>0 && p_towers->at(it)->getBlocks()[ib]->getDlayer()>settings.map_floatPars["th_Layers"]) continue;
      std::vector<const PandoraPlus::CaloBarShower*> tmp_showerX = p_towers->at(it)->getBlocks()[ib]->getShowerXCol();
      std::vector<const PandoraPlus::CaloBarShower*> tmp_showerY = p_towers->at(it)->getBlocks()[ib]->getShowerYCol();

      m_localMaxXCol.insert(m_localMaxXCol.end(), tmp_showerX.begin(), tmp_showerX.end() );
      m_localMaxYCol.insert(m_localMaxYCol.end(), tmp_showerY.begin(), tmp_showerY.end() );
    }

    if(m_localMaxXCol.size()==0 && m_localMaxYCol.size()==0) continue; 
//cout<<"  Local maximum size: barX "<<m_localMaxXCol.size()<<"  barY "<<m_localMaxYCol.size()<<endl;

    std::vector<PandoraPlus::HoughObject> m_HoughObjectsX; m_HoughObjectsX.clear(); 
    std::vector<PandoraPlus::HoughObject> m_HoughObjectsY; m_HoughObjectsY.clear(); 
    for(int il=0; il<m_localMaxXCol.size(); il++){
      PandoraPlus::HoughObject m_obj; m_obj.Clear();
      m_obj.SetLocalMax( m_localMaxXCol[il] );
      m_HoughObjectsX.push_back(m_obj);
    }
    for(int il=0; il<m_localMaxYCol.size(); il++){
      PandoraPlus::HoughObject m_obj; m_obj.Clear();
      m_obj.SetLocalMax( m_localMaxYCol[il] );
      m_HoughObjectsY.push_back(m_obj);
    }

//cout<<"  HoughClusteringAlg: Conformal transformation"<<endl;
    //Conformal transformation  
    ConformalTransformation(m_HoughObjectsX);
    ConformalTransformation(m_HoughObjectsY);

/*
cout<<"    Local max (HoughObjectX): "<<endl;
for(int i=0; i<m_HoughObjectsX.size(); i++) printf("(%.2f, %.2f, %.2f) \t", m_HoughObjectsX[i].getLocalMax()->getPos().x(), m_HoughObjectsX[i].getLocalMax()->getPos().y(), m_HoughObjectsX[i].getLocalMax()->getPos().z() );
cout<<endl;
cout<<"    Conformal Point (HoughObjectX): "<<endl;
for(int i=0; i<m_HoughObjectsX.size(); i++) printf("(%.2f, %.2f) \t", m_HoughObjectsX[i].getConformPointUR().X(), m_HoughObjectsX[i].getConformPointUR().Y());
cout<<endl;

cout<<"    Local max (HoughObjectY): "<<endl;
for(int i=0; i<m_HoughObjectsY.size(); i++) printf("(%.2f, %.2f, %.2f) \t", m_HoughObjectsY[i].getLocalMax()->getPos().x(), m_HoughObjectsY[i].getLocalMax()->getPos().y(), m_HoughObjectsY[i].getLocalMax()->getPos().z() );
cout<<endl;
cout<<"    Conformal Point (HoughObjectY): "<<endl;
for(int i=0; i<m_HoughObjectsY.size(); i++) printf("(%.2f, %.2f) \t", m_HoughObjectsY[i].getConformPointUR().X(), m_HoughObjectsY[i].getConformPointUR().Y());
cout<<endl;
*/


//cout<<"  HoughClusteringAlg: Hough transformation"<<endl;

    //Hough transformation
    HoughTransformation(m_HoughObjectsX);
    HoughTransformation(m_HoughObjectsY);

/*
cout<<"  Check all Hough Objects: "<<endl;

cout<<"  Transformed Hough lines X: "<<endl;
for(int i=0; i<m_HoughObjectsX.size(); i++){
  cout<<"    HoughObj #"<<i<<": Origin localMax position ";
  printf(" (%.3f, %.3f, %.3f) \n", m_HoughObjectsX[i].getLocalMax()->getPos().X(), m_HoughObjectsX[i].getLocalMax()->getPos().Z(), m_HoughObjectsX[i].getLocalMax()->getE());
  printf("    Conformal Point: (%.3f, %.3f) \n", m_HoughObjectsX[i].getConformPointUR().X()-5, m_HoughObjectsX[i].getConformPointUR().Y()-5 );
  cout<<"    LineUR pars: "<<m_HoughObjectsX[i].getHoughLineUR().GetParameter(0)<<"  "<<m_HoughObjectsX[i].getHoughLineUR().GetParameter(1)<<", range: ["<<m_HoughObjectsX[i].getHoughLineUR().GetMaximum()<<", "<<m_HoughObjectsX[i].getHoughLineUR().GetMinimum()<<"]"<<endl;
  cout<<"    LineUL pars: "<<m_HoughObjectsX[i].getHoughLineUL().GetParameter(0)<<"  "<<m_HoughObjectsX[i].getHoughLineUL().GetParameter(1)<<", range: ["<<m_HoughObjectsX[i].getHoughLineUL().GetMaximum()<<", "<<m_HoughObjectsX[i].getHoughLineUL().GetMinimum()<<"]"<<endl;
  cout<<"    LineDR pars: "<<m_HoughObjectsX[i].getHoughLineDR().GetParameter(0)<<"  "<<m_HoughObjectsX[i].getHoughLineDR().GetParameter(1)<<", range: ["<<m_HoughObjectsX[i].getHoughLineDR().GetMaximum()<<", "<<m_HoughObjectsX[i].getHoughLineDR().GetMinimum()<<"]"<<endl;
  cout<<"    LineDL pars: "<<m_HoughObjectsX[i].getHoughLineDL().GetParameter(0)<<"  "<<m_HoughObjectsX[i].getHoughLineDL().GetParameter(1)<<", range: ["<<m_HoughObjectsX[i].getHoughLineDL().GetMaximum()<<", "<<m_HoughObjectsX[i].getHoughLineDL().GetMinimum()<<"]"<<endl;
  cout<<endl;
}


cout<<"  Transformed Hough lines Y: "<<endl;
for(int i=0; i<m_HoughObjectsY.size(); i++){
  cout<<"    HoughObj #"<<i<<": Origin localMax position ";
  printf(" (%.3f, %.3f, %.3f) \n", m_HoughObjectsY[i].getLocalMax()->getPos().X(), m_HoughObjectsY[i].getLocalMax()->getPos().Z(), m_HoughObjectsY[i].getLocalMax()->getE());
  printf("    Conformal Point: (%.3f, %.3f) \n", m_HoughObjectsY[i].getConformPointUR().X()-5, m_HoughObjectsY[i].getConformPointUR().Y()-5 );
  cout<<"    LineUR pars: "<<m_HoughObjectsY[i].getHoughLineUR().GetParameter(0)<<"  "<<m_HoughObjectsY[i].getHoughLineUR().GetParameter(1)<<", range: ["<<m_HoughObjectsY[i].getHoughLineUR().GetMaximum()<<", "<<m_HoughObjectsY[i].getHoughLineUR().GetMinimum()<<"]"<<endl;
  cout<<"    LineUL pars: "<<m_HoughObjectsY[i].getHoughLineUL().GetParameter(0)<<"  "<<m_HoughObjectsY[i].getHoughLineUL().GetParameter(1)<<", range: ["<<m_HoughObjectsY[i].getHoughLineUL().GetMaximum()<<", "<<m_HoughObjectsY[i].getHoughLineUL().GetMinimum()<<"]"<<endl;
  cout<<"    LineDR pars: "<<m_HoughObjectsY[i].getHoughLineDR().GetParameter(0)<<"  "<<m_HoughObjectsY[i].getHoughLineDR().GetParameter(1)<<", range: ["<<m_HoughObjectsY[i].getHoughLineDR().GetMaximum()<<", "<<m_HoughObjectsY[i].getHoughLineDR().GetMinimum()<<"]"<<endl;
  cout<<"    LineDL pars: "<<m_HoughObjectsY[i].getHoughLineDL().GetParameter(0)<<"  "<<m_HoughObjectsY[i].getHoughLineDL().GetParameter(1)<<", range: ["<<m_HoughObjectsY[i].getHoughLineDL().GetMaximum()<<", "<<m_HoughObjectsY[i].getHoughLineDL().GetMinimum()<<"]"<<endl;
  cout<<endl;
}
*/

   
cout<<"  HoughClusteringAlg: Fill bins to get the hills"<<endl;
    //Fill bins to get the hills
    PandoraPlus::HoughSpace m_HoughSpaceX; 
    PandoraPlus::HoughSpace m_HoughSpaceY; 
    FillHoughSpace(m_HoughObjectsX, m_HoughSpaceX);
    FillHoughSpace(m_HoughObjectsY, m_HoughSpaceY);

/*
cout<<"  Save Hough space map"<<endl;
TCanvas *c1 = new TCanvas("c1", "c1", 700,500);
c1->cd();
TH2F *h2_mapXZ = new TH2F();
*h2_mapXZ = m_HoughSpaceX.getSpaceMap();
cout<<" If space map is empty: "<<h2_mapXZ->Integral()<<endl;
h2_mapXZ->Draw("colz");
//c1->Draw();
c1->SaveAs("/cefs/higgs/guofy/CEPCSW_v1.2/run/HoughMapXZ.C");
delete c1, h2_mapXZ;
*/

//cout<<"  HoughClusteringAlg: Find hills"<<endl;
    //Find hills
    FindingHills(m_HoughSpaceX);
    FindingHills(m_HoughSpaceY);

/*
cout<<"  Print Hough hills in HoughSpaceY: "<<endl;
cout<<"  Hill size: "<<m_HoughSpaceY.getHills().size()<<endl;
for(int ih=0; ih<m_HoughSpaceY.getHills().size(); ih++){
  int m_size = m_HoughSpaceY.getHills()[ih].getIndexAlpha().size(); 
  std::vector<int> vec_alpha = m_HoughSpaceY.getHills()[ih].getIndexAlpha();
  std::vector<int> vec_rho = m_HoughSpaceY.getHills()[ih].getIndexRho();

  printf("    In Hill #%d: cell size = %d, alpha range [%d, %d], rho range [%d, %d]. \n", 
          ih, m_size, *min_element(vec_alpha.begin(), vec_alpha.end()), *max_element(vec_alpha.begin(), vec_alpha.end()), 
          *min_element(vec_rho.begin(), vec_rho.end()), *max_element(vec_rho.begin(), vec_rho.end())  );

}
*/

//cout<<"  HoughClusteringAlg: Create output HoughClusters. ";
    //Create output HoughClusters 
    std::vector<const PandoraPlus::LongiCluster*> m_longiClusXCol; m_longiClusXCol.clear();
    std::vector<const PandoraPlus::LongiCluster*> m_longiClusYCol; m_longiClusYCol.clear(); 
    Transform2Clusters(m_HoughSpaceX, m_HoughObjectsX, m_longiClusXCol);
    Transform2Clusters(m_HoughSpaceY, m_HoughObjectsY, m_longiClusYCol);

    for(int ic=0; ic<m_longiClusXCol.size(); ic++) m_datacol.bk_LongiClusCol.push_back( const_cast<PandoraPlus::LongiCluster*>(m_longiClusXCol[ic]) );
    for(int ic=0; ic<m_longiClusYCol.size(); ic++) m_datacol.bk_LongiClusCol.push_back( const_cast<PandoraPlus::LongiCluster*>(m_longiClusYCol[ic]) );



cout<<"  HoughCluster size: "<<m_longiClusXCol.size()<<" / "<<m_longiClusYCol.size()<<endl;
/*
cout<<"  Print LongiClusterX: size = "<<m_longiClusXCol.size()<<endl;
for(int il=0; il<m_longiClusXCol.size(); il++){
  printf("    Clus#%d: Hough Par = (%.3f, %.3f, %.3f), shower layers: \n", il, m_longiClusXCol[il]->getHoughAlpha(), m_longiClusXCol[il]->getHoughRho(), m_longiClusXCol[il]->getHoughIntercept());
  for(int is=0; is<m_longiClusXCol[il]->getBarShowers().size(); is++) 
    printf("      Dlayer = %d, Pos/E = (%.2f, %.2f, %.2f, %.3f) \n", m_longiClusXCol[il]->getBarShowers()[is]->getDlayer(), 
                                                                     m_longiClusXCol[il]->getBarShowers()[is]->getPos().x(),
                                                                     m_longiClusXCol[il]->getBarShowers()[is]->getPos().y(),
                                                                     m_longiClusXCol[il]->getBarShowers()[is]->getPos().z(),
                                                                     m_longiClusXCol[il]->getBarShowers()[is]->getE() );
  cout<<endl;
}

cout<<"  Print LongiClusterY: size = "<<m_longiClusYCol.size()<<endl;
for(int il=0; il<m_longiClusYCol.size(); il++){
  printf("    Clus#%d: Hough Par = (%.3f, %.3f, %.3f), shower layers: \n", il, m_longiClusYCol[il]->getHoughAlpha(), m_longiClusYCol[il]->getHoughRho(), m_longiClusYCol[il]->getHoughIntercept());
  for(int is=0; is<m_longiClusYCol[il]->getBarShowers().size(); is++)
    printf("      Dlayer = %d, Pos/E = (%.2f, %.2f, %.2f, %.3f) \n", m_longiClusYCol[il]->getBarShowers()[is]->getDlayer(),
                                                                     m_longiClusYCol[il]->getBarShowers()[is]->getPos().x(),
                                                                     m_longiClusYCol[il]->getBarShowers()[is]->getPos().y(),
                                                                     m_longiClusYCol[il]->getBarShowers()[is]->getPos().z(),
                                                                     m_longiClusYCol[il]->getBarShowers()[is]->getE() );
  cout<<endl;
}
*/
    p_towers->at(it)->setLongiClusters(m_longiClusXCol, m_longiClusYCol); 
  }//End loop tower

  //m_datacol.TowerCol = p_towers;
//cout<<"End in HoughClusteringAlg"<<endl;
  p_towers = nullptr;
  return StatusCode::SUCCESS;
};



StatusCode HoughClusteringAlg::ClearAlgorithm(){

  return StatusCode::SUCCESS;
};


StatusCode HoughClusteringAlg::ConformalTransformation(std::vector<PandoraPlus::HoughObject>& m_Hobjects){
  if(m_Hobjects.size()==0) return StatusCode::SUCCESS;

  for(int io=0; io<m_Hobjects.size(); io++){
    PandoraPlus::CaloBarShower m_localMax = *(m_Hobjects[io].getLocalMax());

//cout<<"  Object #"<<io<<": local max position: ";
//printf("(%.2f, %.2f, %.2f) \n", m_localMax.getPos().x(), m_localMax.getPos().y(), m_localMax.getPos().z());

    int slayer = m_localMax.getSlayer();
    int module = m_localMax.getModule();
    bool f_isXclus = (slayer==0 ? true : false);
    double rotAngle = -module*PI/4.;
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

    m_Hobjects[io].SetCellSize(10.);
    m_Hobjects[io].SetSlayer(slayer);
    m_Hobjects[io].SetConformalPoint(pos);
  }

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

//cout<<"Origin point in conformal space: (X, Y) = ("<<centerX<<", "<<centerX<<")"<<endl;

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

		m_Hobjects[iobj].SetHoughLine(line1, line2, line3, line4);
  }

  return StatusCode::SUCCESS;
}


StatusCode HoughClusteringAlg::FillHoughSpace(std::vector<PandoraPlus::HoughObject>& m_Hobjects, PandoraPlus::HoughSpace& m_Hspace){
  if(m_Hobjects.size()==0) return StatusCode::SUCCESS;

  m_Hspace.Clear(); 

  TH2F m_houghMap("","",(int)settings.map_floatPars["Nbins_alpha"], -0.1, PI, (int)settings.map_floatPars["Nbins_rho"], -130, 130);
  double width_alpha = (PI+0.1)/(double)settings.map_floatPars["Nbins_alpha"];
  double width_rho = 260./(double)settings.map_floatPars["Nbins_rho"]; 

  for(int io=0; io<m_Hobjects.size(); io++){
//cout<<"  In HoughObj #"<<io<<endl;

    TF1 line_ur = m_Hobjects[io].getHoughLineUR();
    TF1 line_ul = m_Hobjects[io].getHoughLineUL();
    TF1 line_dr = m_Hobjects[io].getHoughLineDR();
    TF1 line_dl = m_Hobjects[io].getHoughLineDL();


    //Loop for alpha bins
		for(int ibin=0; ibin<m_houghMap.GetNbinsX(); ibin++){
      double m_alphaL = ibin*width_alpha-0.1; 
      if(line_ur.Eval(m_alphaL)>130.-2*width_rho
         || line_ur.Eval(m_alphaL)<-130.+2*width_rho
         || line_dl.Eval(m_alphaL)>130.-2*width_rho
         || line_dl.Eval(m_alphaL)<-130.+2*width_rho) continue;


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
      nbin_rho_d = floor( (rho_min+130.)/width_rho );
      nbin_rho_u = floor( (rho_max+130.)/width_rho );

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

  m_Hspace.SetSpaceMap(m_houghMap);

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
      else m_hillCol.push_back(m_hill);
    }
  }}
  m_Hspace.SetHoughHills( m_hillCol );
  return StatusCode::SUCCESS;
}

int HoughClusteringAlg::ExpandingPeak(TH2 &houghMap, int index_a, int index_b, PandoraPlus::HoughSpace::HoughHill& hill){

  if(index_a<0 || index_b<0 || index_a>=settings.map_floatPars["Nbins_alpha"] || index_b>=settings.map_floatPars["Nbins_alpha"]) return 0;
  houghMap.SetBinContent(index_a+1, index_b+1, -1.*houghMap.GetBinContent(index_a+1, index_b+1) );
  int count = 0;
  for(int fl1=-1; fl1<=1; fl1++){
  for(int fl2=-1; fl2<=1; fl2++){
    if(fl1!=0 || fl2!=0){
      if( Abs(houghMap.GetBinContent(index_a+fl1+1, index_b+fl2+1)) > Abs( houghMap.GetBinContent(index_a+1, index_b+1))) return 1;
      if(houghMap.GetBinContent(index_a+fl1+1, index_b+fl2+1) == Abs(houghMap.GetBinContent(index_a+1, index_b+1)) ){
        hill.AddCell(index_a+fl1, index_b+fl2, houghMap.GetBinContent(index_a+fl1+1, index_b+fl2+1) );
        count += ExpandingPeak(houghMap, index_a+fl1, index_b+fl2, hill );
        if(count>0) return 1;
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


//std::vector<const PandoraPlus::LongiCluster*> tmp_clusCol; tmp_clusCol.clear(); 
//cout<<"  Transform2Clusters: input hill size: "<<m_hills.size()<<endl;

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

    if(ave_alpha>alpha_min-alpha_width && ave_alpha<= alpha_min+(alpha_max-alpha_min)/2.) {
    for(int io=0; io<m_Hobjects.size(); io++){
        if( ( m_Hobjects[io].getHoughLineDL().Eval(ave_alpha)<ave_rho+rho_width &&  
              m_Hobjects[io].getHoughLineUR().Eval(ave_alpha)>ave_rho-rho_width) ||
            ( m_Hobjects[io].getHoughLineDL().Eval(ave_alpha)>ave_rho-rho_width &&
              m_Hobjects[io].getHoughLineUR().Eval(ave_alpha)<ave_rho+rho_width) )
          { 
          m_clus->addBarShower( m_Hobjects[io].getLocalMax() ); 
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
          m_clus->addBarShower( m_Hobjects[io].getLocalMax() );
          m_clus->setHoughPars(ave_alpha, ave_rho );
          }
    }}
//tmp_clusCol.push_back(m_clus);

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

//cout<<"  Transform2Cluster: track removal cutflow: "<<endl;
//cout<<"    Initial: "<<tmp_clusCol.size()<<endl;
//cout<<"    After >= 3 hits "<<m_clusCol.size()<<endl;
/*
cout<<"  Print clusters: "<<endl;
for(int ic=0; ic<tmp_clusCol.size(); ic++){
  printf("    Nhit %d, Energy %.2f, HoughPar (%.2f, %.2f).  ", 
          tmp_clusCol[ic].getBarShowers().size(), tmp_clusCol[ic].getE(),
          tmp_clusCol[ic].getHoughAlpha(), tmp_clusCol[ic].getHoughRho() );
  for(int ihit=0; ihit<tmp_clusCol[ic].getBarShowers().size(); ihit++) cout<<tmp_clusCol[ic].getBarShowers()[ihit].getDlayer()<<"  ";
  cout<<endl;
}
cout<<endl;

cout<<"  Print clusters: "<<endl;
for(int ic=0; ic<m_clusCol.size(); ic++)
  printf("    Nhit %d, Energy %.2f, HoughPar (%.2f, %.2f) \n",
          m_clusCol[ic].getBarShowers().size(), m_clusCol[ic].getE(), 
          m_clusCol[ic].getHoughAlpha(), m_clusCol[ic].getHoughRho() );
cout<<endl;
*/

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
//for(int ah = 0; ah<m_nhit; ah++) cout<<m_longiClusCol[ic].getBarShowers()[ah].getDlayer()<<'\t';
//cout<<endl;

    for(int ih=0; ih<m_nhit-1; ih++){
    if(m_longiClusCol[ic]->getBarShowers()[ih+1]->getDlayer() - m_longiClusCol[ic]->getBarShowers()[ih]->getDlayer() > 2){
//cout<<"  Need to depart in hit #"<<ih<<endl;

      PandoraPlus::LongiCluster* clus_head = new PandoraPlus::LongiCluster();
      PandoraPlus::LongiCluster* clus_tail = new PandoraPlus::LongiCluster();

      clus_head->setHoughPars( m_longiClusCol[ic]->getHoughAlpha(), m_longiClusCol[ic]->getHoughRho() );
      clus_head->setIntercept( m_longiClusCol[ic]->getHoughIntercept() );
      clus_tail->setHoughPars( m_longiClusCol[ic]->getHoughAlpha(), m_longiClusCol[ic]->getHoughRho() );
      clus_tail->setIntercept( m_longiClusCol[ic]->getHoughIntercept() );

      for(int jh=0; jh<=ih; jh++) clus_head->addBarShower( m_longiClusCol[ic]->getBarShowers()[jh] );
      for(int jh=ih+1; jh<m_nhit; jh++) clus_tail->addBarShower( m_longiClusCol[ic]->getBarShowers()[jh] );

//cout<<"  Head cluster size: "<<clus_head.getBarShowers().size()<<", isContinueN: "<<clus_head.isContinueN(settings.th_continuetrkN)<<endl;
//cout<<"  Tail cluster size: "<<clus_tail.getBarShowers().size()<<", isContinueN: "<<clus_tail.isContinueN(settings.th_continuetrkN)<<endl;

      if( clus_head->isContinueN(settings.map_floatPars["th_continuetrkN"]) ) m_longiClusCol.push_back(clus_head);
      if( clus_tail->isContinueN(settings.map_floatPars["th_continuetrkN"]) ) m_longiClusCol.push_back(clus_tail);

      delete m_longiClusCol[ic]; m_longiClusCol[ic] = NULL;
      m_longiClusCol.erase(m_longiClusCol.begin()+ic);
      ic--;
      break;
    }}

  }
//cout<<"    After isolated hits removal: "<<m_longiClusCol.size()<<endl;

/*
cout<<"  Print clusters: "<<endl;
for(int ic=0; ic<m_longiClusCol.size(); ic++)
  printf("    Nhit %d, Energy %.2f, HoughPar (%.2f, %.2f) \n",
          m_longiClusCol[ic]->getBarShowers().size(), m_longiClusCol[ic]->getE(), 
          m_longiClusCol[ic]->getHoughAlpha(), m_longiClusCol[ic]->getHoughRho() );
cout<<endl;
*/

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

/*
cout<<"    After subset removal: "<<m_longiClusCol.size()<<endl;
cout<<"  Print clusters: "<<endl;
for(int ic=0; ic<m_longiClusCol.size(); ic++)
  printf("    Nhit %d, Energy %.2f, Layer range [%d, %d], intercept %.2f, energy %.2f \n",
          m_longiClusCol[ic]->getBarShowers().size(), m_longiClusCol[ic]->getE(),
          m_longiClusCol[ic]->getBeginningDlayer(), m_longiClusCol[ic]->getEndDlayer(),
          m_longiClusCol[ic]->getHoughIntercept(), m_longiClusCol[ic]->getE() );
cout<<endl;
*/


  //Cut on energy and intercept
  for(int ic=0; ic<m_longiClusCol.size(); ic++){
    if( fabs(m_longiClusCol[ic]->getHoughIntercept())>=settings.map_floatPars["th_intercept"] || m_longiClusCol[ic]->getE()<settings.map_floatPars["th_AxisE"]){
      delete m_longiClusCol[ic]; m_longiClusCol[ic]=NULL;
      m_longiClusCol.erase(m_longiClusCol.begin()+ic );
      ic--;
    }
  }
//cout<<"    After energy and intercetp cut: "<<m_longiClusCol.size()<<endl;


  //Overlap with other clusters: Iter 1. 
  for(int ic=0; ic<m_longiClusCol.size()-1; ic++){
  for(int jc=ic+1; jc<m_longiClusCol.size(); jc++){
    if(ic>=m_longiClusCol.size()) ic--;

    if( fabs(m_longiClusCol[ic]->getHoughAlpha() - m_longiClusCol[jc]->getHoughAlpha()) > settings.map_floatPars["th_dAlpha1"] || 
        fabs(m_longiClusCol[ic]->getHoughRho() - m_longiClusCol[jc]->getHoughRho()) > settings.map_floatPars["th_dRho1"]  ) 
      continue;

    double m_ratio1 = m_longiClusCol[ic]->OverlapRatioE(m_longiClusCol[jc]);
    double m_ratio2 = m_longiClusCol[jc]->OverlapRatioE(m_longiClusCol[ic]);

//printf("      Tag branch: Clus #%d: (R, E) = (%.2f, %.2f). Clus #%d: (%.2f, %.2f). \n", ic, m_ratio1, m_longiClusCol[ic]->getE(), jc, m_ratio2, m_longiClusCol[jc]->getE() );
    
    if(m_ratio1>settings.map_floatPars["th_overlapE1"] && m_longiClusCol[ic]->getE()<m_longiClusCol[jc]->getE()){
      delete m_longiClusCol[ic]; m_longiClusCol[ic] = NULL;
      m_longiClusCol.erase( m_longiClusCol.begin()+ic );
      ic--;
      break;
    }

    if(m_ratio2>settings.map_floatPars["th_overlapE1"] && m_longiClusCol[jc]->getE()<m_longiClusCol[ic]->getE()){
      delete m_longiClusCol[jc]; m_longiClusCol[jc] = NULL;
      m_longiClusCol.erase( m_longiClusCol.begin()+jc );
      jc--;
    }

  }}

//cout<<"    After 1st iter tagging: "<<m_longiClusCol.size()<<endl;

  //Overlap with other clusters: Iter 2. 
  for(int ic=0; ic<m_longiClusCol.size()-1; ic++){
  for(int jc=ic+1; jc<m_longiClusCol.size(); jc++){
    if(ic>=m_longiClusCol.size()) ic--;

    if( fabs(m_longiClusCol[ic]->getHoughAlpha() - m_longiClusCol[jc]->getHoughAlpha()) > settings.map_floatPars["th_dAlpha2"] ||
        fabs(m_longiClusCol[ic]->getHoughRho() - m_longiClusCol[jc]->getHoughRho()) > settings.map_floatPars["th_dRho2"]  )
      continue;

    double m_ratio1 = m_longiClusCol[ic]->OverlapRatioE(m_longiClusCol[jc]);
    double m_ratio2 = m_longiClusCol[jc]->OverlapRatioE(m_longiClusCol[ic]);

//printf("      Tag branch: Clus #%d: (R, E) = (%.2f, %.2f). Clus #%d: (%.2f, %.2f). \n", ic, m_ratio1, m_longiClusCol[ic]->getE(), jc, m_ratio2, m_longiClusCol[jc]->getE() );

    if(m_ratio1>settings.map_floatPars["th_overlapE2"] && m_longiClusCol[ic]->getE()<m_longiClusCol[jc]->getE()){
      delete m_longiClusCol[ic]; m_longiClusCol[ic]=NULL;
      m_longiClusCol.erase( m_longiClusCol.begin()+ic );
      ic--;
      break;
    }

    if(m_ratio2>settings.map_floatPars["th_overlapE2"] && m_longiClusCol[jc]->getE()<m_longiClusCol[ic]->getE()){
      delete m_longiClusCol[jc]; m_longiClusCol[jc]=NULL;
      m_longiClusCol.erase( m_longiClusCol.begin()+jc );
      jc--;
    }

  }}
//cout<<"    After 2nd iter tagging: "<<m_longiClusCol.size()<<endl;


/*
  for(int ic=0; ic<m_longiClusCol.size()-1; ic++){
  for(int jc=ic+1; jc<m_longiClusCol.size(); jc++){
    if(ic>=m_longiClusCol.size()) ic--;

    double m_ratio1 = m_longiClusCol[ic].OverlapRatio(m_longiClusCol[jc]); // R = Nsame[i|j]/Ntot[i]. 
    double m_ratio2 = m_longiClusCol[jc].OverlapRatio(m_longiClusCol[ic]);

    if(m_ratio1==0 && m_ratio2==0) { continue; }
    else if(m_ratio1==0 || m_ratio2==0) { continue; }

    //Case1: subset
    else if( m_ratio1==1 ){ 
      m_longiClusCol[jc].MergeCluster( m_longiClusCol[ic] );
      m_longiClusCol.erase(m_longiClusCol.begin()+ic );
      ic--;
    }
    else if( m_ratio2==1 ){
      m_longiClusCol[ic].MergeCluster( m_longiClusCol[jc] );
      m_longiClusCol.erase(m_longiClusCol.begin()+jc );
      jc--;
    }

    //Case2: has overlap. Check cluster size first. 
    //  Case 2.1: N_i<Nth, N_j<Nth.
    else if(m_longiClusCol[ic].getBarShowers().size()<settings.th_GoodLayer1 && m_longiClusCol[jc].getBarShowers().size()<settings.th_GoodLayer1){
      m_longiClusCol[ic].MergeCluster( m_longiClusCol[jc] );
      m_longiClusCol.erase(m_longiClusCol.begin()+jc );
      jc--;
    }

    //  Case 2.2: N_i>=Nth, N_j>=Nth
    else if(m_longiClusCol[ic].getBarShowers().size()>=settings.th_GoodLayer1 && m_longiClusCol[jc].getBarShowers().size()>=settings.th_GoodLayer1){
      //Case 2.2.1: R1>Rth && R2>Rth
      if(m_ratio1>settings.th_overlap && m_ratio2>settings.th_overlap){
        m_longiClusCol[ic].MergeCluster( m_longiClusCol[jc] );
        m_longiClusCol.erase(m_longiClusCol.begin()+jc );
        jc--;
      }

      //Case 2.2.2: R1<Rth || R2<Rth
      else{
        std::vector<PandoraPlus::CaloBarShower> m_sharedshs; m_sharedshs.clear();
        for(int is=0; is<m_longiClusCol[ic].getBarShowers().size(); is++){
        for(int js=0; js<m_longiClusCol[jc].getBarShowers().size(); js++){
          if(m_longiClusCol[ic].getBarShowers()[is]==m_longiClusCol[jc].getBarShowers()[js]) m_sharedshs.push_back(m_longiClusCol[ic].getBarShowers()[is]);
        }}
        if(m_longiClusCol[ic].getE()>m_longiClusCol[jc].getE()) m_longiClusCol[jc].RemoveShowers( m_sharedshs );  
        else m_longiClusCol[ic].RemoveShowers( m_sharedshs );
      }
    }

    //  Case 2.3: N_i>=Nth, N_j<Nth
    else{
      //Case 2.3.1: R1<Rth && R2<Rth
      if(m_ratio1<settings.th_overlap && m_ratio2<settings.th_overlap){

        std::vector<PandoraPlus::CaloBarShower> m_sharedshs; m_sharedshs.clear();
        for(int is=0; is<m_longiClusCol[ic].getBarShowers().size(); is++){
        for(int js=0; js<m_longiClusCol[jc].getBarShowers().size(); js++){
          if(m_longiClusCol[ic].getBarShowers()[is]==m_longiClusCol[jc].getBarShowers()[js]) m_sharedshs.push_back(m_longiClusCol[ic].getBarShowers()[is]);
        }}

        if(m_longiClusCol[ic].getE()>m_longiClusCol[jc].getE()) m_longiClusCol[jc].RemoveShowers( m_sharedshs );
        else m_longiClusCol[ic].RemoveShowers( m_sharedshs );
      }

      //Case 2.3.2: R1>=Rth || R2>=Rth
      else{
        m_longiClusCol[ic].MergeCluster( m_longiClusCol[jc] );
        m_longiClusCol.erase(m_longiClusCol.begin()+jc );
        jc--;
      }
    }
*/
/*  for(int ic=0; ic<m_longiClusCol.size(); ic++){
    if( m_longiClusCol[ic].getBarShowers().size()<=settings.th_peak || 
        settings.fl_continuetrk && !m_longiClusCol[ic].isContinue() ||
        !m_longiClusCol[ic].isContinueN(settings.th_continuetrkN) ){

        m_longiClusCol.erase(m_longiClusCol.begin()+ic );
        ic--;
    }
  }
  }}
*/

/*
cout<<"  Print clusters: "<<endl;
for(int ic=0; ic<m_longiClusCol.size(); ic++)
  printf("    Nhit %d, Energy %.2f, HoughPar (%.2f, %.2f, %.2f) \n",
          m_longiClusCol[ic].getBarShowers().size(), m_longiClusCol[ic].getE(), 
          m_longiClusCol[ic].getHoughAlpha(), m_longiClusCol[ic].getHoughRho(), m_longiClusCol[ic].getHoughIntercept() );
cout<<endl;
*/
  return StatusCode::SUCCESS;
};

#endif
