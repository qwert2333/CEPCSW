#ifndef HOUGHCLUSTERINGALG_C
#define HOUGHCLUSTERINGALG_C

#include "Algorithm/HoughClusteringAlg.h"
#include <algorithm>
#include "TCanvas.h"
using namespace std;
using namespace TMath;
void HoughClusteringAlg::Settings::SetInitialValue(){
  Nbins_alpha = 500;
  Nbins_rho = 500;
  th_Layers = 10;
  th_peak = 4;
  fl_continuetrk = true;
  th_continuetrkN = 3;
}

StatusCode HoughClusteringAlg::Initialize(){

  return StatusCode::SUCCESS;
}

StatusCode HoughClusteringAlg::RunAlgorithm( HoughClusteringAlg::Settings& m_settings, PandoraPlusDataCol& m_datacol){
  settings = m_settings; 

  std::vector<CRDEcalEDM::CRDCaloBlock> m_blocks = m_datacol.BlockVec;
  std::vector<CRDEcalEDM::CRDCaloTower> m_towers = m_datacol.TowerCol; 

  std::vector<CRDEcalEDM::CRDCaloHitLongiCluster> m_longiClusXCol; m_longiClusXCol.clear();
  std::vector<CRDEcalEDM::CRDCaloHitLongiCluster> m_longiClusYCol; m_longiClusYCol.clear(); 

  for(int it=0; it<m_towers.size(); it++){
cout<<"  Hough Clustering in tower #"<<it<<endl;
cout<<"  Block number in this tower: "<<m_towers[it].getBlocks().size()<<endl;

    std::vector<CRDEcalEDM::CRDCaloBarShower> m_localMaxXCol; m_localMaxXCol.clear();
    std::vector<CRDEcalEDM::CRDCaloBarShower> m_localMaxYCol; m_localMaxYCol.clear();

    for(int ib=0; ib<m_towers[it].getBlocks().size(); ib++){
      if(m_towers[it].getBlocks()[ib].getDlayer()>settings.th_Layers) continue;
      std::vector<CRDEcalEDM::CRDCaloBarShower> tmp_showerX = m_towers[it].getBlocks()[ib].getShowerXCol();
      std::vector<CRDEcalEDM::CRDCaloBarShower> tmp_showerY = m_towers[it].getBlocks()[ib].getShowerYCol();

      m_localMaxXCol.insert(m_localMaxXCol.end(), tmp_showerX.begin(), tmp_showerX.end() );
      m_localMaxYCol.insert(m_localMaxYCol.end(), tmp_showerY.begin(), tmp_showerY.end() );
    }


    if(m_localMaxXCol.size()==0 && m_localMaxYCol.size()==0) continue; 
cout<<"  Local maximum size: barX "<<m_localMaxXCol.size()<<"  barY "<<m_localMaxYCol.size()<<endl;

    std::vector<CRDEcalEDM::CRDHoughObject> m_HoughObjectsX; m_HoughObjectsX.clear(); 
    std::vector<CRDEcalEDM::CRDHoughObject> m_HoughObjectsY; m_HoughObjectsY.clear(); 
    for(int il=0; il<m_localMaxXCol.size(); il++){
      CRDEcalEDM::CRDHoughObject m_obj; m_obj.Clear();
      m_obj.SetLocalMax( m_localMaxXCol[il] );
      m_HoughObjectsX.push_back(m_obj);
    }
    for(int il=0; il<m_localMaxYCol.size(); il++){
      CRDEcalEDM::CRDHoughObject m_obj; m_obj.Clear();
      m_obj.SetLocalMax( m_localMaxYCol[il] );
      m_HoughObjectsY.push_back(m_obj);
    }

cout<<"  HoughClusteringAlg: Conformal transformation"<<endl;
    //Conformal transformation  
    ConformalTransformation(m_HoughObjectsX);
    ConformalTransformation(m_HoughObjectsY);

/*
cout<<"    Local max (HoughObjectX): "<<endl;
for(int i=0; i<m_HoughObjectsX.size(); i++) printf("(%.2f, %.2f, %.2f) \t", m_HoughObjectsX[i].getLocalMax().getPos().x(), m_HoughObjectsX[i].getLocalMax().getPos().y(), m_HoughObjectsX[i].getLocalMax().getPos().z() );
cout<<endl;
cout<<"    Conformal Point (HoughObjectX): "<<endl;
for(int i=0; i<m_HoughObjectsX.size(); i++) printf("(%.2f, %.2f) \t", m_HoughObjectsX[i].getConformPointUR().X(), m_HoughObjectsX[i].getConformPointUR().Y());
cout<<endl;

cout<<"    Local max (HoughObjectY): "<<endl;
for(int i=0; i<m_HoughObjectsY.size(); i++) printf("(%.2f, %.2f, %.2f) \t", m_HoughObjectsY[i].getLocalMax().getPos().x(), m_HoughObjectsY[i].getLocalMax().getPos().y(), m_HoughObjectsY[i].getLocalMax().getPos().z() );
cout<<endl;
cout<<"    Conformal Point (HoughObjectY): "<<endl;
for(int i=0; i<m_HoughObjectsY.size(); i++) printf("(%.2f, %.2f) \t", m_HoughObjectsY[i].getConformPointUR().X(), m_HoughObjectsY[i].getConformPointUR().Y());
cout<<endl;
*/

cout<<"  HoughClusteringAlg: Hough transformation"<<endl;
    //Hough transformation
    HoughTransformation(m_HoughObjectsX);
    HoughTransformation(m_HoughObjectsY);
/*
cout<<"  Transformed Hough lines X: "<<endl;
for(int i=0; i<m_HoughObjectsX.size(); i++){
  cout<<"    LineUR pars: "<<m_HoughObjectsX[i].getHoughLineUR().GetParameter(0)<<"  "<<m_HoughObjectsX[i].getHoughLineUR().GetParameter(1)<<", range: ["<<m_HoughObjectsX[i].getHoughLineUR().GetMaximum()<<", "<<m_HoughObjectsX[i].getHoughLineUR().GetMinimum()<<"]"<<endl;
  cout<<"    LineUL pars: "<<m_HoughObjectsX[i].getHoughLineUL().GetParameter(0)<<"  "<<m_HoughObjectsX[i].getHoughLineUL().GetParameter(1)<<", range: ["<<m_HoughObjectsX[i].getHoughLineUL().GetMaximum()<<", "<<m_HoughObjectsX[i].getHoughLineUL().GetMinimum()<<"]"<<endl;
  cout<<"    LineDR pars: "<<m_HoughObjectsX[i].getHoughLineDR().GetParameter(0)<<"  "<<m_HoughObjectsX[i].getHoughLineDR().GetParameter(1)<<", range: ["<<m_HoughObjectsX[i].getHoughLineDR().GetMaximum()<<", "<<m_HoughObjectsX[i].getHoughLineDR().GetMinimum()<<"]"<<endl;
  cout<<"    LineDL pars: "<<m_HoughObjectsX[i].getHoughLineDL().GetParameter(0)<<"  "<<m_HoughObjectsX[i].getHoughLineDL().GetParameter(1)<<", range: ["<<m_HoughObjectsX[i].getHoughLineDL().GetMaximum()<<", "<<m_HoughObjectsX[i].getHoughLineDL().GetMinimum()<<"]"<<endl;
}

cout<<"  Transformed Hough lines Y: "<<endl;
for(int i=0; i<m_HoughObjectsY.size(); i++){
  cout<<"    LineUR pars: "<<m_HoughObjectsY[i].getHoughLineUR().GetParameter(0)<<"  "<<m_HoughObjectsY[i].getHoughLineUR().GetParameter(1)<<", range: ["<<m_HoughObjectsY[i].getHoughLineUR().GetMaximum()<<", "<<m_HoughObjectsY[i].getHoughLineUR().GetMinimum()<<"]"<<endl;
  cout<<"    LineUL pars: "<<m_HoughObjectsY[i].getHoughLineUL().GetParameter(0)<<"  "<<m_HoughObjectsY[i].getHoughLineUL().GetParameter(1)<<", range: ["<<m_HoughObjectsY[i].getHoughLineUL().GetMaximum()<<", "<<m_HoughObjectsY[i].getHoughLineUL().GetMinimum()<<"]"<<endl;
  cout<<"    LineDR pars: "<<m_HoughObjectsY[i].getHoughLineDR().GetParameter(0)<<"  "<<m_HoughObjectsY[i].getHoughLineDR().GetParameter(1)<<", range: ["<<m_HoughObjectsY[i].getHoughLineDR().GetMaximum()<<", "<<m_HoughObjectsY[i].getHoughLineDR().GetMinimum()<<"]"<<endl;
  cout<<"    LineDL pars: "<<m_HoughObjectsY[i].getHoughLineDL().GetParameter(0)<<"  "<<m_HoughObjectsY[i].getHoughLineDL().GetParameter(1)<<", range: ["<<m_HoughObjectsY[i].getHoughLineDL().GetMaximum()<<", "<<m_HoughObjectsY[i].getHoughLineDL().GetMinimum()<<"]"<<endl;
}
*/

   
cout<<"  HoughClusteringAlg: Fill bins to get the hills"<<endl;
    //Fill bins to get the hills
    CRDEcalEDM::CRDHoughSpace m_HoughSpaceX; 
    CRDEcalEDM::CRDHoughSpace m_HoughSpaceY; 
    FillHoughSpace(m_HoughObjectsX, m_HoughSpaceX);
    FillHoughSpace(m_HoughObjectsY, m_HoughSpaceY);
/*   
//Print Hough space map:
cout<<"HoughMap X: "<<endl;
for(int ix=0; ix<settings.Nbins_alpha; ix++){
for(int iy=0; iy<settings.Nbins_rho; iy++){
  cout<<m_HoughSpaceX.getSpaceMap().GetBinContent(ix+1, iy+1)<<"  ";
}
cout<<endl;
}
cout<<"HoughMap Y: "<<endl;
for(int ix=0; ix<settings.Nbins_alpha; ix++){
for(int iy=0; iy<settings.Nbins_rho; iy++){
  cout<<m_HoughSpaceY.getSpaceMap().GetBinContent(ix+1, iy+1)<<"  ";
}
cout<<endl;
}
*/

cout<<"  HoughClusteringAlg: Find hills"<<endl;
    //Find hills
    FindingHills(m_HoughSpaceX);
    FindingHills(m_HoughSpaceY);

cout<<"  HoughClusteringAlg: Create output HoughClusters"<<endl;
    //Create output HoughClusters 
    Transform2Clusters(m_HoughSpaceX, m_HoughObjectsX, m_longiClusXCol);
    Transform2Clusters(m_HoughSpaceY, m_HoughObjectsY, m_longiClusYCol);

  }//End loop tower

  //CleanClusters( m_longiClusXCol );
  //CleanClusters( m_longiClusYCol );

  m_datacol.LongiClusXCol = m_longiClusXCol;
  m_datacol.LongiClusYCol = m_longiClusYCol;

cout<<"End in HoughClusteringAlg"<<endl;

  return StatusCode::SUCCESS;
}

StatusCode HoughClusteringAlg::ClearAlgorithm(){

  return StatusCode::SUCCESS;
}


StatusCode HoughClusteringAlg::ConformalTransformation(std::vector<CRDEcalEDM::CRDHoughObject>& m_Hobjects){
  if(m_Hobjects.size()==0) return StatusCode::SUCCESS;

  for(int io=0; io<m_Hobjects.size(); io++){
    CRDEcalEDM::CRDCaloBarShower m_localMax = m_Hobjects[io].getLocalMax();

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


StatusCode HoughClusteringAlg::HoughTransformation(std::vector<CRDEcalEDM::CRDHoughObject>& m_Hobjects){

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


StatusCode HoughClusteringAlg::FillHoughSpace(std::vector<CRDEcalEDM::CRDHoughObject>& m_Hobjects, CRDEcalEDM::CRDHoughSpace& m_Hspace){
  if(m_Hobjects.size()==0) return StatusCode::SUCCESS;

  m_Hspace.Clear(); 

  TH2F m_houghMap("","",settings.Nbins_alpha, -0.1, PI, settings.Nbins_rho, -200, 200);
  double width_alpha = (PI+0.1)/(double)settings.Nbins_alpha;
  double width_rho = 400./(double)settings.Nbins_rho; 

cout<<"Hough object size: "<<m_Hobjects.size()<<endl;  
  for(int io=0; io<m_Hobjects.size(); io++){

    TF1 line_ur = m_Hobjects[io].getHoughLineUR();
    TF1 line_ul = m_Hobjects[io].getHoughLineUL();
    TF1 line_dr = m_Hobjects[io].getHoughLineDR();
    TF1 line_dl = m_Hobjects[io].getHoughLineDL();

    //Loop for alpha bins
		for(int ibin=0; ibin<m_houghMap.GetNbinsX(); ibin++){
      double alpha = ibin*width_alpha-0.1; 
      double rho_u = alpha<PI/2. ? line_ur.Eval(alpha) : line_ul.Eval(alpha);
      double rho_d = alpha<PI/2. ? line_dl.Eval(alpha) : line_dr.Eval(alpha);

      double rho_u_p = (alpha+width_alpha)<PI/2. ? line_ur.Eval(alpha+width_alpha) : line_ul.Eval(alpha+width_alpha);
      double rho_d_p = (alpha+width_alpha)<PI/2. ? line_dl.Eval(alpha+width_alpha) : line_dr.Eval(alpha+width_alpha);

      double rho_min = min(min(min(rho_u,rho_d), rho_u_p), rho_d_p);
      double rho_max = max(max(max(rho_u,rho_d), rho_u_p), rho_d_p);


      double rho_tmp = rho_min; 
      while( floor(rho_tmp/width_rho)<floor(rho_max/width_rho) ){
        m_houghMap.Fill(alpha, rho_tmp);
        rho_tmp += width_rho;
      }
    }
	}//End loop Hough Objects

  m_Hspace.SetSpaceMap(m_houghMap);

  return StatusCode::SUCCESS;
}


StatusCode HoughClusteringAlg::FindingHills(CRDEcalEDM::CRDHoughSpace& m_Hspace){

  std::vector<CRDEcalEDM::CRDHoughSpace::HoughHill> m_hillCol; m_hillCol.clear(); 
  TH2F m_houghMap = m_Hspace.getSpaceMap();
  
  for(int ia=0; ia<settings.Nbins_alpha; ia++){
  for(int ir=0; ir<settings.Nbins_rho; ir++){
/*
cout<<"  Looping for cell: "<<ia<<", "<<ir<<endl;
cout<<"  cell and neighbor: "<<endl;
if (!(ia==0 || ia==settings.Nbins_alpha-1 || ir==0 || ir==settings.Nbins_rho-1)){
cout<<m_houghMap.GetBinContent(ia, ir)<<'\t'<<m_houghMap.GetBinContent(ia+1, ir)<<'\t'<<m_houghMap.GetBinContent(ia+2, ir)<<endl;
cout<<m_houghMap.GetBinContent(ia, ir+1)<<'\t'<<m_houghMap.GetBinContent(ia+1, ir+1)<<'\t'<<m_houghMap.GetBinContent(ia+2, ir+1)<<endl;
cout<<m_houghMap.GetBinContent(ia, ir+2)<<'\t'<<m_houghMap.GetBinContent(ia+1, ir+2)<<'\t'<<m_houghMap.GetBinContent(ia+2, ir+2)<<endl;
}
*/
    if(m_houghMap.GetBinContent(ia+1, ir+1)>settings.th_peak){
      CRDEcalEDM::CRDHoughSpace::HoughHill m_hill; m_hill.Clear();
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

int HoughClusteringAlg::ExpandingPeak(TH2 &houghMap, int index_a, int index_b, CRDEcalEDM::CRDHoughSpace::HoughHill& hill){

  if(index_a<0 || index_b<0 || index_a>=settings.Nbins_alpha || index_b>settings.Nbins_alpha) return 0;
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

StatusCode HoughClusteringAlg::Transform2Clusters( CRDEcalEDM::CRDHoughSpace& m_Hspace, 
                                                   std::vector<CRDEcalEDM::CRDHoughObject>& m_Hobjects,
                                                   std::vector<CRDEcalEDM::CRDCaloHitLongiCluster>& m_longiClusCol )
{
  if(m_Hobjects.size()==0) return StatusCode::SUCCESS;

  std::vector<CRDEcalEDM::CRDCaloHitLongiCluster> m_clusCol; m_clusCol.clear(); 
  std::vector<CRDEcalEDM::CRDHoughSpace::HoughHill> m_hills = m_Hspace.getHills(); 
  double alpha_min = m_Hspace.getAlphaLowEdge();
  double alpha_max = m_Hspace.getAlphaUpEdge();
  double rho_min = m_Hspace.getRhoLowEdge();
  double rho_max = m_Hspace.getRhoUpEdge();
  double alpha_width = m_Hspace.getAlphaBinWidth();
  double rho_width = m_Hspace.getRhoBinWidth();

cout<<"Hill size: "<<m_hills.size()<<endl;
  for(int ih=0; ih<m_hills.size(); ih++){
    CRDEcalEDM::CRDCaloHitLongiCluster m_clus; m_clus.Clear();

    std::vector<int> index_alpha = m_hills[ih].getIndexAlpha();
    std::vector<int> index_rho = m_hills[ih].getIndexRho();
    double sum_alpha=0; double sum_rho=0; 
    for(int i=0; i<index_alpha.size(); i++) sum_alpha += index_alpha[i];
    for(int i=0; i<index_rho.size(); i++) sum_rho += index_rho[i];

    double ave_alpha = alpha_min + (sum_alpha/index_alpha.size()+0.5)*alpha_width;
    double ave_rho = rho_min + (sum_rho/index_rho.size()+0.5)*rho_width;

    if(ave_alpha>alpha_min-alpha_width && ave_alpha<=PI/2.){
    for(int io=0; io<m_Hobjects.size(); io++){
        if( ( m_Hobjects[io].getHoughLineDL().Eval(ave_alpha)<ave_rho+rho_width &&  
              m_Hobjects[io].getHoughLineUR().Eval(ave_alpha)>ave_rho-rho_width) ||
            ( m_Hobjects[io].getHoughLineDL().Eval(ave_alpha)>ave_rho-rho_width &&
              m_Hobjects[io].getHoughLineUR().Eval(ave_alpha)<ave_rho+rho_width) )

           m_clus.AddBarShower( m_Hobjects[io].getLocalMax() );
    }}

    else if(ave_alpha>PI/2. && ave_alpha<alpha_max+alpha_width){
    for(int io=0; io<m_Hobjects.size(); io++){
        if( ( m_Hobjects[io].getHoughLineDR().Eval(ave_alpha)<ave_rho+rho_width &&
              m_Hobjects[io].getHoughLineUL().Eval(ave_alpha)>ave_rho-rho_width) ||
            ( m_Hobjects[io].getHoughLineDR().Eval(ave_alpha)>ave_rho-rho_width &&
              m_Hobjects[io].getHoughLineUL().Eval(ave_alpha)<ave_rho+rho_width) )

           m_clus.AddBarShower( m_Hobjects[io].getLocalMax() );
    }}

    if(m_clus.getBarShowers().size()<=settings.th_peak) continue;
    if(settings.fl_continuetrk && !m_clus.isContinue()) continue;
    if(!m_clus.isContinueN(settings.th_continuetrkN)) continue;
    m_clus.FitAxis();
    m_clusCol.push_back(m_clus);
  }

  CleanClusters(m_clusCol);
  m_longiClusCol.insert(m_longiClusCol.end(), m_clusCol.begin(), m_clusCol.end());

  return StatusCode::SUCCESS;
}


StatusCode HoughClusteringAlg::CleanClusters( std::vector<CRDEcalEDM::CRDCaloHitLongiCluster>& m_longiClusCol ){
  if(m_longiClusCol.size()==0) return StatusCode::SUCCESS;

  for(int ic=0; ic<m_longiClusCol.size(); ic++){
  for(int jc=0; jc<m_longiClusCol.size(); jc++){
    if(ic>=m_longiClusCol.size()) ic--;
    if(ic!=jc && m_longiClusCol[ic].isSubset(m_longiClusCol[jc]) ){
      m_longiClusCol.erase(m_longiClusCol.begin()+jc );   
      jc--;
      if(ic>jc+1) ic--;
    }
  }}

  return StatusCode::SUCCESS;
}

#endif
