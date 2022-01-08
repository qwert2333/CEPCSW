#ifndef HOUGHCLUSTERINGALG_C
#define HOUGHCLUSTERINGALG_C

#include "Algorithm/HoughClusteringAlg.h"
#include <algorithm>
using namespace std;
void HoughClusteringAlg::Settings::SetInitialValue(){
  Nbins_alpha = 800;
  Nbins_rho = 800;
  Layers = 10;
}

StatusCode HoughClusteringAlg::Initialize(){

  return StatusCode::SUCCESS;
}

StatusCode HoughClusteringAlg::RunAlgorithm( HoughClusteringAlg::Settings& m_settings, PandoraPlusDataCol& m_datacol){
  settings = m_settings; 

  std::vector<CRDEcalEDM::CRDCaloLayer> m_layerCol = m_datacol.LayerCol;
  std::vector<CRDEcalEDM::CRDCaloBlock> m_blocks = m_datacol.BlockVec;

cout<<"Block number: "<<m_blocks.size()<<endl;
  for(int ib=0; ib<m_blocks.size(); ib++){
cout<<"  Hough Clustering in block #"<<ib<<endl;
    std::vector<CRDEcalEDM::CRDCaloBarShower> m_localMaxXCol = m_blocks[ib].getShowerXCol();
    std::vector<CRDEcalEDM::CRDCaloBarShower> m_localMaxYCol = m_blocks[ib].getShowerYCol();
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

cout<<"  HoughClusteringAlg: Hough transformation"<<endl;
    //Hough transformation
    HoughTransformation(m_HoughObjectsX);
    HoughTransformation(m_HoughObjectsY);
   
cout<<"  HoughClusteringAlg: Fill bins to get the hills"<<endl;
    //Fill bins to get the hills
    CRDEcalEDM::CRDHoughSpace m_HoughSpaceX; 
    CRDEcalEDM::CRDHoughSpace m_HoughSpaceY; 
    FillHoughSpace(m_HoughObjectsX, m_HoughSpaceX);
    FillHoughSpace(m_HoughObjectsY, m_HoughSpaceY);
   
cout<<"  HoughClusteringAlg: Merge hills"<<endl;
    //Merge hills
    MergingHills(m_HoughSpaceX);
    MergingHills(m_HoughSpaceY);

cout<<"  HoughClusteringAlg: Create output HoughClusters"<<endl;
    //Create output HoughClusters 
    for(int ic=0; ic<m_HoughSpaceX.getHills().size(); ic++) m_datacol.LongiClusXCol.push_back( m_HoughSpaceX.getHills()[ic].TransformToCluster() ); 
    for(int ic=0; ic<m_HoughSpaceY.getHills().size(); ic++) m_datacol.LongiClusYCol.push_back( m_HoughSpaceY.getHills()[ic].TransformToCluster() ); 
  }
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

    int slayer = m_localMax.getSlayer();
    int module = m_localMax.getModule();
    bool f_isXclus = (slayer==0 ? true : false);
    double rotAngle = -module*PI/4.;

    TVector3 p_localMax(0., 0., 0.);
    p_localMax.SetXYZ(m_localMax.getPos().x(), m_localMax.getPos().y(), m_localMax.getPos().z());
    p_localMax.RotateZ(rotAngle);

    TVector2 pos;
    pos.SetX(p_localMax.X());
    if(f_isXclus){
      pos.SetY(p_localMax.Z());
    }
    else{
      pos.SetY(p_localMax.Y());
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
  
	double m_minX = 999;
	double m_maxX = -999;
	double m_minY = 999;
	double m_maxY = -999;
	for(int iobj=0; iobj<m_Hobjects.size(); iobj++){
    if(m_Hobjects[iobj].getConformPointUR().X()<m_minX) m_minX = m_Hobjects[iobj].getConformPointUR().X();
    if(m_Hobjects[iobj].getConformPointUR().X()>m_maxX) m_maxX = m_Hobjects[iobj].getConformPointUR().X();
    if(m_Hobjects[iobj].getConformPointUR().Y()<m_minY) m_minY = m_Hobjects[iobj].getConformPointUR().Y();
    if(m_Hobjects[iobj].getConformPointUR().Y()>m_maxY) m_maxY = m_Hobjects[iobj].getConformPointUR().Y();
	}
  double centreX = (m_minX+m_maxX)/2.; 
	double centerY = (m_minY+m_maxY)/2.; 

  
  for(int iobj=0; iobj<m_Hobjects.size(); iobj++){
    //UR
    TF1 line1("", "[0]*cos(x) + [1]*sin(x)", -0.1, PI/2, 2);  
    line1.SetParameters(m_Hobjects[iobj].getConformPointUR().X()-centreX, m_Hobjects[iobj].getConformPointUR().Y()-centerY);

    //UL
    TF1 line2("", "[0]*cos(x) + [1]*sin(x)", PI/2, PI, 2);
    line2.SetParameters(m_Hobjects[iobj].getConformPointUL().X()-centreX, m_Hobjects[iobj].getConformPointUL().Y()-centerY);

    //DR
    TF1 line3("", "[0]*cos(x) + [1]*sin(x)", PI/2, PI, 2);
    line3.SetParameters(m_Hobjects[iobj].getConformPointDR().X()-centreX, m_Hobjects[iobj].getConformPointDR().Y()-centerY);

    //DL   
    TF1 line4("", "[0]*cos(x) + [1]*sin(x)", -0.1, PI/2, 2);  
    line4.SetParameters(m_Hobjects[iobj].getConformPointDL().X()-centreX, m_Hobjects[iobj].getConformPointDL().Y()-centerY);

		m_Hobjects[iobj].SetHoughLine(line1, line2, line3, line4);
  }

  return StatusCode::SUCCESS;
}


StatusCode HoughClusteringAlg::FillHoughSpace(std::vector<CRDEcalEDM::CRDHoughObject>& m_Hobjects, CRDEcalEDM::CRDHoughSpace& m_Hspace){
  m_Hspace.Clear(); 

  TH2F m_houghMap("","",settings.Nbins_alpha, -0.1, PI, settings.Nbins_rho, -0.13, 0.13);
  double width_alpha = (PI+0.1)/(double)settings.Nbins_alpha;
  double width_rho = 0.26/(double)settings.Nbins_rho; 

  for(int io=0; io<m_Hobjects.size(); io++){
    TF1 line_ur = m_Hobjects[io].getHoughLineUR();
    TF1 line_ul = m_Hobjects[io].getHoughLineUL();
    TF1 line_dr = m_Hobjects[io].getHoughLineDR();
    TF1 line_dl = m_Hobjects[io].getHoughLineDL();

    std::vector<CRDEcalEDM::CRDHoughSpace::HoughCell> m_cells; m_cells.clear(); 
    //Loop for alpha bins
		for(int ibin=0; ibin<m_houghMap.GetNbinsX(); ibin++){
      double alpha = ibin*width_alpha; 
      double rho_u = alpha<PI/2. ? line_ur.Eval(alpha) : line_ul.Eval(alpha);
      double rho_d = alpha<PI/2. ? line_dl.Eval(alpha) : line_dr.Eval(alpha);

      double rho_u_p = (alpha+width_alpha)<PI/2. ? line_ur.Eval(alpha+width_alpha) : line_ul.Eval(alpha+width_alpha);
      double rho_d_p = (alpha+width_alpha)<PI/2. ? line_dl.Eval(alpha+width_alpha) : line_dr.Eval(alpha+width_alpha);

      double rho_min = min(min(min(rho_u,rho_d), rho_u_p), rho_d_p);
      double rho_max = max(max(max(rho_u,rho_d), rho_u_p), rho_d_p);

      double rho_tmp = rho_min; 
      while(floor(rho_tmp/width_rho)!=floor(rho_max/width_rho)){
        m_houghMap.Fill(alpha, rho_tmp);
        CRDEcalEDM::CRDHoughSpace::HoughCell p_cell; p_cell.Clear();
        p_cell.SetIndex(ibin, floor(rho_step/width_rho));
        p_cell.AddObject(m_Hobjects[io]);
        m_cells.push_back(p_cell);
        rho_tmp += width_rho;
      }
    }

    m_Hspace.SetHoughCells(m_cells);
	}//End loop Hough lines. 

  m_Hspace.SetSpaceMap(m_houghMap);

  return StatusCode::SUCCESS;
}


StatusCode HoughClusteringAlg::MergingHills(CRDEcalEDM::CRDHoughSpace& m_Hspace){

  std::vector<CRDEcalEDM::CRDHoughSpace::HoughHill> m_hills; m_hills.clear(); 
  std::vector<CRDEcalEDM::CRDHoughSpace::HoughCell> m_cells = m_Hspace.getCells(); 
  for(int ic=0; ic<m_cells.size(); ic++){
    bool f_inHill=false; 
    for(int ih=0; ih<m_hills.size(); ih++){
      if( m_hills[ih].isNeighbor(m_cells[ic]) ){
        m_hills[ih].AddCell(m_cells[ic]);
        m_cells.erase(m_cells.begin()+ic);
        f_inHill=true;
        break; 
      }
    }
    if(!f_inHill){
      CRDEcalEDM::CRDHoughSpace::HoughHill m_hill; m_hill.Clear(); 
      m_hill.AddCell( m_cells[ic] );
      m_hills.push_back(m_hill);
    }
  }

  return StatusCode::SUCCESS;
}



#endif
