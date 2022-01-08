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

cout<<"  HoughClusteringAlg: Conformal transformation"<<endl;
    //Conformal transformation  
    std::vector<CRDEcalEDM::CRDHoughObject> m_HoughObjectsX; m_HoughObjectsX.clear(); 
    std::vector<CRDEcalEDM::CRDHoughObject> m_HoughObjectsY; m_HoughObjectsY.clear(); 
    ConformalTransformation(m_localMaxXCol, m_HoughObjectsX);
    ConformalTransformation(m_localMaxYCol, m_HoughObjectsY);

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


StatusCode HoughClusteringAlg::ConformalTransformation(std::vector<CRDEcalEDM::CRDCaloBarShower>& m_localMax, std::vector<CRDEcalEDM::CRDHoughObject>& m_Hobjects){

  if(m_localMax.size()==0) return StatusCode::SUCCESS;
  for(int il=0; il<m_localMax.size(); il++){
    CRDEcalEDM::CRDHoughObject m_obj; m_obj.Clear(); 
		m_obj.SetLocalMax( &(m_localMax[il]) );

    int slayer = m_localMax[il].getSlayer();
    int module = m_localMax[il].getModule();
    bool f_isXclus = (slayer==0 ? true : false);
    double rotAngle = -module*PI/4.;
    
		TVector3 p_localMax(0., 0., 0.); 
		p_localMax.SetXYZ(m_localMax[il].getPos().x(), m_localMax[il].getPos().y(), m_localMax[il].getPos().z());
    p_localMax.RotateZ(rotAngle);

    TVector2 pos_u, pos_d; 
		pos_u.SetX(p_localMax.X()/1000);
		pos_d.SetX(p_localMax.X()/1000);
    if(f_isXclus){ 
		  pos_u.SetY((p_localMax.Z()+5)/1000.);
		  pos_d.SetY((p_localMax.Z()+5)/1000.);
    }
		else{
      pos_u.SetY((p_localMax.Y()+5)/1000.);
      pos_d.SetY((p_localMax.Y()+5)/1000.);
		}

		double mod2_u = pos_u.Mod2();
		double mod2_d = pos_d.Mod2();
    pos_u.SetX(2*pos_u.X()/mod2_u); pos_u.SetY(2*pos_u.Y()/mod2_u);
    pos_d.SetX(2*pos_d.X()/mod2_d); pos_d.SetY(2*pos_d.Y()/mod2_d);

    m_obj.SetSlayer(slayer); 
    m_obj.SetConformalPoint(pos_u, pos_d);
    m_Hobjects.push_back(m_obj);
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
		if(m_Hobjects[iobj].getConformPointU().X()<m_minX) m_minX = m_Hobjects[iobj].getConformPointU().X(); 
		if(m_Hobjects[iobj].getConformPointU().X()>m_maxX) m_maxX = m_Hobjects[iobj].getConformPointU().X();
		if(m_Hobjects[iobj].getConformPointU().Y()<m_minY) m_minY = m_Hobjects[iobj].getConformPointU().Y(); 
		if(m_Hobjects[iobj].getConformPointU().Y()>m_maxY) m_maxY = m_Hobjects[iobj].getConformPointU().Y();
	}
  double centreX = (m_minX+m_maxX)/2.; 
	double centerY = (m_minY+m_maxY)/2.; 

  
  for(int iobj=0; iobj<m_Hobjects.size(); iobj++){
    TF1 line_u("", "[0]*cos(x) + [1]*sin(x)", -0.1, PI);
    line_u.SetParameters(m_Hobjects[iobj].getConformPointU().X()-centreX, m_Hobjects[iobj].getConformPointU().Y()-centerY);

    TF1 line_d("", "[0]*cos(x) + [1]*sin(x)", -0.1, PI);
    line_d.SetParameters(m_Hobjects[iobj].getConformPointD().X()-centreX, m_Hobjects[iobj].getConformPointD().Y()-centerY);

		m_Hobjects[iobj].SetHoughLine(line_u, line_d);
  }

  return StatusCode::SUCCESS;
}


StatusCode HoughClusteringAlg::FillHoughSpace(std::vector<CRDEcalEDM::CRDHoughObject>& m_Hobjects, CRDEcalEDM::CRDHoughSpace& m_Hspace){
  m_Hspace.Clear(); 

  TH2F m_houghMap("","",settings.Nbins_alpha, -0.1, PI, settings.Nbins_rho, -0.1, 0.1);
  double width_alpha = (PI+0.1)/(double)settings.Nbins_alpha;
  double width_rho = 0.2/(double)settings.Nbins_rho; 
  for(int io=0; io<m_Hobjects.size(); io++){
    TF1 line_u = m_Hobjects[io].getHoughLineU();
    TF1 line_d = m_Hobjects[io].getHoughLineD();

    std::vector<CRDEcalEDM::CRDHoughSpace::HoughCell> m_cells; m_cells.clear(); 
    //Loop for upper edge
    double rho_tmp = line_u.Eval(-0.1); 
		for(int ibin=0; ibin<m_houghMap.GetNbinsX(); ibin++){
      //Fill Hough space map
      double m_rho = line_u.Eval(-0.1+ibin*width_alpha);
      m_houghMap.Fill(ibin*width_alpha, m_rho);

      //Create HoughCell
      CRDEcalEDM::CRDHoughSpace::HoughCell m_cell; m_cell.Clear(); 
      m_cell.SetIndex(ibin, floor(m_rho/width_rho));
      m_cell.AddObject(m_Hobjects[io]); 
      m_cells.push_back(m_cell);

      //If the line cross several bins
      if( floor(rho_tmp/width_rho)!=floor(m_rho/width_rho) ) { //Not in the same rho bin
      double rho_step = min(rho_tmp, m_rho);
      while( floor(rho_step/width_rho) <= floor(max(rho_tmp, m_rho)/width_rho) ){
        rho_step += width_rho; 
        m_houghMap.Fill(ibin*width_alpha, rho_step);

        CRDEcalEDM::CRDHoughSpace::HoughCell p_cell; p_cell.Clear();
        p_cell.SetIndex(ibin, floor(rho_step/width_rho));
        p_cell.AddObject(m_Hobjects[io]);
        if( std::find( m_cells.begin(), m_cells.end(), p_cell )==m_cells.end() ) m_cells.push_back(p_cell);  
      }}
      rho_tmp = m_rho; 

    }

    //Loop for bottom edge
    rho_tmp = line_d.Eval(-0.1);
    for(int ibin=0; ibin<m_houghMap.GetNbinsX(); ibin++){

      double m_rho = line_d.Eval(-0.1+ibin*width_alpha);
      m_houghMap.Fill(ibin*width_alpha, m_rho);

      CRDEcalEDM::CRDHoughSpace::HoughCell m_cell; m_cell.Clear();
      m_cell.SetIndex(ibin, floor(m_rho/width_rho));
      m_cell.AddObject(m_Hobjects[io]);
      m_cells.push_back(m_cell);

      if( floor(rho_tmp/width_rho)!=floor(m_rho/width_rho) ) { //Not in the same rho bin
      double rho_step = min(rho_tmp, m_rho);
      while( floor(rho_step/width_rho) <= floor(max(rho_tmp, m_rho)/width_rho) ){
        rho_step += width_rho;
        m_houghMap.Fill(ibin*width_alpha, rho_step);

        CRDEcalEDM::CRDHoughSpace::HoughCell p_cell; p_cell.Clear();
        p_cell.SetIndex(ibin, floor(rho_step/width_rho));
        p_cell.AddObject(m_Hobjects[io]);
        if( std::find( m_cells.begin(), m_cells.end(), p_cell )==m_cells.end() ) m_cells.push_back(p_cell);
      }}

      rho_tmp = m_rho;
    } //End loop alpha bins. 
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
