#ifndef CALOHOUGHOBJECT_C
#define CALOHOUGHOBJECT_C

#include "Objects/HoughObject.h"

namespace Cyber{

  HoughObject::HoughObject( const Cyber::Calo1DCluster* _localmax, double _cellSize, double _ecal_inner_radius){
    m_local_max = _localmax;

    setCellSize(_cellSize);
    setCenterPoint(_ecal_inner_radius);
  }


  TVector2 HoughObject::getUpperPoint() const {    
    // Barrel 
    if(m_local_max->getBars()[0]->getSystem()==CaloUnit::System_Barrel){  
      // Slayer 0, or U bars. Bars perpendicular to z axis
      if(m_local_max->getSlayer()==0){
        if(m_center_point.Y()>=0){
          return m_center_point + TVector2(-0.5*m_cell_size, 0.5*m_cell_size);
        }
        else{
          return m_center_point + TVector2(0.5*m_cell_size, 0.5*m_cell_size);
        }
      }
      // Slayer 1, or V bars. Bars parallel to z axis
      else if(m_local_max->getSlayer()==1){
        return m_center_point + TVector2( 0.5*m_cell_size*cos(m_center_point.Phi() + TMath::PiOver2()), 0.5*m_cell_size*sin(m_center_point.Phi() + TMath::PiOver2()) );
      }
      else{
        std::cout<<"Error: Slayer="<<m_local_max->getSlayer()<<", do not use getUpperPoint()!"<<std::endl;
        return TVector2(0,0);
      }
    }

    // Endcap
    else if(m_local_max->getBars()[0]->getSystem()==CaloUnit::System_Endcap){
      if(m_center_point.X()>=0){
        if(m_center_point.Y()>=0){
          return m_center_point + TVector2(-0.5*m_cell_size, 0.5*m_cell_size);
        }
        else{
          return m_center_point + TVector2(0.5*m_cell_size, 0.5*m_cell_size);
        }
      }
      else{
        if(m_center_point.Y()>=0){
          return m_center_point + TVector2(0.5*m_cell_size, 0.5*m_cell_size);
        }
        else{
          return m_center_point + TVector2(-0.5*m_cell_size, 0.5*m_cell_size);
        }
      }
    }

    else{
      std::cout<<"Error: System="<<m_local_max->getBars()[0]->getSystem()<<", do not use getUpperPoint()!"<<std::endl;
      return TVector2(0,0);
    }
    
  }


  TVector2 HoughObject::getLowerPoint() const {
    // Barrel
    if(m_local_max->getBars()[0]->getSystem()==CaloUnit::System_Barrel){  
      // Slayer 0, or U bars. Bars perpendicular to z axis
      if(m_local_max->getSlayer()==0){
        if(m_center_point.Y()>=0){
          return m_center_point + TVector2(0.5*m_cell_size, -0.5*m_cell_size);
        }
        else{
          return m_center_point + TVector2(-0.5*m_cell_size, -0.5*m_cell_size);
        }
      }
      // Slayer 1, or V bars. Bars parallel to z axis
      else if(m_local_max->getSlayer()==1){
        return m_center_point + TVector2( 0.5*m_cell_size*cos(m_center_point.Phi() + 3*TMath::PiOver2()), 0.5*m_cell_size*sin(m_center_point.Phi() + 3*TMath::PiOver2()) );
      }
      else{
        std::cout<<"Error: Slayer="<<m_local_max->getSlayer()<<", do not use getLowerPoint()!"<<std::endl;
        return TVector2(0,0);
      }
    }

    // Endcap
    else if(m_local_max->getBars()[0]->getSystem()==CaloUnit::System_Endcap){
      if(m_center_point.X()>=0){
        if(m_center_point.Y()>=0){
          return m_center_point + TVector2(0.5*m_cell_size, -0.5*m_cell_size);
        }
        else{
          return m_center_point + TVector2(-0.5*m_cell_size, -0.5*m_cell_size);
        }
      }
      else{
        if(m_center_point.Y()>=0){
          return m_center_point + TVector2(-0.5*m_cell_size, -0.5*m_cell_size);
        }
        else{
          return m_center_point + TVector2(0.5*m_cell_size, -0.5*m_cell_size);
        }
      }
    }

    else{
      std::cout<<"Error: System="<<m_local_max->getBars()[0]->getSystem()<<", do not use getLowerPoint()!"<<std::endl;
      return TVector2(0,0);
    }

  }



  void HoughObject::setCenterPoint(double& _ecal_inner_radius){
    if(m_local_max->getBars()[0]->getSystem()==CaloUnit::System_Barrel){  // barrel bars
      if(m_local_max->getSlayer()==0){  // bars perpendicular to z axis
        TVector3 tmp_vec = m_local_max->getPos();
        m_center_point.Set(tmp_vec.Perp(), tmp_vec.z());
      }
      else if(m_local_max->getSlayer()==1){ // bars parallel to z axis
        m_center_point.Set(m_local_max->getPos().x(), m_local_max->getPos().y());
      }
      else{
        std::cout<<"Error: Slayer="<<m_local_max->getSlayer()<<", do not use setCenterPoint()!"<<std::endl;
      }
    }
    else if(m_local_max->getBars()[0]->getSystem()==CaloUnit::System_Endcap){ // endcap bars
      if(m_local_max->getSlayer()==0){ // bars parrallel to y axis
        m_center_point.Set(m_local_max->getPos().z(), m_local_max->getPos().x());
      }
      else if(m_local_max->getSlayer()==1){ // bars parallel to x axis
        m_center_point.Set(m_local_max->getPos().z(), m_local_max->getPos().y());
      }
      else{
        std::cout<<"Error: Slayer="<<m_local_max->getSlayer()<<", do not use setCenterPoint()!"<<std::endl;
      }
    }
    else{
      std::cout<<"Error: System="<<m_local_max->getBars()[0]->getSystem()<<", do not use setCenterPoint()!"<<std::endl;
    }
    
    
  }


  void HoughObject::setHoughLine(TF1& _line1, TF1& _line2){
    m_Hough_line_1 = _line1;  
    m_Hough_line_2 = _line2;
    //m_Hough_line_3 = _line3;  
    //m_Hough_line_4 = _line4;
  }



};
#endif
