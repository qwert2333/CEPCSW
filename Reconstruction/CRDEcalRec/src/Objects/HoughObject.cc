#ifndef CALOHOUGHOBJECT_C
#define CALOHOUGHOBJECT_C

#include "Objects/HoughObject.h"

namespace PandoraPlus{

  HoughObject::HoughObject( const PandoraPlus::Calo1DCluster* _localmax, double _cellSize, double _ecal_inner_radius){
    m_local_max = _localmax;

    setCellSize(_cellSize);
    setCenterPoint(_ecal_inner_radius);
  }


  void HoughObject::setCenterPoint(double& _ecal_inner_radius){
    if(m_local_max->getSlayer()==0){
        m_center_point.Set( (m_local_max->getDlayer()-1)*20. + _ecal_inner_radius + m_cell_size*0.5 , 
                          m_local_max->getPos().z());
    }
    else if(m_local_max->getSlayer()==1){
        m_center_point.Set(m_local_max->getPos().x(), m_local_max->getPos().y());
    }
    else{
        std::cout<<"Error: Slayer="<<m_local_max->getSlayer()<<", do not use setCenterPoint()!"<<std::endl;
    }
  }


  void HoughObject::setHoughLine(int _module, TF1& _line1, TF1& _line2, TF1& _line3, TF1& _line4){
    m_Hough_line_1 = _line1;  
    m_Hough_line_2 = _line2;
    m_Hough_line_3 = _line3;  
    m_Hough_line_4 = _line4;
  }



};
#endif