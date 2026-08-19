#ifndef CALOHOUGHOBJECT_C
#define CALOHOUGHOBJECT_C

#include "Objects/HoughObject.h"

namespace PandoraPlus{

  HoughObject::HoughObject( const PandoraPlus::Calo1DCluster* _localmax, double _cellSize, double _ecal_inner_radius, double _phi){
    m_local_max = _localmax;

    setCellSize(_cellSize);
    setCenterPoint(_ecal_inner_radius, _phi);
  }


  void HoughObject::setCenterPoint(double& _ecal_inner_radius, double _phi){
    if(m_local_max->getSlayer()==0){
      if(_phi==0)
        m_center_point.Set( (m_local_max->getDlayer()-1)*20. + _ecal_inner_radius + m_cell_size*0.5 , 
                          m_local_max->getPos().z());
      else{
        double intPart, fracPart;
        fracPart = modf((_phi+TMath::Pi())/(TMath::Pi()/4.), &intPart);  // yyy: _phi + TMath::Pi() ranges from 0 to 2pi
        if(fracPart<0.489 || fracPart>0.711)  //Not in crack region.
          m_center_point.Set( (m_local_max->getDlayer()-1)*20. + _ecal_inner_radius + m_cell_size*0.5 ,
                          m_local_max->getPos().z());
        else{
          int iCrack = intPart+2;
          if(iCrack>=8) iCrack = iCrack-8;

          double tmp_phi = _phi;
          while(tmp_phi<0.) tmp_phi += TMath::Pi() / 4.;  
          while(tmp_phi>=TMath::Pi() / 4.) tmp_phi -= TMath::Pi() / 4.;

          int imodule = m_local_max->getTowerID()[0][0];
          if(imodule==iCrack){
            double Rref = _ecal_inner_radius/cos( tmp_phi );
            m_center_point.Set( (m_local_max->getDlayer()-1)*20. + Rref + m_cell_size*0.5 ,
                          m_local_max->getPos().z());
          }
          else{
            double Rref = _ecal_inner_radius/cos( TMath::Pi() / 4. - tmp_phi );
            m_center_point.Set( (m_local_max->getDlayer()-1)*20. + Rref + m_cell_size*0.5 ,
                          m_local_max->getPos().z());
          }

        }
      }
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
