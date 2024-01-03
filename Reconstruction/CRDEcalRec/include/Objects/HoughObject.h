#ifndef CALOHOUGHOBJECT_H
#define CALOHOUGHOBJECT_H

#include "Objects/Calo1DCluster.h"
#include "TVector2.h"
#include "TF1.h"


namespace PandoraPlus {

  class HoughObject{
  public:
    HoughObject() {};
    HoughObject( const PandoraPlus::Calo1DCluster* _localmax, double _cellSize, double _ecal_inner_radius);
    ~HoughObject() { };


    TVector2 getCenterPoint() const { return m_center_point; }
    TVector2 getPointU()  const { return m_center_point + TVector2(0,  m_cell_size/TMath::Sqrt(2)); }
    TVector2 getPointD()  const { return m_center_point + TVector2(0, -m_cell_size/TMath::Sqrt(2)); }
    TVector2 getPointL()  const { return m_center_point + TVector2(-m_cell_size/TMath::Sqrt(2), 0); }
    TVector2 getPointR()  const { return m_center_point + TVector2( m_cell_size/TMath::Sqrt(2), 0); }
    TVector2 getPointUR() const { return m_center_point + TVector2( m_cell_size/2,  m_cell_size/2); }
    TVector2 getPointDL() const { return m_center_point + TVector2(-m_cell_size/2, -m_cell_size/2); }
    TVector2 getPointUL() const { return m_center_point + TVector2(-m_cell_size/2,  m_cell_size/2); }
    TVector2 getPointDR() const { return m_center_point + TVector2( m_cell_size/2, -m_cell_size/2); }

    TF1 getHoughLine1() const { return m_Hough_line_1; }
    TF1 getHoughLine2() const { return m_Hough_line_2; }
    TF1 getHoughLine3() const { return m_Hough_line_3; }
    TF1 getHoughLine4() const { return m_Hough_line_4; }

    int getModule() const { return (m_local_max->getTowerID())[0][0]; }
    int getSlayer() const { return m_local_max->getSlayer(); }
    double getE() const { return m_local_max->getEnergy(); }
    double getCellSize() const { return m_cell_size; }
    const PandoraPlus::Calo1DCluster* getLocalMax() const { return m_local_max; }

    void setCellSize(double _cs) { m_cell_size=_cs; }
    void setCenterPoint(double& _ecal_inner_radius);
    void setHoughLine(int module, TF1& line1, TF1& line2, TF1& line3, TF1& line4);


  private:
    const PandoraPlus::Calo1DCluster* m_local_max;  //Local max
    double m_cell_size;
    TVector2 m_center_point;  // center position

    TF1 m_Hough_line_1;    // ur or u
    TF1 m_Hough_line_2;    // dl or d
    TF1 m_Hough_line_3;    // ul or l
    TF1 m_Hough_line_4;    // dr or r
    // The above conversion is only for the octagon barrel ECAL

};

};
#endif
