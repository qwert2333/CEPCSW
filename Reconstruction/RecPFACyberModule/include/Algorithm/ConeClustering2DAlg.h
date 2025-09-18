#ifndef _CONECLUSTERING2D_ALG_H
#define _CONECLUSTERING2D_ALG_H

#include "CyberDataCol.h"
#include "Tools/Algorithm.h"
#include "TMath.h"

using namespace Cyber;
class ConeClustering2DAlg: public Cyber::Algorithm{
public: 

  ConeClustering2DAlg(){};
  ~ConeClustering2DAlg(){};

  class Factory : public Cyber::AlgorithmFactory
  {
  public: 
    Cyber::Algorithm* CreateAlgorithm() const{ return new ConeClustering2DAlg(); } 

  };

  StatusCode ReadSettings(Cyber::Settings& m_settings);
  StatusCode Initialize( CyberDataCol& m_datacol );
  StatusCode RunAlgorithm( CyberDataCol& m_datacol );
  StatusCode ClearAlgorithm();

  //Self defined algorithms
  StatusCode LongiConeLinking( std::map<int, std::vector<const Cyber::Calo1DCluster*> >& orderedShower, 
                               std::vector<Cyber::CaloHalfCluster*>& ClusterCol, 
                               std::vector<std::shared_ptr<Cyber::CaloHalfCluster>>& bk_HFclus );
  TVector2 GetProjectedAxis( const Cyber::CaloHalfCluster* m_shower );
  TVector2 GetProjectedRelR( const Cyber::Calo1DCluster* m_shower1, const Cyber::Calo1DCluster* m_shower2 );

private: 

  std::vector<Cyber::CaloHalfCluster*> p_HalfClusterV;
  std::vector<Cyber::CaloHalfCluster*> p_HalfClusterU;

  std::vector<const Cyber::Calo1DCluster*> m_localMaxVCol;
  std::vector<const Cyber::Calo1DCluster*> m_localMaxUCol;
  std::vector<const Cyber::CaloHalfCluster*> const_longiClusVCol;
  std::vector<const Cyber::CaloHalfCluster*> const_longiClusUCol;

};
#endif
