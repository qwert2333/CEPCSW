#ifndef _CONECLUSTERING2D_ALG_H
#define _CONECLUSTERING2D_ALG_H

#include "PandoraPlusDataCol.h"
#include "Tools/Algorithm.h"
#include "TMath.h"

using namespace PandoraPlus;
class ConeClustering2DAlg: public PandoraPlus::Algorithm{
public: 

  ConeClustering2DAlg(){};
  ~ConeClustering2DAlg(){};

  class Factory : public PandoraPlus::AlgorithmFactory
  {
  public: 
    PandoraPlus::Algorithm* CreateAlgorithm() const{ return new ConeClustering2DAlg(); } 

  };

  StatusCode ReadSettings(PandoraPlus::Settings& m_settings);
  StatusCode Initialize( PandoraPlusDataCol& m_datacol );
  StatusCode RunAlgorithm( PandoraPlusDataCol& m_datacol );
  StatusCode ClearAlgorithm();

  //Self defined algorithms
  StatusCode LongiConeLinking( std::map<int, std::vector<const PandoraPlus::Calo1DCluster*> >& orderedShower, 
                               std::vector<PandoraPlus::CaloHalfCluster*>& ClusterCol, 
                               std::vector<std::shared_ptr<PandoraPlus::CaloHalfCluster>>& bk_HFclus );
  TVector2 GetProjectedAxis( const PandoraPlus::CaloHalfCluster* m_shower );
  TVector2 GetProjectedRelR( const PandoraPlus::Calo1DCluster* m_shower1, const PandoraPlus::Calo1DCluster* m_shower2 );

private: 

  std::vector<PandoraPlus::CaloHalfCluster*> p_HalfClusterV;
  std::vector<PandoraPlus::CaloHalfCluster*> p_HalfClusterU;

  std::vector<const PandoraPlus::Calo1DCluster*> m_localMaxVCol;
  std::vector<const PandoraPlus::Calo1DCluster*> m_localMaxUCol;
  std::vector<const PandoraPlus::CaloHalfCluster*> const_longiClusVCol;
  std::vector<const PandoraPlus::CaloHalfCluster*> const_longiClusUCol;

};
#endif
