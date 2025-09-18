#ifndef _TRACKCLUSTERCONNECTING_ALG_H
#define _TRACKCLUSTERCONNECTING_ALG_H

#include "Tools/Algorithm.h"

using namespace Cyber;
class TrackClusterConnectingAlg: public Cyber::Algorithm{
public: 

  TrackClusterConnectingAlg(){};
  ~TrackClusterConnectingAlg(){};

  class Factory : public Cyber::AlgorithmFactory
  {
  public: 
    Cyber::Algorithm* CreateAlgorithm() const{ return new TrackClusterConnectingAlg(); } 

  };

  StatusCode ReadSettings(Cyber::Settings& m_settings);
  StatusCode Initialize( CyberDataCol& m_datacol );
  StatusCode RunAlgorithm( CyberDataCol& m_datacol );
  StatusCode ClearAlgorithm();

  //Self defined algorithms
  double GetMinR2Trk( const Cyber::Calo3DCluster* p_clus, const Cyber::Track* m_trk);
  StatusCode PFOCreating( std::vector<const Cyber::Calo3DCluster*>& m_clusters,
                          std::vector<const Cyber::Track*>& m_trks,
                          std::vector<std::shared_ptr<Cyber::PFObject>>& m_PFOs );

  StatusCode EcalChFragAbsorption( std::vector<const Cyber::Calo3DCluster*>& m_clusters,
                                   std::vector<const Cyber::Track*>& m_trks,
                                   std::vector<std::shared_ptr<Cyber::Calo3DCluster>>& m_newclusCol ); 

  StatusCode HcalExtrapolatingMatch(std::vector<const Cyber::Calo3DCluster*>& m_clusters, std::vector<std::shared_ptr<Cyber::PFObject>>& m_PFOs);

private: 

  std::vector<const Cyber::Calo3DCluster*> m_EcalClusters;
  std::vector<const Cyber::Calo3DCluster*> m_HcalClusters;
  std::vector<const Cyber::Track*> m_tracks; 

  std::vector<std::shared_ptr<Cyber::Calo3DCluster>> m_absorbedEcal; 
  std::vector<std::shared_ptr<Cyber::PFObject>> m_PFObjects; 

  CyberDataCol m_bkCol;

  static bool compTrkP( std::shared_ptr<Cyber::PFObject> pfo1, std::shared_ptr<Cyber::PFObject> pfo2 )
    { return pfo1->getTrackMomentum() > pfo2->getTrackMomentum(); }

};
#endif
