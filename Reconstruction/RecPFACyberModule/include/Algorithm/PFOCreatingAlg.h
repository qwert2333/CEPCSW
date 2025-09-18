#ifndef _PFOCREATING_ALG_H
#define _PFOCREATING_ALG_H

#include "Tools/Algorithm.h"

using namespace Cyber;
class PFOCreatingAlg: public Cyber::Algorithm{
public: 

  PFOCreatingAlg(){};
  ~PFOCreatingAlg(){};

  class Factory : public Cyber::AlgorithmFactory
  {
  public: 
    Cyber::Algorithm* CreateAlgorithm() const{ return new PFOCreatingAlg(); } 

  };

  StatusCode ReadSettings( Cyber::Settings& m_settings );
  StatusCode Initialize( CyberDataCol& m_datacol );
  StatusCode RunAlgorithm( CyberDataCol& m_datacol );
  StatusCode ClearAlgorithm();

  std::vector<Cyber::Track*> getTracks() const { return m_tracks; }
  std::vector<Cyber::Calo3DCluster*> getECALClusters() const { return m_ecal_clusters; }
  std::vector<Cyber::Calo3DCluster*> getHCALClusters() const { return m_hcal_clusters; }

  //Self defined algorithms
  // Get canditate clusters in HCAL for charged particles
  StatusCode GetChargedHCALCandidates(const Cyber::Track* _track,
                                      std::vector<Cyber::Calo3DCluster*>& _hcal_clusters,
                                      std::vector<Cyber::Calo3DCluster*>& _hcal_clus_candidate);
  // Get nearby HCAL clusters
  StatusCode GetNearbyHCALCandidates( Cyber::Calo3DCluster* _ecal_cluster,
                                  std::vector<Cyber::Calo3DCluster*>& _hcal_clusters,
                                  std::vector<Cyber::Calo3DCluster*>& _hcal_clus_candidate);
  // If a neutral cluster in ECAL reach the outermost ECAL boundary
  bool isReachOuterMostECAL(Cyber::Calo3DCluster* _ecal_cluster);
  // erase the used_elements in the left_elements
  template<typename T1, typename T2> StatusCode CleanUsedElements(std::vector<T1>& _used_elements,
                                                    std::vector<T2>& _left_elements);
  template<typename T1, typename T2> StatusCode CleanUsedElement(T1 _used_elements,
                                                    std::vector<T2>& _left_elements);
  // Create PFO with:
  //   1. tracks with no clusters in ECAL and HCAL
  //   2. HCAL clusters
  StatusCode CreateLeftPFO( std::vector<Cyber::Track*>& _tracks,
                            std::vector<Cyber::Calo3DCluster*>& _hcal_clusters,
                            std::vector<std::shared_ptr<Cyber::PFObject>>& _pfobjects);
  
  

private: 
  std::vector<Cyber::Track*> m_tracks;
  std::vector<Cyber::Calo3DCluster*> m_ecal_clusters;
  std::vector<Cyber::Calo3DCluster*> m_hcal_clusters;

  std::vector<std::shared_ptr<Cyber::PFObject>> m_pfobjects;
  

};

#endif
