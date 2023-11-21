#ifndef _PFOCREATING_ALG_H
#define _PFOCREATING_ALG_H

#include "Tools/Algorithm.h"

using namespace PandoraPlus;
class PFOCreatingAlg: public PandoraPlus::Algorithm{
public: 

  PFOCreatingAlg(){};
  ~PFOCreatingAlg(){};

  class Factory : public PandoraPlus::AlgorithmFactory
  {
  public: 
    PandoraPlus::Algorithm* CreateAlgorithm() const{ return new PFOCreatingAlg(); } 

  };

  StatusCode ReadSettings( PandoraPlus::Settings& m_settings );
  StatusCode Initialize( PandoraPlusDataCol& m_datacol );
  StatusCode RunAlgorithm( PandoraPlusDataCol& m_datacol );
  StatusCode ClearAlgorithm();

  std::vector<PandoraPlus::Track*> getTracks() const { return m_tracks; }
  std::vector<PandoraPlus::Calo3DCluster*> getECALClusters() const { return m_ecal_clusters; }
  std::vector<PandoraPlus::Calo3DCluster*> getHCALClusters() const { return m_hcal_clusters; }

  //Self defined algorithms
  // Get canditate clusters in HCAL for charged particles
  StatusCode GetChargedHCALCandidates(const PandoraPlus::Track* _track,
                                      std::vector<PandoraPlus::Calo3DCluster*>& _hcal_clusters,
                                      std::vector<PandoraPlus::Calo3DCluster*>& _hcal_clus_candidate);
  // Get nearby HCAL clusters
  StatusCode GetNearbyHCALCandidates( PandoraPlus::Calo3DCluster* _ecal_cluster,
                                  std::vector<PandoraPlus::Calo3DCluster*>& _hcal_clusters,
                                  std::vector<PandoraPlus::Calo3DCluster*>& _hcal_clus_candidate);
  // If a neutral cluster in ECAL reach the outermost ECAL boundary
  bool isReachOuterMostECAL(PandoraPlus::Calo3DCluster* _ecal_cluster);
  // erase the used_elements in the left_elements
  template<typename T1, typename T2> StatusCode CleanUsedElements(std::vector<T1>& _used_elements,
                                                    std::vector<T2>& _left_elements);
  template<typename T1, typename T2> StatusCode CleanUsedElement(T1 _used_elements,
                                                    std::vector<T2>& _left_elements);
  // Create PFO with:
  //   1. tracks with no clusters in ECAL and HCAL
  //   2. HCAL clusters
  StatusCode CreateLeftPFO( std::vector<PandoraPlus::Track*>& _tracks,
                            std::vector<PandoraPlus::Calo3DCluster*>& _hcal_clusters,
                            std::vector<std::shared_ptr<PandoraPlus::PFObject>>& _pfobjects);
  
  

private: 
  std::vector<PandoraPlus::Track*> m_tracks;
  std::vector<PandoraPlus::Calo3DCluster*> m_ecal_clusters;
  std::vector<PandoraPlus::Calo3DCluster*> m_hcal_clusters;

  std::vector<std::shared_ptr<PandoraPlus::PFObject>> m_pfobjects;
  

};

#endif
