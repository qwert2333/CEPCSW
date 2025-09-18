#ifndef _PFORECLUSTERING_ALG_H
#define _PFORECLUSTERING_ALG_H

#include "Tools/Algorithm.h"

using namespace Cyber;
class PFOReclusteringAlg: public Cyber::Algorithm{
public: 

  PFOReclusteringAlg(){};
  ~PFOReclusteringAlg(){};

  class Factory : public Cyber::AlgorithmFactory
  {
  public: 
    Cyber::Algorithm* CreateAlgorithm() const{ return new PFOReclusteringAlg(); } 

  };

  StatusCode ReadSettings(Cyber::Settings& m_settings);
  StatusCode Initialize( CyberDataCol& m_datacol );
  StatusCode RunAlgorithm( CyberDataCol& m_datacol );
  StatusCode ClearAlgorithm();

  StatusCode ReCluster_MergeToChg( std::vector< std::shared_ptr<Cyber::PFObject> >& m_chargedPFOs,
                                   std::vector< std::shared_ptr<Cyber::PFObject> >& m_neutralPFOs );

  StatusCode ReCluster_SplitFromChg( std::vector< std::shared_ptr<Cyber::PFObject> >& m_chargedPFOs,
                                     std::vector< std::shared_ptr<Cyber::PFObject> >& m_neutralPFOs );

private: 
  
  std::vector<std::shared_ptr<Cyber::PFObject>>* p_PFObjects; 

  CyberDataCol m_bkCol;

  static bool compTrkP( std::shared_ptr<Cyber::PFObject> pfo1, std::shared_ptr<Cyber::PFObject> pfo2 )
    { return pfo1->getTrackMomentum() > pfo2->getTrackMomentum(); }
};
#endif
