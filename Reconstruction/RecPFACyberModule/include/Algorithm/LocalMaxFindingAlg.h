#ifndef _LOCALMAXFINDING_ALG_H
#define _LOCALMAXFINDING_ALG_H

#include <set>
#include "Tools/Algorithm.h"
using namespace Cyber;

class LocalMaxFindingAlg: public Cyber::Algorithm{
public: 

  LocalMaxFindingAlg(){};
  ~LocalMaxFindingAlg(){};

  class Factory : public Cyber::AlgorithmFactory
  {
  public: 
    Cyber::Algorithm* CreateAlgorithm() const { return new LocalMaxFindingAlg(); }

  };

  StatusCode ReadSettings(Cyber::Settings& m_settings);
  StatusCode Initialize( CyberDataCol& m_datacol );
  StatusCode RunAlgorithm( CyberDataCol& m_datacol );
  StatusCode ClearAlgorithm();

  StatusCode GetLocalMax( const Cyber::Calo1DCluster* m_1dClus, std::vector<std::shared_ptr<Cyber::Calo1DCluster>>& m_output);
  StatusCode GetLocalMaxBar( std::vector<const Cyber::CaloUnit*>& barCol, std::vector<const Cyber::CaloUnit*>& localMaxCol );
  std::vector<const Cyber::CaloUnit*>  getNeighbors(const Cyber::CaloUnit* seed, std::vector<const Cyber::CaloUnit*>& barCol);

private: 

  std::vector<std::shared_ptr<Cyber::CaloHalfCluster>>* p_HalfClusU = nullptr;
  std::vector<std::shared_ptr<Cyber::CaloHalfCluster>>* p_HalfClusV = nullptr;

  static bool compBar( const Cyber::CaloUnit* bar1, const Cyber::CaloUnit* bar2 )
      { return bar1->getBar() < bar2->getBar(); }

};
#endif
