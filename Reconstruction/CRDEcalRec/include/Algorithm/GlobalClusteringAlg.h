#ifndef GLOBALCLUSTERING_ALG_H
#define GLOBALCLUSTERING_ALG_H

#include "Tools/Algorithm.h"

#include "time.h"
#include <TTimeStamp.h> 
#include <ctime>

#include <cstdlib>

using namespace PandoraPlus;

class GlobalClusteringAlg : public PandoraPlus::Algorithm{
public: 

  GlobalClusteringAlg(){};
  ~GlobalClusteringAlg(){};
 
  class Factory : public PandoraPlus::AlgorithmFactory
  {
  public:
    PandoraPlus::Algorithm* CreateAlgorithm() const{ return new GlobalClusteringAlg(); }
 
  };
 
  StatusCode ReadSettings(PandoraPlus::Settings& m_settings); 
  StatusCode Initialize( PandoraPlusDataCol& m_datacol );
  StatusCode RunAlgorithm( PandoraPlusDataCol& m_datacol );
  StatusCode ClearAlgorithm();

  template<typename T1, typename T2> StatusCode Clustering(std::vector<T1*>& m_input, std::vector<T2*>& m_output);
  
  //geometry construction
/*   int m_module = settings.map_intPars["m_module"];
  int m_modulestart = settings.map_intPars["m_modulestart"];
  int m_part = settings.map_intPars["m_part"];
  int m_stave = settings.map_intPars["m_stave"];
  int m_superlayer = settings.map_intPars["m_superlayer"];
  int m_startnumber = settings.map_intPars["m_startnumber"];
  int m_phibarnumber = settings.map_intPars["m_phibarnumber"];
  int m_zbarnumber = settings.map_intPars["m_zbarnumber"]; */

private:

  std::vector<PandoraPlus::CaloUnit*> m_bars; 
  std::vector<PandoraPlus::CaloUnit*> m_processbars;        
  std::vector<PandoraPlus::CaloUnit*> m_restbars;           
  std::vector<PandoraPlus::Calo1DCluster*> m_1dclusters;    
  std::vector<PandoraPlus::CaloHalfCluster*> m_halfclusters; 
  //std::vector<PandoraPlus::Calo2DCluster*> m_2dclusters;     m_2dclusters.clear();
  //std::vector<PandoraPlus::Calo3DCluster*> m_3dclusters;     m_3dclusters.clear();
  
};
#endif
