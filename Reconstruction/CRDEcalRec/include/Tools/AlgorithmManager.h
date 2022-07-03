#ifndef ALGORITHM_MANAGER_H
#define ALGORITHM_MANAGER_H

#include "Tools/Algorithm.h"

namespace PandoraPlus{

  class Algorithm;
  class AlgorithmFactory;

  class AlgorithmManager{
  public: 

    AlgorithmManager(){};
    ~AlgorithmManager(){ Clean(); };
   
    void Clean(){
      m_algorithmNames.clear();
      for(auto iter : m_algorithmFactoryMap) delete iter.second;
      for(auto iter : m_algorithmMap) delete iter.second;
   
      m_algorithmFactoryMap.clear();
      m_algorithmMap.clear();
    }
   
    StatusCode RegisterAlgorithmFactory(const std::string& algorithmName, AlgorithmFactory *const algorithmFactory)
    {
      if(!m_algorithmFactoryMap.insert(AlgorithmFactoryMap::value_type(algorithmName, algorithmFactory)).second )
        return StatusCode::FAILURE; 
   
      return StatusCode::SUCCESS;
    }
   
    StatusCode RegisterAlgorithm( const std::string& algorithmName, Settings& m_algorithmSetting ){
      m_algorithmNames.push_back(algorithmName);
      m_algorithmMap.insert(AlgorithmMap::value_type(algorithmName, m_algorithmFactoryMap[algorithmName]->CreateAlgorithm()) );
      m_algorithmMap[algorithmName]->ReadSettings(m_algorithmSetting);
   
      return StatusCode::SUCCESS;
    }
   
    //StatusCode CreateAlgorithm(const std::map<std::string, Settings>& m_algorithmSettings); 
    StatusCode RunAlgorithm( PandoraPlusDataCol& m_datacol ){
      for(auto iter : m_algorithmNames){
cout<<"Processing Algorithm: "<<iter<<endl;
        m_algorithmMap[iter]->Initialize();
        m_algorithmMap[iter]->RunAlgorithm(m_datacol);
        m_algorithmMap[iter]->ClearAlgorithm();
      }
      return StatusCode::SUCCESS;
    }

  private: 
    typedef std::map<const std::string, Algorithm *const> AlgorithmMap;
    typedef std::map< std::string, AlgorithmFactory *const> AlgorithmFactoryMap;

    AlgorithmMap              m_algorithmMap;
    AlgorithmFactoryMap       m_algorithmFactoryMap;
    std::vector<std::string>  m_algorithmNames;

  };

};

#endif
