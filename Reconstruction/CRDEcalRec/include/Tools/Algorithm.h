#ifndef ALGORITHM_TEMP_H
#define ALGORITHM_TEMP_H

#include "GaudiAlg/GaudiAlgorithm.h"

#include "PandoraPlusDataCol.h"

namespace PandoraPlus{

  class Settings{
  public: 
    Settings(){};
    ~Settings(){ Clear(); }
   
    void Clear() { map_intPars.clear();  map_floatPars.clear(); map_boolPars.clear(); map_stringPars.clear(); map_stringVecPars.clear(); }
    std::map<std::string, int>    map_intPars; 
    std::map<std::string, double> map_floatPars; 
    std::map<std::string, bool>   map_boolPars; 
    std::map<std::string, std::string> map_stringPars; 
    std::map<std::string, std::vector<std::string> > map_stringVecPars;

  };

  class Algorithm{
  public:
    Algorithm() {};
    virtual ~Algorithm() {};

    virtual StatusCode ReadSettings(Settings& m_settings) = 0;
    virtual StatusCode Initialize() = 0;
    virtual StatusCode RunAlgorithm( PandoraPlusDataCol& m_datacol ) = 0;
    virtual StatusCode ClearAlgorithm() = 0;


    Settings settings; 
  };

  class AlgorithmFactory{
  public:

    virtual ~AlgorithmFactory(){}; 
    virtual Algorithm* CreateAlgorithm() const = 0;
  };

};
#endif
