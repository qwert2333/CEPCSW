#ifndef _EXAMPLE_ALG_C
#define _EXAMPLE_ALG_C

#include "Algorithm/ExampleAlg.h"

StatusCode ExampleAlg::ReadSettings(Settings& m_settings){
  settings = m_settings;

  //Initialize parameters
  if(settings.map_floatPars.find("Par1")==settings.map_floatPars.end()) settings.map_floatPars["Par1"] = 0;
  return StatusCode::SUCCESS;
};

StatusCode ExampleAlg::Initialize(){
  std::cout<<"Initialize ExampleAlg"<<std::endl;

  return StatusCode::SUCCESS;
};

StatusCode ExampleAlg::RunAlgorithm( PandoraPlusDataCol& m_datacol ){
  std::cout<<"Excute ExampleAlg"<<std::endl;
  std::cout<<"  Run SelfAlg1"<<std::endl;
  SelfAlg1();

  return StatusCode::SUCCESS;
};

StatusCode ExampleAlg::ClearAlgorithm(){
  std::cout<<"End run ExampleAlg. Clean it."<<std::endl;

  return StatusCode::SUCCESS;
};


StatusCode ExampleAlg::SelfAlg1(){
  std::cout<<"  Processing SelfAlg1: print Par1 = "<<settings.map_floatPars["Par1"]<<std::endl;

  return StatusCode::SUCCESS;
};


#endif
