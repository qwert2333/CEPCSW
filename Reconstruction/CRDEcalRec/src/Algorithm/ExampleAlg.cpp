#ifndef _EXAMPLE_ALG_C
#define _EXAMPLE_ALG_C

#include "Algorithm/ExampleAlg.h"

StatusCode ExampleAlg::ReadSettings(Settings& m_settings){
  settings = m_settings;

  //Initialize parameters
  if(settings.map_stringPars.find("Par1")==settings.map_stringPars.end()) settings.map_stringPars["Par1"] = "HERE";
  if(settings.map_floatPars.find("Par2")==settings.map_floatPars.end()) settings.map_floatPars["Par2"] = 0;
  if(settings.map_boolPars.find("Par3")==settings.map_boolPars.end()) settings.map_boolPars["Par3"] = 0;

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
  std::cout<<"  Processing SelfAlg1: print Setting parameters"<<std::endl;

  for(auto &iter : settings.map_stringPars) cout<<iter.first<<": "<<iter.second<<endl;
  for(auto &iter : settings.map_floatPars) cout<<iter.first<<": "<<iter.second<<endl;
  for(auto &iter : settings.map_boolPars) cout<<iter.first<<": "<<iter.second<<endl;


  return StatusCode::SUCCESS;
};


#endif
