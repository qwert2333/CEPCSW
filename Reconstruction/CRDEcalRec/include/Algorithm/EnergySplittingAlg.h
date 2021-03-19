#ifndef _ENERGYSPLITTING_ALG_H
#define _ENERGYSPLITTING_ALG_H

#include "PandoraPlusDataCol.h"

#include "TVector3.h"
#include "TVector.h"
#include "TMatrix.h"

using namespace CRDEcalEDM;

class EnergySplittingAlg {
public: 

  class Settings {
  public: 
    Settings() {};
    void SetInitialValue();

    double Sth_split;
    double Eth_SeedAbs;
    double Eth_ShowerAbs;
    double Eth_ClusAbs;
    double Eth_SeedWithNeigh;
    double Eth_SeedWithTot;
    double Eth_ShowerWithTot;
    double Eth_ClusterWithTot;
    double Etot;

  };

  EnergySplittingAlg();

  StatusCode Initialize();

  StatusCode RunAlgorithm( EnergySplittingAlg::Settings& m_settings, PandoraPlusDataCol& m_datacol);

  CRDEcalEDM::CRDCaloLayer ClusterinLayer( CRDEcalEDM::DigiBlock& m_blocks);

  std::vector<CRDEcalEDM::CRDCaloBarCluster> Clustering(std::vector<CRDEcalEDM::CRDCaloBar>& barCol);

  std::vector<CRDEcalEDM::CRDCaloBar>  findSeeds( CRDEcalEDM::CRDCaloBarCluster& m_cluster );

  std::vector<CRDEcalEDM::CRDCaloBar>  getNeighbors(CRDEcalEDM::CRDCaloBarCluster& m_cluster, CRDEcalEDM::CRDCaloBar& seed);

  std::vector<CRDEcalEDM::CRDCaloBarShower>  ClusterSplitting( CRDEcalEDM::CRDCaloBarCluster& m_cluster);

  void CalculateInitialEseed( const std::vector<CRDEcalEDM::CRDCaloBar>& Seeds, const dd4hep::Position* pos, double* Eseed);

  double GetShowerProfile(const dd4hep::Position& p_bar, const dd4hep::Position& p_seed );

  Settings settings;

private: 

};

#endif
