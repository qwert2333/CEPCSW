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
    void SetValues( double _Sth_split=200, double _Eth_seedAbs=0.02, double _Eth_showerAbs=0.03, double _Eth_clusAbs=0.03,
                    double _Eth_seedwnei=0.4, double _Eth_seedwtot=0., double _Eth_showerwtot=0.01, double _Eth_cluswtot=0.01, bool _isforced=true );

    double Sth_split;
    double Eth_SeedAbs;
    double Eth_ShowerAbs;
    double Eth_ClusAbs;
    double Eth_SeedWithNeigh;
    double Eth_SeedWithTot;
    double Eth_ShowerWithTot;
    double Eth_ClusterWithTot;
    double Etot;

    bool isForced;
    std::map<int, int> map_mode; // map<Dlayer, ForcedMode> 

  };

  EnergySplittingAlg();
  ~EnergySplittingAlg(){};

  StatusCode Initialize();

  StatusCode RunAlgorithm( EnergySplittingAlg::Settings& m_settings, PandoraPlusDataCol& m_datacol);

  CRDEcalEDM::CRDCaloLayer ClusterinLayer( CRDEcalEDM::CRDCaloBlock& m_blocks);

  std::vector<CRDEcalEDM::CRDCaloBarCluster> Clustering(std::vector<CRDEcalEDM::CRDCaloBar>& barCol);

  bool ClusteringWithExpShower( std::vector<CRDEcalEDM::CRDCaloBar>& barCol,  std::vector<CRDEcalEDM::CRDCaloBarCluster>& outClus, const std::vector<CRDEcalEDM::CRDExpEMShower>& ExpEMShowers );

  std::vector<CRDEcalEDM::CRDExpEMShower> ExpShowerInCluster( CRDEcalEDM::CRDCaloBarCluster& cluster, const std::vector<CRDEcalEDM::CRDExpEMShower>& ExpEMShowers);

  int NMatchedSeedShower( std::vector<CRDEcalEDM::CRDCaloBar>& seeds,  const std::vector<CRDEcalEDM::CRDExpEMShower>& ExpEMShowers );

  std::vector<CRDEcalEDM::CRDCaloBar>  findSeeds( CRDEcalEDM::CRDCaloBarCluster& m_cluster, int NforcedSeed = -1 );

  std::vector<CRDEcalEDM::CRDCaloBar>  getNeighbors(CRDEcalEDM::CRDCaloBarCluster& m_cluster, CRDEcalEDM::CRDCaloBar& seed);

  std::vector<CRDEcalEDM::CRDCaloBarShower>  ClusterSplitting( CRDEcalEDM::CRDCaloBarCluster& m_cluster);
  std::vector<CRDEcalEDM::CRDCaloBarShower>  ClusterSplittingWithExpShower( CRDEcalEDM::CRDCaloBarCluster& m_cluster, const std::vector<CRDEcalEDM::CRDExpEMShower>& m_expShowers );

  void CalculateInitialEseed( const std::vector<CRDEcalEDM::CRDCaloBar>& Seeds, const dd4hep::Position* pos, double* Eseed);

  double GetShowerProfile(const dd4hep::Position& p_bar, const dd4hep::Position& p_seed );

  static bool compE( CRDEcalEDM::CRDCaloBar& bar1, CRDEcalEDM::CRDCaloBar& bar2 ) { return bar1.getEnergy()>bar2.getEnergy(); }




  Settings settings;

private: 

};


#endif
