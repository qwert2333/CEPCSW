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
    void SetValues( double _Sth_split=-1, double _Eth_seedAbs=0.005, double _Eth_showerAbs=0.0, double _Eth_clusAbs=0.0,
                    double _Eth_seedwnei=0.5, double _Eth_seedwtot=0., double _Eth_showerwtot=0., double _Eth_cluswtot=0., 
                    bool _usecandi=true, int _debug=0 );
    void SetScndMomentThres( double _th ){ Sth_split = _th; }
    void SetEseedThres( double _th ) { Eth_SeedAbs = _th; }
    void SetEseedRel( double _th ) { Eth_SeedWithNeigh = _th; }
    void SetUseCandidate( bool _flag ) { UseCandidate = _flag; }
    void SetDebug(int _lv) { Debug=_lv; }

    double Sth_split;
    double Eth_SeedAbs;
    double Eth_ShowerAbs;
    double Eth_ClusAbs;
    double Eth_SeedWithNeigh;
    double Eth_SeedWithTot;
    double Eth_ShowerWithTot;
    double Eth_ClusterWithTot;
    double Etot;
    int    Debug;

    bool UseCandidate;
    std::map<int, int> map_mode; // map<Dlayer, ForcedMode> 

  };

  EnergySplittingAlg();
  ~EnergySplittingAlg(){};

  StatusCode Initialize();

  StatusCode RunAlgorithm( EnergySplittingAlg::Settings& m_settings, PandoraPlusDataCol& m_datacol);

  CRDEcalEDM::CRDCaloLayer ClusterinLayer( CRDEcalEDM::CRDCaloBlock& m_blocks);

  bool Clustering( std::vector<CRDEcalEDM::CRDCaloBar>& barCol,  std::vector<CRDEcalEDM::CRDCaloBarCluster>& outClus);

  bool Clustering( std::vector<CRDEcalEDM::CRDCaloBar>& barCol,  std::vector<CRDEcalEDM::CRDCaloBarCluster>& outClus, const std::vector<CRDEcalEDM::CRDShowerCandidate>& m_CandiCol );

  std::vector<CRDEcalEDM::CRDShowerCandidate> CandidateInCluster( CRDEcalEDM::CRDCaloBarCluster& cluster, const std::vector<CRDEcalEDM::CRDShowerCandidate>& m_CandiCol);

  std::vector<CRDEcalEDM::CRDCaloBar> GetMatchedSeeds( std::vector<CRDEcalEDM::CRDCaloBar>& seeds, const std::vector<CRDEcalEDM::CRDShowerCandidate>& m_CandiCol );

  std::vector<CRDEcalEDM::CRDCaloBar>  findSeeds( CRDEcalEDM::CRDCaloBarCluster& m_cluster);

  std::vector<CRDEcalEDM::CRDCaloBar>  findSeeds( CRDEcalEDM::CRDCaloBarCluster& m_cluster, const std::vector<CRDEcalEDM::CRDShowerCandidate>& m_CandiCol, bool f_force );

  std::vector<CRDEcalEDM::CRDCaloBar>  getNeighbors(CRDEcalEDM::CRDCaloBarCluster& m_cluster, CRDEcalEDM::CRDCaloBar& seed);

  bool  ClusterSplitting( CRDEcalEDM::CRDCaloBarCluster& m_cluster, std::vector<CRDEcalEDM::CRDCaloBarShower>& outshCol  );
  bool  ClusterSplitting( CRDEcalEDM::CRDCaloBarCluster& m_cluster, std::vector<CRDEcalEDM::CRDCaloBarShower>& outshCol, const std::vector<CRDEcalEDM::CRDShowerCandidate>& m_EMcandidates, const std::vector<CRDEcalEDM::CRDShowerCandidate>& m_Trkcandidates );

  bool MergeToClosestCluster( CRDEcalEDM::CRDCaloBarCluster& iclus, std::vector<CRDEcalEDM::CRDCaloBarCluster>& clusvec ); 

  void CalculateInitialEseed( const std::vector<CRDEcalEDM::CRDCaloBar>& Seeds, const dd4hep::Position* pos, double* Eseed);

  double GetShowerProfile(const dd4hep::Position& p_bar, const dd4hep::Position& p_seed );

  static bool compE( CRDEcalEDM::CRDCaloBar& bar1, CRDEcalEDM::CRDCaloBar& bar2 ) { return bar1.getEnergy()>bar2.getEnergy(); }




  Settings settings;

private: 

};


#endif
