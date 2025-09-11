#ifndef FINALPID_SVC_H
#define FINALPID_SVC_H

#include "FinalPIDSvc/IFinalPIDSvc.h"
#include "DetInterface/IGeomSvc.h"
#include <GaudiKernel/Service.h>

#include <DDRec/DetectorData.h>
#include <DDRec/CellIDPositionConverter.h>
#include <DD4hep/Segmentations.h>

#include "FinalPIDSvc/WorkingPoint.h"
#include "FinalPIDSvc/MuonExtrapolator.h"

#include "TFile.h"
#include "TTree.h"
#include "TGraph2D.h"
#include "TH2F.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TMath.h"
#include <TMatrixD.h>
#include <TDecompSVD.h>
#include <xgboost/c_api.h>

/**
 * @class FinalPIDSvc
 * @brief PID implementation service
 * @author Geliang Liu (glliu@ihep.ac.cn)
*/

class FinalPIDSvc : public extends<Service, IFinalPIDSvc>
{
    public:

        FinalPIDSvc(const std::string& name, ISvcLocator* svc);
        virtual ~FinalPIDSvc();

        void SetCollections( const edm4hep::TrackerHitCollection* barrelhits, const edm4hep::TrackerHitCollection* endcaphits, const edm4hep::RecTofCollection* tofcol, const edm4hep::RecDqdxCollection* dqdxcol, const edm4hep::ReconstructedParticleCollection* PFO) override;

        void MatchMuonHitsToTracks() override;

        void SetWP_mu(int muWP) override;
        void SetWP_ele(int eleWP) override;
        void SetWP_pho(double phoWP) override;

        void Set_dR_max(double dR_max) override;

        // void Set_muHits_settings(bool use_dd, int n_muHits, bool forceMuHits) override;
        // void Set_muHits_parameters_mindR(double* sigma, double *lambda) override;
        // void Set_muHits_parameters_dd(double* sigma, double *lambda)  override;

        // void AddVar(int var) override;
        // void RemoveVar(int var) override;

        bool LoadPFO(const edm4hep::ReconstructedParticle pfo) override;
        
        // void ComputeChi2(edm4hep::ReconstructedParticle& pfo) override;

        void ApplyModel() override;

        // void ComputeChi2Total() override;

        int GetType() override;

        double GetChi2Total(int i_pdg) override; //total chi2
        // double GetChi2_NDF(int i_pdg) override; //total chi2 divided by NDF

        double GetChi2TPC(int i_pdg) override; //chi2 from TPC
        double GetChi2TOF(int i_pdg) override; //chi2 from TOF
        // double GetChi2Eecalp(int i_pdg) override; //chi2 from Eecal/p for electron
        // double GetChi2Eecal(int i_pdg) override; //chi2 from Eecal for muon
        // double GetChi2Ehcal(int i_pdg) override; //chi2 from Ehcal for lepton
        // double GetChi2Rhad(int i_pdg) override;//chi2 from Rhad for lepton
        // double GetChi2MindR(int i_pdg) override;//chi2 from min dR for muon
        // double GetChi2MindR1(int i_pdg) override;//chi2 from min dR for muon
        // double GetChi2MindR2(int i_pdg) override;//chi2 from min dR for muon
        // double GetChi2dd(int i_pdg) override;//chi2 from min dR for muon
        // double GetChi2dd1(int i_pdg) override;//chi2 from min dR for muon
        // double GetChi2dd2(int i_pdg) override;//chi2 from min dR for muon

        // double GetNDF(int i_pdg) override;
        float GetProb(int i_pdg) override;

        double GetP() override;
        double GetTheta() override;
        double GetPhi() override;

        double GetTof() override;
        double GetDndx() override;
        double GetE(bool isHCAL) override;
        double GetEp(bool isHCAL) override;
        double GetL(bool isHCAL) override;
        double GetR90(bool isHCAL) override;
        double GetWeta2(bool isHCAL) override;
        double GetWphi2(bool isHCAL) override;
        double GetTimeFirst(bool isHCAL) override;
        double GetTimeLast(bool isHCAL) override;

        double GetMindR(int i) override;
        double GetMindR_last() override;
        double Getdd(int i) override;
        double Getdd_last() override;

        int GetNhcal() override;
        int GetNmuon() override;

        int GetCharge() override;

        StatusCode initialize() override;
        StatusCode finalize() override;

    private:
        // Gaudi::Property<std::string> m_Model_Folder{this, "Model_Folder", ""};
        Gaudi::Property<std::string> m_input_eleID_WP{this, "input_eleID_WP", ""};
        Gaudi::Property<std::string> m_input_muID_WP{this, "input_muID_WP", ""};
        Gaudi::Property<std::string> m_input_lepID_model{this, "input_lepID_model", ""};
        Gaudi::Property<std::string> m_input_hadID_model{this, "input_hadID_model", ""};
        Gaudi::Property<std::string> m_input_phoID_model{this, "input_phoID_model", ""};
        // Gaudi::Property<std::string> m_WP_Folder{this, "WP_Folder", ""};
        // Gaudi::Property<std::string> m_LepID_Folder{this, "LepID_Folder", ""};
        // Gaudi::Property<std::string> m_HadID_Folder{this, "HadID_Folder", ""};
        // Gaudi::Property<std::string> m_PhoID_Folder{this, "PhoID_Folder", ""};
        Gaudi::Property<bool> m_readFromData{this, "readFromData", 0};
        Gaudi::Property<bool> m_computeVar{this, "computeVar", 1};

        const edm4hep::TrackerHitCollection* _barrelhits;
        const edm4hep::TrackerHitCollection* _endcaphits;
        const edm4hep::RecTofCollection* _tofcol;
        const edm4hep::RecDqdxCollection* _dqdxcol;
        const edm4hep::ReconstructedParticleCollection* _PFO;

        static constexpr int N_p_bins_lepton=16;
        static constexpr int N_p_bins_photon=14;
        static constexpr int N_theta_bins=8;
        const float p_bins_lepton[N_p_bins_lepton+1]={0,0.5,1,2,3,4,5,6,8,10,20,30,40,50,60,70,80};
        const float p_bins_photon[N_p_bins_photon+1]={1,2,3,4,5,6,8,10,20,30,40,50,60,70,80};
        const float theta_bins[N_theta_bins+1]={8,20,35,45,55,65,75,85,90};

        MuonExtrapolator *_muonExtrapolator;

        SmartIF<IGeomSvc> m_geosvc;
        dd4hep::Detector* m_dd4hep;

        dd4hep::DDSegmentation::BitFieldCoder* m_decoder_barrel;
        dd4hep::DDSegmentation::BitFieldCoder* m_decoder_endcap;
        
        int _muWP=WP::is90;
        int _eleWP=WP::is90;
        double _phoWP=0.1;
        double _dR_max=1e10+1;
        
        std::map<int,TrackState> _extrap_TS;
        std::map<int,std::vector<int>> _superlayer_PFO_to_MuonHits;
        std::map<int,std::vector<int>> _sector_PFO_to_MuonHits;
        std::map<int,std::vector<int>> _idx_PFO_to_MuonHits;
        std::map<int,std::vector<double>> _dR_PFO_to_MuonHits; //delta R
        std::map<int,std::vector<double>> _dd_PFO_to_MuonHits; //distance in the layer
        std::map<int,std::vector<double>> _mindR_PFO_to_MuonHits;
        std::map<int,std::vector<double>> _mindd_PFO_to_MuonHits; //Notice that this is not the minimum dd, but the dd corresponding to the minimum dR
        TGraph2D *_mu_dPhi_corr, *_mu_dEta_corr;

        TFile *_ele_Input, *_mu_Input;
        TH2F *_ele_WP[WP::num_of_WPs-2], *_mu_WP[WP::num_of_WPs-2];

        BoosterHandle _model_lepton[N_p_bins_lepton][N_theta_bins];
        BoosterHandle _model_hadron[N_p_bins_lepton][N_theta_bins];
        BoosterHandle _model_photon[N_p_bins_photon][N_theta_bins];
        DMatrixHandle dmatrix;

        double _chi2_TPC[5], _chi2_TOF[5], _chi2_Total[5];
        float _prob_lepton[3], _prob_hadron[3], _prob_photon[2];
        
        double _p, _theta, _theta_bkp, _phi;
        double _tof, _dndx;
        double _Eecalp, _Eecal, _Lecal, _R90ecal, _Weta2ecal, _Wphi2ecal;
        double _Ehcalp, _Ehcal, _Lhcal, _R90hcal, _Weta2hcal, _Wphi2hcal;
        double _Tecal_first, _Tecal_last, _Thcal_first, _Thcal_last;

        double _mindR[3], _dd[3], _mindR_last, _dd_last;
        int _Nhcal, _Nmuon;
        int _charge;

        void FillTPCTOF(bool doPID, const edm4hep::ReconstructedParticle pfo);

        double safe_interpolation(TGraph2D *graph, double p, double theta);
        double find_content(TH2F* h, double p, double theta);
        int get_bin_index(float x, const float *x_bins, const int N_bins);

        TVector3 solve_barrel_plane(TVector3 x0, TVector3 dx, double r);
        TVector3 solve_endcap_plane(TVector3 x0, TVector3 dx, double z);

        void computeShowerShapes(double E, TVector3 x, std::vector<double> E_hit, std::vector<double> Time_hit, std::vector<TVector3> x_hit, double& R90, double& Weta2, double& Wphi2, double& T_first, double& T_last);

        const std::map<int, int> PDGIDs = {
            {0, -11},
            {1, -13},
            {2, 211},
            {3, 321},
            {4, 2212},
        };

        const std::map<int, double> ParticleMass = {
            {11, 0.000511},
            {13, 0.105658},
            {211, 0.139570},
            {321, 0.493677},
            {2212, 0.938272},
            {22, 0.},
            //{130, 0.497611},
            {130, 0.} //Temp: set K_0 mass as 0, due to imperfect neutral hadron reconstruction in CyberPFA. 
        };

};

#endif

