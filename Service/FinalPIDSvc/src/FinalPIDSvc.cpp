#include "FinalPIDSvc.h"
#include "DataHelper/HelixClass.h"

#include "DD4hep/Detector.h"
#include <DD4hep/Objects.h>
#include <DDRec/CellIDPositionConverter.h>

#include <vector>
#include <math.h>
#include <stdlib.h>
#include <iostream>
#include <cassert>
#include <nlohmann/json.hpp>
using json = nlohmann::json;

DECLARE_COMPONENT(FinalPIDSvc)

using namespace std;
using namespace edm4hep;

FinalPIDSvc::FinalPIDSvc(const std::string& name, ISvcLocator* svc)
    : base_class(name, svc)
{
    _ele_Input=nullptr;
    _mu_Input=nullptr;
    _muonExtrapolator = nullptr;
}

FinalPIDSvc::~FinalPIDSvc()
{
    if (_ele_Input != nullptr) delete _ele_Input;
    if (_mu_Input != nullptr) delete _mu_Input;
    if (_muonExtrapolator != nullptr) delete _muonExtrapolator;
}

void FinalPIDSvc::SetCollections( const edm4hep::TrackerHitCollection* barrelhits, const edm4hep::TrackerHitCollection* endcaphits, const edm4hep::RecTofCollection* tofcol, const edm4hep::RecDqdxCollection* dqdxcol, const edm4hep::ReconstructedParticleCollection* PFO) 
{
    _barrelhits=barrelhits;
    _endcaphits=endcaphits;
    _tofcol=tofcol;
    _dqdxcol=dqdxcol;
    _PFO=PFO;
}

void FinalPIDSvc::MatchMuonHitsToTracks()
{
    _superlayer_PFO_to_MuonHits.clear();
    _sector_PFO_to_MuonHits.clear();
    _idx_PFO_to_MuonHits.clear();
    _dR_PFO_to_MuonHits.clear();
    _dd_PFO_to_MuonHits.clear();
    _mindR_PFO_to_MuonHits.clear();
    _mindd_PFO_to_MuonHits.clear();

    _extrap_TS.clear();

    for(const auto& pfo : *_PFO) {
        if ( isnan(pfo.getMomentum()[0]) || isnan(pfo.getMomentum()[1]) || isnan(pfo.getMomentum()[2]) || pfo.getEnergy() == 0) {
            continue;
        }
        if (pfo.tracks_size()==0) continue;

        TrackState newst;
        int idx;

        for (auto trk : pfo.getTracks()){
            idx=trk.getObjectID().index;
            newst=_muonExtrapolator->extrap_Simple(trk);
        }
        _extrap_TS[idx]=newst;
    }

    for (const auto &hit: *_barrelhits) {
        TVector3 x_hit=TVector3(hit.getPosition().x,hit.getPosition().y,hit.getPosition().z);
        int idx_mu=hit.getObjectID().index;

        double mindR=1e10;
        int idx_pfo=-1;
        TVector3 x0;
        TVector3 dx;

        for(const auto& pfo : *_PFO) {
            if ( isnan(pfo.getMomentum()[0]) || isnan(pfo.getMomentum()[1]) || isnan(pfo.getMomentum()[2]) || pfo.getEnergy() == 0) {
                continue;
            }
            if (pfo.tracks_size()==0) continue;

            TLorentzVector p_pfo=TLorentzVector(pfo.getMomentum()[0],pfo.getMomentum()[1],pfo.getMomentum()[2],pfo.getEnergy());
            
            double p=p_pfo.P();
            double theta=p_pfo.Theta()*TMath::RadToDeg();
            if (theta>90) theta=180-theta;

            int idx;
            for (auto trk : pfo.getTracks()){
                idx=trk.getObjectID().index;
            }
            TrackState newst=_extrap_TS[idx];

            HelixClass helix;
            helix.Initialize_Canonical( newst.phi, newst.D0, newst.Z0, newst.omega, newst.tanLambda, 3.0 );

            TVector3 p_trk(helix.getMomentum()[0],helix.getMomentum()[1],helix.getMomentum()[2]);
            //track extrapolation correction: notice that it changes momentum. While since it is not actually used, it doens't matter.
            p_trk.SetPtEtaPhi(p_trk.Pt(),p_trk.Eta()+safe_interpolation(_mu_dEta_corr,p,theta),p_trk.Phi()+safe_interpolation(_mu_dPhi_corr,p,theta));

            TVector3 x_ref=TVector3(newst.referencePoint.x,newst.referencePoint.y,newst.referencePoint.z);

            TVector3 x_rel=x_hit-x_ref;
            double dPhi=x_rel.DeltaPhi(p_trk);
            double dEta=x_rel.Eta()-p_trk.Eta();
            double dR=TMath::Sqrt(dPhi*dPhi+dEta*dEta);

            if (dR<mindR) {
                idx_pfo=idx;
                mindR=dR;
                x0=x_ref;
                dx=p_trk;
            }
        }

        unsigned long long cellid = hit.getCellID();
        _superlayer_PFO_to_MuonHits[idx_pfo].push_back(m_decoder_barrel->get(cellid, "Superlayer"));
        _sector_PFO_to_MuonHits[idx_pfo].push_back(m_decoder_barrel->get(cellid, "Fe"));

        _idx_PFO_to_MuonHits[idx_pfo].push_back(idx_mu);
        _dR_PFO_to_MuonHits[idx_pfo].push_back(mindR);

        double r_hit=x_hit.Pt();
        TVector3 x_extrap=solve_barrel_plane(x0,dx,r_hit);
        double dd=TMath::Sqrt(r_hit*r_hit*x_extrap.DeltaPhi(x_hit)*x_extrap.DeltaPhi(x_hit)+(x_extrap[2]-x_hit[2])*(x_extrap[2]-x_hit[2])); //For barrel, the distance at this r plane is computed as sqrt[(rdφ)^2+(dz)^2]
        _dd_PFO_to_MuonHits[idx_pfo].push_back(dd);
    }

    for (const auto &hit: *_endcaphits) {
        
        TVector3 x_hit=TVector3(hit.getPosition().x,hit.getPosition().y,hit.getPosition().z);
        int idx_mu=hit.getObjectID().index;

        double mindR=1e10;
        int idx_pfo=-1;
        TVector3 x0;
        TVector3 dx;

        for(const auto& pfo : *_PFO) {
            if ( isnan(pfo.getMomentum()[0]) || isnan(pfo.getMomentum()[1]) || isnan(pfo.getMomentum()[2]) || pfo.getEnergy() == 0) {
                continue;
            }
            if (pfo.tracks_size()==0) continue;

            TLorentzVector p_pfo=TLorentzVector(pfo.getMomentum()[0],pfo.getMomentum()[1],pfo.getMomentum()[2],pfo.getEnergy());
            
            double p=p_pfo.P();
            double theta=p_pfo.Theta()*TMath::RadToDeg();
            if (theta>90) theta=180-theta;
            int charge=pfo.getCharge();

            TrackState newst;
            int idx;

            for (auto trk : pfo.getTracks()){
                idx=trk.getObjectID().index;
                newst=_muonExtrapolator->extrap_Simple(trk);
            }

            HelixClass helix;
            helix.Initialize_Canonical( newst.phi, newst.D0, newst.Z0, newst.omega, newst.tanLambda, 3.0 );

            TVector3 p_trk(helix.getMomentum()[0],helix.getMomentum()[1],helix.getMomentum()[2]);
            //track extrapolation correction: notice that it changes momentum. While since it is not actually used, it doens't matter.
            p_trk.SetPtEtaPhi(p_trk.Pt(),p_trk.Eta()+safe_interpolation(_mu_dEta_corr,p,theta),p_trk.Phi()-charge*safe_interpolation(_mu_dPhi_corr,p,theta));

            TVector3 x_ref=TVector3(newst.referencePoint.x,newst.referencePoint.y,newst.referencePoint.z);

            TVector3 x_rel=x_hit-x_ref;
            double dPhi=x_rel.DeltaPhi(p_trk);
            double dEta=x_rel.Eta()-p_trk.Eta();
            double dR=TMath::Sqrt(dPhi*dPhi+dEta*dEta);

            if (dR<mindR) {
                idx_pfo=idx;
                mindR=dR;
                x0=x_ref;
                dx=p_trk;
            }
        }

        unsigned long long cellid = hit.getCellID();
        _superlayer_PFO_to_MuonHits[idx_pfo].push_back(m_decoder_endcap->get(cellid, "Superlayer"));
        _sector_PFO_to_MuonHits[idx_pfo].push_back(m_decoder_endcap->get(cellid, "Endcap"));

        _idx_PFO_to_MuonHits[idx_pfo].push_back(idx_mu);
        _dR_PFO_to_MuonHits[idx_pfo].push_back(mindR);

        double z_hit=x_hit[2];
        TVector3 x_extrap=solve_endcap_plane(x0,dx,z_hit);
        double dd=TMath::Sqrt((x_extrap[0]-x_hit[0])*(x_extrap[0]-x_hit[0])+(x_extrap[1]-x_hit[1])*(x_extrap[1]-x_hit[1]));//For endcap, the distance at this z plane is computed as sqrt[(dx)^2+(dy)^2]
        _dd_PFO_to_MuonHits[idx_pfo].push_back(dd);
    }
    
    for (auto it = _dR_PFO_to_MuonHits.begin(); it != _dR_PFO_to_MuonHits.end(); ++it) {

        int i_pfo=it->first;
        std::vector<double> vec_dR=it->second;
        std::vector<double> vec_dd=_dd_PFO_to_MuonHits.at(i_pfo);
        std::vector<int> vec_layer=_superlayer_PFO_to_MuonHits.at(i_pfo);
        std::vector<int> vec_sector=_sector_PFO_to_MuonHits.at(i_pfo);

        double mindR=1e10;
        double mindd=1e10;
        int layer1=-1;
        int sector1=-1;
        if (!vec_dR.empty()) {
            for (size_t i=0;i<vec_dR.size();i++) {
                if (vec_dR.at(i)<_dR_max && vec_dR.at(i)<mindR) {
                    mindR=vec_dR.at(i);
                    mindd=vec_dd.at(i);
                    layer1=vec_layer.at(i);
                    sector1=vec_sector.at(i);
                }
            }
        }
        _mindR_PFO_to_MuonHits[i_pfo].push_back(mindR);
        _mindd_PFO_to_MuonHits[i_pfo].push_back(mindd);

        mindR=1e10;
        mindd=1e10;
        int layer2=-1;
        int sector2=-1;
        if (!vec_dR.empty()) {
            for (size_t i=0;i<vec_dR.size();i++) {
                if ((vec_layer.at(i)!=layer1 || vec_sector.at(i)!=sector1) && vec_dR.at(i)<_dR_max && vec_dR.at(i)<mindR) {
                    mindR=vec_dR.at(i);
                    mindd=vec_dd.at(i);
                    layer2=vec_layer.at(i);
                    sector2=vec_sector.at(i);
                }
            }
        }
        _mindR_PFO_to_MuonHits[i_pfo].push_back(mindR);
        _mindd_PFO_to_MuonHits[i_pfo].push_back(mindd);

        mindR=1e10;
        mindd=1e10;
        if (!vec_dR.empty()) {
            for (size_t i=0;i<vec_dR.size();i++) {
                if ((vec_layer.at(i)!=layer1 || vec_sector.at(i)!=sector1) && (vec_layer.at(i)!=layer2 || vec_sector.at(i)!=sector2) && vec_dR.at(i)<_dR_max && vec_dR.at(i)<mindR) {
                    mindR=vec_dR.at(i);
                    mindd=vec_dd.at(i);
                }
            }
        }
        _mindR_PFO_to_MuonHits[i_pfo].push_back(mindR);
        _mindd_PFO_to_MuonHits[i_pfo].push_back(mindd);

        mindR=1e10;
        mindd=1e10;
        if (!vec_dR.empty()) {
            for (size_t i=0;i<vec_dR.size();i++) {
                if ((vec_layer.at(i)==5 || vec_layer.at(i)==6) && vec_dR.at(i)<_dR_max && vec_dR.at(i)<mindR) {
                    mindR=vec_dR.at(i);
                    mindd=vec_dd.at(i);
                }
            }
        }
        _mindR_PFO_to_MuonHits[i_pfo].push_back(mindR);
        _mindd_PFO_to_MuonHits[i_pfo].push_back(mindd);
    }
}

void FinalPIDSvc::SetWP_mu(int muWP) 
{
    _muWP=muWP;
}

void FinalPIDSvc::SetWP_ele(int eleWP)
{
    _eleWP=eleWP;
}

void FinalPIDSvc::SetWP_pho(double phoWP)
{
    _phoWP=phoWP;
}

void FinalPIDSvc::Set_dR_max(double dR_max)
{
    _dR_max=dR_max;
}

bool FinalPIDSvc::LoadPFO(const edm4hep::ReconstructedParticle pfo)
{
    if ( isnan(pfo.getMomentum()[0]) || isnan(pfo.getMomentum()[1]) || isnan(pfo.getMomentum()[2]) || pfo.getEnergy() == 0) {
        return false;
    }

    if (m_computeVar) {

        TLorentzVector p_pfo=TLorentzVector(pfo.getMomentum()[0],pfo.getMomentum()[1],pfo.getMomentum()[2],pfo.getEnergy());

        _charge=pfo.getCharge();

        _p=p_pfo.P();
        if (_p==0) return false;

        _theta=p_pfo.Theta()*TMath::RadToDeg();
        _phi=p_pfo.Phi();
        if (_theta>90) _theta_bkp=180-_theta;
        else _theta_bkp=_theta;

        _Eecal=0;
        _Ehcal=0;
        _Nhcal=0;
        TVector3 x3ecal=TVector3(0,0,0);
        TVector3 x3hcal=TVector3(0,0,0);

        std::vector<TVector3> x3ecal_hit;
        std::vector<TVector3> x3hcal_hit;
        std::vector<double> Eecal_hit;
        std::vector<double> Ehcal_hit;
        std::vector<double> Eecal_time;
        std::vector<double> Ehcal_time;

        for (unsigned i=0;i<pfo.clusters_size();i++) {
            edm4hep::Cluster cluster=pfo.getClusters(i);
            double x=cluster.getPosition().x;
            double y=cluster.getPosition().y;
            double z=cluster.getPosition().z;
            double e=cluster.getEnergy();
            if (isnan(x) || isnan(y) || isnan(z) || isnan(e)) continue;

            // positions.emplace_back(x, y, z);
            // energies.push_back(e);

            if (TMath::Sqrt(x * x + y * y) < 2130 && fabs(z) < 3230) {
                _Eecal += e;
                x3ecal += e*TVector3(x,y,z);
                // info() << "ECAL cluster size:" << cluster.hits_size() << endmsg;
                for (unsigned j=0;j<cluster.hits_size();j++) {
                    double x_hit=cluster.getHits(j).getPosition().x;
                    double y_hit=cluster.getHits(j).getPosition().y;
                    double z_hit=cluster.getHits(j).getPosition().z;
                    double e_hit=cluster.getHits(j).getEnergy();
                    x3ecal_hit.emplace_back(x_hit, y_hit, z_hit);
                    Eecal_hit.push_back(e_hit);
                    Eecal_time.push_back(cluster.getHits(j).getTime());
                }
            }
            else {
                _Ehcal += e;
                x3hcal += e*TVector3(x,y,z);
                _Nhcal += 1;
                // info() << "HCAL cluster size:" << cluster.hits_size() << endmsg;
                for (unsigned j=0;j<cluster.hits_size();j++) {
                    double x_hit=cluster.getHits(j).getPosition().x;
                    double y_hit=cluster.getHits(j).getPosition().y;
                    double z_hit=cluster.getHits(j).getPosition().z;
                    double e_hit=cluster.getHits(j).getEnergy();
                    x3hcal_hit.emplace_back(x_hit, y_hit, z_hit);
                    Ehcal_hit.push_back(e_hit);
                    Ehcal_time.push_back(cluster.getHits(j).getTime());
                }
            }
        }

        _Eecalp=_Eecal/_p;
        _Ehcalp=_Ehcal/_p;

        if (_Eecal>0) x3ecal=x3ecal*(1./_Eecal);
        _Lecal=x3ecal.Mag();
        if (_Ehcal>0) x3hcal=x3hcal*(1./_Ehcal);
        _Lhcal=x3hcal.Mag();

        computeShowerShapes(std::accumulate(Eecal_hit.begin(), Eecal_hit.end(), 0.0), x3ecal, Eecal_hit, Eecal_time, x3ecal_hit, _R90ecal, _Weta2ecal, _Wphi2ecal, _Tecal_first, _Tecal_last);
        computeShowerShapes(std::accumulate(Ehcal_hit.begin(), Ehcal_hit.end(), 0.0), x3hcal, Ehcal_hit, Ehcal_time, x3hcal_hit, _R90hcal, _Weta2hcal, _Wphi2hcal, _Thcal_first, _Thcal_last);

        int idx;
        for (auto trk : pfo.getTracks()){
            idx=trk.getObjectID().index;
        }

        if (_mindR_PFO_to_MuonHits.find(idx) != _mindR_PFO_to_MuonHits.end()) {
            for (size_t i=0;i<3;i++) {
                _mindR[i]=_mindR_PFO_to_MuonHits[idx].at(i);
                _dd[i]=_mindd_PFO_to_MuonHits[idx].at(i);
            }
            _mindR_last=_mindR_PFO_to_MuonHits[idx].at(3);
            _dd_last=_mindd_PFO_to_MuonHits[idx].at(3);
        }
        else {
            for (size_t i=0;i<3;i++) {
                _mindR[i]=1e10;
                _dd[i]=1e10;
            }
            _mindR_last=1e10;
            _dd_last=1e10;
        }

        if (_dR_PFO_to_MuonHits.find(idx) != _dR_PFO_to_MuonHits.end()) {
            _Nmuon=_dR_PFO_to_MuonHits.at(idx).size();
        }
        else {
            _Nmuon=0;
        }
    
        bool doPID=true;
        if ( _tofcol->size() == 0 && _dqdxcol->size() == 0 ) doPID=false;
        FillTPCTOF(doPID,pfo);
    }

    if (m_readFromData) {
        if (pfo.getCharge()==0) {
            if (pfo.particleIDs_size()!=2) {
                error()<<"Neutral particles should have 2 ParticleIDs instead of "<<pfo.particleIDs_size()<<endmsg;
                return false;
            }
            for (unsigned int i_fl=0;i_fl<2;i_fl++) {
                auto pid=pfo.getParticleIDs(i_fl);
                _prob_photon[i_fl]=pid.getLikelihood();
            }
        }
        else {
            if (pfo.particleIDs_size()!=5) {
                error()<<"Charged particles should have 5 ParticleIDs instead of "<<pfo.particleIDs_size()<<endmsg;
                return false;
            }
            for (unsigned int i_fl=0;i_fl<2;i_fl++) {
                auto pid=pfo.getParticleIDs(i_fl);
                _prob_lepton[i_fl]=pid.getLikelihood();
            }
            _prob_lepton[2]=0;
            for (unsigned int i_fl=2;i_fl<5;i_fl++) {
                auto pid=pfo.getParticleIDs(i_fl);
                _prob_lepton[2]+=pid.getLikelihood();
                _prob_hadron[i_fl-2]=pid.getLikelihood();
            }
            if (_prob_lepton[2]>0) {
                for (unsigned int i_fl=2;i_fl<5;i_fl++) {
                    _prob_hadron[i_fl-2]/=_prob_lepton[2];
                }
            }
        }
    }
    else ApplyModel();

    return true;
}

int FinalPIDSvc::GetType()
{
    if (_charge==0) {
        // if (_Eecal>0 && _Ehcal/_Eecal<_phoWP)
        if (_prob_photon[0]>=0.5)
            return 22;
        else return 130;
    }
    else {
        if (_muWP!=WP::noLep) {
            if (_muWP==WP::Best) {
                float* maxElement = std::max_element(_prob_lepton, _prob_lepton + 3);
                int index = static_cast<int>(maxElement - _prob_lepton);
                if (index==1) return 13;
            }
            else {
                if (_prob_lepton[1]>=find_content(_mu_WP[_muWP],_p,_theta_bkp)) return 13;
            }
        }
        if (_eleWP!=WP::noLep) {
            if (_eleWP==WP::Best) {
                float* maxElement = std::max_element(_prob_lepton, _prob_lepton + 3);
                int index = static_cast<int>(maxElement - _prob_lepton);
                if (index==0) return 11;
            }
            else {
                if (_prob_lepton[0]>=find_content(_ele_WP[_eleWP],_p,_theta_bkp)) return 11;
            }
        }
        float* maxElement = std::max_element(_prob_hadron, _prob_hadron + 3);
        int index = static_cast<int>(maxElement - _prob_hadron);
        return PDGIDs.at(index+2);
    }
}

double FinalPIDSvc::GetChi2Total(int i_pdg) { return _chi2_Total[i_pdg]; }

double FinalPIDSvc::GetChi2TPC(int i_pdg) { return _chi2_TPC[i_pdg]; }

double FinalPIDSvc::GetChi2TOF(int i_pdg) { return _chi2_TOF[i_pdg]; }

float FinalPIDSvc::GetProb(int i_pdg) {
    if (_charge==0) return _prob_photon[i_pdg]; 
    else {
        if (i_pdg<2) return _prob_lepton[i_pdg];
        else return _prob_lepton[2]*_prob_hadron[i_pdg-2];
    }
}

double FinalPIDSvc::GetP() { return _p; }

double FinalPIDSvc::GetTheta() { return _theta; }

double FinalPIDSvc::GetPhi() { return _phi; }

double FinalPIDSvc::GetTof() { return _tof; }

double FinalPIDSvc::GetDndx() { return _dndx; }

double FinalPIDSvc::GetE(bool isHCAL) { return isHCAL?_Ehcal:_Eecal; }

double FinalPIDSvc::GetEp(bool isHCAL) { return isHCAL?_Ehcalp:_Eecalp; }

double FinalPIDSvc::GetL(bool isHCAL) { return isHCAL?_Lhcal:_Lecal; }

double FinalPIDSvc::GetR90(bool isHCAL) { return isHCAL?_R90hcal:_R90ecal; }

double FinalPIDSvc::GetWeta2(bool isHCAL) { return isHCAL?_Weta2hcal:_Weta2ecal; }

double FinalPIDSvc::GetWphi2(bool isHCAL) { return isHCAL?_Wphi2hcal:_Wphi2ecal; }

double FinalPIDSvc::GetTimeFirst(bool isHCAL) { return isHCAL?_Thcal_first:_Tecal_first; }

double FinalPIDSvc::GetTimeLast(bool isHCAL) { return isHCAL?_Thcal_last:_Tecal_last; }

double FinalPIDSvc::GetMindR(int i) { return _mindR[i]; }

double FinalPIDSvc::Getdd(int i) { return _dd[i]; }

double FinalPIDSvc::GetMindR_last() { return _mindR_last; }

double FinalPIDSvc::Getdd_last() { return _dd_last; }

int FinalPIDSvc::GetNhcal() {return _Nhcal; }

int FinalPIDSvc::GetNmuon() {return _Nmuon; }

int FinalPIDSvc::GetCharge() { return _charge; }

StatusCode FinalPIDSvc::initialize()
{
    if (!m_readFromData && !m_computeVar) {
        error()<<"You cannot compute PID without computing variables !!!!"<<endmsg;
        return StatusCode::FAILURE;
    }

    _ele_Input=new TFile(m_input_eleID_WP.value().c_str(),"read");

    for (int i=0;i<WP::num_of_WPs-2;i++) {
        _ele_WP[i]=(TH2F *)_ele_Input->Get(WP::WPs.at(i));
    }

    _mu_Input=new TFile(m_input_muID_WP.value().c_str(),"read");

    _mu_dPhi_corr=(TGraph2D *)_mu_Input->Get("dPhi corr_new");
    _mu_dEta_corr=(TGraph2D *)_mu_Input->Get("dEta corr_new");
    for (int i=0;i<WP::num_of_WPs-2;i++) {
        _mu_WP[i]=(TH2F *)_mu_Input->Get(WP::WPs.at(i));
    }

    if (m_computeVar) {
        double BfieldRBound = (3535.0 + 4235.0)/2.0;
        double BfieldZBound = 4075.0;

        _muonExtrapolator = new MuonExtrapolator( BfieldRBound, BfieldZBound );

        m_geosvc = service<IGeomSvc>("GeomSvc");
        if (!m_geosvc) {
            error()<<"Failed to get the GeomSvc." << endmsg;
            return StatusCode::FAILURE;
        }
        m_dd4hep = m_geosvc->lcdd();
        if ( !m_dd4hep ) {
            error()<<"Failed to get the lcdd()." << endmsg;
            return StatusCode::FAILURE;
        }
        m_decoder_barrel = m_geosvc->getDecoder("MuonBarrelCollection");
        m_decoder_endcap = m_geosvc->getDecoder("MuonEndcapCollection");
        if(!m_decoder_barrel || !m_decoder_endcap){
            error()<<"Failed to get the decoder."<< endmsg;
            return StatusCode::FAILURE;
        }
    }

    //lepton ID models
    if (!m_readFromData) {
        std::ifstream in_lepton(m_input_lepID_model.value().c_str());
        if (!in_lepton) {
            std::cerr << "Could not open big JSON file\n";
            return StatusCode::FAILURE;
        }
        std::stringstream buffer_lepton;
        buffer_lepton << in_lepton.rdbuf();
        std::string cont_lepton = buffer_lepton.str();
        json j_lepton = json::parse(cont_lepton);
        
        for (int ip = 0; ip < N_p_bins_lepton; ++ip) {
            for (int itheta = 0; itheta < N_theta_bins; ++itheta) {
                if (XGBoosterCreate(nullptr, 0, &_model_lepton[ip][itheta]) != 0) {
                    std::cerr << "Failed to create model[" << ip << "][" << itheta << "]\n";
                    return StatusCode::FAILURE;
                }
                XGBoosterSetParam(_model_lepton[ip][itheta], "nthread", "4");
                std::stringstream ss;
                ss.str("");
                if (ip<2) ss << Form("model_p%.1f_theta%.0f",p_bins_lepton[ip],theta_bins[itheta]);
                else ss << Form("model_p%.0f_theta%.0f",p_bins_lepton[ip],theta_bins[itheta]);
                std::string filename = ss.str();

                if (!j_lepton.contains(filename)) {
                    std::cerr << "Lepton model not found: " << filename << "\n";
                    return StatusCode::FAILURE;
                }
                std::string model_str = j_lepton[filename].dump();
                if (XGBoosterLoadModelFromBuffer(_model_lepton[ip][itheta], model_str.data(), model_str.size()) != 0) {
                    std::cerr << "Failed to load lepton model: " << filename << "\n";
                    return StatusCode::FAILURE;
                }
            }
        }
    }

    //hadron ID models
    if (!m_readFromData) {
        std::ifstream in_hadron(m_input_hadID_model.value().c_str());
        if (!in_hadron) {
            std::cerr << "Could not open big JSON file\n";
            return StatusCode::FAILURE;
        }
        std::stringstream buffer_hadron;
        buffer_hadron << in_hadron.rdbuf();
        std::string cont_hadron = buffer_hadron.str();
        json j_hadron = json::parse(cont_hadron);

        for (int ip = 0; ip < N_p_bins_lepton; ++ip) {
            for (int itheta = 0; itheta < N_theta_bins; ++itheta) {
                if (XGBoosterCreate(nullptr, 0, &_model_hadron[ip][itheta]) != 0) {
                    std::cerr << "Failed to create model[" << ip << "][" << itheta << "]\n";
                    return StatusCode::FAILURE;
                }
                XGBoosterSetParam(_model_hadron[ip][itheta], "nthread", "4");
                std::stringstream ss;
                ss.str("");
                if (ip<2) ss << Form("model_p%.1f_theta%.0f",p_bins_lepton[ip],theta_bins[itheta]);
                else ss << Form("model_p%.0f_theta%.0f",p_bins_lepton[ip],theta_bins[itheta]);
                std::string filename = ss.str();

                if (!j_hadron.contains(filename)) {
                    std::cerr << "Hadron model not found: " << filename << "\n";
                    return StatusCode::FAILURE;
                }
                std::string model_str = j_hadron[filename].dump();
                if (XGBoosterLoadModelFromBuffer(_model_hadron[ip][itheta], model_str.data(), model_str.size()) != 0) {
                    std::cerr << "Failed to load hadron model: " << filename << "\n";
                    return StatusCode::FAILURE;
                }
            }
        }
    }

    //photon ID models
    if (!m_readFromData) {
        std::ifstream in_photon(m_input_phoID_model.value().c_str());
        if (!in_photon) {
            std::cerr << "Could not open big JSON file\n";
            return StatusCode::FAILURE;
        }
        std::stringstream buffer_photon;
        buffer_photon << in_photon.rdbuf();
        std::string cont_photon = buffer_photon.str();
        json j_photon = json::parse(cont_photon);

        for (int ip = 0; ip < N_p_bins_photon; ++ip) {
            for (int itheta = 0; itheta < N_theta_bins; ++itheta) {
                if (XGBoosterCreate(nullptr, 0, &_model_photon[ip][itheta]) != 0) {
                    std::cerr << "Failed to create model[" << ip << "][" << itheta << "]\n";
                    return StatusCode::FAILURE;
                }
                XGBoosterSetParam(_model_photon[ip][itheta], "nthread", "4");
                std::stringstream ss;
                ss.str("");
                ss << Form("model_p%.1f_theta%.0f",p_bins_photon[ip],theta_bins[itheta]);
                std::string filename = ss.str();

                if (!j_photon.contains(filename)) {
                    std::cerr << "Hadron model not found: " << filename << "\n";
                    return StatusCode::FAILURE;
                }
                std::string model_str = j_photon[filename].dump();
                if (XGBoosterLoadModelFromBuffer(_model_photon[ip][itheta], model_str.data(), model_str.size()) != 0) {
                    std::cerr << "Failed to load hadron model: " << filename << "\n";
                    return StatusCode::FAILURE;
                }
            }
        }
    }

    return StatusCode::SUCCESS;
}

StatusCode FinalPIDSvc::finalize()
{
    if (_ele_Input != nullptr) delete _ele_Input;
    if (_mu_Input != nullptr) delete _mu_Input;
    if (_muonExtrapolator != nullptr) delete _muonExtrapolator;
    if (!m_readFromData) {
        for (int ip = 0; ip < N_p_bins_lepton; ++ip) {
            for (int itheta = 0; itheta < N_theta_bins; ++itheta) {
                if (_model_lepton[ip][itheta] != nullptr) {
                    XGBoosterFree(_model_lepton[ip][itheta]);
                    _model_lepton[ip][itheta] = nullptr;
                }
            }
        }
        for (int ip = 0; ip < N_p_bins_lepton; ++ip) {
            for (int itheta = 0; itheta < N_theta_bins; ++itheta) {
                if (_model_hadron[ip][itheta] != nullptr) {
                    XGBoosterFree(_model_hadron[ip][itheta]);
                    _model_hadron[ip][itheta] = nullptr;
                }
            }
        }
        for (int ip = 0; ip < N_p_bins_photon; ++ip) {
            for (int itheta = 0; itheta < N_theta_bins; ++itheta) {
                if (_model_photon[ip][itheta] != nullptr) {
                    XGBoosterFree(_model_photon[ip][itheta]);
                    _model_photon[ip][itheta] = nullptr;
                }
            }
        }
    }
    return StatusCode::SUCCESS;
}

void FinalPIDSvc::FillTPCTOF(bool doPID, const edm4hep::ReconstructedParticle pfo)
{
    _tof=-1;
    _dndx=-1;
    for (int i_pdg = 0; i_pdg < 5; i_pdg++){
        _chi2_TPC[i_pdg]=_chi2_TOF[i_pdg]=-1;
    }
    if (doPID) {
        for (auto trk : pfo.getTracks()){

            for (auto dqdx : *_dqdxcol){
                if (dqdx.getTrack() == trk){
                    _dndx=dqdx.getDQdx().value;
                    for (int i_pdg = 0; i_pdg < 5; i_pdg++){
                        _chi2_TPC[i_pdg] = dqdx.getHypotheses(i_pdg).chi2;
                    }
                }
            }

            for (auto tof : *_tofcol){
                if (tof.getTrack() == trk){
                    _tof = tof.getTime();
                    std::array<float, 5> tofexpts = tof.getTimeExp();
                    double tofexpterr = tof.getSigma();
                    for (int i_pdg = 0; i_pdg < 5; i_pdg++){
                        _chi2_TOF[i_pdg] = std::pow( (tofexpts[i_pdg] - _tof) / tofexpterr, 2);
                    }
                }
            }
        }
    }
    for (int i_pdg=0;i_pdg<5;i_pdg++) {
        _chi2_Total[i_pdg]=0;
        if (_chi2_TPC[i_pdg]>=0) {
            _chi2_Total[i_pdg]+=_chi2_TPC[i_pdg];
        }
        if (_chi2_TOF[i_pdg]>=0) {
            _chi2_Total[i_pdg]+=_chi2_TOF[i_pdg];
        }
    }
}

void FinalPIDSvc::ApplyModel()
{
    if (_charge!=0) {
        int ip=get_bin_index(_p,p_bins_lepton,N_p_bins_lepton);
        int itheta=get_bin_index(_theta_bkp,theta_bins,N_theta_bins);

        if (p_bins_lepton[ip]>25 && theta_bins[itheta]>25) {
            const int ncol = 15;
            const float data[ncol] = {(float)_tof, (float)_dndx, (float)_Eecalp, (float)_Ehcalp, (float)_Lecal, (float)_Lhcal, (float)_R90ecal, (float)_R90hcal, (float)_Weta2ecal, (float)_Weta2hcal, (float)_Wphi2ecal, (float)_Wphi2hcal, (float)_mindR[0], (float)_Nhcal, (float)_Nmuon};

            XGDMatrixCreateFromMat(data, 1, ncol, std::numeric_limits<float>::quiet_NaN(), &dmatrix);

            bst_ulong out_len;
            const float* out_result = nullptr;
            XGBoosterPredict(_model_lepton[ip][itheta], dmatrix, 0, 0, 0, &out_len, &out_result);
            for (int i=0;i<3;i++) _prob_lepton[i]=out_result[i];
        }
        else {
            const int ncol = 18;
            const float data[ncol] = {(float)_tof, (float)_dndx, (float)_Eecalp, (float)_Ehcalp, (float)_Lecal, (float)_Lhcal, (float)_R90ecal, (float)_R90hcal, (float)_Weta2ecal, (float)_Weta2hcal, (float)_Wphi2ecal, (float)_Wphi2hcal, (float)_mindR[0], (float)_mindR[1], (float)_mindR[2], (float)_mindR_last, (float)_Nhcal, (float)_Nmuon};

            XGDMatrixCreateFromMat(data, 1, ncol, std::numeric_limits<float>::quiet_NaN(), &dmatrix);

            bst_ulong out_len;
            const float* out_result = nullptr;
            XGBoosterPredict(_model_lepton[ip][itheta], dmatrix, 0, 0, 0, &out_len, &out_result);
            for (int i=0;i<3;i++) _prob_lepton[i]=out_result[i];
        }

        const int ncol = 7;
        const float data[ncol] = {(float)_Eecalp, (float)_chi2_TOF[2], (float)_chi2_TPC[2], (float)_chi2_TOF[3], (float)_chi2_TPC[3], (float)_chi2_TOF[4], (float)_chi2_TPC[4]};

        XGDMatrixCreateFromMat(data, 1, ncol, std::numeric_limits<float>::quiet_NaN(), &dmatrix);

        bst_ulong out_len;
        const float* out_result = nullptr;
        XGBoosterPredict(_model_hadron[ip][itheta], dmatrix, 0, 0, 0, &out_len, &out_result);
        for (int i=0;i<3;i++) _prob_hadron[i]=out_result[i];
    }
    else {
        int ip=get_bin_index(_p,p_bins_photon,N_p_bins_photon);
        int itheta=get_bin_index(_theta_bkp,theta_bins,N_theta_bins);

        const int ncol = 13;
        const float data[ncol] = {(float)_tof, (float)_dndx, (float)_Eecalp, (float)_Ehcalp, (float)_Lecal, (float)_Lhcal, (float)_R90ecal, (float)_R90hcal, (float)_Weta2ecal, (float)_Weta2hcal, (float)_Wphi2ecal, (float)_Wphi2hcal, (float)_Nhcal};

        XGDMatrixCreateFromMat(data, 1, ncol, std::numeric_limits<float>::quiet_NaN(), &dmatrix);

        bst_ulong out_len;
        const float* out_result = nullptr;
        XGBoosterPredict(_model_photon[ip][itheta], dmatrix, 0, 0, 0, &out_len, &out_result);
        _prob_photon[0]=1.-out_result[0];
        _prob_photon[1]=out_result[0];
    }
}

double FinalPIDSvc::safe_interpolation(TGraph2D *graph, double p, double theta)
{
    assert(graph && "Error: graph is a nullptr!");

    double z = graph->Interpolate(p, theta);

    if (z != 0) {
        return z;
    } else {
        int n_points = graph->GetN();
        double min_dist = DBL_MAX;
        double best_z = 0.0;

        double p_i, theta_i, z_i;

        for (int i = 0; i < n_points; ++i) {
            graph->GetPoint(i, p_i, theta_i, z_i);
            // if (p_i<4) continue;
            double dist = std::sqrt((p - p_i) * (p - p_i) + (theta - theta_i) * (theta - theta_i));
            if (dist < min_dist) {
                min_dist = dist;
                best_z = z_i;
            }
        }

        return best_z;
    }
}

double FinalPIDSvc::find_content(TH2F* h, double p, double theta) {
    TAxis* xAxis = h->GetXaxis();
    TAxis* yAxis = h->GetYaxis();

    double dx = xAxis->GetBinWidth(1);
    double dy = yAxis->GetBinWidth(1);

    double p_min = xAxis->GetXmin();
    double p_max = xAxis->GetXmax();
    double theta_min = yAxis->GetXmin();
    double theta_max = yAxis->GetXmax();

    double epsilon_x = 1e-6 * dx;
    double epsilon_y = 1e-6 * dy;

    double p_clamped = std::min(std::max(p, p_min), p_max - epsilon_x);
    double theta_clamped = std::min(std::max(theta, theta_min), theta_max - epsilon_y);

    int globalBin = h->FindBin(p_clamped, theta_clamped);
    return h->GetBinContent(globalBin);
}

int FinalPIDSvc::get_bin_index(float x, const float *x_bins, const int N_bins) 
{
    if (x <= x_bins[0]) {
        return 0;
    }

    if (x >= x_bins[N_bins]) {
        return N_bins - 1;
    }

    int left = 0, right = N_bins - 1;
    while (left <= right) {
        int mid = left + (right - left) / 2;

        if (x >= x_bins[mid] && x < x_bins[mid + 1]) {
            return mid;
        } else if (x < x_bins[mid]) {
            right = mid - 1;
        } else {
            left = mid + 1;
        }
    }

    return -1;
}

TVector3 FinalPIDSvc::solve_barrel_plane(TVector3 x0, TVector3 dx, double r)
{
    double A=dx[0]*dx[0]+dx[1]*dx[1];
    double B=2*(x0[0]*dx[0]+x0[1]*dx[1]);
    double C=x0[0]*x0[0]+x0[1]*x0[1]-r*r;
    if (A!=0 && B*B-4*A*C>=0) {
        double s=(TMath::Sqrt(B*B-4*A*C)-B)/2/A;
        return x0+s*dx;
    }
    else {
        cout << "No solution in barrel!" <<endl;
        return x0;
    }
}

TVector3 FinalPIDSvc::solve_endcap_plane(TVector3 x0, TVector3 dx, double z)
{
    if (dx[2]!=0) {
        double s=(z-x0[2])/dx[2];
        return x0+s*dx;
    }
    else {
        cout << "No solution in endcap!" <<endl;
        return x0;
    }
}

void FinalPIDSvc::computeShowerShapes(double E, TVector3 x, std::vector<double> E_hit, std::vector<double> Time_hit, std::vector<TVector3> x_hit, double& R90, double& Weta2, double& Wphi2, double& T_first, double& T_last)
{
    TMatrixD cov(3, 3);
    cov.Zero();
    R90=Weta2=Wphi2=0;

    if (!Time_hit.empty()) {
        T_first = *std::min_element(Time_hit.begin(), Time_hit.end());
        T_last  = *std::max_element(Time_hit.begin(), Time_hit.end());
    }
    else {
        T_first=T_last=-1;
    }

    if (E>0 && E_hit.size()>2) {
        double Eta_avg=0;
        double Phi_avg=0;
        for (size_t i = 0; i < x_hit.size(); ++i) {
            TVector3 r = x_hit[i] - x;
            for (int a = 0; a < 3; ++a) {
                for (int b = 0; b < 3; ++b) {
                    cov[a][b] += E_hit[i] * r[a] * r[b];
                }
            }
            Eta_avg += x_hit[i].Eta()*E_hit[i];
            Phi_avg += x_hit[i].Phi()*E_hit[i];
        }
        cov *= (1.0 / E);
        Eta_avg *= (1.0 / E);
        Phi_avg *= (1.0 / E);

        TDecompSVD svd(cov);
        if (svd.Decompose()) {
            TMatrixD eigenVecs = svd.GetU();  // Principal components
            TVector3 axis(
                eigenVecs(0, 0),
                eigenVecs(1, 0),
                eigenVecs(2, 0)
            );
            axis = axis.Unit();

            std::vector<std::pair<double, float>> radialEnergies;
            for (size_t i = 0; i < x_hit.size(); ++i) {
                TVector3 diff = x_hit[i] - x;
                double proj = diff.Dot(axis);
                TVector3 closest = x + proj * axis;
                double r = (x_hit[i] - closest).Mag();

                radialEnergies.emplace_back(r, E_hit[i]);

                Weta2+=E_hit[i]*(x_hit[i].Eta()-Eta_avg)*(x_hit[i].Eta()-Eta_avg);
                Wphi2+=E_hit[i]*(x_hit[i].Phi()-Phi_avg)*(x_hit[i].Phi()-Phi_avg);
            }
            Weta2 *= (1.0 / E);
            Wphi2 *= (1.0 / E);

            std::sort(radialEnergies.begin(), radialEnergies.end());
            double cumE = 0;
            for (const auto& [r, e] : radialEnergies) {
                cumE += e;
                if (cumE >= 0.9 * E) {
                    R90 = r;
                    break;
                }
            }
        }
    }
}