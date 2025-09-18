#ifndef TRACK_CREATOR_C
#define TRACK_CREATOR_C

#include "Tools/TrackCreator.h"

namespace Cyber{

  TrackCreator::TrackCreator(const Settings& m_settings) : settings( m_settings ){

  };


  StatusCode TrackCreator::CreateTracks( CyberDataCol& m_DataCol, 
                                         std::vector<DataHandle<edm4hep::TrackCollection>*>& r_TrackCols, 
                                         DataHandle<edm4hep::MCRecoTrackParticleAssociationCollection>* r_MCParticleTrkCol ){

    if(r_TrackCols.size()==0 || settings.map_stringVecPars.at("trackCollections").size()==0) return StatusCode::SUCCESS;

    //Save readin collections
    m_DataCol.collectionMap_Track.clear(); 
    for(int icol=0; icol<r_TrackCols.size(); icol++){
      const edm4hep::TrackCollection* const_TrkCol = r_TrackCols[icol]->get(); 

      std::vector<edm4hep::Track> m_TrkCol; m_TrkCol.clear();
      for(unsigned int icol=0; icol<const_TrkCol->size(); icol++){
        edm4hep::Track m_trk = const_TrkCol->at(icol);
        if(!m_trk.isAvailable() || m_trk.getType()==0) continue;
        m_TrkCol.push_back(m_trk);
      }

      m_DataCol.collectionMap_Track[ settings.map_stringVecPars.at("trackCollections")[icol] ] = m_TrkCol; 
    }


    //Convert to local objects
    std::vector<std::shared_ptr<Cyber::Track>> m_trkCol; m_trkCol.clear();
    const edm4hep::MCRecoTrackParticleAssociationCollection* const_MCPTrkAssoCol = r_MCParticleTrkCol->get();

    for(auto iter : m_DataCol.collectionMap_Track){
      auto const_TrkCol = iter.second; 
      for(int itrk=0; itrk<const_TrkCol.size(); itrk++){

        //Cyber::Track* m_trk = new Cyber::Track();
        std::shared_ptr<Cyber::Track> m_trk = std::make_shared<Cyber::Track>();
        std::vector<Cyber::TrackState> m_trkstates;

        for(int its=0; its<const_TrkCol[itrk].trackStates_size(); its++){
          Cyber::TrackState m_trkst;
          m_trkst.D0 = const_TrkCol[itrk].getTrackStates(its).D0;
          m_trkst.Z0 = const_TrkCol[itrk].getTrackStates(its).Z0;
          m_trkst.phi0 = const_TrkCol[itrk].getTrackStates(its).phi;
          m_trkst.tanLambda = const_TrkCol[itrk].getTrackStates(its).tanLambda;
          m_trkst.Omega = const_TrkCol[itrk].getTrackStates(its).omega;
          m_trkst.Kappa = m_trkst.Omega*1000./(0.299792458*settings.map_floatPars.at("BField"));   
          m_trkst.location = const_TrkCol[itrk].getTrackStates(its).location;
          m_trkst.referencePoint.SetXYZ( const_TrkCol[itrk].getTrackStates(its).referencePoint[0],
                                         const_TrkCol[itrk].getTrackStates(its).referencePoint[1],
                                         const_TrkCol[itrk].getTrackStates(its).referencePoint[2] );

          m_trkstates.push_back(m_trkst);        
        }
        m_trk->setTrackStates("Input", m_trkstates);
        m_trk->setType(const_TrkCol[itrk].getType());
        m_trk->setOriginTrack( const_TrkCol[itrk] );

        for(int ilink=0; ilink<const_MCPTrkAssoCol->size(); ilink++){
          if( const_TrkCol[itrk] == const_MCPTrkAssoCol->at(ilink).getRec() ) {
            m_trk->addLinkedMCP( std::make_pair(const_MCPTrkAssoCol->at(ilink).getSim(), const_MCPTrkAssoCol->at(ilink).getWeight()) );
            break;
          }
        }


        //Assign PID from dNdx and TOF
        if( m_DataCol.dNdxCol && m_DataCol.tofCol ){
          std::array<double, 5> chi2s; chi2s.fill(0);

          for (auto dqdx : *m_DataCol.dNdxCol){
            if (dqdx.getTrack() == const_TrkCol[itrk]){
              for (int i = 0; i < 5; i++){
                chi2s[i] = dqdx.getHypotheses(i).chi2;
              }
            }
          }

          for (auto tof : *m_DataCol.tofCol){
            if (tof.getTrack() == const_TrkCol[itrk]){
              double toft = tof.getTime();
              std::array<float, 5> tofexpts = tof.getTimeExp();
              double tofexpterr = tof.getSigma();
              std::array<double, 5> tofchi2s;
              for (int i = 0; i < 5; i++){
                tofchi2s[i] = std::pow( (tofexpts[i] - toft) / tofexpterr, 2);
                chi2s[i] += tofchi2s[i];
              }
            }
          }

          int minchi2idx = std::distance(chi2s.begin(), std::min_element(chi2s.begin(), chi2s.end()));
          int pdgid = m_trk->getCharge() * PDGIDs.at(minchi2idx);

          m_trk->setPID(pdgid);    
//cout<<"  Readin track #"<<itrk<<": pid "<<pdgid<<", truth MC pid "<<m_trk->getLeadingMCP().getPDG()<<endl;
        }

        m_trkCol.push_back(m_trk);
      }
    }

    if(settings.map_boolPars.at("DoCleanTrack")) SelectGoodTrack(m_trkCol);
    m_DataCol.TrackCol = m_trkCol;

    //Track extrapolation
    //  Write settings: geometry description
    //  m_TrkExtraSettings.map_floatPar["Nlayers"] = 28;

    if(settings.map_boolPars.at("UseTruthMatchTrk")) m_TrkExtraSettings.map_intPars["Input_track"] = 1;
    m_TrkExtraAlg = new TrackExtrapolatingAlg();
    m_TrkExtraAlg->ReadSettings(m_TrkExtraSettings);
    m_TrkExtraAlg->Initialize( m_DataCol );
    m_TrkExtraAlg->RunAlgorithm( m_DataCol );
    m_TrkExtraAlg->ClearAlgorithm();
    delete m_TrkExtraAlg;
    m_TrkExtraAlg = nullptr;

    return StatusCode::SUCCESS;
  }


  StatusCode TrackCreator::CreateTracksFromMCParticle(CyberDataCol& m_DataCol, 
                                          DataHandle<edm4hep::MCParticleCollection>& r_MCParticleCol){


    // Convert MC charged particles to local Track objects.
    // Assuming the tracks are ideal helixes.

    // Get charged MC particles with generatorStatus==1
    const edm4hep::MCParticleCollection* const_MCPCol = r_MCParticleCol.get();
    std::vector<edm4hep::MCParticle> m_MCPvec; m_MCPvec.clear(); 
    for(int i=0; i<const_MCPCol->size(); i++){
      edm4hep::MCParticle m_MCp = const_MCPCol->at(i);
      if(m_MCp.getGeneratorStatus()==1 && m_MCp.getCharge()!=0)
        m_MCPvec.push_back(m_MCp);
    }

    // Convert to local Track objects
    std::vector<std::shared_ptr<Cyber::Track>> m_trkCol; m_trkCol.clear();
    for(auto mcp : m_MCPvec){
      std::shared_ptr<Cyber::Track> m_trk = std::make_shared<Cyber::Track>();
      std::vector<Cyber::TrackState> m_trkstates;

      TVector3 mcp_vertex(mcp.getVertex().x, mcp.getVertex().y, mcp.getVertex().z);
      TVector3 mcp_p(mcp.getMomentum().x, mcp.getMomentum().y, mcp.getMomentum().z);
      double mcp_pT = TMath::Sqrt(mcp_p.X()*mcp_p.X() + mcp_p.Y()*mcp_p.Y());
      double charge = mcp.getCharge();

      // Evaluate track state at vertex
      Cyber::TrackState m_trkst;  
      m_trkst.location = 1;  // At IP
      m_trkst.D0 = 0;
      m_trkst.Z0 = 0;
      m_trkst.phi0 = TMath::ATan2(mcp_p.Y(), mcp_p.X());
      m_trkst.Kappa = 1 / mcp_pT;
      if(charge<0) m_trkst.Kappa = -m_trkst.Kappa;
      m_trkst.tanLambda = mcp_p.Z() / mcp_pT;
      m_trkst.Omega = 0.299792458 * settings.map_floatPars.at("BField") / 1000. / mcp_pT;
      m_trkst.referencePoint = mcp_vertex;
      m_trkstates.push_back(m_trkst);

      m_trk->setTrackStates("Input", m_trkstates);
      m_trk->setType(0);  // It is a "MC track" and not detected by any tracker system
      m_trk->addLinkedMCP( std::make_pair(mcp, 1.) );
      m_trkCol.push_back(m_trk);
    }
    m_DataCol.TrackCol = m_trkCol;

    m_TrkExtraSettings.map_intPars["Input_track"] = 1;

    m_TrkExtraAlg = new TrackExtrapolatingAlg();
    m_TrkExtraAlg->ReadSettings(m_TrkExtraSettings);
    m_TrkExtraAlg->Initialize( m_DataCol );
    m_TrkExtraAlg->RunAlgorithm( m_DataCol );
    m_TrkExtraAlg->ClearAlgorithm();
    delete m_TrkExtraAlg;
    m_TrkExtraAlg = nullptr;

    return StatusCode::SUCCESS;
  }


  StatusCode TrackCreator::SelectGoodTrack(std::vector<std::shared_ptr<Cyber::Track>>& trkCol){
    //Use truth matched tracks
    if(settings.map_boolPars.at("UseTruthMatchTrk")){

      std::map<edm4hep::MCParticle, std::vector<std::shared_ptr<Cyber::Track>>> truthMap;
      for(int itrk=0; itrk<trkCol.size(); itrk++){
        truthMap[ trkCol[itrk]->getLeadingMCP() ].push_back( trkCol[itrk] );
      }
      std::vector<std::shared_ptr<Cyber::Track>> tmp_selTrk; tmp_selTrk.clear();
      for(auto iter: truthMap){
        if(iter.second.size()<=0) continue;
        int m_index = -1;
        double minIP = 999999;
        for(int itrk=0; itrk<iter.second.size(); itrk++){
          double tmp_IP = sqrt( iter.second[itrk]->getD0()*iter.second[itrk]->getD0() + iter.second[itrk]->getZ0()*iter.second[itrk]->getZ0() );
          if(tmp_IP<minIP){
            minIP = tmp_IP;
            m_index = itrk;
          }
        }

        if(m_index>=0) tmp_selTrk.push_back( iter.second[m_index] );
      }
      for(int itrk=0; itrk<trkCol.size(); itrk++){
        if( find(tmp_selTrk.begin(), tmp_selTrk.end(), trkCol[itrk])==tmp_selTrk.end() ){
          trkCol.erase(trkCol.begin() + itrk);
          itrk--;
        }
      }
//cout<<"Matched track size: "<<trkCol.size()<<endl;

      //Create truth track from linked MC particle
      std::vector<std::shared_ptr<Cyber::Track>> truthtrk; 
      for(int itrk=0; itrk<trkCol.size(); itrk++){
        edm4hep::MCParticle mcp = trkCol[itrk]->getLeadingMCP();

        std::shared_ptr<Cyber::Track> m_trk = std::make_shared<Cyber::Track>();
        std::vector<Cyber::TrackState> m_trkstates;        

        TVector3 mcp_vertex(mcp.getVertex().x, mcp.getVertex().y, mcp.getVertex().z);
        TVector3 mcp_p(mcp.getMomentum().x, mcp.getMomentum().y, mcp.getMomentum().z);
        double mcp_pT = TMath::Sqrt(mcp_p.X()*mcp_p.X() + mcp_p.Y()*mcp_p.Y());
        double charge = mcp.getCharge();

        // Evaluate track state at vertex
        Cyber::TrackState m_trkst;
        m_trkst.location = 1;  // At IP
        m_trkst.D0 = 0;
        m_trkst.Z0 = 0;
        m_trkst.phi0 = TMath::ATan2(mcp_p.Y(), mcp_p.X());
        m_trkst.Kappa = 1 / mcp_pT;
        if(charge<0) m_trkst.Kappa = -m_trkst.Kappa;
        m_trkst.tanLambda = mcp_p.Z() / mcp_pT;
        m_trkst.Omega = 0.3 * settings.map_floatPars.at("BField") / 1000. / mcp_pT;
        m_trkst.referencePoint = mcp_vertex;
        m_trkstates.push_back(m_trkst);

        m_trk->setTrackStates("Input", m_trkstates);
        m_trk->setType(0);  // It is a "MC track" and not detected by any tracker system
        m_trk->addLinkedMCP( std::make_pair(mcp, 1.) );
        truthtrk.push_back(m_trk);       
      }
//cout<<"Created truth track size: "<<truthtrk.size()<<endl;
      trkCol = truthtrk;

      std::sort(trkCol.begin(), trkCol.end(),  compTrkP);
      return StatusCode::SUCCESS;
    }


    //Select tracks with some custom criteria
    float trk_P, trk_Pt, trk_Nhit, trk_cosT, trk_startR, trk_startZ, trk_endR, trk_endZ, trk_length, trk_IP;
    std::ifstream file(settings.map_stringPars.at("TrackIDWeightFile"));
    if(!file.good()){
      std::cout << "ERROR: Did not find BDT weight file. Will skip the track cleaning. " << std::endl;
      std::sort(trkCol.begin(), trkCol.end(),  compTrkP);
      return StatusCode::SUCCESS;
    }

    TMVA::Reader *mva_rdr = new TMVA::Reader("Silent");
    mva_rdr->AddVariable("trk_P", &trk_P);
    mva_rdr->AddVariable("trk_Pt", &trk_Pt);
    mva_rdr->AddVariable("trk_Nhit", &trk_Nhit);
    mva_rdr->AddVariable("trk_IP", &trk_IP);
    mva_rdr->AddVariable("trk_cosT", &trk_cosT);
    mva_rdr->AddVariable("trk_startR", &trk_startR);
    mva_rdr->AddVariable("trk_startZ", &trk_startZ);
    mva_rdr->AddVariable("trk_endR", &trk_endR);
    mva_rdr->AddVariable("trk_endZ", &trk_endZ);
    mva_rdr->AddVariable("trk_length", &trk_length);
    
    mva_rdr->BookMVA(settings.map_stringPars.at("TrackIDMethod"), settings.map_stringPars.at("TrackIDWeightFile"));    
//cout<<"Use BDT track cleaning. BDT cut "<<settings.map_floatPars.at("BDTCut")<<endl;
    for(int itrk=0; itrk<trkCol.size(); itrk++){
      //float tmp_IP = sqrt( trkCol[itrk]->getD0()*trkCol[itrk]->getD0() + trkCol[itrk]->getZ0()*trkCol[itrk]->getZ0() );
      trk_P = trkCol[itrk]->getMomentum();
      trk_Pt = trkCol[itrk]->getPt();
      trk_Nhit = trkCol[itrk]->getTrackerHits();
      trk_cosT = trkCol[itrk]->getP3().CosTheta();
      trk_startR = trkCol[itrk]->getStartPoint().Perp();
      trk_startZ = trkCol[itrk]->getStartPoint().z();
      trk_endR = trkCol[itrk]->getEndPoint().Perp();
      trk_endZ = trkCol[itrk]->getEndPoint().z();
      trk_length = trkCol[itrk]->getEndPoint().Perp() - trkCol[itrk]->getStartPoint().Perp();
      trk_IP = sqrt(trkCol[itrk]->getD0()*trkCol[itrk]->getD0() + trkCol[itrk]->getZ0()*trkCol[itrk]->getZ0());
      float BDTout = mva_rdr->EvaluateMVA(settings.map_stringPars.at("TrackIDMethod"));
//printf("  In track #%d: P %.3f, Pt %.3f, cosT %.3f, start (%.3f, %.3f), end (%.3f, %.3f), Nhit %.0f, length %.3f, IP %.3f, BDT %.3f \n ", itrk, trk_P, trk_Pt, trk_cosT, trk_startR, trk_startZ, trk_endR, trk_endZ, trk_Nhit, trk_length, trk_IP, BDTout);
      if(BDTout<settings.map_floatPars.at("BDTCut")){
        trkCol.erase(trkCol.begin() + itrk);
        itrk--;
      }
    }
    std::sort(trkCol.begin(), trkCol.end(),  compTrkP);
    delete mva_rdr;

    return StatusCode::SUCCESS;
  }

};
#endif
