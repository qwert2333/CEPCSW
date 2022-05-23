#ifndef PANDORAPLUS_ALG_C
#define PANDORAPLUS_ALG_C

#include "PandoraPlusPFAlg.h"


#define C 299.79  // unit: mm/ns
#define PI 3.141592653
using namespace std;
using namespace dd4hep;


DECLARE_COMPONENT( PandoraPlusPFAlg )

PandoraPlusPFAlg::PandoraPlusPFAlg(const std::string& name, ISvcLocator* svcLoc)
  : GaudiAlgorithm(name, svcLoc),
    _nEvt(0)
{
 
   
}

StatusCode PandoraPlusPFAlg::initialize()
{
  //Initialize Creator settings

  //Initialize Creators
  m_pTrackCreator      = new TrackCreator( m_pTrackCreatorSettings );
  m_pEcalHitsCreator   = new EcalHitsCreator( m_EcalHitsCreatorSettings );

  //Initialize services
  m_edmsvc = service<ICRDEcalSvc>("CRDEcalSvc");
  if ( !m_edmsvc )  throw "CRDEcalDigiAlg :Failed to find CRDEcalSvc ...";
	
  rndm.SetSeed(m_seed);
  std::cout<<"PandoraPlusPFAlg::initialize"<<std::endl;

  //Readin collections
  for(auto m_trkname : name_TrackCol) r_TrkCol.push_back(new TrkType(m_trkname, Gaudi::DataHandle::Reader, this)); 
  for(auto m_caloname : name_EcalHits) r_EcalHitCol.push_back(new CaloType(m_caloname, Gaudi::DataHandle::Reader, this));

  //Output file for code testing. 

  return GaudiAlgorithm::initialize();
}

StatusCode PandoraPlusPFAlg::execute()
{
  if(_nEvt==0) std::cout<<"PandoraPlusPFAlg::execute Start"<<std::endl;
  std::cout<<"Processing event: "<<_nEvt<<std::endl;
  if(_nEvt<m_Nskip){ _nEvt++;  return GaudiAlgorithm::initialize(); }


  //InitializeForNewEvent(); 
  m_DataCol.Clear();

  //Readin collections 
  //MCParticle
  const edm4hep::MCParticleCollection* const_MCPCol = r_MCParticleCol.get();
  for(int imc=0; imc<const_MCPCol->size(); imc++){
    edm4hep::MCParticle m_MCp = const_MCPCol->at(imc);
    m_DataCol.collectionMap_MC[name_MCParticleCol].push_back(m_MCp);
  }


  m_pTrackCreator->CreateTracks( m_DataCol, r_TrkCol );
  m_pEcalHitsCreator->CreateEcalHits( m_DataCol, r_EcalHitCol );

  //Perform PFA algorithm



  //Write Ana tuples


  //Reset
  m_edmsvc->ClearSystem();

  std::cout<<"Event: "<<_nEvt<<" is done"<<std::endl;
  _nEvt ++ ;
  return StatusCode::SUCCESS;
}

StatusCode PandoraPlusPFAlg::finalize()
{

  info() << "Processed " << _nEvt << " events " << endmsg;
  return GaudiAlgorithm::finalize();
}

#endif
