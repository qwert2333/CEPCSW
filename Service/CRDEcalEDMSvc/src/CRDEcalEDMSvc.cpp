#include "GearSvc/IGearSvc.h"
#include "gear/GearMgr.h"
#include "CRDEcalEDMSvc.h"

DECLARE_COMPONENT(CRDEcalEDMSvc)

CRDEcalEDMSvc::CRDEcalEDMSvc(const std::string& name, ISvcLocator* svc)
  : base_class(name, svc){
}

CRDEcalEDMSvc::~CRDEcalEDMSvc(){
}

StatusCode CRDEcalEDMSvc::initialize(){
  return StatusCode::SUCCESS;
}
StatusCode CRDEcalEDMSvc::finalize(){
  return StatusCode::SUCCESS;
}
