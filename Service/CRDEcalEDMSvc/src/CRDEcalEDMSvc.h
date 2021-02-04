#ifndef CRDEcalEDMSvc_h
#define CRDEcalEDMSvc_h

#include <GaudiKernel/Service.h>

class CRDEcalEDMSvc : public extends<Service, ICRDEcalEDMSvc>{
 public:
  CRDEcalEDMSvc(const std::string& name, ISvcLocator* svc);
  ~CRDEcalEDMSvc();

  StatusCode initialize() override;
  StatusCode finalize() override;

};

#endif
