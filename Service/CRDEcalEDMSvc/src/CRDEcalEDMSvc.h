#ifndef CRDECALEDMSvc_h
#define CRDECALEDMSvc_h

#include <GaudiKernel/Service.h>
#include "CRDEcalEDMSvc/ICRDEcalEDMSvc.h"
#include "CRDEcalEDMSvc/CRDCaloBar.h"

class CRDEcalEDMSvc: public extends<Service, ICRDEcalEDMSvc> {
public:
    using extends::extends;

    StatusCode initialize() override;
    StatusCode finalize() override;

    void setDigiSystem( std::vector<CRDEcalEDM::DigiBlock> blockVec ) { m_DigiBlockVec = blockVec; }
    std::vector<CRDEcalEDM::DigiBlock> getDigiSystem() { return m_DigiBlockVec; }

private:
    std::vector<CRDEcalEDM::DigiBlock> m_DigiBlockVec;

};

#endif
