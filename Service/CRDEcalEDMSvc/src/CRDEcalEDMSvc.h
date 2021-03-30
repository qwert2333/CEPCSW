#ifndef CRDECALEDMSvc_h
#define CRDECALEDMSvc_h

#include <GaudiKernel/Service.h>
#include "CRDEcalEDMSvc/ICRDEcalEDMSvc.h"

class CRDEcalEDMSvc: public extends<Service, ICRDEcalEDMSvc> {
public:
    using extends::extends;

    StatusCode initialize() override;
    StatusCode finalize() override;

    void setDigiSystem( std::vector<CRDEcalEDM::CRDCaloBlock>& blockVec ) { m_CaloBlockVec = blockVec; }
    void getDigiSystem( std::vector<CRDEcalEDM::CRDCaloBlock>& blockVec ) { blockVec = m_CaloBlockVec; }
    void ClearSystem() { m_CaloBlockVec.clear(); }

private:
    std::vector<CRDEcalEDM::CRDCaloBlock> m_CaloBlockVec;

};

#endif
