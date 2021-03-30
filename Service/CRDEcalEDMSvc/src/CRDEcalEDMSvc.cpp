#include "CRDEcalEDMSvc.h"

DECLARE_COMPONENT(CRDEcalEDMSvc)

StatusCode
CRDEcalEDMSvc::initialize() {
    m_CaloBlockVec.clear();
    StatusCode sc = Service::initialize();

    return sc;
}

StatusCode
CRDEcalEDMSvc::finalize() {
    // clear or reset
    m_CaloBlockVec.clear();

    StatusCode sc = Service::finalize();
    return sc;
}

