#ifndef CRDECALSvc_C
#define CRDECALSvc_C

#include "CRDEcalSvc.h"

DECLARE_COMPONENT(CRDEcalSvc)

StatusCode
CRDEcalSvc::initialize() {
    m_digiHits.clear();
    StatusCode sc = Service::initialize();

    return sc;
}

StatusCode
CRDEcalSvc::finalize() {
    // clear or reset
    m_digiHits.clear();

    StatusCode sc = Service::finalize();
    return sc;
}
#endif
