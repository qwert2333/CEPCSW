//====================================================================
//  Detector description implementation for Chunxiu Liu's EcalMatrix
//--------------------------------------------------------------------
//
//  Author     : Tao Lin
//               Examples from lcgeo
//                   lcgeo/detector/calorimeter/
//
//====================================================================
#include "DD4hep/DetFactoryHelper.h"
#include "XML/Layering.h"
#include "XML/Utilities.h"
#include "DDRec/DetectorData.h"
#include "DDSegmentation/Segmentation.h"

#define MYDEBUG(x) std::cout << __FILE__ << ":" << __LINE__ << ": " << x << std::endl;
#define MYDEBUGVAL(x) std::cout << __FILE__ << ":" << __LINE__ << ": " << #x << ": " << x << std::endl;

using dd4hep::rec::LayeredCalorimeterData;
static dd4hep::Ref_t create_detector(dd4hep::Detector& theDetector,
                                     xml_h e,
                                     dd4hep::SensitiveDetector sens) {

    xml_det_t x_det = e;

    std::string det_name = x_det.nameStr();
    std::string det_type = x_det.typeStr();
    MYDEBUGVAL(det_name);
    MYDEBUGVAL(det_type);
    int detid = x_det.id();

    double posx = 180*dd4hep::cm;
    double sl_Length = 28.*dd4hep::cm;  //x. longitude
    double sl_height = 40.*dd4hep::cm;  //y. height
    double sl_width  = 1.*dd4hep::cm;   //z. 
    //int Nlayers = 28;
    int Nslices = sl_height/sl_width;    //z direction. 


    //Define detector and motherVolume(world)
    dd4hep::DetElement sdet(det_name, x_det.id());
    dd4hep::Volume motherVol = theDetector.pickMotherVolume(sdet);

    //Envelope volume: 80*80*28 cm^3 box
    dd4hep::Material  air(theDetector.material("Air"));
    dd4hep::Volume envelopeVol(det_name, dd4hep::Box(sl_Length/2, sl_height/2, Nslices*sl_width/2),  air);
    envelopeVol.setVisAttributes(theDetector, "VisibleRed");
    dd4hep::PlacedVolume envelopePlv = motherVol.placeVolume(envelopeVol, dd4hep::Position(posx+sl_Length/2.,0,0));
    envelopePlv.addPhysVolID("system",x_det.id());
    sdet.setPlacement(envelopePlv);

    dd4hep::Material mat_BGO(theDetector.material("G4_BGO"));

    dd4hep::Volume sliceVol("slice", dd4hep::Box(sl_Length/2, sl_height/2, sl_width/2), mat_BGO);
    sliceVol.setVisAttributes(theDetector, "VisibleRed");
    sliceVol.setSensitiveDetector(sens);

    for(int is=0;is<Nslices;is++){
      dd4hep::PlacedVolume plv_s = envelopeVol.placeVolume(sliceVol, dd4hep::Position(0,0,(2*is+1)*sl_width/2-sl_height/2));
		plv_s.addPhysVolID("slice",is);
      std::string sname = "Slice_"+std::to_string(is);
      dd4hep::DetElement slDet(sdet, sname, detid);	
      slDet.setPlacement(plv_s);
    }

    sens.setType("calorimeter");
    MYDEBUG("create_detector DONE. ");
    return sdet;
}

DECLARE_DETELEMENT(CRDEcalBlock, create_detector)
