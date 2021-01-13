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
    double bar_Length = 40.*dd4hep::cm; 
    double bar_width  = 1.*dd4hep::cm;
    int Ngroups = 7;
    int Nbars = bar_Length/bar_width;


    //Define detector and motherVolume(world)
    dd4hep::DetElement sdet(det_name, x_det.id());
    dd4hep::Volume motherVol = theDetector.pickMotherVolume(sdet);

    //Envelope volume: 80*80*28 cm^3 box
    dd4hep::Material  air(theDetector.material("Air"));
    dd4hep::Volume envelopeVol(det_name, dd4hep::Box(Ngroups*4*bar_width, bar_Length*2, bar_Length*2),  air);
    envelopeVol.setVisAttributes(theDetector, "InvisibleWithChildren");
    dd4hep::PlacedVolume envelopePlv = motherVol.placeVolume(envelopeVol, dd4hep::Position(posx,0,0));
    envelopePlv.addPhysVolID("system",x_det.id());
    sdet.setPlacement(envelopePlv);

    //Loop to place crystal bars
		
    dd4hep::Material mat_BGO(theDetector.material("G4_BGO"));

    dd4hep::Volume bar_s0("bar_s0", dd4hep::Box(bar_width/2, bar_Length/2, bar_width/2), mat_BGO);
    bar_s0.setVisAttributes(theDetector, "Invisible");
    bar_s0.setSensitiveDetector(sens);
    dd4hep::Volume bar_s1("bar_s1", dd4hep::Box(bar_width/2, bar_width/2, bar_Length/2), mat_BGO);
    bar_s1.setVisAttributes(theDetector, "Invisible");
    bar_s1.setSensitiveDetector(sens);


    for(int igr=1; igr<=Ngroups;igr++){

      dd4hep::Volume group("group", dd4hep::Box(bar_width*2, bar_Length, bar_Length), mat_BGO);
      //group.setVisAttributes(theDetector, "InvisibleWithChildren");
      std::string groupname = "Group_"+std::to_string(igr);
      dd4hep::DetElement groupdet(sdet, groupname, detid);

      dd4hep::Volume block("block", dd4hep::Box(bar_width, bar_Length/2, bar_Length/2), mat_BGO);
      block.setVisAttributes(theDetector, "VisibleRed");
      dd4hep::DetElement sd(groupdet, "Block0", detid);

      //layer 0: bars along phi. length=barz_s0. Bar num=Nbars
      for(int ibar0=1;ibar0<=Nbars;ibar0++){
        dd4hep::PlacedVolume plv_bar0 = block.placeVolume(bar_s0, dd4hep::Position(-bar_width/2, 0, (2*ibar0-1)*bar_width/2-bar_Length/2));
        plv_bar0.addPhysVolID("slayer",0).addPhysVolID("bar",ibar0);
        std::string barname0 = "CrystalBar_s0_"+std::to_string(ibar0);
        dd4hep::DetElement bardet0(sd, barname0, detid);
        bardet0.setPlacement(plv_bar0);
      }

      //layer1 
      for(int ibar1=1;ibar1<=Nbars;ibar1++){
         dd4hep::PlacedVolume plv_bar1 = block.placeVolume(bar_s1, dd4hep::Position(bar_width/2, (2*ibar1-1)*bar_width/2-bar_Length/2, 0));
         plv_bar1.addPhysVolID("slayer",1).addPhysVolID("bar",ibar1);
         std::string barname1 = "CrystalBar_s1_"+std::to_string(ibar1);
         dd4hep::DetElement bardet1(sd, barname1, detid);
         bardet1.setPlacement(plv_bar1);
      }

      dd4hep::PlacedVolume plv0 = group.placeVolume(block, dd4hep::Position(-bar_width,0,0));
      plv0.addPhysVolID("block",1);
      sd.setPlacement(plv0);
      dd4hep::PlacedVolume plv1 = group.placeVolume(block, dd4hep::Position(bar_width,bar_Length/2,bar_Length/2));
      plv1.addPhysVolID("block",2);
      sd.setPlacement(plv1);
      dd4hep::PlacedVolume plv2 = group.placeVolume(block, dd4hep::Position(bar_width,-bar_Length/2,bar_Length/2));
      plv2.addPhysVolID("block",3);
      sd.setPlacement(plv2);
      dd4hep::PlacedVolume plv3 = group.placeVolume(block, dd4hep::Position(bar_width,-bar_Length/2,-bar_Length/2));
      plv3.addPhysVolID("block",4);
      sd.setPlacement(plv3);
      dd4hep::PlacedVolume plv4 = group.placeVolume(block, dd4hep::Position(bar_width,bar_Length/2,-bar_Length/2));
      plv4.addPhysVolID("block",5);
      sd.setPlacement(plv4);

      dd4hep::PlacedVolume groupVol = envelopeVol.placeVolume(group, dd4hep::Position((2*igr-1)*2*bar_width,0,0));
      groupVol.addPhysVolID("group",igr);
      groupdet.setPlacement(groupVol);
    }

    sens.setType("calorimeter");
    MYDEBUG("create_detector DONE. ");
    return sdet;
}

DECLARE_DETELEMENT(CRDEcalModule, create_detector)
