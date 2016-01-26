//----------------------------------------------------------------------------
//
// "Copyright Centre National d'Etudes Spatiales"
//
// License:  LGPL-2
//
// See LICENSE.txt file in the top level directory for more details.
//
//----------------------------------------------------------------------------
// $Id$

#include <ossimTerraSarXSarSensorModel.h>
#include <ossim/base/ossimXmlDocument.h>

namespace ossimplugins
{

ossimTerraSarXSarSensorModel::ossimTerraSarXSarSensorModel()
{}
  

ossimTerraSarXSarSensorModel::ossimTerraSarXSarSensorModel(const ossimTerraSarXSarSensorModel& m)
  : ossimSarSensorModel(m)
{

}

/** Destructor */
ossimTerraSarXSarSensorModel::~ossimTerraSarXSarSensorModel()
{}

  void ossimTerraSarXSarSensorModel::readAnnotationFile(const std::string & annotationXml, const std::string & geoXml)
  {
  ossimRefPtr<ossimXmlDocument> xmlDoc = new ossimXmlDocument(annotationXml);
  
  //isGRD parse variant?
  std::string product_type = xmlDoc->getRoot()->findFirstNode("productInfo/productVariantInfo/productVariant")->getText();

  std::cout << "type " <<  product_type << std::endl;

  isGRD = (product_type == "MGD" || product_type == "GEC" || product_type == "EEC");

  // First, lookup position/velocity records
  std::vector<ossimRefPtr<ossimXmlNode> > xnodes;
  xmlDoc->findNodes("/level1Product/platform/orbit/stateVec",xnodes);
  
  std::cout << "Number of states " << xnodes.size() << std::endl;

  for(std::vector<ossimRefPtr<ossimXmlNode> >::iterator itNode = xnodes.begin(); itNode!=xnodes.end();++itNode)
    {
      OrbitRecordType orbitRecord;

    // Retrieve acquisition time
    ossimString att1 = "timeUTC";
    ossimString s;
    s = (*itNode)->findFirstNode(att1)->getText();
    s = s.replaceAllThatMatch("T"," ");
    orbitRecord.azimuthTime = boost::posix_time::time_from_string(s);

    // Retrieve ECEF position
    att1 = "posX";
    orbitRecord.position[0] = (*itNode)->findFirstNode(att1)->getText().toDouble();
    att1 = "posY";
    orbitRecord.position[1] = (*itNode)->findFirstNode(att1)->getText().toDouble();
    att1 = "posZ";
    orbitRecord.position[2] = (*itNode)->findFirstNode(att1)->getText().toDouble();

    // Retrieve ECEF velocity
     ossimEcefVector vel;
    att1 = "velX";
    orbitRecord.velocity[0] = (*itNode)->findFirstNode(att1)->getText().toDouble();
    att1 = "velY";
    orbitRecord.velocity[1] = (*itNode)->findFirstNode(att1)->getText().toDouble();
    att1 = "velZ";
    orbitRecord.velocity[2] = (*itNode)->findFirstNode(att1)->getText().toDouble();

    //Add one orbits record
    std::cout << "Add theOrbitRecords" << std::endl;
    theOrbitRecords.push_back(orbitRecord);
    }

  //Parse the near range time (in seconds)
  theNearRangeTime = xmlDoc->getRoot()->findFirstNode("productInfo/sceneInfo/rangeTime/firstPixel")->getText().toDouble();

  std::cout << "theNearRangeTime " << theNearRangeTime << std::endl;

  //Parse the range sampling rate
  theRangeSamplingRate = xmlDoc->getRoot()->findFirstNode("instrument/settings/RSF")->getText().toDouble();
  
  std::cout << "theRangeSamplingRate " << theRangeSamplingRate << std::endl;
  
  //Parse the range resolution
  theRangeResolution = xmlDoc->getRoot()->findFirstNode("productSpecific/complexImageInfo/slantRangeResolution")->getText().toDouble();

  std::cout << "theRangeResolution " << theRangeResolution << std::endl;
  
  //Parse the radar frequency 
  theRadarFrequency = xmlDoc->getRoot()->findFirstNode("instrument/settings/settingRecord/PRF")->getText().toDouble();

  std::cout << "theRadarFrequency " << theRadarFrequency << std::endl;

  //Manage only strip map product for now (one burst)

  //Parse azimuth time start/stop
  ossimString s1 = xmlDoc->getRoot()->findFirstNode("productInfo/sceneInfo/start/timeUTC")->getText();
  s1 = s1.replaceAllThatMatch("T"," ");
  TimeType azimuthTimeStart = boost::posix_time::time_from_string(s1);

  std::cout << "azimuthTimeStart " << azimuthTimeStart << std::endl;

  ossimString s2 = xmlDoc->getRoot()->findFirstNode("productInfo/sceneInfo/stop/timeUTC")->getText();
  s2 = s2.replaceAllThatMatch("T"," ");
  TimeType azimuthTimeStop = boost::posix_time::time_from_string(s2);

  std::cout << "azimuthTimeStop " << azimuthTimeStop << std::endl;

  double td = (azimuthTimeStop - azimuthTimeStart).total_microseconds(); 

  // numberOfRows  
  unsigned int numberOfRows = xmlDoc->getRoot()->findFirstNode("productInfo/imageDataInfo/imageRaster/numberOfRows")->getText().toUInt16();

  std::cout << "numberOfRows " << numberOfRows << std::endl;

  //Compute azimuth time interval
  theAzimuthTimeInterval = td / static_cast<double> (numberOfRows);

  std::cout << "theAzimuthTimeInterval " << theAzimuthTimeInterval  << " and 1/prf: " << (1 / theRadarFrequency) * 1000000 << std::endl;

  //For Terrasar-X only 1 burst is supported for now
  BurstRecordType burstRecord;

  burstRecord.startLine = 0;  
  burstRecord.azimuthStartTime = azimuthTimeStart;
  burstRecord.azimuthStopTime = azimuthTimeStop;

  burstRecord.endLine = numberOfRows - 1;

  theBurstRecords.push_back(burstRecord);

  //GRD (detected product)
  if(isGRD)
    {
    std::cerr << "Detected product not handle yet." << std::endl;  
    }

  //Parse GCPs
  ossimRefPtr<ossimXmlDocument> xmlGeo = new ossimXmlDocument(geoXml);

  xnodes.clear();
  xmlGeo->findNodes("/geoReference/geolocationGrid/gridPoint",xnodes);
  
  std::cout<<"Found "<<xnodes.size()<<" GCPs"<<std::endl;

  for(std::vector<ossimRefPtr<ossimXmlNode> >::iterator itNode = xnodes.begin(); itNode!=xnodes.end();++itNode)
    {
      GCPRecordType gcpRecord;
    
      // Get delta acquisition time
      ossimString att1 = "t";
      double deltaAzimuth = (*itNode)->findFirstNode(att1)->getText().toDouble();
      gcpRecord.azimuthTime = azimuthTimeStart + boost::posix_time::microseconds(deltaAzimuth * 1000000);

      //Get delta range time
      att1 = "tau";
      gcpRecord.slantRangeTime = theNearRangeTime + (*itNode)->findFirstNode(att1)->getText().toDouble();
    
      att1 = "col";
      gcpRecord.imPt.x = (*itNode)->findFirstNode(att1)->getText().toDouble();
    
      att1 = "row";
      gcpRecord.imPt.y = (*itNode)->findFirstNode(att1)->getText().toDouble();
      
      ossimGpt geoPoint;
      att1 = "lat";
      gcpRecord.worldPt.lat = (*itNode)->findFirstNode(att1)->getText().toDouble();
      att1 = "lon";
      gcpRecord.worldPt.lon = (*itNode)->findFirstNode(att1)->getText().toDouble();
      att1 = "height";
      gcpRecord.worldPt.hgt = (*itNode)->findFirstNode(att1)->getText().toDouble();

      theGCPRecords.push_back(gcpRecord);
    }

  this->optimizeTimeOffsetsFromGcps();
}

} // namespace ossimplugins
