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

#include <ossimSentinel1SarSensorModel.h>
#include <ossim/base/ossimXmlDocument.h>

namespace ossimplugins
{

ossimSentinel1SarSensorModel::ossimSentinel1SarSensorModel()
{}
  

ossimSentinel1SarSensorModel::ossimSentinel1SarSensorModel(const ossimSentinel1SarSensorModel& m)
  : ossimSarSensorModel(m)
{}

/** Destructor */
ossimSentinel1SarSensorModel::~ossimSentinel1SarSensorModel()
{}

void ossimSentinel1SarSensorModel::readAnnotationFile(const std::string & annotationXml)
{
  ossimRefPtr<ossimXmlDocument> xmlDoc = new ossimXmlDocument(annotationXml);

  //Parse specific metadata for Sentinel1
  //TODO add as members to the Sentinel1SarSensorModel
  std::string product_type = xmlDoc->getRoot()->findFirstNode("adsHeader/productType")->getText();
  std::string mode = xmlDoc->getRoot()->findFirstNode("adsHeader/mode")->getText();
  std::string swath = xmlDoc->getRoot()->findFirstNode("adsHeader/swath")->getText();
  std::string polarisation = xmlDoc->getRoot()->findFirstNode("adsHeader/polarisation")->getText();
  
  //TODO add as member of the base class
  isGRD = (product_type == "GRD");

  // First, lookup position/velocity records
  std::vector<ossimRefPtr<ossimXmlNode> > xnodes;
  xmlDoc->findNodes("/product/generalAnnotation/orbitList/orbit",xnodes);

  //TODO uncomment and adapt following code from s1_inverse to fill
  //SarSensorModel structure

  for(std::vector<ossimRefPtr<ossimXmlNode> >::iterator itNode = xnodes.begin(); itNode!=xnodes.end();++itNode)
    {
    OrbitRecordType orbitRecord;

    // Retrieve acquisition time
    ossimString att1 = "time";
    ossimString s;
    s = (*itNode)->findFirstNode(att1)->getText();
    s = s.replaceAllThatMatch("T"," ");
    orbitRecord.azimuthTime = boost::posix_time::time_from_string(s);
   
    // Retrieve ECEF position
    att1 = "position";
    ossimString att2 = "x";
    orbitRecord.position[0] = (*itNode)->findFirstNode(att1)->findFirstNode(att2)->getText().toDouble();
    att2 = "y";
    orbitRecord.position[1] = (*itNode)->findFirstNode(att1)->findFirstNode(att2)->getText().toDouble();
    att2 = "z";
    orbitRecord.position[2] = (*itNode)->findFirstNode(att1)->findFirstNode(att2)->getText().toDouble();

    // Retrieve ECEF velocity
     ossimEcefVector vel;
    att1 = "velocity";
    att2 = "x";
    orbitRecord.velocity[0] = (*itNode)->findFirstNode(att1)->findFirstNode(att2)->getText().toDouble();
    att2 = "y";
    orbitRecord.velocity[1] = (*itNode)->findFirstNode(att1)->findFirstNode(att2)->getText().toDouble();
    att2 = "z";
    orbitRecord.velocity[2] = (*itNode)->findFirstNode(att1)->findFirstNode(att2)->getText().toDouble();

    //Add one orbits record
    theOrbitRecords.push_back(orbitRecord);
    }

  //Parse the near range time (in seconds)
  theNearRangeTime = xmlDoc->getRoot()->findFirstNode("imageAnnotation/imageInformation/slantRangeTime")->getText().toDouble();

  //Parse the range sampling rate
  theRangeSamplingRate = xmlDoc->getRoot()->findFirstNode("generalAnnotation/productInformation/rangeSamplingRate")->getText().toDouble();

  //Parse the range resolution
  theRangeResolution = xmlDoc->getRoot()->findFirstNode("imageAnnotation/imageInformation/rangePixelSpacing")->getText().toDouble();

  //Parse the radar frequency 
  theRadarFrequency = xmlDoc->getRoot()->findFirstNode("generalAnnotation/productInformation/radarFrequency")->getText().toDouble();

  //Parse azimuth time interval
  theAzimuthTimeInterval = xmlDoc->getRoot()->findFirstNode("imageAnnotation/imageInformation/azimuthTimeInterval")->getText().toDouble()*1000000;

  
  // Now read burst records as well
  xnodes.clear();  
  xmlDoc->findNodes("/product/swathTiming/burstList/burst",xnodes);

  if(xnodes.empty())
    {  
    BurstRecordType burstRecord;

    burstRecord.startLine = 0;
    
    ossimString s = xmlDoc->getRoot()->findFirstNode("imageAnnotation/imageInformation/productFirstLineUtcTime")->getText();

    s = s.replaceAllThatMatch("T"," ");
    
    burstRecord.azimuthStartTime = boost::posix_time::time_from_string(s);

    std::cout<< burstRecord.azimuthStartTime<<std::endl;

    s = xmlDoc->getRoot()->findFirstNode("imageAnnotation/imageInformation/productLastLineUtcTime")->getText();
    s = s.replaceAllThatMatch("T"," ");
    
    burstRecord.azimuthStopTime = boost::posix_time::time_from_string(s);
  
    burstRecord.endLine = xmlDoc->getRoot()->findFirstNode("imageAnnotation/imageInformation/numberOfLines")->getText().toUInt16()-1;

    theBurstRecords.push_back(burstRecord);
    }
  else
    {
    unsigned int linesPerBurst = xmlDoc->getRoot()->findFirstNode("swathTiming/linesPerBurst")->getText().toUInt16();
    unsigned int burstId(0);
    
    for(std::vector<ossimRefPtr<ossimXmlNode> >::iterator itNode = xnodes.begin(); itNode!=xnodes.end();++itNode,++burstId)
      { 
      BurstRecordType burstRecord;

      ossimString att1 = "azimuthTime";
      ossimString s;
      s = (*itNode)->findFirstNode(att1)->getText();
      s = s.replaceAllThatMatch("T"," ");
      ossimSarSensorModel::TimeType azTime(boost::posix_time::time_from_string(s));

      att1 = "firstValidSample";
      s = (*itNode)->findFirstNode(att1)->getText();

      long first_valid(0), last_valid(0);
      bool begin_found(false), end_found(false);
      
      std::vector<ossimString> ssp = s.split(" ");

      for(auto sIt = ssp.begin(); sIt != ssp.end() && !end_found;++sIt)
        {
        if(!begin_found)
          {
          if(*sIt!="-1")
            {
            begin_found = true;
            }
          else
            {
            ++first_valid;
            }
          ++last_valid;
          }
        else
          {
          if(!end_found && *sIt=="-1")
            {
            end_found = true;
            }
            else
            {
            ++last_valid;
          }
         }
        }

      burstRecord.startLine = burstId*linesPerBurst + first_valid;
      burstRecord.endLine = burstId*linesPerBurst + last_valid;
      
      burstRecord.azimuthStartTime = azTime + boost::posix_time::microseconds(first_valid*theAzimuthTimeInterval);
      burstRecord.azimuthStopTime = azTime + boost::posix_time::microseconds(last_valid*theAzimuthTimeInterval);
      
      theBurstRecords.push_back(burstRecord);
      }
   }

  if(isGRD)
    {
    xnodes.clear();  
    xmlDoc->findNodes("/product/coordinateConversion/coordinateConversionList/coordinateConversion",xnodes);

    for(std::vector<ossimRefPtr<ossimXmlNode> >::iterator itNode = xnodes.begin(); itNode!=xnodes.end();++itNode)
      {
      CoordinateConversionRecordType coordRecord;
      
      ossimString att1 = "azimuthTime";
      ossimString s;
      s = (*itNode)->findFirstNode(att1)->getText();
      s = s.replaceAllThatMatch("T"," ");
      coordRecord.azimuthTime = boost::posix_time::time_from_string(s);
    
      att1 = "sr0";
      coordRecord.rg0 = (*itNode)->findFirstNode(att1)->getText().toDouble();

      att1 = "srgrCoefficients";
      s = (*itNode)->findFirstNode(att1)->getText();
      std::vector<ossimString> ssplit = s.split(" ");

      for(auto cIt = ssplit.begin();cIt !=ssplit.end();++cIt)
        {
        coordRecord.coefs.push_back(cIt->toDouble());
        }
      assert(!coordRecord.coefs.empty()&&"The srgr record has empty coefs vector.");

      theSlantRangeToGroundRangeRecords.push_back(coordRecord);      
      }
    }
  
  xnodes.clear();
  xmlDoc->findNodes("/product/geolocationGrid/geolocationGridPointList/geolocationGridPoint",xnodes);

  for(std::vector<ossimRefPtr<ossimXmlNode> >::iterator itNode = xnodes.begin(); itNode!=xnodes.end();++itNode)
    {
    GCPRecordType gcpRecord;
    
    // Retrieve acquisition time
    ossimString att1 = "azimuthTime";
    ossimString s;
    s = (*itNode)->findFirstNode(att1)->getText();
    s = s.replaceAllThatMatch("T"," ");
    gcpRecord.azimuthTime = boost::posix_time::time_from_string(s);

    att1 = "slantRangeTime";
    gcpRecord.slantRangeTime = (*itNode)->findFirstNode(att1)->getText().toDouble();
    
    att1 = "pixel";
    gcpRecord.imPt.x = (*itNode)->findFirstNode(att1)->getText().toDouble();


    // In TOPSAR products, GCPs are weird (they fall in black lines
    // between burst. This code allows to move them to a valid area of
    // the image.
    if(theBurstRecords.size()>2)
      {
      ossimSarSensorModel::TimeType acqStart;
      bool burstFound(false);
      unsigned long acqStartLine(0);
    
      for(std::vector<BurstRecordType>::reverse_iterator bIt = theBurstRecords.rbegin();bIt!=theBurstRecords.rend() && !burstFound;++bIt)
      {
      if(gcpRecord.azimuthTime >= bIt->azimuthStartTime && gcpRecord.azimuthTime < bIt->azimuthStopTime)
        {
        burstFound = true;
        acqStart = bIt->azimuthStartTime;
        acqStartLine = bIt->startLine;        
        }
      }

    if(!burstFound)
      {
      if(gcpRecord.azimuthTime < theBurstRecords.front().azimuthStartTime)
        {
        acqStart = theBurstRecords.front().azimuthStartTime;
        acqStartLine = theBurstRecords.front().startLine;
        } 
      else if (gcpRecord.azimuthTime >=  theBurstRecords.front().azimuthStopTime)
        {
        acqStart = theBurstRecords.back().azimuthStartTime;
        acqStartLine = theBurstRecords.back().startLine;
        }
      }
    boost::posix_time::time_duration timeSinceStart = (gcpRecord.azimuthTime-acqStart);
      
    double timeSinceStartInMicroSeconds = timeSinceStart.total_microseconds();
    gcpRecord.imPt.y= timeSinceStartInMicroSeconds/theAzimuthTimeInterval + acqStartLine;
    }
    
    else
      {
      att1 = "line";
      gcpRecord.imPt.y = (*itNode)->findFirstNode(att1)->getText().toDouble();
       }
    ossimGpt geoPoint;
    att1 = "latitude";
    gcpRecord.worldPt.lat = (*itNode)->findFirstNode(att1)->getText().toDouble();
    att1 = "longitude";
    gcpRecord.worldPt.lon = (*itNode)->findFirstNode(att1)->getText().toDouble();
    att1 = "height";
    gcpRecord.worldPt.hgt = (*itNode)->findFirstNode(att1)->getText().toDouble();

    theGCPRecords.push_back(gcpRecord);
    }

}

} // namespace ossimplugins
