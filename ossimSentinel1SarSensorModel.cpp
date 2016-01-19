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

  std::cout<<"Product type: "<<product_type<<std::endl;
  std::cout<<"Mode: "<<mode<<", swath: "<<swath<<", polarisation: "<<polarisation<<std::endl<<std::endl;
  
  // First, lookup position/velocity records
  std::vector<ossimRefPtr<ossimXmlNode> > xnodes;
  xmlDoc->findNodes("/product/generalAnnotation/orbitList/orbit",xnodes);

  std::cout<<"Reading orbit records ..."<<std::endl;
  std::cout<<"Number of orbit records found: "<<xnodes.size()<<std::endl;

  //TODO uncomment and adapt following code from s1_inverse to fill
  //SarSensorModel structure
/*
  for(auto itNode = xnodes.begin(); itNode!=xnodes.end();++itNode)
    {

    // Retrieve acquisition time
    ossimString att1 = "time";
    ossimString s;
    s = (*itNode)->findFirstNode(att1)->getText();
    s = s.replaceAllThatMatch("T"," ");
    std::cout<<s;
    boost::posix_time::ptime acqTime(boost::posix_time::time_from_string(s));
   
    // Retrieve ECEF position
    ossimEcefPoint pos;
    att1 = "position";
    ossimString att2 = "x";
    pos[0] = atof((*itNode)->findFirstNode(att1)->findFirstNode(att2)->getText().c_str());
    att2 = "y";
    pos[1] = atof((*itNode)->findFirstNode(att1)->findFirstNode(att2)->getText().c_str());
    att2 = "z";
    pos[2] = atof((*itNode)->findFirstNode(att1)->findFirstNode(att2)->getText().c_str());

    std::cout<<", ECEF pos: "<<pos[0]<<", "<<pos[1]<<", "<<pos[2];

    // Retrieve ECEF velocity
     ossimEcefVector vel;
    att1 = "velocity";
    att2 = "x";
    vel[0] = atof((*itNode)->findFirstNode(att1)->findFirstNode(att2)->getText().c_str());
    att2 = "y";
    vel[1] = atof((*itNode)->findFirstNode(att1)->findFirstNode(att2)->getText().c_str());
    att2 = "z";
    vel[2] = atof((*itNode)->findFirstNode(att1)->findFirstNode(att2)->getText().c_str());

    std::cout<<", ECEF vel: "<<vel[0]<<", "<<vel[1]<<", "<<vel[2]<<std::endl;

    //Add one orbits record
    theOrbitRecords.push_back(make_tuple(acqTime,pos,vel));
    }
  
  
  std::cout<<"Reading other useful values ..."<<std::endl;
  ossimString s = xmlDoc->getRoot()->findFirstNode("imageAnnotation/imageInformation/productFirstLineUtcTime")->getText();

    s = s.replaceAllThatMatch("T"," ");
  std::cout<<"Acquisition start time: "<<s<<std::endl;
  acqStartTime = boost::posix_time::time_from_string(s);

  s = xmlDoc->getRoot()->findFirstNode("imageAnnotation/imageInformation/productLastLineUtcTime")->getText();
  s = s.replaceAllThatMatch("T"," ");
  std::cout<<"Acquisition stop time: "<<s<<std::endl;
  acqStopTime = boost::posix_time::time_from_string(s);

  nbLines = xmlDoc->getRoot()->findFirstNode("imageAnnotation/imageInformation/numberOfLines")->getText().toUInt16();
  nbSamples = xmlDoc->getRoot()->findFirstNode("imageAnnotation/imageInformation/numberOfSamples")->getText().toUInt16();

  std::cout<<"Image size: "<<nbSamples<<" x "<<nbLines<<std::endl;


  //Parse the near range time (in seconds)

  nearRangeTime = xmlDoc->getRoot()->findFirstNode("imageAnnotation/imageInformation/slantRangeTime")->getText().toDouble();
  std::cout<<"Near range time: "<<nearRangeTime<<" s"<<std::endl;

  theNearRangeDistance = nearRangeTime * C /2;

  std::cout<<"Near range distance: "<<nearRangeDistance<< "m"<<std::endl;

  //Parse the range sampling rate

  rangeSamplingRate = xmlDoc->getRoot()->findFirstNode("generalAnnotation/productInformation/rangeSamplingRate")->getText().toDouble();
  double estimatedRangeRes = (1/rangeSamplingRate)*C/2;
  std::cout<<"Range sampling rate: "<<rangeSamplingRate<<" Hz (estimated range res: "<<estimatedRangeRes<<" m)"<<std::endl;

  //Parse the range resolution

  theRangeResolution = xmlDoc->getRoot()->findFirstNode("imageAnnotation/imageInformation/rangePixelSpacing")->getText().toDouble();
  std::cout<<"Range res from product: "<<rangeRes<<" m"<<std::endl;

  //Parse the radar frequency
  
  double theRadarFrequency = xmlDoc->getRoot()->findFirstNode("generalAnnotation/productInformation/radarFrequency")->getText().toDouble();

  std::cout<<"Radar frequency: "<<theRadarFrequency<<" Hz"<<std::endl;
  
  boost::posix_time::time_duration td = (acqStopTime - acqStartTime);
  
  double acquisition_duration = td.total_microseconds();

  std::cout<<"Acquisition duration: "<<boost::posix_time::to_simple_string(td)<<" ("<<acquisition_duration/1000000<<" s)"<<std::endl;

  double estimatedAzimuthTimeIntervalInMicroSeconds = acquisition_duration/nbLines;
  double prf = 1000000/estimatedAzimuthTimeIntervalInMicroSeconds;

  //Parse azimuth time interval

  theAzimuthTimeInterval = xmlDoc->getRoot()->findFirstNode("imageAnnotation/imageInformation/azimuthTimeInterval")->getText().toDouble()*1000000;

  std::cout<<"Estimated prf: "<<prf<<" Hz ("<<estimatedAzimuthTimeIntervalInMicroSeconds/1000000<<" s)"<<std::endl;
  std::cout<<"Azimuth time interval from product: "<<azimuthTimeIntervalInMicroSeconds/1000000<<" s"<<std::endl;
  
  std::cout<<"Done."<<std::endl<<std::endl;
  
  // Now read burst records as well
  BurstRecordVectorType burstRecords;

  std::cout<<"Reading burst records ..."<<std::endl;
  xnodes.clear();  
  xmlDoc->findNodes("/product/swathTiming/burstList/burst",xnodes);

  if(xnodes.empty())
    {
  std::cout<<"No burst records found, skipping"<<std::endl;
    }
  else
    {
    std::cout<<"Number of burst records found: "<<xnodes.size()<<std::endl; 

    unsigned int linesPerBurst = xmlDoc->getRoot()->findFirstNode("swathTiming/linesPerBurst")->getText().toUInt16();

    std::cout<<"LinesPerBurst: "<<linesPerBurst<<std::endl;

    unsigned int burstId(0);
    
    for(auto itNode = xnodes.begin(); itNode!=xnodes.end();++itNode,++burstId)
      { 
      ossimString att1 = "azimuthTime";
      ossimString s;
      s = (*itNode)->findFirstNode(att1)->getText();
      s = s.replaceAllThatMatch("T"," ");
      boost::posix_time::ptime azTime(boost::posix_time::time_from_string(s));

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

      unsigned long burstFirstValidLine = burstId*linesPerBurst + first_valid;
      unsigned long burstLastValidLine  = burstId*linesPerBurst + last_valid;
      
      
      boost::posix_time::ptime burstFirstValidTime = azTime + boost::posix_time::microseconds(first_valid*azimuthTimeIntervalInMicroSeconds);
      boost::posix_time::ptime burstLastValidTime = azTime + boost::posix_time::microseconds(last_valid*azimuthTimeIntervalInMicroSeconds);

      std::cout<<"Burst #"<<burstId<<std::endl;
      std::cout<<"FirstValidSample: "<<burstFirstValidLine<<" ("<<burstFirstValidTime<<")"<<std::endl;
      std::cout<<"LastValidSample: "<<burstLastValidLine<<" ("<<burstLastValidTime<<")"<<std::endl;
      
      burstRecords.push_back(std::make_tuple(burstFirstValidTime, burstFirstValidLine,burstLastValidTime,burstLastValidLine));
      }
   }
  std::cout<<"Done."<<std::endl<<std::endl;

  if(isGrd)
    {
    std::cout<<"Reading Slant range to Ground range coefficients ..."<<std::endl;

    xnodes.clear();  
    xmlDoc->findNodes("/product/coordinateConversion/coordinateConversionList/coordinateConversion",xnodes);
    std::cout<<"Number of records found: "<<xnodes.size()<<std::endl;

    for(auto itNode = xnodes.begin(); itNode!=xnodes.end();++itNode)
      {
      ossimString att1 = "azimuthTime";
      ossimString s;
      s = (*itNode)->findFirstNode(att1)->getText();
      s = s.replaceAllThatMatch("T"," ");
      boost::posix_time::ptime azTime(boost::posix_time::time_from_string(s));
    
      att1 = "sr0";
      double sr0 = (*itNode)->findFirstNode(att1)->getText().toDouble();

      att1 = "srgrCoefficients";
      s = (*itNode)->findFirstNode(att1)->getText();
      std::vector<ossimString> ssplit = s.split(" ");

      std::vector<double> coefs;

      for(auto cIt = ssplit.begin();cIt !=ssplit.end();++cIt)
        {
        coefs.push_back(cIt->toDouble());
        }

      srgrRecords.push_back(std::make_tuple(azTime,sr0,coefs));      
      }

    std::cout<<"Done."<<std::endl<<std::endl;
    }
  
  
  std::cout<<"Reading ground control points ..."<<std::endl;
  GCPVectorType gcps;
  xnodes.clear();
  xmlDoc->findNodes("/product/geolocationGrid/geolocationGridPointList/geolocationGridPoint",xnodes);
  std::cout<<"Number of GCPs found: "<<xnodes.size()<<std::endl;
  
  for(auto itNode = xnodes.begin(); itNode!=xnodes.end();++itNode)
    {
    // Retrieve acquisition time
    ossimString att1 = "azimuthTime";
    ossimString s;
    s = (*itNode)->findFirstNode(att1)->getText();
    s = s.replaceAllThatMatch("T"," ");
    boost::posix_time::ptime azTime(boost::posix_time::time_from_string(s));

    att1 = "slantRangeTime";
    double slantRangeTime = (*itNode)->findFirstNode(att1)->getText().toDouble();
    
    ossimDpt imPoint;
    att1 = "pixel";
    imPoint.x = (*itNode)->findFirstNode(att1)->getText().toDouble();
    // att1 = "line";
    // imPoint.y = (*itNode)->findFirstNode(att1)->getText().toDouble();

    if(!burstRecords.empty())
      {
    boost::posix_time::ptime acqStart;
    bool burstFound(false);
    unsigned long acqStartLine(0);
    
    for(auto bIt = burstRecords.rbegin();bIt!=burstRecords.rend() && !burstFound;++bIt)
      {
      if(azTime > std::get<0>(*bIt) && azTime < std::get<2>(*bIt))
        {
        burstFound = true;
        acqStart = std::get<0>(*bIt);
        acqStartLine = std::get<1>(*bIt);        
        }
      }

    if(!burstFound)
      {
      if(azTime < std::get<0>(burstRecords.front()))
        {
        acqStart = std::get<0>(burstRecords.front());
        acqStartLine = std::get<1>(burstRecords.front());
        } 
      else if (azTime > std::get<0>(burstRecords.back()))
        {
        acqStart = std::get<0>(burstRecords.back());
        acqStartLine = std::get<1>(burstRecords.back());
        }
      }
    boost::posix_time::time_duration timeSinceStart = (azTime-acqStart);
      
    double timeSinceStartInMicroSeconds = timeSinceStart.total_microseconds();
    imPoint.y = timeSinceStartInMicroSeconds/azimuthTimeIntervalInMicroSeconds + acqStartLine;
    }
    
    else
      {
       att1 = "line";
       imPoint.y = (*itNode)->findFirstNode(att1)->getText().toDouble();
       }
    ossimGpt geoPoint;
    att1 = "latitude";
    geoPoint.lat = (*itNode)->findFirstNode(att1)->getText().toDouble();
    att1 = "longitude";
    geoPoint.lon = (*itNode)->findFirstNode(att1)->getText().toDouble();
    att1 = "height";
    geoPoint.hgt = (*itNode)->findFirstNode(att1)->getText().toDouble();

    gcps.push_back(make_tuple(azTime,slantRangeTime,imPoint,geoPoint));
    }
  std::cout<<"Done."<<std::endl<<std::endl;
  */
}

} // namespace ossimplugins
