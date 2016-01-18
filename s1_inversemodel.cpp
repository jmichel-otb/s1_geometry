#include <tuple>
#include <vector>

#include "boost/date_time/posix_time/posix_time.hpp"

#include <ossim/base/ossimXmlDocument.h>
#include <ossim/base/ossimRefPtr.h>
#include <ossim/base/ossimXmlNode.h>
#include <ossim/base/ossimEcefPoint.h>
#include <ossim/base/ossimEcefVector.h>
#include <ossim/base/ossimDpt.h>
#include <ossim/base/ossimGpt.h>


typedef std::tuple<boost::posix_time::ptime,ossimEcefPoint,ossimEcefVector> RecordType;
typedef std::tuple<boost::posix_time::ptime,double,ossimDpt,ossimGpt> GCPType;
typedef std::vector<RecordType> RecordVectorType;
typedef std::vector<GCPType> GCPVectorType;
typedef std::tuple<boost::posix_time::ptime,unsigned long,boost::posix_time::ptime,unsigned long> BurstRecordType;
typedef std::vector<BurstRecordType> BurstRecordVectorType;


const double C = 299792458;

double compute_doppler(const double& radarFreq,const ossimEcefVector& vel, const ossimEcefPoint& sensorPos, const ossimEcefPoint& inputPos)
{
  // Eq. 19, p. 25
  ossimEcefVector s2gVec = inputPos - sensorPos;

  return 2 * radarFreq * vel.dot(s2gVec);
}


ossimDpt inverse_loc(const double & radarFreq, const boost::posix_time::ptime acqStartTime, const double & azimythTimeIntervalInMicroSeconds, const double & nearRangeTime, const double & rangeSamplingRate, const RecordVectorType & records, const BurstRecordVectorType& burstRecords, const ossimGpt& worldPoint, boost::posix_time::ptime& estimatedTime, double & estimatedSlantRangeTime)
{
  // First convert lat/lon to ECEF
  ossimEcefPoint inputPt(worldPoint);

  auto it = records.begin();

  double lastDoppler = compute_doppler(radarFreq, std::get<2>(*it),std::get<1>(*it),inputPt);
  auto lastRecord = it;
  
  bool lastDopplerSign = lastDoppler < 0;
  bool found = false;
  
  ++it;

  auto currentRecord = it;
  double currentDoppler = lastDoppler;

  // Look for the consecutive records where doppler freq changes sign
  while(it!=records.end() && !found)
    {
    currentRecord = it;
    currentDoppler = compute_doppler(radarFreq, std::get<2>(*it),std::get<1>(*it),inputPt);

    bool currentDopplerSign = currentDoppler< 0;

    //std::cout<<"last doppler: "<<lastDoppler<<", current doppler: "<<currentDoppler<<std::endl;
    
    if(currentDopplerSign != lastDopplerSign)
      {
      found = true;
      }
    else
      {
      lastDoppler = currentDoppler;
      lastRecord = currentRecord;
      ++it;
      }
    }

  // TODO: Handle not found case here

  // now interpolate time and sensor position
  double interp = std::abs(lastDoppler)/(std::abs(lastDoppler)+std::abs(currentDoppler));

  //std::cout<<"Interpolation coef: "<<interp<<std::endl;
  
  // Here I do not now if microsecond is sufficient
  boost::posix_time::time_duration delta_td = std::get<0>(*currentRecord) - std::get<0>(*lastRecord);
  double deltat = delta_td.total_microseconds();
  
  //std::cout<<"Corresponding deltat: "<<deltat<<" µs"<<std::endl;
  
  boost::posix_time::time_duration td = boost::posix_time::microseconds(static_cast<unsigned long>(floor(interp * deltat+0.5)));
  
  estimatedTime = std::get<0>(*lastRecord)+td;
  
  //std::cout<<"Times: "<<boost::posix_time::to_simple_string(std::get<0>(*lastRecord))<<" < "<<boost::posix_time::to_simple_string(estimatedTime)<<" < "<<boost::posix_time::to_simple_string(std::get<0>(*currentRecord))<<std::endl;

  ossimEcefVector deltaPos = (std::get<1>(*currentRecord)-std::get<1>(*lastRecord));
  deltaPos = deltaPos * interp;  
  ossimEcefPoint interSensorPos = std::get<1>(*lastRecord) + deltaPos;
  
  //std::cout<<"Sensor positions: "<<std::get<1>(*lastRecord)<<" <
  //"<<interSensorPos<<" < "<<std::get<1>(*currentRecord)<<std::endl;
 
  
  // TODO: Handle possible bistatic bias correction here (not needed
  // for S1)

  // Bistatic correction (eq 25, p 28)


  // TODO: Look this in metadata (if it should be applied)
  bool bistatic_correction = false;

  if(bistatic_correction)
    {
    double halftrange = 1000000 * (interSensorPos-inputPt).magnitude()/C;
    boost::posix_time::time_duration bistatic_td = boost::posix_time::microseconds(static_cast<unsigned long>(floor(halftrange+0.5)));
    estimatedTime = estimatedTime + bistatic_td;

    std::cout<<"Bistatic td: "<<boost::posix_time::to_simple_string(bistatic_td)<<std::endl;

    interp = static_cast<double>((td+bistatic_td).total_microseconds())/static_cast<double>(delta_td.total_microseconds());

    deltaPos = (std::get<1>(*currentRecord)-std::get<1>(*lastRecord));
    deltaPos = deltaPos * interp;  
    interSensorPos = std::get<1>(*lastRecord) + deltaPos;
    }
  
  ossimDpt resp;

  if(burstRecords.empty())
    {
    // Now compute fractional line index
    // Eq 22 p 27

    boost::posix_time::time_duration timeSinceStart = (estimatedTime-acqStartTime);
  
    double timeSinceStartInMicroSeconds = timeSinceStart.total_microseconds();
  
    //std::cout<<"timeSinceStartInMicroSeconds: "<<timeSinceStartInMicroSeconds<<" µs"<<std::endl;
  
    resp.y = timeSinceStartInMicroSeconds/azimythTimeIntervalInMicroSeconds;
    }
  else
    {
  boost::posix_time::ptime acqStart;
    bool burstFound(false);
    unsigned long acqStartLine(0);

   
    
    for(auto bIt = burstRecords.rbegin();bIt!=burstRecords.rend() && !burstFound;++bIt)
      {
      if(estimatedTime > std::get<0>(*bIt) && estimatedTime < std::get<2>(*bIt))
        {
        burstFound = true;
        acqStart = std::get<0>(*bIt);
        acqStartLine = std::get<1>(*bIt);        
        }
      }

    if(!burstFound)
      {
      if(estimatedTime < std::get<0>(burstRecords.front()))
        {
        acqStart = std::get<0>(burstRecords.front());
        acqStartLine = std::get<1>(burstRecords.front());
        } 
      else if (estimatedTime > std::get<0>(burstRecords.back()))
        {
        acqStart = std::get<0>(burstRecords.back());
        acqStartLine = std::get<1>(burstRecords.back());
        }
      }

    std::cout<<boost::posix_time::to_simple_string(acqStart)<<" "<<acqStartLine<<std::endl;
    

      boost::posix_time::time_duration timeSinceStart = (estimatedTime-acqStart);
      
      double timeSinceStartInMicroSeconds = timeSinceStart.total_microseconds();
      
      std::cout<<"timeSinceStartInMicroSeconds: "<<timeSinceStartInMicroSeconds<<" µs "<< timeSinceStartInMicroSeconds/azimythTimeIntervalInMicroSeconds <<std::endl;
  
      resp.y = timeSinceStartInMicroSeconds/azimythTimeIntervalInMicroSeconds + acqStartLine;


      
    }
  
  // TODO: Here this is different for products with burst
  
  //std::cout<<"Estimated line position: "<<resp.y<<std::endl;

  // Now compute fractional sample index
  // Eq 23 and 24 p. 28

  double range_distance = (interSensorPos-inputPt).magnitude();
  estimatedSlantRangeTime = 2*range_distance/C;

  std::cout<<estimatedSlantRangeTime<<" "<<nearRangeTime<<std::endl;
  
  resp.x = (estimatedSlantRangeTime - nearRangeTime)*rangeSamplingRate; 

  // TODO: Handle GRD products here
  
  //std::cout<<"Estimated sample position: "<<resp.x<<std::endl;

  return resp;
}



int main(int argc, char * argv[])
{
  std::cout.precision(9);
  
  if(argc != 2)
    return EXIT_FAILURE;
  
  std::string annotationXml = argv[1];

  unsigned int nbLines = 0;
  unsigned int nbSamples = 0;
  boost::posix_time::ptime acqStartTime;
  boost::posix_time::ptime acqStopTime;
  double nearRangeTime(0.);
  double rangeSamplingRate(0.);
  
  RecordVectorType records;

  ossimRefPtr<ossimXmlDocument> xmlDoc = new ossimXmlDocument(annotationXml);
  

  std::string product_type = xmlDoc->getRoot()->findFirstNode("adsHeader/productType")->getText();
  std::string mode = xmlDoc->getRoot()->findFirstNode("adsHeader/mode")->getText();
  std::string swath = xmlDoc->getRoot()->findFirstNode("adsHeader/swath")->getText();
  std::string polarisation = xmlDoc->getRoot()->findFirstNode("adsHeader/polarisation")->getText();

  bool isGrd = (product_type == "GRD");


  std::cout<<"Product type: "<<product_type<<std::endl;
  std::cout<<"Mode: "<<mode<<", swath: "<<swath<<", polarisation: "<<polarisation<<std::endl<<std::endl;
  
  // First, lookup position/velocity records
  std::vector<ossimRefPtr<ossimXmlNode> > xnodes;
  xmlDoc->findNodes("/product/generalAnnotation/orbitList/orbit",xnodes);

  std::cout<<"Reading orbit records ..."<<std::endl;
  std::cout<<"Number of orbit records found: "<<xnodes.size()<<std::endl;

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

    records.push_back(make_tuple(acqTime,pos,vel));
    }
  std::cout<<"done."<<std::endl<<std::endl;
  

  std::cout<<"Reading other useful values ..."<<std::endl;
  ossimString s = xmlDoc->getRoot()->findFirstNode("imageAnnotation/imageInformation/productFirstLineUtcTime")->getText();
  //ossimString s = xmlDoc->getRoot()->findFirstNode("generalAnnotation/downlinkInformationList/downlinkInformation/firstLineSensingTime")->getText();
    s = s.replaceAllThatMatch("T"," ");
  std::cout<<"Acquisition start time: "<<s<<std::endl;
  acqStartTime = boost::posix_time::time_from_string(s);

  s = xmlDoc->getRoot()->findFirstNode("imageAnnotation/imageInformation/productLastLineUtcTime")->getText();
  //s = xmlDoc->getRoot()->findFirstNode("generalAnnotation/downlinkInformationList/downlinkInformation/lastLineSensingTime")->getText();
  s = s.replaceAllThatMatch("T"," ");
  std::cout<<"Acquisition stop time: "<<s<<std::endl;
  acqStopTime = boost::posix_time::time_from_string(s);

  nbLines = xmlDoc->getRoot()->findFirstNode("imageAnnotation/imageInformation/numberOfLines")->getText().toUInt16();
  nbSamples = xmlDoc->getRoot()->findFirstNode("imageAnnotation/imageInformation/numberOfSamples")->getText().toUInt16();

  std::cout<<"Image size: "<<nbSamples<<" x "<<nbLines<<std::endl;

  nearRangeTime = xmlDoc->getRoot()->findFirstNode("imageAnnotation/imageInformation/slantRangeTime")->getText().toDouble();
  std::cout<<"Near range time: "<<nearRangeTime<<" s"<<std::endl;

  double nearRangeDistance = nearRangeTime * C /2;

  std::cout<<"Near range distance: "<<nearRangeDistance<< "m"<<std::endl;

  rangeSamplingRate = xmlDoc->getRoot()->findFirstNode("generalAnnotation/productInformation/rangeSamplingRate")->getText().toDouble();
  double rangeRes = (1/rangeSamplingRate)*C/2;
  std::cout<<"Range sampling rate: "<<rangeSamplingRate<<" Hz (estimated range res: "<<rangeRes<<" m)"<<std::endl;

  double radarFrequency = xmlDoc->getRoot()->findFirstNode("generalAnnotation/productInformation/radarFrequency")->getText().toDouble();
  std::cout<<"Radar frequency: "<<radarFrequency<<" Hz"<<std::endl;
  
  boost::posix_time::time_duration td = (acqStopTime - acqStartTime);
  
  double acquisition_duration = td.total_microseconds();

  std::cout<<"Acquisition duration: "<<boost::posix_time::to_simple_string(td)<<" ("<<acquisition_duration/1000000<<" s)"<<std::endl;

  double estimatedAzimuthTimeIntervalInMicroSeconds = acquisition_duration/nbLines;
  double prf = 1000000/estimatedAzimuthTimeIntervalInMicroSeconds;

  double azimuthTimeIntervalInMicroSeconds = xmlDoc->getRoot()->findFirstNode("imageAnnotation/imageInformation/azimuthTimeInterval")->getText().toDouble()*1000000;

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

  unsigned int count = 1;

  for(auto itGcp = gcps.begin();itGcp!=gcps.end();++itGcp,++count)
    {
  boost::posix_time::ptime estimatedTime;
  double estimatedSlantRangeTime;

  ossimDpt estimatedPos = inverse_loc(radarFrequency,acqStartTime,azimuthTimeIntervalInMicroSeconds, nearRangeTime, rangeSamplingRate,records, burstRecords,std::get<3>(*itGcp),estimatedTime,estimatedSlantRangeTime);
  
  std::cout<<"ProcessingGCP #"<<count<<":"<<std::endl;
  std::cout<<"Position: "<<get<2>(*itGcp)<<", estimated: "<<estimatedPos<<", residual: "<<estimatedPos-get<2>(*itGcp)<<std::endl;
  std::cout<<"Azimuth time: "<<boost::posix_time::to_simple_string(get<0>(*itGcp))<<", estimated: "<<boost::posix_time::to_simple_string(estimatedTime)<<", residual: "<<boost::posix_time::to_simple_string(estimatedTime-get<0>(*itGcp))<<std::endl;
  std::cout<<"Slant range time: "<<get<1>(*itGcp)<<", estimated: "<<estimatedSlantRangeTime<<", residual: "<<get<1>(*itGcp)-estimatedSlantRangeTime<<std::endl;
  std::cout<<std::endl;
}
  

  return EXIT_SUCCESS;
}
