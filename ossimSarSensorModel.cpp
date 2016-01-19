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

#include <ossimSarSensorModel.h>

namespace ossimplugins
{

ossimSarSensorModel::ossimSarSensorModel()
  : theOrbitRecords(),
    theGCPRecords(),
    theBurstRecords(),
    theSlantRangeToGroundRangeRecords(),
    theRadarFrequency(0.),
    theAzimuthTimeInterval(0.),
    theNearRangeTime(0.),
    theRangeSamplingRate(0.),
    theRangeResolution(0.),
    theBistaticCorrectionNeeded(false),
    isGRD(false)
{}
  

ossimSarSensorModel::ossimSarSensorModel(const ossimSarSensorModel& m)
{
  this->theOrbitRecords = m.theOrbitRecords;
  this->theGCPRecords = m.theGCPRecords;
  this->theBurstRecords = m.theBurstRecords;
  this->theSlantRangeToGroundRangeRecords = m.theSlantRangeToGroundRangeRecords;
  this->theRadarFrequency = m.theRadarFrequency;
  this->theAzimuthTimeInterval = m.theAzimuthTimeInterval;
  this->theNearRangeTime = m.theNearRangeTime;
  this->theRangeSamplingRate = m.theRangeSamplingRate;
  this->theRangeResolution = m.theRangeResolution;
  this->theBistaticCorrectionNeeded = m.theBistaticCorrectionNeeded;
  this->isGRD = m.isGRD;
}

/** Destructor */
ossimSarSensorModel::~ossimSarSensorModel()
{}

void ossimSarSensorModel::lineSampleHeightToWorld(const ossimDpt& imPt, const double & heightEllipsoid, ossimGpt& worldPt) const
{
// Not implemented yet
}

void ossimSarSensorModel::lineSampleToWorld(const ossimDpt& imPt, ossimGpt& worldPt) const
{
// Not implemented yet
}

void ossimSarSensorModel::worldToLineSample(const ossimGpt& worldPt, ossimDpt & imPt) const
{
  // First compute azimuth and range time
  TimeType azimuthTime;
  double rangeTime;
  
  bool success = worldToAzimuthRangeTime(worldPt, azimuthTime, rangeTime);

  if(!success)
    {
    imPt.makeNan();
    return;
    }

  // Convert azimuth time to line
  azimuthTimeToLine(azimuthTime,imPt.x);
  
  if(isGRD)
    {
    // GRD case
    double groundRange(0);
    double nearGroundRange(0);
    slantRangeToGroundRange(rangeTime,azimuthTime,groundRange);
    slantRangeToGroundRange(theNearRangeTime,azimuthTime,groundRange);

    // Eq 32 p. 31
    imPt.x = (groundRange - nearGroundRange)/theRangeResolution;
    }
  else
    {
    // SLC case
    // Eq 23 and 24 p. 28
    imPt.x = (rangeTime - theNearRangeTime)*theRangeSamplingRate;
    }
  
}

bool ossimSarSensorModel::worldToAzimuthRangeTime(const ossimGpt& worldPt, TimeType & azimuthTime, double & rangeTime) const
{
  // First convert lat/lon to ECEF
  ossimEcefPoint inputPt(worldPt);

  // Compute zero doppler time
  TimeType interpTime;
  ossimEcefPoint interpSensorPos;
  ossimEcefVector interpSensorVel;
  
  bool success = zeroDopplerLookup(inputPt,azimuthTime,interpSensorPos,interpSensorVel);

  if(!success)
    {
    return false;
    }
  
  if(theBistaticCorrectionNeeded)
    {
    // Compute bistatic correction if needed
    DurationType bistaticCorrection;
    computeBistaticCorrection(inputPt,interpSensorPos,bistaticCorrection);

    // Update interpolated azimuth time
    azimuthTime += bistaticCorrection;

    // Update sensor position and velocity
    interpolateSensorPosVel(interpTime,interpSensorPos,interpSensorVel);
    }

  // rangeTime is the round-tripping time to target
  double range_distance = (interpSensorPos-inputPt).magnitude();
  rangeTime = 2*range_distance/C;

  return true;
}
  
void ossimSarSensorModel::computeRangeDoppler(const ossimEcefPoint & inputPt, const ossimEcefPoint & sensorPos, const ossimEcefVector sensorVel, double & range, double & doppler) const
{
  // eq. 19, p. 25
  ossimEcefVector s2gVec = inputPt - sensorPos;

  doppler = 2 * theRadarFrequency * sensorVel.dot(s2gVec);
  range = s2gVec.magnitude();
}

void ossimSarSensorModel::interpolateSensorPosVel(const TimeType & azimuthTime, ossimEcefPoint& sensorPos, ossimEcefVector& sensorVel, unsigned int deg) const
{
  assert(!theOrbitRecords.empty()&&"The orbit records vector is empty");
  
  // Lagrangian interpolation of sensor position and velocity
  
  unsigned int nBegin(0), nEnd(0);

  sensorPos[0] = 0;
  sensorPos[1] = 0;
  sensorPos[2] = 0;

  sensorVel[0] = 0;
  sensorVel[1] = 0;
  sensorVel[2] = 0;

  // First, we search for the correct set of record to use during
  // interpolation

  // If there are less records than degrees, use them all
  if(theOrbitRecords.size()<deg)
    {
    nEnd = theOrbitRecords.size()-1;
    }
  else
    {
    // Search for the deg number of records around the azimuth time
    unsigned int t_min_idx = 0;
    DurationType t_min = azimuthTime - theOrbitRecords.front().azimuthTime;

    if(t_min.is_negative())
      t_min = t_min.invert_sign();
    
    unsigned int count = 0;
    
    for(std::vector<OrbitRecordType>::const_iterator it = theOrbitRecords.begin();it!=theOrbitRecords.end();++it,++count)
      {
      DurationType current_time = azimuthTime-it->azimuthTime;

      if(current_time.is_negative())
        current_time = current_time.invert_sign();
      
      if(t_min > current_time)
        {
        t_min_idx = count;
        t_min = current_time;
        }
      }  
    nBegin = std::max((int)t_min_idx-(int)deg/2+1,(int)0);
    nEnd = std::min(nBegin+deg-1,(unsigned int)theOrbitRecords.size());
    nBegin = nEnd<theOrbitRecords.size()-1 ? nBegin : nEnd-deg+1;
    }

  // Compute lagrangian interpolation using records from nBegin to nEnd
  for(unsigned int i = nBegin; i < nEnd; ++i)
    {
    double w = 1.;
    
    for(unsigned int j = nBegin; j < nEnd; ++j)
      {

      if(j!=i)
        {
        double td1 = (azimuthTime - theOrbitRecords[j].azimuthTime).total_microseconds();
        double td2 = (theOrbitRecords[i].azimuthTime - theOrbitRecords[j].azimuthTime).total_microseconds();
        w*=td1/td2;
        }
      }

    sensorPos[0]+=w*theOrbitRecords[i].position[0];
    sensorPos[1]+=w*theOrbitRecords[i].position[1];
    sensorPos[2]+=w*theOrbitRecords[i].position[2];

    sensorVel[0]+=w*theOrbitRecords[i].velocity[0];
    sensorVel[1]+=w*theOrbitRecords[i].velocity[1];
    sensorVel[2]+=w*theOrbitRecords[i].velocity[2];
    }
}

void ossimSarSensorModel::slantRangeToGroundRange(const double & slantRange, const TimeType & azimuthTime, double & groundRange) const
{
  assert(theSlantRangeToGroundRangeRecords.empty()&&"The slant range to ground range records vector is empty.");

  // First, we need to find the correct pair of records for interpolation
  std::vector<CoordinateConversionRecordType>::const_iterator it = theSlantRangeToGroundRangeRecords.begin();

  CoordinateConversionRecordType srgrRecord  = *it;

  if(azimuthTime > srgrRecord.azimuthTime)
    {
    std::vector<CoordinateConversionRecordType>::const_iterator  previousRecord = it;
    
    ++it;

    std::vector<CoordinateConversionRecordType>::const_iterator nextRecord = it;

    bool found = false;

    // Look for the correct record
    while(it!=theSlantRangeToGroundRangeRecords.end() && !found)
      {
      if(azimuthTime >= previousRecord->azimuthTime
         && azimuthTime < nextRecord->azimuthTime)
        {
        found = true;
        }
      else
        {
        previousRecord = nextRecord;
        ++it;
        nextRecord = it;
        }      
      }
    if(!found)
      {
      srgrRecord = *previousRecord;
      }
    else
      {
      // If azimuth time is between 2 records, interpolate
      double interp = (azimuthTime-(previousRecord->azimuthTime)).total_microseconds()/static_cast<double>((nextRecord->azimuthTime-previousRecord->azimuthTime).total_microseconds());

      double sr0 = (1-interp) * previousRecord->rg0 + interp*nextRecord->rg0;
      
      std::vector<double> coefs;
      std::vector<double>::const_iterator pIt = previousRecord->coefs.begin();
      std::vector<double>::const_iterator nIt = nextRecord->coefs.begin();
      
      for(;pIt != previousRecord->coefs.end() && nIt != nextRecord->coefs.end();++pIt,++nIt)
        {
        coefs.push_back(interp*(*nIt)+(1-interp)*(*pIt));
        }

      srgrRecord.rg0 = sr0;
      srgrRecord.coefs = coefs;
      }
    }

// Now that we have the interpolated coefs, compute ground range
// from slant range
  double sr_minus_sr0 = slantRange -srgrRecord.rg0;

  double ground_range = 0;
  
  for(std::vector<double>::const_reverse_iterator cIt = srgrRecord.coefs.rbegin();cIt!=srgrRecord.coefs.rend();++cIt)
    {
    ground_range = *cIt + sr_minus_sr0*ground_range;
    }
}


bool ossimSarSensorModel::zeroDopplerLookup(const ossimEcefPoint & inputPt, TimeType & interpAzimuthTime, ossimEcefPoint & interpSensorPos, ossimEcefVector & interpSensorVel) const
{
  assert(!theOrbitRecords.empty()&&"Orbit records vector is empty()");

  
  std::vector<OrbitRecordType>::const_iterator it = theOrbitRecords.begin();

  double range1(0.), doppler1(0.), range2(0.), doppler2(0.);

  // Compute range and doppler of first record
  computeRangeDoppler(inputPt,it->position, it->velocity, range1, doppler1);

  std::vector<OrbitRecordType>::const_iterator record1 = it;
  
  bool dopplerSign1 = doppler1 < 0;
  bool found = false;
  
  ++it;

  std::vector<OrbitRecordType>::const_iterator record2;

  // Look for the consecutive records where doppler freq changes sign
  // Note: implementing a bisection algorithm here might be faster
  while(it!=theOrbitRecords.end() && !found)
    {
    record2 = it;

    // compute range and doppler of current record
    computeRangeDoppler(inputPt,it->position, it->velocity, range2, doppler2);
   
    bool dopplerSign2 = doppler2 <0;

    // If a change of sign is detected
    if(dopplerSign1 != dopplerSign2)
      {
      found = true;
      }
    else
      {
      doppler1 = doppler2;
      record1 = record2;
      ++it;
      }
    }

  // If not found, pass error to caller (not a programming error, but
  // eronous input parameters
  if(!found)
    {
    return false;
    }

  // now interpolate time and sensor position
  double interpDenom = std::abs(doppler1)+std::abs(doppler2);

  assert(interpDenom>0&&"Both doppler frequency are null in interpolation weight computation");
  
  double interp = std::abs(doppler2)/interpDenom;
  
  // Note that microsecond precision is used here
  DurationType delta_td = record2->azimuthTime - record1->azimuthTime;
  double deltat = static_cast<double>(delta_td.total_microseconds());

  // Compute interpolated time offset wrt record1
  DurationType td = boost::posix_time::microseconds(static_cast<unsigned long>(floor(interp * deltat+0.5)));

  // Compute interpolated azimuth time
  interpAzimuthTime = record1->azimuthTime+td;
  
  // Interpolate sensor position and velocity
  interpolateSensorPosVel(interpAzimuthTime,interpSensorPos, interpSensorVel);

  return true;
}

void ossimSarSensorModel::computeBistaticCorrection(const ossimEcefPoint & inputPt, const ossimEcefPoint & sensorPos, DurationType & bistaticCorrection) const
{
  // Bistatic correction (eq 25, p 28)
  double halftrange = 1000000 * (sensorPos-inputPt).magnitude()/C;
  bistaticCorrection= boost::posix_time::microseconds(static_cast<unsigned long>(floor(halftrange+0.5)));
}


void ossimSarSensorModel::azimuthTimeToLine(const TimeType & azimuthTime, double & line) const
{
  assert(!theBurstRecords.empty()&&"Burst records are empty (at least one burst should be available)");

  std::vector<BurstRecordType>::const_iterator currentBurst = theBurstRecords.begin();

  bool burstFound(false);
  
  // Look for the correct burst. In most cases the number of burst
  // records will be 1 (except for TOPSAR Sentinel1 products)
  for(; currentBurst!= theBurstRecords.end() && !burstFound; ++currentBurst)
    {
    
    if(azimuthTime >= currentBurst->azimuthStartTime
       && azimuthTime < currentBurst->azimuthStopTime)
      {
      burstFound = true;
      }
    }

  // If no burst is found, we will use the first (resp. last burst to
  // extrapolate line
  if(!burstFound && theBurstRecords.size()>1)
    {
    if(azimuthTime < theBurstRecords.front().azimuthStartTime)
      {
      currentBurst = theBurstRecords.begin();
      } 
    else if (azimuthTime > theBurstRecords.back().azimuthStopTime)
      {
      currentBurst = theBurstRecords.end()-1;
      }
    }

  DurationType timeSinceStart = azimuthTime - currentBurst->azimuthStartTime;
  double timeSinceStartInMicroSeconds = static_cast<double>(timeSinceStart.total_microseconds());

  // Eq 22 p 27
  line = (timeSinceStartInMicroSeconds/theAzimuthTimeInterval) + currentBurst->startLine;
}

}
