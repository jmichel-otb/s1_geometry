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
#include <ossim/base/ossimLsrSpace.h>

namespace ossimplugins
{

ossimSarSensorModel::ossimSarSensorModel()
  : theOrbitRecords(),
    theGCPRecords(),
    theBurstRecords(),
    theSlantRangeToGroundRangeRecords(),
    theGroundRangeToSlantRangeRecords(),
    theRadarFrequency(0.),
    theAzimuthTimeInterval(0.),
    theNearRangeTime(0.),
    theRangeSamplingRate(0.),
    theRangeResolution(0.),
    theBistaticCorrectionNeeded(false),
    isGRD(false)
{}
  

ossimSarSensorModel::ossimSarSensorModel(const ossimSarSensorModel& m)
  : ossimSensorModel(m)
{
  this->theOrbitRecords = m.theOrbitRecords;
  this->theGCPRecords = m.theGCPRecords;
  this->theBurstRecords = m.theBurstRecords;
  this->theSlantRangeToGroundRangeRecords = m.theSlantRangeToGroundRangeRecords;
  this->theGroundRangeToSlantRangeRecords = m.theGroundRangeToSlantRangeRecords;
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

void ossimSarSensorModel::lineSampleHeightToWorld(const ossimDpt& imPt, const double & heightAboveEllipsoid, ossimGpt& worldPt) const
{
  assert(!theGCPRecords.empty()&&"theGCPRecords is empty.");
  
  // Not implemented yet
  double rangeTime;
  TimeType azimuthTime;

  bool success = lineSampleToAzimuthRangeTime(imPt,azimuthTime,rangeTime);

  if(!success)
    {
    worldPt.makeNan();
    return;
    }

  // Find the closest GCP
  double distance = (imPt-theGCPRecords.front().imPt).length();  

  std::vector<GCPRecordType>::const_iterator refGcp = theGCPRecords.begin();

  for(std::vector<GCPRecordType>::const_iterator gcpIt = theGCPRecords.begin();
      gcpIt!=theGCPRecords.end();++gcpIt)
    {
    if((imPt-gcpIt->imPt).length() < distance)
      {
      distance = (imPt-gcpIt->imPt).length();
      refGcp = gcpIt;
      }
    }

  ossimGpt refPt = refGcp->worldPt;

  // Set the height reference
  ossim_float64 hgtSet;
  if ( ossim::isnan(heightAboveEllipsoid) )
    {
    hgtSet = refPt.height();
    }
  else
    {
    hgtSet = heightAboveEllipsoid;
    }
  ossimHgtRef hgtRef(AT_HGT, hgtSet);

  ossimEcefPoint sensorPos;
  ossimEcefVector sensorVel;
  
  interpolateSensorPosVel(azimuthTime,sensorPos,sensorVel);
  
  //double range, doppler;
  ossimEcefPoint ellPt;
  projToSurface(refPt,azimuthTime,rangeTime,&hgtRef,ellPt);

  ossimGpt gpt(ellPt);

  
  worldPt = ossimGpt(ellPt);
}

void ossimSarSensorModel::lineSampleToWorld(const ossimDpt& imPt, ossimGpt& worldPt) const
{
// Not implemented yet
}

void ossimSarSensorModel::worldToLineSample(const ossimGpt& worldPt, ossimDpt & imPt) const
{
  assert(theRangeResolution>0&&"theRangeResolution is null.");
  
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
  azimuthTimeToLine(azimuthTime,imPt.y);
  
  if(isGRD)
    {  
    // GRD case
    double groundRange(0);
    double nearGroundRange(0);
    slantRangeToGroundRange(rangeTime*C/2,azimuthTime,groundRange);
    
    // Eq 32 p. 31
    imPt.x = groundRange/theRangeResolution;
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
  double rangeDistance = (interpSensorPos-inputPt).magnitude();
  rangeTime = 2*rangeDistance/C;

  return true;
}

bool ossimSarSensorModel::lineSampleToAzimuthRangeTime(const ossimDpt & imPt, TimeType & azimuthTime, double & rangeTime) const
{
  // First compute azimuth time here
  bool success = lineToAzimuthTime(imPt.y,azimuthTime);

  if(!success)
    {
    return false;
    }
  
  
  // Then compute range time
  if(isGRD)
    {
    // Handle grd case here
    double slantRange;
    groundRangeToSlantRange(imPt.x*theRangeResolution,azimuthTime, slantRange);
    rangeTime = 2*slantRange/C;
    }
  else
    {
    rangeTime = theNearRangeTime + imPt.x*(1/theRangeSamplingRate);
    }

  return true;
}
  
void ossimSarSensorModel::computeRangeDoppler(const ossimEcefPoint & inputPt, const ossimEcefPoint & sensorPos, const ossimEcefVector sensorVel, double & range, double & doppler) const
{
  assert(theRadarFrequency>0&&"theRadarFrequency is null");
  
  // eq. 19, p. 25
  ossimEcefVector s2gVec = inputPt - sensorPos;

  range = s2gVec.magnitude();

  double coef = -2*C/(theRadarFrequency*range);
 
  doppler = coef * sensorVel.dot(s2gVec);
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
  applyCoordinateConversion(slantRange,azimuthTime,theSlantRangeToGroundRangeRecords,groundRange);
}

void ossimSarSensorModel::groundRangeToSlantRange(const double & groundRange, const TimeType & azimuthTime, double & slantRange) const
{
  applyCoordinateConversion(groundRange,azimuthTime,theGroundRangeToSlantRangeRecords,slantRange);
}

void ossimSarSensorModel::applyCoordinateConversion(const double & in, const TimeType& azimuthTime, const std::vector<CoordinateConversionRecordType> & records, double & out) const
{
  assert(!records.empty()&&"The records vector is empty.");

  // First, we need to find the correct pair of records for interpolation
  std::vector<CoordinateConversionRecordType>::const_iterator it = records.begin();

  CoordinateConversionRecordType srgrRecord;
  
  std::vector<CoordinateConversionRecordType>::const_iterator  previousRecord = it;
  ++it;
  
  std::vector<CoordinateConversionRecordType>::const_iterator nextRecord = it;
  
  bool found = false;
  
  // Look for the correct record
  while(it!=records.end() && !found)
    {
    nextRecord = it;
    
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
    if(azimuthTime < records.front().azimuthTime)
      {
      srgrRecord = records.front();
      }
    else if(azimuthTime >= records.back().azimuthTime)
      {
      srgrRecord = records.back();
      }
    }
  else
    {  
    assert(!previousRecord->coefs.empty()&&"previousRecord coefficients vector is empty.");
    assert(!nextRecord->coefs.empty()&&"nextRecord coefficients vector is empty.");
    
    // If azimuth time is between 2 records, interpolate
    double interp = (azimuthTime-(previousRecord->azimuthTime)).total_microseconds()/static_cast<double>((nextRecord->azimuthTime-previousRecord->azimuthTime).total_microseconds());
    
    srgrRecord.rg0 = (1-interp) * previousRecord->rg0 + interp*nextRecord->rg0;
    
    srgrRecord.coefs.clear();
    std::vector<double>::const_iterator pIt = previousRecord->coefs.begin();
    std::vector<double>::const_iterator nIt = nextRecord->coefs.begin();
    
    for(;pIt != previousRecord->coefs.end() && nIt != nextRecord->coefs.end();++pIt,++nIt)
      {
      srgrRecord.coefs.push_back(interp*(*nIt)+(1-interp)*(*pIt));
      }
    
    assert(!srgrRecord.coefs.empty()&&"Slant range to ground range interpolated coefficients vector is empty.");
    }

  // Now that we have the interpolated coefs, compute ground range
  // from slant range
  double sr_minus_sr0 =  in-srgrRecord.rg0;

  assert(!srgrRecord.coefs.empty()&&"Slant range to ground range coefficients vector is empty.");

  out = 0;
  
  for(std::vector<double>::const_reverse_iterator cIt = srgrRecord.coefs.rbegin();cIt!=srgrRecord.coefs.rend();++cIt)
    {
    out = *cIt + sr_minus_sr0*out;
    }



}


bool ossimSarSensorModel::zeroDopplerLookup(const ossimEcefPoint & inputPt, TimeType & interpAzimuthTime, ossimEcefPoint & interpSensorPos, ossimEcefVector & interpSensorVel) const
{
  assert(!theOrbitRecords.empty()&&"Orbit records vector is empty()");

  
  std::vector<OrbitRecordType>::const_iterator it = theOrbitRecords.begin();

  double doppler1(0.), doppler2(0.);

  // Compute range and doppler of first record
  // NOTE: here we only use the scalar product with vel and discard
  // the constant coef as it has no impact on doppler sign
  
  doppler1 = (inputPt-it->position).dot(it->velocity);
  
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
    doppler2 = (inputPt-it->position).dot(it->velocity);
 
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
  
  double interp = std::abs(doppler1)/interpDenom;
  
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
  for(std::vector<BurstRecordType>::const_iterator it = theBurstRecords.begin(); it!= theBurstRecords.end() && !burstFound; ++it)
    {
    
    if(azimuthTime >= it->azimuthStartTime
       && azimuthTime < it->azimuthStopTime)
      {
      burstFound = true;
      currentBurst = it;
      }
    }

  // If no burst is found, we will use the first (resp. last burst to
  // extrapolate line
  if(!burstFound)
    {
  
    if(theBurstRecords.size()>1)
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
    else
      {
      // Fall back to the only record
      currentBurst = theBurstRecords.begin();
      }
    }

  DurationType timeSinceStart = azimuthTime - currentBurst->azimuthStartTime;
  double timeSinceStartInMicroSeconds = static_cast<double>(timeSinceStart.total_microseconds());

  // Eq 22 p 27
  line = (timeSinceStartInMicroSeconds/theAzimuthTimeInterval) + currentBurst->startLine;
}

bool ossimSarSensorModel::lineToAzimuthTime(const double & line, TimeType & azimuthTime) const
{
  assert(!theBurstRecords.empty()&&"Burst records are empty (at least one burst should be available)");

  std::vector<BurstRecordType>::const_iterator currentBurst = theBurstRecords.begin();

  bool burstFound(false);

  if(theBurstRecords.size() == 1)
    {
    burstFound = true;
    }
  else
    {
    // Look for the correct burst. In most cases the number of burst
    // records will be 1 (except for TOPSAR Sentinel1 products)
    for(std::vector<BurstRecordType>::const_iterator it = theBurstRecords.begin(); it!= theBurstRecords.end() && !burstFound; ++it)
      {
    
      if(line >= it->startLine
         && line < it->endLine)
        {
        burstFound = true;
        currentBurst = it;
        }
      }
    }

  // If no burst is found, we will use the first (resp. last burst to
  // extrapolate line
  if(!burstFound)
    {
    return false;
    }

  double timeSinceStartInMicroSeconds = (line - currentBurst->startLine)*theAzimuthTimeInterval;
  
  DurationType timeSinceStart = boost::posix_time::microseconds(timeSinceStartInMicroSeconds);

  // Eq 22 p 27
  azimuthTime = currentBurst->azimuthStartTime + timeSinceStart;

  return true;
}



bool ossimSarSensorModel::projToSurface(const ossimEcefPoint& initPt, const TimeType & azimuthTime, const double & rangeTime, const ossimHgtRef * hgtRef, ossimEcefPoint & ellPt) const
{ 
  // Set slopes for tangent plane
  ossim_float64 sx  = 0.0;
  ossim_float64 sy  = 0.0;

  // Set tangent plane normal vector in ENU
  ossimColumnVector3d tpn(sx, sy, 1.0);
  ossimColumnVector3d tpnn(-sx, -sy, 1.0);
  
  // Initialize at OP point
   ossimEcefPoint rg(initPt);

   // Matrices
   NEWMAT::SymmetricMatrix BtB(3);
   NEWMAT::ColumnVector BtF(3);
   NEWMAT::ColumnVector F(3);
   NEWMAT::ColumnVector dR(3);

   ossimEcefPoint sensorPos;
   ossimEcefVector sensorVel;

   interpolateSensorPosVel(azimuthTime,sensorPos,sensorVel);

   double range = rangeTime * C/2;
   
   // Ensure we make at least one loop
   bool init = true;
   
   ossim_int32 iter = 0;
   // Implement simple Levenberg-Marquart algorithm
   // Convergence thresholds:
   // F(1) Thresh on range diff = 0.01 m
   // F(2) Thresh on doppler diff = 0.001 Hz
   // F(3) Thresh on elevation diff = 0.1 m
   while ((init || abs(F(1))>=0.1 || abs(F(2))>=0.001 || abs(F(3))>=0.1) && iter<10)
     {
     if(init)
       init =false;
     
     // Compute current latitude/longitude estimate
     ossimGpt pg(rg);
   
     // Set reference point @ desired elevation
     ossim_float64 atHgt = hgtRef->getRefHeight(pg);
     pg.height(atHgt);
     ossimEcefPoint rt(pg);
      
      // Define ENU space at reference point
      ossimLsrSpace enu(pg);
      
      // Rotate normal vector to ECF
      ossimEcefVector st = enu.lsrToEcefRotMatrix()*tpn;
      
      // Compute current range & Doppler estimate
      ossim_float64 rngComp;
      ossim_float64 dopComp;
      computeRangeDoppler(rg, sensorPos, sensorVel, rngComp, dopComp);
      
      // Compute current height estimate
      ossim_float64 diffHgt = st.dot(rg-rt);
      
      // Compute current fr, fd, ft
      F(1) = rngComp-range;
      F(2) = dopComp;
      F(3) = diffHgt;
      
      // Delta use for partial derivatives estimation (in meters)
      double d = 10.;

      // Compute partial derivative
      ossimEcefVector p_fr, p_fd, p_ft,dx(d,0,0),dy(0,d,0),dz(0,0,d);

      ossim_float64 rdx,rdy,rdz, fdx,fdy,fdz;

      computeRangeDoppler(rg+dx,sensorPos,sensorVel,rdx,fdx);
      p_fr[0] = (rngComp - rdx)/d;
      p_fd[0] = (dopComp - fdx)/d;
      p_fr[0]= (st.dot(-dx))/d;

      computeRangeDoppler(rg+dy,sensorPos,sensorVel,rdy,fdy);
      p_fr[1] = (rngComp - rdy)/d;
      p_fd[1] = (dopComp - fdy)/d;
      p_ft[1] = (st.dot(-dy))/d;

      computeRangeDoppler(rg+dz,sensorPos,sensorVel,rdz,fdz);
      p_fr[2] = (rngComp - rdz)/d;
      p_fd[2] = (dopComp - fdz)/d;
      p_ft[2] = (st.dot(-dz))/d;
      
      // Form B-matrix
      NEWMAT::Matrix B = ossimMatrix3x3::create(p_fr[0], p_fr[1], p_fr[2],
                                                p_fd[0], p_fd[1], p_fd[2],
                                                p_ft[0], p_ft[1], p_ft[2]);
      
      // Form coefficient matrix & discrepancy vector
      BtF << B.t()*F;
      BtB << B.t()*B;
      // Solve system
      dR = solveLeastSquares(BtB, BtF);
      // Update estimate
      for (ossim_int32 k=0; k<3; k++)
         rg[k] += dR(k+1);
      iter++;
   }
   // Set intersection for return
   ellPt = rg;
   return true;
}



//*************************************************************************************************
// Infamous DUP
//*************************************************************************************************
ossimObject* ossimSarSensorModel::dup() const
{
  return new ossimSarSensorModel(*this);
}
bool ossimSarSensorModel::useForward() const
{
  return false;
}

bool ossimSarSensorModel::autovalidateInverseModelFromGCPs(const double & xtol, const double & ytol, const double azTimeTol, const double & rangeTimeTol) const
{
  bool success = true;

  unsigned int gcpId = 1;
  
  for(std::vector<GCPRecordType>::const_iterator gcpIt = theGCPRecords.begin(); gcpIt!=theGCPRecords.end();++gcpIt,++gcpId)
    {
    ossimDpt estimatedImPt;
    TimeType estimatedAzimuthTime;
    double   estimatedRangeTime;

    bool thisSuccess = true;
    
    // Estimate times
    bool s1 = this->worldToAzimuthRangeTime(gcpIt->worldPt,estimatedAzimuthTime,estimatedRangeTime);
    this->worldToLineSample(gcpIt->worldPt,estimatedImPt);

    if(!s1)
      {
      thisSuccess = false;
      }

    if(std::abs(estimatedImPt.x - gcpIt->imPt.x) > xtol)
      {
      thisSuccess = false;
      }

    if(std::abs(estimatedImPt.y - gcpIt->imPt.y) > ytol)
      {
      thisSuccess = false;
      }

    if(std::abs((estimatedAzimuthTime-gcpIt->azimuthTime).total_microseconds()>azTimeTol))
      {
      thisSuccess = false;
      }
       
    if(std::abs(estimatedRangeTime - gcpIt->slantRangeTime)>rangeTimeTol)
      {
      thisSuccess = false;
      }

    bool verbose = true;

    success = success && thisSuccess;
    
    if(verbose && !thisSuccess)
      {
    
      std::cout<<"GCP #"<<gcpId<<std::endl;
      std::cout<<"Azimuth time: ref="<<gcpIt->azimuthTime<<", predicted: "<<estimatedAzimuthTime<<", res="<<boost::posix_time::to_simple_string(estimatedAzimuthTime-gcpIt->azimuthTime)<<std::endl;
      std::cout<<"Slant range time: ref="<<gcpIt->slantRangeTime<<", predicted: "<<estimatedRangeTime<<", res="<<std::abs(estimatedRangeTime - gcpIt->slantRangeTime)<<std::endl;
      std::cout<<"Image point: ref="<<gcpIt->imPt<<", predicted="<<estimatedImPt<<", res="<<estimatedImPt-gcpIt->imPt<<std::endl;
      std::cout<<std::endl;
      }
    }

  if(success)
    {
    std::cout<<"All GCPs within "<<ytol <<" azimuth pixel, "<<xtol<<" range pixel, "<<azTimeTol<<" s of azimuth time, "<<rangeTimeTol<<" of range time"<<std::endl;

    }
    

  return success;
}


bool ossimSarSensorModel::autovalidateForwardModelFromGCPs(const double& resTol) const
{
  bool success = true;
  bool verbose = true;
  
  unsigned int gcpId = 1;
  
  for(std::vector<GCPRecordType>::const_iterator gcpIt = theGCPRecords.begin(); gcpIt!=theGCPRecords.end();++gcpIt,++gcpId)
    {
    ossimGpt estimatedWorldPt;
    ossimGpt refPt = gcpIt->worldPt;

    double estimatedRangeTime;
    TimeType estimatedAzimuthTime;

    lineSampleToAzimuthRangeTime(gcpIt->imPt,estimatedAzimuthTime,estimatedRangeTime);
    
    lineSampleHeightToWorld(gcpIt->imPt,gcpIt->worldPt.height(),estimatedWorldPt);

    double res = refPt.distanceTo(estimatedWorldPt);

    if(res>resTol)
      {
      success = false;

      if(verbose)
        {
        std::cout<<"GCP #"<<gcpId<<std::endl;
        std::cout<<"Azimuth time: ref="<<gcpIt->azimuthTime<<", predicted: "<<estimatedAzimuthTime<<", res="<<boost::posix_time::to_simple_string(estimatedAzimuthTime-gcpIt->azimuthTime)<<std::endl;
        std::cout<<"Slant range time: ref="<<gcpIt->slantRangeTime<<", predicted: "<<estimatedRangeTime<<", res="<<std::abs(estimatedRangeTime - gcpIt->slantRangeTime)<<std::endl;
        std::cout<<"Im point: "<<gcpIt->imPt<<std::endl;
        std::cout<<"World point: ref="<<refPt<<", predicted="<<estimatedWorldPt<<", res="<<res<<" m"<<std::endl;
        std::cout<<std::endl;
        }
      }
    }

  return success;
}
}
