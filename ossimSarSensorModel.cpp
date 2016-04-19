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

namespace {// Anonymous namespace
    const bool k_verbose = false; // global verbose constant; TODO: use an option
}// Anonymous namespace

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
    isGRD(false),
    theAzimuthTimeOffset(0),
    theRangeTimeOffset(0)
{}

void ossimSarSensorModel::lineSampleHeightToWorld(const ossimDpt& imPt, const double & heightAboveEllipsoid, ossimGpt& worldPt) const
{
  assert(!theGCPRecords.empty()&&"theGCPRecords is empty.");

  // Find the closest GCP
  double distance = (imPt-theGCPRecords.front().imPt).length();

  std::vector<GCPRecordType>::const_iterator refGcp = theGCPRecords.begin();

  for(std::vector<GCPRecordType>::const_iterator gcpIt = theGCPRecords.begin();
      gcpIt!=theGCPRecords.end();++gcpIt)
    {
    const double currentDistance = (imPt-gcpIt->imPt).length();

    if(currentDistance < distance)
      {
      distance = currentDistance;
      refGcp = gcpIt;
      }
    }

  // Set the height reference
  ossim_float64 hgtSet;
  if ( ossim::isnan(heightAboveEllipsoid) )
    {
    hgtSet = refGcp->worldPt.height();
    }
  else
    {
    hgtSet = heightAboveEllipsoid;
    }
  const ossimHgtRef hgtRef(AT_HGT, hgtSet);

  ossimEcefPoint ellPt;

  // Simple iterative inversion of inverse model starting at closest gcp
  projToSurface(*refGcp,imPt,&hgtRef,ellPt);

  worldPt = ossimGpt(ellPt);
}

void ossimSarSensorModel::lineSampleToWorld(const ossimDpt& imPt, ossimGpt& worldPt) const
{
  assert(!theGCPRecords.empty()&&"theGCPRecords is empty.");

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

  const ossimGpt refPt = refGcp->worldPt;

  const ossimHgtRef hgtRef(AT_DEM);

  ossimEcefPoint ellPt;

  // Simple iterative inversion of inverse model starting at closest gcp
  projToSurface(*refGcp,imPt,&hgtRef,ellPt);

  worldPt = ossimGpt(ellPt);
}

void ossimSarSensorModel::worldToLineSample(const ossimGpt& worldPt, ossimDpt & imPt) const
{
  assert(theRangeResolution>0&&"theRangeResolution is null.");

  // First compute azimuth and range time
  TimeType azimuthTime;
  double rangeTime;

  const bool success = worldToAzimuthRangeTime(worldPt, azimuthTime, rangeTime);

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

  const bool success = zeroDopplerLookup(inputPt,azimuthTime,interpSensorPos,interpSensorVel);

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
  const double rangeDistance = (interpSensorPos-inputPt).magnitude();
  rangeTime = theRangeTimeOffset + 2*rangeDistance/C;

  return true;
}

void ossimSarSensorModel::lineSampleToAzimuthRangeTime(const ossimDpt & imPt, TimeType & azimuthTime, double & rangeTime) const
{
  // First compute azimuth time here
  lineToAzimuthTime(imPt.y,azimuthTime);

  // Then compute range time
  if(isGRD)
    {
    // Handle grd case here
    double slantRange;
    groundRangeToSlantRange(imPt.x*theRangeResolution,azimuthTime, slantRange);
    rangeTime = theRangeTimeOffset + 2*slantRange/C;
    }
  else
    {
    rangeTime = theRangeTimeOffset + theNearRangeTime + imPt.x*(1/theRangeSamplingRate);
    }
}

void ossimSarSensorModel::computeRangeDoppler(const ossimEcefPoint & inputPt, const ossimEcefPoint & sensorPos, const ossimEcefVector sensorVel, double & range, double & doppler) const
{
  assert(theRadarFrequency>0&&"theRadarFrequency is null");

  // eq. 19, p. 25
  const ossimEcefVector s2gVec = inputPt - sensorPos;

  range = s2gVec.magnitude();

  const double coef = -2*C/(theRadarFrequency*range);

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
        const double td1 = (azimuthTime - theOrbitRecords[j].azimuthTime).total_microseconds();
        const double td2 = (theOrbitRecords[i].azimuthTime - theOrbitRecords[j].azimuthTime).total_microseconds();
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
        const double interp = (azimuthTime-(previousRecord->azimuthTime)).total_microseconds()/static_cast<double>((nextRecord->azimuthTime-previousRecord->azimuthTime).total_microseconds());

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
    const double sr_minus_sr0 =  in-srgrRecord.rg0;

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

  double doppler2(0.);

  // Compute range and doppler of first record
  // NOTE: here we only use the scalar product with vel and discard
  // the constant coef as it has no impact on doppler sign

  double doppler1 = (inputPt-it->position).dot(it->velocity);

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

      const bool dopplerSign2 = doppler2 <0;

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
  const double interpDenom = std::abs(doppler1)+std::abs(doppler2);

  assert(interpDenom>0&&"Both doppler frequency are null in interpolation weight computation");

  const double interp = std::abs(doppler1)/interpDenom;

  // Note that microsecond precision is used here
  const DurationType delta_td = record2->azimuthTime - record1->azimuthTime;
  const double deltat = static_cast<double>(delta_td.total_microseconds());

  // Compute interpolated time offset wrt record1
  const DurationType td = boost::posix_time::microseconds(static_cast<unsigned long>(floor(interp * deltat+0.5)));
  const DurationType offset = boost::posix_time::microseconds(static_cast<unsigned long>(floor(theAzimuthTimeOffset+0.5)));

  // Compute interpolated azimuth time
  interpAzimuthTime = record1->azimuthTime + td + offset;

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

  const DurationType timeSinceStart = azimuthTime - currentBurst->azimuthStartTime;
  const double timeSinceStartInMicroSeconds = static_cast<double>(timeSinceStart.total_microseconds());

  // Eq 22 p 27
  line = (timeSinceStartInMicroSeconds/theAzimuthTimeInterval) + currentBurst->startLine;
}

void ossimSarSensorModel::lineToAzimuthTime(const double & line, TimeType & azimuthTime) const
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

    if(!burstFound)
      {
      if(line < theBurstRecords.front().startLine)
        {
        currentBurst = theBurstRecords.begin();
        }
      else if (line >= theBurstRecords.back().endLine)
        {
        currentBurst = theBurstRecords.end()-1;
        }
      }

    }

    const double timeSinceStartInMicroSeconds = (line - currentBurst->startLine)*theAzimuthTimeInterval;

    const DurationType timeSinceStart = boost::posix_time::microseconds(timeSinceStartInMicroSeconds);
    const DurationType offset = boost::posix_time::microseconds(static_cast<unsigned long>(floor(theAzimuthTimeOffset+0.5)));
  // Eq 22 p 27
  azimuthTime = currentBurst->azimuthStartTime + timeSinceStart + offset;
}



bool ossimSarSensorModel::projToSurface(const GCPRecordType & initGcp, const ossimDpt & target, const ossimHgtRef * hgtRef, ossimEcefPoint & ellPt) const
{
  // Initialize current estimation
  ossimEcefPoint currentEstimation(initGcp.worldPt);

  // Compute corresponding image position
  ossimDpt currentImPoint(initGcp.imPt);

  ossim_float64 currentImResidual = (target-currentImPoint).length();
  double currentHeightResidual = initGcp.worldPt.height() - hgtRef->getRefHeight(initGcp.worldPt);

  bool init = true;

  unsigned int iter = 0;

  NEWMAT::SymmetricMatrix BtB(3);
  NEWMAT::ColumnVector BtF(3);
  NEWMAT::ColumnVector F(3);
  NEWMAT::ColumnVector dR(3);

  // Stop condition: img residual < 1e-2 pixels, height residual <
  // 0.01 m, nb iter < 50. init ensure that loop runs at least once.
  while((init || (currentImResidual > 0.01 || std::abs(currentHeightResidual) > 0.01))  && iter < 50)
    {
    if(init)
      init =false;

    // std::cout<<"Iter: "<<iter<<", Res: im="<<currentImResidual<<", hgt="<<currentHeightResidual<<'\n';

    // compute residuals
    F(1) = target.x - currentImPoint.x;
    F(2) = target.y - currentImPoint.y;
    F(3) = currentHeightResidual;

    // std::cout<<"F("<<iter<<")="<<F<<'\n';

    // Delta use for partial derivatives estimation (in meters)
    const double d = 10.;

    // Compute partial derivatives
    ossimEcefVector p_fx, p_fy, p_fh,dx(d,0,0),dy(0,d,0),dz(0,0,d);
    ossimDpt tmpImPt;

    ossim_float64 rdx,rdy,rdz, fdx,fdy,fdz;

    ossimGpt currentEstimationWorld(currentEstimation);
    ossimGpt tmpGpt = ossimGpt(currentEstimation+dx);
    worldToLineSample(tmpGpt,tmpImPt);
    p_fx[0] = (currentImPoint.x-tmpImPt.x)/d;
    p_fy[0] = (currentImPoint.y-tmpImPt.y)/d;
    p_fh[0] = (currentEstimationWorld.height()-tmpGpt.height())/d;

    tmpGpt = ossimGpt(currentEstimation+dy);
    worldToLineSample(tmpGpt,tmpImPt);
    p_fx[1] = (currentImPoint.x-tmpImPt.x)/d;
    p_fy[1] = (currentImPoint.y-tmpImPt.y)/d;
    p_fh[1] = (currentEstimationWorld.height()-tmpGpt.height())/d;

    tmpGpt = ossimGpt(currentEstimation+dz);
    worldToLineSample(tmpGpt,tmpImPt);
    p_fx[2] = (currentImPoint.x-tmpImPt.x)/d;
    p_fy[2] = (currentImPoint.y-tmpImPt.y)/d;
    p_fh[2] = (currentEstimationWorld.height()-tmpGpt.height())/d;

    // Form B-matrix
    NEWMAT::Matrix B = ossimMatrix3x3::create(p_fx[0], p_fx[1], p_fx[2],
                                              p_fy[0], p_fy[1], p_fy[2],
                                              p_fh[0], p_fh[1], p_fh[2]);

    // std::cout<<"B: "<<B<<'\n';

    // Invert system
    dR = B.i() * F;

    // Update estimate
    for (ossim_int32 k=0; k<3; k++)
      {
      currentEstimation[k] -= dR(k+1);
      }

    // std::cout<<"dR: "<<dR<<'\n';

    currentEstimationWorld=ossimGpt(currentEstimation);

    // Update residuals
    const ossim_float64 atHgt = hgtRef->getRefHeight(currentEstimationWorld);
    currentHeightResidual = atHgt - currentEstimationWorld.height();

    ossimDpt newImPoint;
    worldToLineSample(currentEstimationWorld,currentImPoint);

    // std::cout<<currentImPoint<<'\n';

    currentImResidual = (currentImPoint-target).length();

    ++iter;
    }

  // std::cout<<"Iter: "<<iter<<", Res: im="<<currentImResidual<<", hgt="<<currentHeightResidual<<'\n';

  ellPt = currentEstimation;
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
  if(theGCPRecords.empty())
    {
      return false;
    }

  bool success = true;

  unsigned int gcpId = 1;

  for(std::vector<GCPRecordType>::const_iterator gcpIt = theGCPRecords.begin(); gcpIt!=theGCPRecords.end();++gcpIt,++gcpId)
  {
      ossimDpt estimatedImPt;
      TimeType estimatedAzimuthTime;
      double   estimatedRangeTime;

      // Estimate times
      const bool s1 = this->worldToAzimuthRangeTime(gcpIt->worldPt,estimatedAzimuthTime,estimatedRangeTime);
      this->worldToLineSample(gcpIt->worldPt,estimatedImPt);

      const bool thisSuccess
          =  s1
          && (std::abs(estimatedImPt.x - gcpIt->imPt.x) <= xtol)
          && (std::abs(estimatedImPt.y - gcpIt->imPt.y) <= ytol)
          && (std::abs((estimatedAzimuthTime-gcpIt->azimuthTime).total_microseconds() <= azTimeTol))
          && (std::abs(estimatedRangeTime - gcpIt->slantRangeTime) <= rangeTimeTol)
          ;

      const bool verbose = k_verbose;

      success = success && thisSuccess;

      if(verbose)
      {
          std::cout<<"GCP #"<<gcpId<<'\n';
          std::cout<<"Azimuth time: ref="<<gcpIt->azimuthTime<<", predicted: "<<estimatedAzimuthTime<<", res="<<boost::posix_time::to_simple_string(estimatedAzimuthTime-gcpIt->azimuthTime)<<'\n';
          std::cout<<"Slant range time: ref="<<gcpIt->slantRangeTime<<", predicted: "<<estimatedRangeTime<<", res="<<std::abs(estimatedRangeTime - gcpIt->slantRangeTime)<<'\n';
          std::cout<<"Image point: ref="<<gcpIt->imPt<<", predicted="<<estimatedImPt<<", res="<<estimatedImPt-gcpIt->imPt<<'\n';
          std::cout<<'\n';
      }
  }

  if(success)
    {
    std::cout<<"All GCPs within "<<ytol <<" azimuth pixel, "<<xtol<<" range pixel, "<<azTimeTol<<" s of azimuth time, "<<rangeTimeTol<<" of range time\n";
    }

  return success;
}


bool ossimSarSensorModel::autovalidateForwardModelFromGCPs(const double& resTol)
{

  // First, split half of the gcps to serve as tests, and remove them
  // temporarily from theGCPRecord.
  std::vector<GCPRecordType> gcpRecordSave, testGcps,refGcps;

  gcpRecordSave = theGCPRecords;

  unsigned int count = 0;

  for(std::vector<GCPRecordType>::const_iterator gcpIt = theGCPRecords.begin(); gcpIt!=theGCPRecords.end();++gcpIt,++count)
    {
    if(count%2 == 0)
      {
      refGcps.push_back(*gcpIt);
      }
    else
      {
      testGcps.push_back(*gcpIt);
      }

    }

  theGCPRecords.swap(refGcps);

  bool success = true;
  const bool verbose = k_verbose;

  unsigned int gcpId = 1;

  for(std::vector<GCPRecordType>::const_iterator gcpIt = testGcps.begin(); gcpIt!=testGcps.end();++gcpIt,++gcpId)
    {
    ossimGpt estimatedWorldPt;
    const ossimGpt refPt = gcpIt->worldPt;

    double estimatedRangeTime;
    TimeType estimatedAzimuthTime;

    lineSampleToAzimuthRangeTime(gcpIt->imPt,estimatedAzimuthTime,estimatedRangeTime);

    lineSampleHeightToWorld(gcpIt->imPt,gcpIt->worldPt.height(),estimatedWorldPt);

    const double res = refPt.distanceTo(estimatedWorldPt);

    if(res>resTol || estimatedWorldPt.hasNans())
      {
      success = false;

      if(verbose)
        {
        std::cout<<"GCP #"<<gcpId<<'\n';
        std::cout<<"Azimuth time: ref="<<gcpIt->azimuthTime<<", predicted: "<<estimatedAzimuthTime<<", res="<<boost::posix_time::to_simple_string(estimatedAzimuthTime-gcpIt->azimuthTime)<<'\n';
        std::cout<<"Slant range time: ref="<<gcpIt->slantRangeTime<<", predicted: "<<estimatedRangeTime<<", res="<<std::abs(estimatedRangeTime - gcpIt->slantRangeTime)<<'\n';
        std::cout<<"Im point: "<<gcpIt->imPt<<'\n';
        std::cout<<"World point: ref="<<refPt<<", predicted="<<estimatedWorldPt<<", res="<<res<<" m\n";
        std::cout<<'\n';
        }
      }
    }

  theGCPRecords = gcpRecordSave;

  return success;
}

void ossimSarSensorModel::optimizeTimeOffsetsFromGcps()
{

  double cumulRangeTime(0), cumulAzimuthTime(0);
  unsigned int count=0;

  // First, fix the azimuth time
  for(std::vector<GCPRecordType>::const_iterator gcpIt = theGCPRecords.begin(); gcpIt!=theGCPRecords.end();++gcpIt)
    {
    ossimDpt estimatedImPt;
    TimeType estimatedAzimuthTime;
    double   estimatedRangeTime;

    bool thisSuccess = true;

    // Estimate times
    const bool s1 = this->worldToAzimuthRangeTime(gcpIt->worldPt,estimatedAzimuthTime,estimatedRangeTime);

    if(s1)
      {
      cumulAzimuthTime+=-(estimatedAzimuthTime-gcpIt->azimuthTime).total_microseconds();
      ++count;
      }
    }

  theAzimuthTimeOffset= cumulAzimuthTime/=count;

  // Then, fix the range time
  count=0;

  for(std::vector<GCPRecordType>::const_iterator gcpIt = theGCPRecords.begin(); gcpIt!=theGCPRecords.end();++gcpIt)
    {
    ossimDpt estimatedImPt;
    TimeType estimatedAzimuthTime;
    double   estimatedRangeTime;

    bool thisSuccess = true;

    // Estimate times
    const bool s1 = this->worldToAzimuthRangeTime(gcpIt->worldPt,estimatedAzimuthTime,estimatedRangeTime);

    if(s1)
      {
      cumulRangeTime+=-estimatedRangeTime+gcpIt->slantRangeTime;
      ++count;
      }
    }

  theRangeTimeOffset = cumulRangeTime/=count;
}
}
