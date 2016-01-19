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

#ifndef ossimSarSensorModel_HEADER
#define ossimSarSensorModel_HEADER

#if defined(__GNUC__) || defined(__clang__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-parameter"
#pragma GCC diagnostic ignored "-Woverloaded-virtual"
#pragma GCC diagnostic ignored "-Wshadow"
#include <ossim/projection/ossimSensorModel.h>
#pragma GCC diagnostic pop
#else
#include <ossimPluginConstants.h>
#include <ossim/projection/ossimSensorModel.h>
#endif

#include "boost/date_time/posix_time/posix_time.hpp"

namespace ossimplugins
{


class ossimSarSensorModel : public ossimSensorModel
{
public:

  typedef boost::posix_time::ptime         TimeType;
  typedef boost::posix_time::time_duration DurationType;

  typedef struct
  {
    TimeType azimuthTime;
    ossimEcefPoint position;
    ossimEcefVector velocity;
  } OrbitRecordType;

  typedef struct
  {
    TimeType azimuthTime;
    double   slantRangeTime;
    ossimDpt imPt;
    ossimGpt worldPt;
  } GCPRecordType;

  typedef struct
  {
    TimeType azimuthStartTime;
    unsigned long startLine;
    TimeType azimuthStopTime;
    unsigned long endLine;
  } BurstRecordType;

  typedef struct
  {
    TimeType azimuthTime;
    double rg0;
    std::vector<double> coefs;
  } CoordinateConversionRecordType;
  
  /** Constructor */
  ossimSarSensorModel();
  
  /** Copy constructor */
  ossimSarSensorModel(const ossimSarSensorModel& m);

  /** Destructor */
  virtual ~ossimSarSensorModel();

  virtual void lineSampleHeightToWorld(const ossimDpt& imPt, const double & heightEllipsoid, ossimGpt& worldPt) const;

  virtual void lineSampleToWorld(const ossimDpt& imPt, ossimGpt& worldPt) const;

  virtual void worldToLineSample(const ossimGpt& worldPt, ossimDpt & imPt) const;
  
  virtual void worldToLineSampleAzimuthRangeTime(const ossimGpt& worldPt, ossimDpt & imPt, double & azimuthTime, double & rangeTime) const;
  
protected:

  virtual void computeRangeDoppler(const ossimEcefPoint & inputPt, double & range, double & doppler) const;

  virtual void interpolateSensorPosVel(const TimeType & azimuthTime, ossimEcefPoint& sensorPos, ossimEcefVector& sensorVel) const;

  virtual void slantRangeToGroundRange(const double & slantRange, const TimeType & azimuthTime, double & groundRange) const;

  std::vector<OrbitRecordType> theOrbitRecords;

  std::vector<GCPRecordType>   theGCPRecords;

  std::vector<BurstRecordType> theBurstRecords;

  std::vector<CoordinateConversionRecordType> theSlantRangeToGroundRangeRecords;

  double theRadarFrequency; // in Hz

  double theAzimuthTimeInterval; // in microseconds

  double theNearRangeTime; // in seconds

  double theRangeSamplingRate; // in Hz

  double theRangeResolution; // in meters

  bool isGRD;

  const double C = 299792458;

};

}

#endif
