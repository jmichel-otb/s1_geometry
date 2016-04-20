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

#include <boost/config.hpp>

#if defined(__GNUC__) || defined(__clang__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-parameter"
// #pragma GCC diagnostic ignored "-Woverloaded-virtual"
// #pragma GCC diagnostic ignored "-Wshadow"
#include <ossim/projection/ossimSensorModel.h>
#include <ossim/elevation/ossimHgtRef.h>

#pragma GCC diagnostic pop
#else
#include <ossim/projection/ossimSensorModel.h>
#include <ossim/elevation/ossimHgtRef.h>
#endif

#include <boost/date_time/posix_time/posix_time.hpp>

namespace ossimplugins
{

class ossimSarSensorModel : public ossimSensorModel
{
public:

  typedef boost::posix_time::ptime         TimeType;
  typedef boost::posix_time::time_duration DurationType;

  struct OrbitRecordType
  {
    TimeType        azimuthTime;
    ossimEcefPoint  position;
    ossimEcefVector velocity;
  };

  struct GCPRecordType
  {
    TimeType azimuthTime;
    double   slantRangeTime;
    ossimDpt imPt;
    ossimGpt worldPt;
  };

  struct BurstRecordType
  {
    TimeType      azimuthStartTime;
    unsigned long startLine;
    TimeType      azimuthStopTime;
    unsigned long endLine;
  };

  struct CoordinateConversionRecordType
  {
    TimeType            azimuthTime;
    double              rg0;
    std::vector<double> coefs;
  };

  /** Constructor */
  ossimSarSensorModel();

#if ! (defined(BOOST_NO_DEFAULTED_FUNCTIONS) || defined(BOOST_NO_CXX1_DEFAULTED_FUNCTIONS))
  /** Copy constructor */
  ossimSarSensorModel(ossimSarSensorModel const& m) = default;
  /** Move constructor */
  ossimSarSensorModel(ossimSarSensorModel && m) = default;

  /** Destructor */
  virtual ~ossimSarSensorModel() = default;
#endif

  virtual void lineSampleHeightToWorld(const ossimDpt& imPt, const double & heightEllipsoid, ossimGpt& worldPt) const;

  virtual void lineSampleToWorld(const ossimDpt& imPt, ossimGpt& worldPt) const;


  /** This method implement inverse sar geolocation using method found
   *  in ESA document "Guide to ASAR geocoding" (ref
   *  RSL-ASAR-GC-AD). Equation numbers can be found in source code
   *  comments.
   *
   * \param[in] worldPt World point to geocode
   * \param[out] imPt Corresponding estimated image point
   */
  virtual void worldToLineSample(const ossimGpt& worldPt, ossimDpt & imPt) const;

  /**
   * Sub-routine of lineSampleToWorld that computes azimuthTime and
   * slant range time from worldPoint
   *
   * \param[in] worldPoint World point to geocode
   * \param[out] azimuthTime Estimated zero-doppler azimuth time
   * \param[out] rangeTime Estimated range time
   * \return True if sucess, false otherwise. In this case,
   * azimuthTime and rangeTime will not be modified.
   */
  /*virtual*/ bool worldToAzimuthRangeTime(const ossimGpt& worldPt, TimeType & azimuthTime, double & rangeTime) const;

  // TODO: document me
  /*virtual*/ void lineSampleToAzimuthRangeTime(const ossimDpt & imPt, TimeType & azimuthTime, double & rangeTime) const;

  // TODO: document me
  bool autovalidateInverseModelFromGCPs(const double & xtol = 1, const double & ytol = 1, const double azTimeTol = 500, const double &rangeTimeTo=0.0000000001) const;

  // TODO: document me
  bool autovalidateForwardModelFromGCPs(double resTol = 10);

  //Pure virtual in base class
  bool useForward() const;

  void optimizeTimeOffsetsFromGcps();

  /**
   * Returns pointer to a new instance, copy of this.
   */
  virtual ossimObject* dup() const;

  //TODO: Add virtual method readAnnotationFile?

protected:

  /**
   * Compute range and doppler frequency from an input point, sensor
   * position and velocity.
   *
   * \param[in] inputPt The target point
   * \param[in] sensorPos The sensor position
   * \param[in] sensorvel The sensor velocity
   * \param[out] range Estimated range
   * \param[out] doppler Estimated doppler frequency
   */
  virtual void computeRangeDoppler(const ossimEcefPoint & inputPt, const ossimEcefPoint & sensorPos, const ossimEcefVector sensorVel, double & range, double & doppler) const;

  /**
   * Interpolate sensor position and velocity at given azimuth time
   * using lagragian interpolation of orbital records.
   *
   * \param[in] azimuthTime The time at which to interpolate
   * \param[out] sensorPos Interpolated sensor position
   * \param[out] sensorvel Interpolated sensor velocity
   * \param[in] deg Degree of lagragian interpolation
   */
  /*virtual*/ void interpolateSensorPosVel(const TimeType & azimuthTime, ossimEcefPoint& sensorPos, ossimEcefVector& sensorVel, unsigned int deg = 8) const;

  /**
   * Convert slant range to ground range by interpolating slant range
   * to ground range coefficients.
   *
   * \param[in] slantRangeTime The slantRange to convert (meters)
   * \param[in] azimuthTime The corresponding azimuth time
   * \param[out] groundRange The estimated ground range (meters)
   */
  /*virtual*/ void slantRangeToGroundRange(const double & slantRange, const TimeType & azimuthTime, double & groundRange) const;

  // TODO: Document me
  /*virtual*/ void groundRangeToSlantRange(const double & groundRange, const TimeType & azimuthTime, double & slantRange) const;

  // TODO: Document me
  /*virtual*/ void applyCoordinateConversion(const double & in, const TimeType& azimuthTime, const std::vector<CoordinateConversionRecordType> & records, double & out) const;
  /**
   * Estimate the zero-doppler azimuth time and corresponding sensor
   * position and velocity from the inputPt.
   *
   * \param[in] inputPt The point to estimated zero-doppler time on
   * \param[out] interpAzimuthTime Interpolated azimuth time
   * \param[out] interpSensorPos Interpolated sensor position
   * \param[out] interpSensorVel Interpolated sensor velocity
   * \return True if success, false otherwise. In this case, output
   * parameters are left untouched.
   */
  /*virtual*/ bool zeroDopplerLookup(const ossimEcefPoint & inputPt, TimeType & interpAzimuthTime, ossimEcefPoint & interpSensorPos, ossimEcefVector & interpSensorVel) const;

  /**
   * Compute the bistatic correction to apply.
   *
   * \param[in] inputPt The point to compute bistatic correction on
   * \param[in] sensorPos The corresponding sensor position
   * \param[out] bistaticCorrection The estimated bistatic correction
   */
  /*virtual*/ void computeBistaticCorrection(const ossimEcefPoint & inputPt, const ossimEcefPoint & sensorPos, DurationType & bistaticCorrection) const;

  /**
   * Convert azimuth time to fractional line.
   *
   * \param[in] azimuthTime The azimuth time to convert
   * \param[out] The estimated fractional line
   */
  /*virtual*/ // TODO: check why virtual
  void azimuthTimeToLine(const TimeType & azimuthTime, double & line) const;

  // TODO: document me
  /*virtual*/ // TODO: check why virtual
  void lineToAzimuthTime(const double & line, TimeType & azimuthTime) const;

  // TODO: document me
  /*virtual*/ // TODO: check why virtual
  bool projToSurface(const GCPRecordType & initGcp, const ossimDpt & target, const ossimHgtRef & hgtRef, ossimEcefPoint & ellpt) const;

  /**
   * Finds closest GCP.
   *
   * \param[in] imPt  «imPt-explanations»
   * \return the closest GCP record to \c imPt.
   * \throw None
   * \pre `theGCPRecords` shall not be empty.
   */
  GCPRecordType const& findClosestGCP(ossimDpt const& imPt) const;


  std::vector<OrbitRecordType>                theOrbitRecords;
  std::vector<GCPRecordType>                  theGCPRecords;
  std::vector<BurstRecordType>                theBurstRecords;
  std::vector<CoordinateConversionRecordType> theSlantRangeToGroundRangeRecords;
  std::vector<CoordinateConversionRecordType> theGroundRangeToSlantRangeRecords;

  double                                      theRadarFrequency; // in Hz
  double                                      theAzimuthTimeInterval; // in microseconds
  double                                      theNearRangeTime; // in seconds
  double                                      theRangeSamplingRate; // in Hz
  double                                      theRangeResolution; // in meters
  bool                                        theBistaticCorrectionNeeded; // Do we need to compute
                                                                           // bistatic correction ?
  bool                                        isGRD; // True if the product is GRD. False if it is SLC
  double                                      theAzimuthTimeOffset; // Offset in microseconds
  double                                      theRangeTimeOffset; // Offset in seconds;

  static const double C;
private:
  /** Disabled assignment operator.  */
  ossimSarSensorModel& operator=(ossimSarSensorModel const& rhs);
};

}

#endif
