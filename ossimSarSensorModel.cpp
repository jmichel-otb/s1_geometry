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
{}
  

ossimSarSensorModel::ossimSarSensorModel(const ossimSarSensorModel& m)
{}

/** Destructor */
ossimSarSensorModel::~ossimSarSensorModel()
{}

void ossimSarSensorModel::lineSampleHeightToWorld(const ossimDpt& imPt, const double & heightEllipsoid, ossimGpt& worldPt) const
{}

void ossimSarSensorModel::lineSampleToWorld(const ossimDpt& imPt, ossimGpt& worldPt) const
{}

void ossimSarSensorModel::worldToLineSample(const ossimGpt& worldPt, ossimDpt & imPt) const
{}

void ossimSarSensorModel::worldToLineSampleAzimuthRangeTime(const ossimGpt& worldPt, ossimDpt & imPt, double & azimuthTime, double & rangeTime) const
{}
  
void ossimSarSensorModel::computeRangeDoppler(const ossimEcefPoint & inputPt, double & range, double & doppler) const
{}

void ossimSarSensorModel::interpolateSensorPosVel(const TimeType & azimuthTime, ossimEcefPoint& sensorPos, ossimEcefVector& sensorVel) const
{}

void ossimSarSensorModel::slantRangeToGroundRange(const double & slantRange, const TimeType & azimuthTime, double & groundRange) const
{}

}
