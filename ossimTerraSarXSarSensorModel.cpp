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
{}

/** Destructor */
ossimTerraSarXSarSensorModel::~ossimTerraSarXSarSensorModel()
{}

void ossimTerraSarXSarSensorModel::readAnnotationFile(const std::string & annotationXml)
{
  //Parse specific metadata for TerraSarX
  std::cerr << "Not implemented yet." << std::endl;
}

} // namespace ossimplugins
