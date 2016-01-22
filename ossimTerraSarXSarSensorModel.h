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

#ifndef ossimTerraSarXSarSensorModel_HEADER
#define ossimTerraSarXSarSensorModel_HEADER

#include "ossimSarSensorModel.h"

namespace ossimplugins
{


class ossimTerraSarXSarSensorModel : public ossimSarSensorModel
{
public:

  /** Constructor */
  ossimTerraSarXSarSensorModel();
  
  /** Copy constructor */
  ossimTerraSarXSarSensorModel(const ossimTerraSarXSarSensorModel& m);

  /** Destructor */
  virtual ~ossimTerraSarXSarSensorModel();

  //Not implemented yet
  void readAnnotationFile(const std::string & annotationXml);

protected:
  /*
  std::string theProductType;
  std::string theMode;
  std::string theSwath;
  std::string thePolarisation; 
  */
};

} // end namespace

#endif
