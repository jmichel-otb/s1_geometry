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

#ifndef ossimSentinel1SarSensorModel_HEADER
#define ossimSentinel1SarSensorModel_HEADER

#include "ossimSarSensorModel.h"

namespace ossimplugins
{


class ossimSentinel1SarSensorModel : public ossimSarSensorModel
{
public:

  /** Constructor */
  ossimSentinel1SarSensorModel();
  
  /** Copy constructor */
  ossimSentinel1SarSensorModel(const ossimSentinel1SarSensorModel& m);

  /** Destructor */
  virtual ~ossimSentinel1SarSensorModel();

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
