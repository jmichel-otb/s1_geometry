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

#include "ossimSentinel1SarSensorModel.h"


int main(int argc, char * argv[])
{
  std::cout.precision(9);
  
  if(argc != 2)
    return EXIT_FAILURE;
  
  std::string annotationXml = argv[1];

  ossimplugins::ossimSentinel1SarSensorModel * sensor = new ossimplugins::ossimSentinel1SarSensorModel();

  sensor->readAnnotationFile(annotationXml);

  delete sensor;
  return EXIT_SUCCESS;
}
