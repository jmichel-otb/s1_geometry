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

#include "ossimTerraSarXSarSensorModel.h"


int main(int argc, char * argv[])
{
  std::cout.precision(9);
  
  if(argc != 2)
    return EXIT_FAILURE;
  
  std::string annotationXml = argv[1];

  ossimplugins::ossimTerraSarXSarSensorModel * sensor = new ossimplugins::ossimTerraSarXSarSensorModel();

  sensor->readAnnotationFile(annotationXml);

  bool validate = sensor->autovalidateInverseModelFromGCPs();

  delete sensor;
  
  return validate;
}
