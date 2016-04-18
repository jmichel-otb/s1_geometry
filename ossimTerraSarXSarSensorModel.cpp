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

void ossimplugins::ossimTerraSarXSarSensorModel::readAnnotationFile(const std::string & annotationXml, const std::string & geoXml)
{
    ossimRefPtr<ossimXmlDocument> xmlDoc = new ossimXmlDocument(annotationXml);

    //isGRD parse variant?
    std::string product_type = xmlDoc->getRoot()->findFirstNode("productInfo/productVariantInfo/productVariant")->getText();

    std::cout << "type " <<  product_type << std::endl;

    isGRD = (product_type == "MGD" || product_type == "GEC" || product_type == "EEC");

    // First, lookup position/velocity records
    std::vector<ossimRefPtr<ossimXmlNode> > xnodes;
    xmlDoc->findNodes("/level1Product/platform/orbit/stateVec",xnodes);

    std::cout << "Number of states " << xnodes.size() << std::endl;

    for(std::vector<ossimRefPtr<ossimXmlNode> >::iterator itNode = xnodes.begin(); itNode!=xnodes.end();++itNode)
    {
        OrbitRecordType orbitRecord;

        // Retrieve acquisition time
        ossimString att1 = "timeUTC";
        ossimString s = (*itNode)->findFirstNode(att1)->getText();
        std::replace(s.begin(), s.end(), 'T', ' ');
        orbitRecord.azimuthTime = boost::posix_time::time_from_string(s);

        // Retrieve ECEF position
        att1 = "posX";
        orbitRecord.position[0] = (*itNode)->findFirstNode(att1)->getText().toDouble();
        att1 = "posY";
        orbitRecord.position[1] = (*itNode)->findFirstNode(att1)->getText().toDouble();
        att1 = "posZ";
        orbitRecord.position[2] = (*itNode)->findFirstNode(att1)->getText().toDouble();

        // Retrieve ECEF velocity
        att1 = "velX";
        orbitRecord.velocity[0] = (*itNode)->findFirstNode(att1)->getText().toDouble();
        att1 = "velY";
        orbitRecord.velocity[1] = (*itNode)->findFirstNode(att1)->getText().toDouble();
        att1 = "velZ";
        orbitRecord.velocity[2] = (*itNode)->findFirstNode(att1)->getText().toDouble();

        //Add one orbits record
        std::cout << "Add theOrbitRecords\n";
        theOrbitRecords.push_back(orbitRecord);
    }

    //Parse the near range time (in seconds)
    theNearRangeTime = xmlDoc->getRoot()->findFirstNode("productInfo/sceneInfo/rangeTime/firstPixel")->getText().toDouble();

    std::cout << "theNearRangeTime " << theNearRangeTime << std::endl;

    //Parse the range sampling rate
    theRangeSamplingRate = xmlDoc->getRoot()->findFirstNode("instrument/settings/RSF")->getText().toDouble();

    std::cout << "theRangeSamplingRate " << theRangeSamplingRate << std::endl;

    //Parse the range resolution
    theRangeResolution = xmlDoc->getRoot()->findFirstNode("productSpecific/complexImageInfo/slantRangeResolution")->getText().toDouble();

    std::cout << "theRangeResolution " << theRangeResolution << std::endl;

    //Parse the radar frequency
    theRadarFrequency = xmlDoc->getRoot()->findFirstNode("instrument/settings/settingRecord/PRF")->getText().toDouble();

    std::cout << "theRadarFrequency " << theRadarFrequency << std::endl;

    //Manage only strip map product for now (one burst)

    //Parse azimuth time start/stop
    ossimString s1 = xmlDoc->getRoot()->findFirstNode("productInfo/sceneInfo/start/timeUTC")->getText();
    std::replace(s1.begin(), s1.end(), 'T', ' ');
    TimeType azimuthTimeStart = boost::posix_time::time_from_string(s1);

    std::cout << "azimuthTimeStart " << azimuthTimeStart << std::endl;

    ossimString s2 = xmlDoc->getRoot()->findFirstNode("productInfo/sceneInfo/stop/timeUTC")->getText();
    std::replace(s2.begin(), s2.end(), 'T', ' ');
    TimeType azimuthTimeStop = boost::posix_time::time_from_string(s2);

    std::cout << "azimuthTimeStop " << azimuthTimeStop << std::endl;

    const double td = (azimuthTimeStop - azimuthTimeStart).total_microseconds();

    // numberOfRows
    unsigned int numberOfRows = xmlDoc->getRoot()->findFirstNode("productInfo/imageDataInfo/imageRaster/numberOfRows")->getText().toUInt16();

    std::cout << "numberOfRows " << numberOfRows << std::endl;

    //Compute azimuth time interval
    theAzimuthTimeInterval = td / static_cast<double> (numberOfRows);

    std::cout << "theAzimuthTimeInterval " << theAzimuthTimeInterval  << " and 1/prf: " << (1 / theRadarFrequency) * 1000000 << std::endl;

    //For Terrasar-X only 1 burst is supported for now
    BurstRecordType burstRecord;

    burstRecord.startLine = 0;
    burstRecord.azimuthStartTime = azimuthTimeStart;
    burstRecord.azimuthStopTime = azimuthTimeStop;

    burstRecord.endLine = numberOfRows - 1;

    theBurstRecords.push_back(burstRecord);

    //GRD (detected product)
    if(isGRD)
    {
        //Retrieve Slant Range to Ground range coeddifcients
        CoordinateConversionRecordType coordRecord;

        //Get azimuth time start (again)
        coordRecord.azimuthTime = azimuthTimeStart;

        //Set ground range origin to 0 (FIXME?)
        coordRecord.rg0 = 0.;

        //Read  coefficients
        xnodes.clear();

        const unsigned int polynomialDegree = xmlDoc->getRoot()->findFirstNode("productSpecific/projectedImageInfo/slantToGroundRangeProjection/polynomialDegree")->getText().toUInt16();

        std::cout << "Number of coefficients " << polynomialDegree << std::endl;

        ossimString path = "/level1Product/productSpecific/projectedImageInfo/slantToGroundRangeProjection/coefficient";
        const ossimString EXP = "exponent";
        ossimString s;

        xmlDoc->findNodes(path, xnodes);

        if ( xnodes.size() )
        {
            for (unsigned int i = 0; i < xnodes.size(); ++i)
            {
                if (xnodes[i].valid())
                {
                    xnodes[i]->getAttributeValue(s, EXP);
                    coordRecord.coefs.push_back(xnodes[i]->getText().toDouble());
                    std::cout << "Coef number " << i << " value: " << xnodes[i]->getText().toDouble() << std::endl;
                }
            }
        }
        assert(!coordRecord.coefs.empty()&&"The srgr record has empty coefs vector.");

        theSlantRangeToGroundRangeRecords.push_back(coordRecord);
    }

    //Parse GCPs
    ossimRefPtr<ossimXmlDocument> xmlGeo = new ossimXmlDocument(geoXml);

    xnodes.clear();
    xmlGeo->findNodes("/geoReference/geolocationGrid/gridPoint",xnodes);

    std::cout<<"Found "<<xnodes.size()<<" GCPs\n";

    for(std::vector<ossimRefPtr<ossimXmlNode> >::iterator itNode = xnodes.begin(); itNode!=xnodes.end();++itNode)
    {
        GCPRecordType gcpRecord;

        // Get delta acquisition time
        ossimString att1 = "t";
        const double deltaAzimuth = (*itNode)->findFirstNode(att1)->getText().toDouble();
        gcpRecord.azimuthTime = azimuthTimeStart + boost::posix_time::microseconds(deltaAzimuth * 1000000);

        //Get delta range time
        att1 = "tau";
        gcpRecord.slantRangeTime = theNearRangeTime + (*itNode)->findFirstNode(att1)->getText().toDouble();

        att1 = "col";
        gcpRecord.imPt.x = (*itNode)->findFirstNode(att1)->getText().toDouble() - 1.;

        att1 = "row";
        gcpRecord.imPt.y = (*itNode)->findFirstNode(att1)->getText().toDouble() - 1.;

        ossimGpt geoPoint;
        att1 = "lat";
        gcpRecord.worldPt.lat = (*itNode)->findFirstNode(att1)->getText().toDouble();
        att1 = "lon";
        gcpRecord.worldPt.lon = (*itNode)->findFirstNode(att1)->getText().toDouble();
        att1 = "height";
        gcpRecord.worldPt.hgt = (*itNode)->findFirstNode(att1)->getText().toDouble();

        theGCPRecords.push_back(gcpRecord);
    }

    this->optimizeTimeOffsetsFromGcps();
}
