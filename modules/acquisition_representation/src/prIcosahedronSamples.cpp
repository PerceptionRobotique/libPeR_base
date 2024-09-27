#include <per/prIcosahedronSamples.h>

prIcosahedronSamples::prIcosahedronSamples(std::string coordsPath)
{
    this->coordsPath = coordsPath;
    this->maxSubdivLevels = 7;
}

int prIcosahedronSamples::loadDelaunaySamples(unsigned int subdivLevels)
{
    this->subdivLevels = subdivLevels;
    std::string filename = coordsPath + "/delaunayIcosahedronVertices/delaunayIcosahedron3DVertices_" + std::to_string(subdivLevels) + "subdiv.xml";
    xmlDocPtr doc = xmlParseFile(filename.c_str());
    xmlNodePtr node = xmlDocGetRootElement(doc);

    xmlNodePtr nodeMain, node3DPoint, nodeXYZ;

    double x, y, z;

    nodeMain = node->xmlChildrenNode;
    while (nodeMain->type != XML_ELEMENT_NODE)
        nodeMain = nodeMain->next;

    if (xmlC2S(nodeMain->name) == LABEL_NB_SAMPLES)
    {
        xmlReadIntChild(doc, nodeMain, nbVertices);
    }

    // prCartesian3DPointVec samplesCoordinates[nbVertices];
    samplesCoords = new prCartesian3DPointVec[nbVertices];

    nodeMain = nodeMain->next;
    while (nodeMain->type != XML_ELEMENT_NODE)
        nodeMain = nodeMain->next;

    if (xmlC2S(nodeMain->name) == LABEL_COORDS)
    {
        int i = 0;
        for (node3DPoint = nodeMain->xmlChildrenNode; node3DPoint != NULL; node3DPoint = node3DPoint->next)
        {
            if (node3DPoint->type != XML_ELEMENT_NODE)
                continue;

            if (xmlC2S(node3DPoint->name) == LABEL_POINT)
            {
                for (nodeXYZ = node3DPoint->xmlChildrenNode; nodeXYZ != NULL; nodeXYZ = nodeXYZ->next)
                {

                    if (nodeXYZ->type != XML_ELEMENT_NODE)
                        continue;

                    if (xmlC2S(nodeXYZ->name) == LABEL_X)
                    {
                        xmlReadDoubleChild(doc, nodeXYZ, x);
                    }
                    if (xmlC2S(nodeXYZ->name) == LABEL_Y)
                    {
                        xmlReadDoubleChild(doc, nodeXYZ, y);
                    }
                    if (xmlC2S(nodeXYZ->name) == LABEL_Z)
                    {
                        xmlReadDoubleChild(doc, nodeXYZ, z);
                    }
                }
                if (i <= nbVertices)
                {
                    // samplesCoordinates[i] = prCartesian3DPointVec(x, y, z);
                    samplesCoords[i] = prCartesian3DPointVec(x, y, z);
                    i++;
                }
            }
        }
    }
    xmlFreeDoc(doc);

    return 0;
}

int prIcosahedronSamples::loadVoronoiSamples(unsigned int subdivLevels)
{
    this->subdivLevels = subdivLevels;
    std::string filename = coordsPath + "/voronoiIcosahedronVertices/voronoiIcosahedron3DVertices_" + std::to_string(subdivLevels) + "subdiv.xml";
    xmlDocPtr doc = xmlParseFile(filename.c_str());
    xmlNodePtr node = xmlDocGetRootElement(doc);

    xmlNodePtr nodeFaces, nodeFace, nodeDelaunay3DPoint, node3DPoint, nodeXYZ;

    double xD, yD, zD, xV, yV, zV;

    nodeFaces = node->xmlChildrenNode;
    while (nodeFaces->type != XML_ELEMENT_NODE)
        nodeFaces = nodeFaces->next;

    if (xmlC2S(nodeFaces->name) == LABEL_NB_SAMPLES)
    {
        xmlReadIntChild(doc, nodeFaces, nbVertices);
    }

    // prCartesian3DPointVec samplesCoordinates[nbVertices];
    samplesCoords = new prCartesian3DPointVec[nbVertices];

    if (!samplesVoronoiMap.empty())
    {
        samplesVoronoiMap.clear();
    }

    nodeFaces = nodeFaces->next;
    while (nodeFaces->type != XML_ELEMENT_NODE)
        nodeFaces = nodeFaces->next;

    if (xmlC2S(nodeFaces->name) == LABEL_FACES)
    {
        int i = 0;
        for (nodeFace = nodeFaces->xmlChildrenNode; nodeFace != NULL; nodeFace = nodeFace->next)
        {
            if (nodeFace->type != XML_ELEMENT_NODE)
                continue;

            if (xmlC2S(nodeFace->name) == LABEL_FACE)
            {
                // std::cout << "should be in the face " << i << std::endl;
                for (nodeDelaunay3DPoint = nodeFace->xmlChildrenNode; nodeDelaunay3DPoint != NULL; nodeDelaunay3DPoint = nodeDelaunay3DPoint->next)
                {

                    if (nodeDelaunay3DPoint->type != XML_ELEMENT_NODE)
                        continue;

                    if (xmlC2S(nodeDelaunay3DPoint->name) == LABEL_DELAU_POINT)
                    {
                        // std::cout << "should be in the delaunay coordinates" << std::endl;
                        for (nodeXYZ = nodeDelaunay3DPoint->xmlChildrenNode; nodeXYZ != NULL; nodeXYZ = nodeXYZ->next)
                        {

                            if (nodeXYZ->type != XML_ELEMENT_NODE)
                                continue;

                            if (xmlC2S(nodeXYZ->name) == LABEL_X)
                            {
                                xmlReadDoubleChild(doc, nodeXYZ, xD);
                            }
                            if (xmlC2S(nodeXYZ->name) == LABEL_Y)
                            {
                                xmlReadDoubleChild(doc, nodeXYZ, yD);
                            }
                            if (xmlC2S(nodeXYZ->name) == LABEL_Z)
                            {
                                xmlReadDoubleChild(doc, nodeXYZ, zD);
                            }
                        }
                        if (i <= nbVertices)
                        {
                            // samplesCoordinates[i] = prCartesian3DPointVec(x, y, z);
                            samplesCoords[i] = prCartesian3DPointVec(xD, yD, zD);
                            // i++;
                        }
                    }

                    if (xmlC2S(nodeDelaunay3DPoint->name) == LABEL_COORDS)
                    {
                        // std::cout << "should be in the delaunay coordinates" << std::endl;
                        std::vector<prCartesian3DPointVec> localVerticesVec;
                        for (node3DPoint = nodeDelaunay3DPoint->xmlChildrenNode; node3DPoint != NULL; node3DPoint = node3DPoint->next)
                        {

                            if (node3DPoint->type != XML_ELEMENT_NODE)
                                continue;

                            if (xmlC2S(node3DPoint->name) == LABEL_POINT)
                            {
                                for (nodeXYZ = node3DPoint->xmlChildrenNode; nodeXYZ != NULL; nodeXYZ = nodeXYZ->next)
                                {

                                    if (nodeXYZ->type != XML_ELEMENT_NODE)
                                        continue;

                                    if (xmlC2S(nodeXYZ->name) == LABEL_X)
                                    {
                                        xmlReadDoubleChild(doc, nodeXYZ, xV);
                                    }
                                    if (xmlC2S(nodeXYZ->name) == LABEL_Y)
                                    {
                                        xmlReadDoubleChild(doc, nodeXYZ, yV);
                                    }
                                    if (xmlC2S(nodeXYZ->name) == LABEL_Z)
                                    {
                                        xmlReadDoubleChild(doc, nodeXYZ, zV);
                                    }
                                }
                                localVerticesVec.push_back(prCartesian3DPointVec(xV, yV, zV));
                                // std::cout << xV << ", " << yV << ", " << zV << std::endl;
                            }
                        }
                        if (i < nbVertices)
                        {
                            samplesVoronoiMap.insert({i, localVerticesVec});
                            i++;
                        }
                    }
                }
            }
        }
    }
    xmlFreeDoc(doc);

    return 0;
}

void prIcosahedronSamples::createNeighborMap(double limitDist)
{
    if (!neighMap.empty())
    {
        neighMap.clear();
    }

    for (int i = 0; i < nbVertices; i++)
    {
        std::vector<unsigned int> localNeigh;
        prCartesian3DPointVec XS = samplesCoords[i];
        for (int j = 0; j < nbVertices; j++)
        {
            // if (i != j)
            // {
            prCartesian3DPointVec XSNeigh = samplesCoords[j];
            double dot = XS.get_X() * XSNeigh.get_X() + XS.get_Y() * XSNeigh.get_Y() + XS.get_Z() * XSNeigh.get_Z();
            dot = (dot < -1.) ? -1. : (dot > 1.) ? 1.
                                                 : dot;
            dist = acos(dot);

            if (dist < limitDist)
            {
                localNeigh.push_back(j);
            }
            // }
        }
        neighMap.insert({i, localNeigh});
    }
}