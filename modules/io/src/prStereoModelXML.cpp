#include <per/prStereoModelXML.h>

//a eliminer
#include <per/prOmniXML.h>
#include <per/prPerspectiveXML.h>

prStereoModelXML::prStereoModelXML(std::string fic):prSensorModelXML(fic)
{
    
}

void prStereoModelXML::operator<<(prStereoModel & system)
{
    xmlOpenToWrite();
    writeSystem(system);
    xmlWriteToFile();
    xmlEndAccessFile();
}

void prStereoModelXML::operator>>(prStereoModel & system)
{
    xmlOpenToParse();
    readSystem(system);
    xmlEndAccessFile();
}

void prStereoModelXML::operator>>(std::vector<CameraModelType> & cameras)
{
    xmlOpenToParse();
    readCameraTypes(cameras);
    xmlEndAccessFile();
}

void prStereoModelXML::operator>>(int & nbCams)
{
    xmlOpenToParse();
    readNbCams(nbCams);
    xmlEndAccessFile();
}

xmlNodePtr prStereoModelXML::putSystem(prStereoModel & system)
{
   	xmlNodePtr node_system;
    
    xmlNodePtr node_camera;
   	xmlNodePtr node_pose;
    xmlNodePtr node_tmp;
    char str[21];
    
    //System - top parent
    node_system = xmlNewNode(NULL,(xmlChar*)LABEL_XML_SYSTEM);
    
    //number of cameras in the system
    node_tmp = xmlNewComment((xmlChar*)"Number of cameras in the system");
    xmlAddChild(node_system,node_tmp);
    sprintf(str,"%d",system.get_nbsens());
    xmlNewTextChild(node_system,NULL,(xmlChar*)LABEL_XML_NBCAMS,(xmlChar*)str);
    
    for(int icam = 0 ; icam < system.get_nbsens() ; icam++)
    {
        switch (((prCameraModel *)system.sen[icam])->getType())
        {
            case Omni :
            case Fisheye :
                node_camera = prOmniXML::putModel((prOmni *)(system.sen[icam]));
                break;
            case Persp :
                node_camera = prPerspectiveXML::putModel((prPerspective *)(system.sen[icam]));
                break;
            case Paraboloid :
                //                    node_camera = CParaboloidXml::putModel((CParaboloid *)(system.cam[icam]));
                break;
            default:
                std::cout << "prStereoModelXML::putSystem : Type de camera inconnu" << std::endl;
                continue;
                break;
        }
        //Camera number
        node_tmp = xmlNewComment((xmlChar*)"Camera number on the system");
        xmlAddChild(node_camera,node_tmp);
        sprintf(str,"%d",icam);
        xmlNewTextChild(node_camera,NULL,(xmlChar*)LABEL_XML_NOCAM,(xmlChar*)str);
        
        
        node_tmp = xmlNewComment((xmlChar*)"Relative pose to the system origin");
        xmlAddChild(node_camera,node_tmp);
        
        node_pose = xmlNewNode(NULL,(xmlChar*)LABEL_XML_POSE);
        
        vpPoseVector r(system.sjMr[icam]);
        //Pose tX
        sprintf(str,"%.10f",r[0]);
        xmlNewTextChild(node_pose,NULL,(xmlChar*)LABEL_XML_TX,(xmlChar*)str);
        //Pose tY
        sprintf(str,"%.10f",r[1]);
        xmlNewTextChild(node_pose,NULL,(xmlChar*)LABEL_XML_TY,(xmlChar*)str);
        //Pose tZ
        sprintf(str,"%.10f",r[2]);
        xmlNewTextChild(node_pose,NULL,(xmlChar*)LABEL_XML_TZ,(xmlChar*)str);
        //Pose theta uX
        sprintf(str,"%.10f",r[3]);
        xmlNewTextChild(node_pose,NULL,(xmlChar*)LABEL_XML_THETAUX,(xmlChar*)str);
        //Pose theta uY
        sprintf(str,"%.10f",r[4]);
        xmlNewTextChild(node_pose,NULL,(xmlChar*)LABEL_XML_THETAUY,(xmlChar*)str);
        //Pose theta uZ
        sprintf(str,"%.10f",r[5]);
        xmlNewTextChild(node_pose,NULL,(xmlChar*)LABEL_XML_THETAUZ,(xmlChar*)str);
        
        xmlAddChild(node_camera, node_pose);
        
        xmlAddChild(node_system, node_camera);
    }
    
    return node_system;
}

int prStereoModelXML::getSystem(xmlNodePtr node_system,prStereoModel & system)
{
    /*std::string camera_name_tmp = "";
     double u0,v0,au,av,vald;
     int vali;
     ModelType model_type;
     char *val_char;
     
     for (node_model = node_model->xmlChildrenNode; node_model != NULL;  node_model = node_model->next)
     {
     if(node_model->type != XML_ELEMENT_NODE) continue;
     
     if(xmlC2S(node_model->name) == LABEL_XML_CAMERA_NAME)
     {
     xmlReadCharChild (node_model, &val_char);
     camera_name_tmp = val_char;
     }
     
     if(xmlC2S(node_model->name) == LABEL_XML_CAMERA_TYPE)
     {
     xmlReadIntChild (node_model, vali);
     model_type = (ModelType)vali;
     }
     
     if(xmlC2S(node_model->name) == LABEL_XML_U0)
     {
     xmlReadDoubleChild (node_model, vald);
     u0=vald;
     }
     
     if(xmlC2S(node_model->name) == LABEL_XML_V0)
     {
     xmlReadDoubleChild (node_model, vald);
     v0=vald;
     }
     
     if(xmlC2S(node_model->name) == LABEL_XML_AU)
     {
     xmlReadDoubleChild (node_model, vald);
     au=vald;
     }
     
     if(xmlC2S(node_model->name) == LABEL_XML_AV)
     {
     xmlReadDoubleChild (node_model, vald);
     av=vald;
     }
     
     }
     cam->init(au,av,u0,v0);
     cam->setName(camera_name_tmp);
     cam->setType(model_type);*/
    
    return 0;
}


void prStereoModelXML::writeSystem(prStereoModel & system)
{
    xmlAddChild(node,putSystem(system));
}

int prStereoModelXML::readNbCams(int &nbCams)
{
    int vali;
    
    xmlNodePtr node_system;
    
    xmlNodePtr node_tmp;
    
    node_system = node->xmlChildrenNode;
    while (node_system->type != XML_ELEMENT_NODE)
        node_system = node_system->next;
    
    if(xmlC2S(node_system->name) != LABEL_XML_SYSTEM)
        return -1;
    
    //The XML file has a LABEL_XML_SYSTEM node
    node_tmp = node_system->xmlChildrenNode;
    while (node_tmp->type != XML_ELEMENT_NODE)
        node_tmp = node_tmp->next;
    
    if(xmlC2S(node_tmp->name) != LABEL_XML_NBCAMS)
        return -2;
    xmlReadIntChild (doc, node_tmp, vali);
    nbCams = (CameraModelType)vali;
    
    return 0;
}

int prStereoModelXML::readCameraTypes(std::vector<CameraModelType> &cameras)
{
    int vali, nbCams, icam;
    CameraModelType type_model;
    
    xmlNodePtr node_system;
    xmlNodePtr node_tmp;
    xmlNodePtr node_camera;
    
    node_system = node->xmlChildrenNode;
    while (node_system->type != XML_ELEMENT_NODE)
        node_system = node_system->next;
    
    if(xmlC2S(node_system->name) != LABEL_XML_SYSTEM)
        return -1;
    
    //The XML file has a LABEL_XML_SYSTEM node
    node_tmp = node_system->xmlChildrenNode;
    while (node_tmp->type != XML_ELEMENT_NODE)
        node_tmp = node_tmp->next;
    
    if(xmlC2S(node_tmp->name) != LABEL_XML_NBCAMS)
        return -2;
    xmlReadIntChild (doc, node_tmp, vali);
    nbCams = (CameraModelType)vali;
    
    cameras.reserve(nbCams);
    
    for (node_tmp = node_system->xmlChildrenNode; node_tmp != NULL;  node_tmp = node_tmp->next)
    {
        if(node_tmp->type != XML_ELEMENT_NODE) continue;
        
        if(xmlC2S(node_tmp->name) == LABEL_XML_CAMERA)
        {
            for (node_camera = node_tmp->xmlChildrenNode; node_camera != NULL;  node_camera = node_camera->next)
            {
                if(node_camera->type != XML_ELEMENT_NODE) continue;
                
                if(xmlC2S(node_camera->name) == LABEL_XML_CAMERA_TYPE)
                {
                    xmlReadIntChild (doc, node_camera, vali);
                    type_model = (CameraModelType)vali;
                }
                
                if(xmlC2S(node_camera->name) == LABEL_XML_NOCAM)
                {
                    xmlReadIntChild (doc, node_camera, vali);
                    icam = (CameraModelType)vali;
                }
            }
            cameras[icam] = type_model;
        }
    }
    
    return 0;
}

int prStereoModelXML::readSystem(prStereoModel & system)
{
    int vali, nbCams, icam;
    double vald;
    CameraModelType model_type;
    
    xmlNodePtr node_system;
    xmlNodePtr node_tmp;
    xmlNodePtr node_camera;
    xmlNodePtr node_pose;
    
    node_system = node->xmlChildrenNode;
    while (node_system->type != XML_ELEMENT_NODE)
        node_system = node_system->next;
    
    if(xmlC2S(node_system->name) != LABEL_XML_SYSTEM)
        return -1;
    
    //The XML file has a LABEL_XML_SYSTEM node
    node_tmp = node_system->xmlChildrenNode;
    while (node_tmp->type != XML_ELEMENT_NODE)
        node_tmp = node_tmp->next;
    
    if(xmlC2S(node_tmp->name) != LABEL_XML_NBCAMS)
        return -2;
    xmlReadIntChild (doc, node_tmp, vali);
    nbCams = (CameraModelType)vali;
    
    for (node_tmp = node_system->xmlChildrenNode; node_tmp != NULL;  node_tmp = node_tmp->next)
    {
        if(node_tmp->type != XML_ELEMENT_NODE) continue;
        
        if(xmlC2S(node_tmp->name) == LABEL_XML_CAMERA)
        {
            node_camera = node_tmp->xmlChildrenNode;
            while ( (node_camera != NULL) && (xmlC2S(node_camera->name) != LABEL_XML_NOCAM) )
                node_camera = node_camera->next;
            
            xmlReadIntChild (doc, node_camera, vali);
            icam = vali;
            
            model_type = (CameraModelType)prCameraModelXML::getModelType(doc, node_tmp);
            
            switch(model_type)
            {
                case Omni :
                case Fisheye :
                    system.sen[icam] = new prOmni();
                    prOmniXML::getModel(doc, node_tmp, (prOmni *)(system.sen[icam]));
                    break;
                case Persp :
                    system.sen[icam] = new prPerspective();
                    prPerspectiveXML::getModel(doc, node_tmp, (prPerspective *)(system.sen[icam]));
                    break;
                case Paraboloid :
                    //                    system.cam[icam] = new CParaboloid();
                    //                    node_camera = CParaboloidXml::putModel((CParaboloid *)(system.cam[icam]));
                    break;
                default:
                    std::cout << "prStereoModelXML::putSystem : Type de camera inconnu" << std::endl;
                    continue;
                    break;
            }
            
            node_camera = node_tmp->xmlChildrenNode;
            while ( (node_camera != NULL) && (xmlC2S(node_camera->name) != LABEL_XML_POSE) )
                node_camera = node_camera->next;
            
            vpPoseVector r;
            for (node_pose = node_camera->xmlChildrenNode; node_pose != NULL;  node_pose = node_pose->next)
            {
                if(node_pose->type != XML_ELEMENT_NODE) continue;
                
                if(xmlC2S(node_pose->name) == LABEL_XML_TX)
                {
                    xmlReadDoubleChild (doc, node_pose, vald);
                    r[0] = vald;
                }
                
                if(xmlC2S(node_pose->name) == LABEL_XML_TY)
                {
                    xmlReadDoubleChild (doc, node_pose, vald);
                    r[1] = vald;
                }
                
                if(xmlC2S(node_pose->name) == LABEL_XML_TZ)
                {
                    xmlReadDoubleChild (doc, node_pose, vald);
                    r[2] = vald;
                }
                
                if(xmlC2S(node_pose->name) == LABEL_XML_THETAUX)
                {
                    xmlReadDoubleChild (doc, node_pose, vald);
                    r[3] = vald;
                }
                
                if(xmlC2S(node_pose->name) == LABEL_XML_THETAUY)
                {
                    xmlReadDoubleChild (doc, node_pose, vald);
                    r[4] = vald;
                }
                
                if(xmlC2S(node_pose->name) == LABEL_XML_THETAUZ)
                {
                    xmlReadDoubleChild (doc, node_pose, vald);
                    r[5] = vald;
                }
            }
            system.sjMr[icam].buildFrom(r);
        }
    }
    
    /*
     for (node = node->xmlChildrenNode; node != NULL;  node = node->next)
     {
     if(node->type != XML_ELEMENT_NODE) continue;
     if(xmlC2S(node->name) == LABEL_XML_CAMERA)
     {
     nbCamera++;
     getModel(node,cam);			
     }
     }
     
     //renvoie la derni\E8re camera par defaut
     */
    return 0;
}

