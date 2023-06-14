#include <per/prCameraModelXML.h>

prCameraModelXML::prCameraModelXML(std::string fic):prSensorModelXML(fic)
{
    
}

void prCameraModelXML::operator<<(prCameraModel* cam)
{
    xmlOpenToWrite();
    writeModel(cam);
    xmlWriteToFile();
    xmlEndAccessFile();
}

void prCameraModelXML::operator>>(prCameraModel* cam)
{
    xmlOpenToParse();
    readModel(cam);
    xmlEndAccessFile();
}

xmlNodePtr prCameraModelXML::putModel(prCameraModel* cam)
{
    xmlNodePtr node_camera;
    xmlNodePtr node_tmp;
    char str[21];
    
    //Camera - parent
    node_camera = xmlNewNode(NULL,(xmlChar*)LABEL_XML_CAMERA);
    
    //Camera name
    node_tmp = xmlNewComment((xmlChar*)"Name of the camera");
    xmlAddChild(node_camera,node_tmp);
    xmlNewTextChild(node_camera,NULL,(xmlChar*)LABEL_XML_CAMERA_NAME,(xmlChar*)cam->getName().c_str());
    
    //Camera type
    node_tmp = xmlNewComment((xmlChar*)"Type of the camera");
    xmlAddChild(node_camera,node_tmp);
    sprintf(str,"%d",cam->getType());
    xmlNewTextChild(node_camera,NULL,(xmlChar*)LABEL_XML_CAMERA_TYPE,(xmlChar*)str);
    
    //Camera au av
    node_tmp = xmlNewComment((xmlChar*)"Pixel ratio");
    xmlAddChild(node_camera,node_tmp);
    sprintf(str,"%.10f",cam->getau());
    xmlNewTextChild(node_camera,NULL,(xmlChar*)LABEL_XML_AU,(xmlChar*)str);
    sprintf(str,"%.10f",cam->getav());
    xmlNewTextChild(node_camera,NULL,(xmlChar*)LABEL_XML_AV,(xmlChar*)str);
    
    //Camera u0 v0
    node_tmp = xmlNewComment((xmlChar*)"Principal point");
    xmlAddChild(node_camera,node_tmp);
    sprintf(str,"%.10f",cam->getu0());
    xmlNewTextChild(node_camera,NULL,(xmlChar*)LABEL_XML_U0,(xmlChar*)str);
    sprintf(str,"%.10f",cam->getv0());
    xmlNewTextChild(node_camera,NULL,(xmlChar*)LABEL_XML_V0,(xmlChar*)str);
    
    return node_camera;
}

int prCameraModelXML::getModel(xmlDocPtr &doc, xmlNodePtr node_model,prCameraModel* cam)
{
    std::string camera_name_tmp = "";
    double u0,v0,au,av,vald;
    int vali;
    CameraModelType model_type;
    char *val_char;
    
    for (node_model = node_model->xmlChildrenNode; node_model != NULL;  node_model = node_model->next)
    {
        if(node_model->type != XML_ELEMENT_NODE) continue;
        
        if(xmlC2S(node_model->name) == LABEL_XML_CAMERA_NAME)
        {
            xmlReadCharChild (doc, node_model, &val_char);
            camera_name_tmp = val_char;
        }
        
        if(xmlC2S(node_model->name) == LABEL_XML_CAMERA_TYPE)
        {
            xmlReadIntChild (doc, node_model, vali);
            model_type = (CameraModelType)vali;
        }
        
        if(xmlC2S(node_model->name) == LABEL_XML_U0)
        {
            xmlReadDoubleChild (doc, node_model, vald);
            u0=vald;
        }
        
        if(xmlC2S(node_model->name) == LABEL_XML_V0)
        {
            xmlReadDoubleChild (doc, node_model, vald);
            v0=vald;
        }
        
        if(xmlC2S(node_model->name) == LABEL_XML_AU)
        {
            xmlReadDoubleChild (doc, node_model, vald);
            au=vald;
        }
        
        if(xmlC2S(node_model->name) == LABEL_XML_AV)
        {
            xmlReadDoubleChild (doc, node_model, vald);
            av=vald;
        }
        
    }
    cam->init(au,av,u0,v0);
    //cam->setName(camera_name_tmp);
    //cam->setType(model_type);
    
    return 0;
}


void prCameraModelXML::writeModel(prCameraModel* cam)
{
    xmlAddChild(node,putModel(cam));
}

int prCameraModelXML::readModel(prCameraModel* cam)
{
    int nbCamera = 0;
    
    for (node = node->xmlChildrenNode; node != NULL;  node = node->next)
    {
        if(node->type != XML_ELEMENT_NODE) continue;
        if(xmlC2S(node->name) == LABEL_XML_CAMERA)
        {
            nbCamera++;
            getModel(doc, node,cam);
        }
    }
    
    //renvoie la derniËre camera par defaut
    return nbCamera;
}

int prCameraModelXML::getModelType(xmlDocPtr& doc, xmlNodePtr node_model)
{
    int vali;
    CameraModelType model_type;
    
    for (node_model = node_model->xmlChildrenNode; node_model != NULL;  node_model = node_model->next)
    {
        if(node_model->type != XML_ELEMENT_NODE) continue;
        
        if(xmlC2S(node_model->name) == LABEL_XML_CAMERA_TYPE)
        {
            xmlReadIntChild (doc,node_model, vali);
            return vali;
        }
    }
    
    return -1;
}
