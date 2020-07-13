#include <per/prOmniXML.h>

prOmniXML::prOmniXML(std::string fic):prCameraModelXML(fic)
{
    
}

xmlNodePtr prOmniXML::putModel(prOmni* cam)
{
    xmlNodePtr node_tmp;
    xmlNodePtr node_model = prCameraModelXML::putModel((prCameraModel*)cam);
    char str[21];
    
    node_tmp = xmlNewComment((xmlChar*)"Xi");
    xmlAddChild(node_model,node_tmp);
    sprintf(str,"%.10f",cam->getXi());
    xmlNewTextChild(node_model,NULL,(xmlChar*)LABEL_XML_XI,(xmlChar*)str);
    
    return node_model;
}

int prOmniXML::getModel(xmlDocPtr& doc, xmlNodePtr node_model,prOmni* cam)
{
    std::string camera_name_tmp = "";
    double u0,v0,au,av,vald,xi;
    int vali;
    CameraModelType model_type;
    char *val_char;
    
    for (node_model = node_model->xmlChildrenNode; node_model != NULL;  node_model = node_model->next)
    {
        if(node_model->type != XML_ELEMENT_NODE) continue;
        
        if(xmlC2S(node_model->name) == LABEL_XML_CAMERA_NAME)
        {
            xmlReadCharChild (doc,node_model, &val_char);
            camera_name_tmp = val_char;
        }
        
        if(xmlC2S(node_model->name) == LABEL_XML_CAMERA_TYPE)
        {
            xmlReadIntChild (doc,node_model, vali);
            model_type = (CameraModelType)vali;
        }
        
        if(xmlC2S(node_model->name) == LABEL_XML_U0)
        {
            xmlReadDoubleChild (doc,node_model, vald);
            u0=vald;
        }
        
        if(xmlC2S(node_model->name) == LABEL_XML_V0)
        {
            xmlReadDoubleChild (doc,node_model, vald);
            v0=vald;
        }
        
        if(xmlC2S(node_model->name) == LABEL_XML_AU)
        {
            xmlReadDoubleChild (doc,node_model, vald);
            au=vald;
        }
        
        if(xmlC2S(node_model->name) == LABEL_XML_AV)
        {
            xmlReadDoubleChild (doc,node_model, vald);
            av=vald;
        }
        
        if(xmlC2S(node_model->name) == LABEL_XML_XI)
        {
            xmlReadDoubleChild (doc,node_model, vald);
            xi=vald;
        }
        
    }
    cam->init(au,av,u0,v0,xi);
    //cam->setName(camera_name_tmp);
    //cam->setType(model_type);
    
    return 0;
}

void prOmniXML::writeModel(prOmni* cam)
{
    xmlAddChild(node,putModel(cam));
}

int prOmniXML::readModel(prOmni* cam)
{
    int nbCamera = 0;
    
    for (node = node->xmlChildrenNode; node != NULL;  node = node->next)
    {
        if(node->type != XML_ELEMENT_NODE) continue;
        if(xmlC2S(node->name) == LABEL_XML_CAMERA)
        {
            nbCamera++;
            getModel(doc,node,cam);
        }
    }
    
    //renvoie la derniËre camera par defaut
    return nbCamera;
}


void prOmniXML::operator<<(prOmni &cam)
{
    xmlOpenToWrite();
    writeModel(&cam);
    xmlWriteToFile();
    xmlEndAccessFile();
}

void prOmniXML::operator>>(prOmni &cam)
{
    xmlOpenToParse();
    readModel(&cam);	
    xmlEndAccessFile();
}

