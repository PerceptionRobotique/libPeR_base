#include <per/prXML.h>

prXML::prXML()
{
    
}

prXML::prXML(std::string fic)
{
    setPath(fic);
    
}

void prXML::xmlOpenToParse()
{
    doc = xmlParseFile(file.c_str());
    
    node = xmlDocGetRootElement(doc);
    
}

void prXML::xmlOpenToWrite()
{
    doc = xmlReadFile(file.c_str(),NULL,XML_PARSE_NOWARNING + XML_PARSE_NOERROR + XML_PARSE_NOBLANKS);
    doc = xmlNewDoc ((xmlChar*)"1.0");
    node = xmlNewNode(NULL,(xmlChar*)LABEL_XML_ROOT);
    xmlDocSetRootElement(doc,node);
    xmlNodePtr node_tmp = xmlNewComment((xmlChar*)"PR team of the MIS lab - UPJV - France");
    xmlAddChild(node,node_tmp);
    node = xmlDocGetRootElement(doc);
}

void prXML::setPath(std::string nPath)
{
    file = nPath;
}

std::string prXML::getPath()
{
    return file;
}

void prXML::xmlWriteToFile()
{
    xmlSaveFormatFile(file.c_str(),doc,1);
    
}

void prXML::xmlEndAccessFile()
{
    xmlFreeDoc(doc);
}
void prXML::xmlReadCharChild(xmlDocPtr &doc, xmlNodePtr &nodein,char **res)
{
    xmlNodePtr cur;
    
    cur = nodein->xmlChildrenNode;
    *res = (char *) xmlNodeListGetString(doc, cur, 1);
}

void prXML::xmlReadIntChild(xmlDocPtr &doc, xmlNodePtr &nodein,int &res)
{
    char * val_char;
    char * control_convert;
    int val_int;
    xmlNodePtr cur;
    
    cur = nodein ->xmlChildrenNode;
    val_char = (char *) xmlNodeListGetString(doc, cur, 1);
    val_int = strtol ((char *)val_char, &control_convert, 10);
    if (val_char == control_convert) val_int = 0;
    res = val_int;
}

void prXML::xmlReadDoubleChild(xmlDocPtr &doc, xmlNodePtr nodein, double &res)
{
    char * val_char;
    char * control_convert;
    double val_double;
    xmlNodePtr cur;
    
    cur = nodein ->xmlChildrenNode;
    val_char = (char *) xmlNodeListGetString(doc, cur, 1);
    val_double = strtod ((char *)val_char, &control_convert);
    if (val_char == control_convert) val_double = 0;
    res = val_double;
}

std::string prXML::xmlC2S(const xmlChar *str)
{
    std::string tmp=(char*)str;
    return tmp;
}
