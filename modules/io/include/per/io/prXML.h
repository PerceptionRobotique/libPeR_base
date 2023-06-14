/*!
 \file prXML.h
 \brief Header file for the prXML class
 \author Guillaume CARON
 \version 0.1
 \date january 2017
 */

#if !defined(_PRXML_H)
#define _PRXML_H

 /*!< XML tag definition*/
#define LABEL_XML_ROOT                               "root"

#include <per/prcommon.h>
#include <libxml/xmlmemory.h>
#include <string>

#include <per/prIo.h>

/*!
 \class PER_EXPORT prXML prXML.h <per/prXML.h>
 \brief Base class for writing to XML files and reading from XML files.
 */
class PER_EXPORT prXML : public prIo {
protected:
    xmlNodePtr node; /*!< Smart pointer to a XML tree node*/
    xmlDocPtr doc; /*!< Smart pointer to a XML document*/
    std::string file; /*!< XML filename*/
    
    /*!
     * \fn prXML(std::string ficIn)
     * \brief Constructor of a prXML object
     * \param ficIn the XML filename string
     */
    prXML(std::string ficIn);

    /*!
     * \fn prXML()
     * \brief Default constructor of a prXML object
     */
    prXML();
    
    
    /*!
     * \fn void xmlOpenToParse()
     * \brief Open the XML file for reading
     * \return Nothing
     */
    void xmlOpenToParse();
    
    /*!
     * \fn void xmlOpenToWrite()
     * \brief Open the XML file for writing
     * \return Nothing
     */
    void xmlOpenToWrite();
    
    /*!
     * \fn void xmlEndAccessFile()
     * \brief Close the XML file
     * \return Nothing
     */
    void xmlEndAccessFile();
    
    /*!
     * \fn void setPath(std::string nPath)
     * \brief Sets the XML filename
     * \param nPath the XML filename string
     * \return Nothing
     */
    void setPath(std::string nPath);
    
    /*!
     * \fn std::string getPath()
     * \brief Gets the XML filename
     * \return the XML filename string
     */
    std::string getPath();
    
    /*!
     * \fn static void xmlReadCharChild(xmlDocPtr&, xmlNodePtr &node,char **res)
     * \brief Reads the next child of a XML tree node of type char
     * \param doc the XML document smart pointer
     * \param node the XML node smart pointer
     * \param res a pointer to the array of char that will contain the reading as output
     * \return Nothing
     */
    static void xmlReadCharChild(xmlDocPtr & doc, xmlNodePtr &node,char **res);
    
    /*!
     * \fn static void xmlReadDoubleChild(xmlDocPtr&, xmlNodePtr &node,double &res)
     * \brief Reads the next child of a XML tree node of type double
     * \param doc the XML document smart pointer
     * \param node the XML node smart pointer
     * \param res the double variable that will contain the reading as output
     * \return Nothing
     */
    static void xmlReadDoubleChild(xmlDocPtr&, xmlNodePtr node, double &res);
    
    /*!
     * \fn static void xmlReadIntChild(xmlDocPtr&, xmlNodePtr &node,int &res)
     * \brief Reads the next child of a XML tree node of type int
     * \param doc the XML document smart pointer
     * \param node the XML node smart pointer
     * \param res the int variable that will contain the reading as output
     * \return Nothing
     */
    static void xmlReadIntChild(xmlDocPtr&,xmlNodePtr &node,int &res);
    
    /*!
     * \fn static std::string xmlC2S(const xmlChar *str)
     * \brief Converts a characters string from the xmlChar type to the std::sting type
     * \param str the xmlChar type characters string
     * \return the converted string
     */
    static std::string xmlC2S(const xmlChar *str);
    
public:
    /*!
     * \fn void xmlWriteToFile()
     * \brief Effectively write the XML tree to a file
     * \return the converted string
     */
    void xmlWriteToFile();

    virtual ~prXML() override = default;

};

#endif  //_PRXML_H
