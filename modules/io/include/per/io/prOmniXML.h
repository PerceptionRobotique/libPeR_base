/*!
 \file prOmniXML.h
 \brief Header file for the prOmniXML class
 \author Guillaume CARON
 \version 0.1
 \date january 2017
 */

#if !defined(_PROMNIXML_H)
#define _PROMNIXML_H

#include <per/prCameraModelXML.h>
#include <per/prOmni.h>

 /*!< XML tag definition */
#define LABEL_XML_XI                                "xi"

/*!
 \class prOmniXML prOmniXML.h <per/prOmniXML.h>
 \brief Class for writing omnidirectional camera models to XML files and reading them from XML files.
 */
class prOmniXML : public prCameraModelXML
{
private:
    /*!
     * \fn int readModel(prOmni* cam)
     * \brief Reads a camera model from the XML tree and saves parameters in a prOmni camera model
     * \param cam the pointer to the last prOmni camera model (must be allocated before calling readModel) of the XML tree
     * \return the number of cameras in the XML tree
     */
    int readModel(prOmni* cam);
    
    /*!
     * \fn int writeModel(prOmni* cam)
     * \brief Writes a prOmni camera model to the XML tree
     * \param cam the pointer to the prOmni camera model to save
     * \return Nothing
     */
    void writeModel(prOmni* cam);
    
public:
    /*!
     * \fn static xmlNodePtr putModel(prOmni* cam)
     * \brief Insert the prOmni camera model parameters in the XML tree
     * \param cam the pointer to the prOmni camera model that holds the parameters to save
     * \return the smart pointer to the created XML node
     */
    static xmlNodePtr putModel(prOmni* cam);
    
    /*!
     * \fn static int getModel(xmlDocPtr &doc, xmlNodePtr node,prOmni* cam)
     * \brief Gets parameters of the prOmni camera model cam from the XML tree node
     * \param doc the XML document smart pointer
     * \param node the XML node at which to read the prOmni camera model parameters
     * \param cam the prOmni camera model that holds the parameters to save
     * \return 0 if runned fine
     */
    static int getModel(xmlDocPtr &doc, xmlNodePtr node,prOmni* cam);
    
    /*!
     * \fn prOmniXML(std::string fic)
     * \brief Constructor of a prOmniXML object
     * \param fic the XML filename string
     */
    prOmniXML(std::string fic);
    
    /*!
     * \fn void operator<<(prOmni & cam)
     * \brief << operator override for writing a prOmni camera model to a XML file
     *
     * \param cam the prOmni camera model to write
     * \return Nothing.
     */
    void operator<<(prOmni & cam);
    
    /*!
     * \fn void operator>>(prOmni & cam)
     * \brief >> operator override for reading a prOmni camera model from a XML file
     *
     * \param cam the prOmni camera model that will hold the loaded parameters
     * \return Nothing.
     */
    void operator>>(prOmni & cam);
};

#endif  //_PROMNIXML_H
