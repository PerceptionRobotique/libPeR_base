/*!
 \file prSensorModelXML.h
 \brief Header file for the prSensorModelXML class
 \author Guillaume CARON
 \version 0.1
 \date january 2017
 */

#if !defined(_PRSENSORMODELXML_H)
#define _PRSENSORMODELXML_H

#include <per/prXML.h>

/*!
 \class prSensorModelXML prSensorModelXML.h <per/prSensorModelXML.h>
 \brief Base class for writing sensor models to XML files and reading them from XML files.
 */
class prSensorModelXML : public prXML
{
public:
    /*!
     * \fn prSensorModelXML(std::string ficIn)
     * \brief Constructor of a prSensorModelXML object
     * \param ficIn the XML filename string
     */
    prSensorModelXML(std::string fic);
};

#endif  //_PRSENSORMODELXML_H
