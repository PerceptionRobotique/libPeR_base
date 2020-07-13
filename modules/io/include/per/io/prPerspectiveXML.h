/*!
 \file prPerspectiveXML.h
 \brief Header file for the prPerspectiveXML class
 \author Guillaume CARON
 \version 0.1
 \date january 2017
 */

#if !defined(_PRPERSPECTIVEXML_H)
#define _PRPERSPECTIVEXML_H

#include <per/prCameraModelXML.h>
#include <per/prPerspective.h>

/*!
 \class prPerspectiveXML prPerspectiveXML.h <per/prPerspectiveXML.h>
 \brief Class for writing perspective camera models to XML files and reading them from XML files.
 */
class prPerspectiveXML : public prCameraModelXML
{
public:
    /*!
     * \fn prPerspectiveXML(std::string fic)
     * \brief Constructor of a prPerspectiveXML object
     * \param fic the XML filename string
     */
    prPerspectiveXML(std::string fic);
    
    /*!
     * \fn void operator<<(prPerspective & cam)
     * \brief << operator override for writing a prPerspective camera model to a XML file
     *
     * \param cam the prPerspective camera model to write
     * \return Nothing.
     */
    void operator<<(prPerspective & cam);
    
    /*!
     * \fn void operator>>(prPerspective & cam)
     * \brief >> operator override for reading a prPerspective camera model from a XML file
     *
     * \param cam the prPerspective camera model that will hold the loaded parameters
     * \return Nothing.
     */
    void operator>>(prPerspective & cam);
};

#endif  //_PRPERSPECTIVEXML_H
