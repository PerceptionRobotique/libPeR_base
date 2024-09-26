/*!
 \file prPerspectiveXML.h
 \brief Header file for the prPerspectiveXML class
 \author Guillaume CARON
 \version 0.1
 \date january 2017
 */

#if !defined(_PREQUIRECTANGULARXML_H)
#define _PREQUIRECTANGULARXML_H

#include <per/prcommon.h>
#include <per/prCameraModelXML.h>
#include <per/prEquirectangular.h>

/*!
 \class PER_EXPORT prPerspectiveXML prPerspectiveXML.h <per/prPerspectiveXML.h>
 \brief Class for writing perspective camera models to XML files and reading them from XML files.
 */
class PER_EXPORT prEquirectangularXML : public prCameraModelXML
{
public:
    /*!
     * \fn prEquirectangularXML(std::string fic)
     * \brief Constructor of a prEquirectangularXML object
     * \param fic the XML filename string
     */
    prEquirectangularXML(std::string fic);

    virtual ~prEquirectangularXML() override = default;
    
    /*!
     * \fn void operator<<(prEquirectangular & cam)
     * \brief << operator override for writing a prEquirectangular camera model to a XML file
     *
     * \param cam the prEquirectangular camera model to write
     * \return Nothing.
     */
    void operator<<(prEquirectangular & cam);
    
    /*!
     * \fn void operator>>(prEquirectangular & cam)
     * \brief >> operator override for reading a prEquirectangular camera model from a XML file
     *
     * \param cam the prEquirectangular camera model that will hold the loaded parameters
     * \return Nothing.
     */
    void operator>>(prEquirectangular & cam);
};

#endif  //_PREQUIRECTANGULARXML_H
