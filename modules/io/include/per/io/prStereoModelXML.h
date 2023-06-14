/*!
 \file prStereoModelXML.h
 \brief Header file for the prStereoModelXML class
 \author Guillaume CARON
 \version 0.1
 \date january 2017
 */

#if !defined(_PRSTEREOMODELXML_H)
#define _PRSTEREOMODELXML_H

#include <per/prcommon.h>
#include <per/prSensorModelXML.h>
#include <per/prStereoModel.h>
#include <per/prCameraModel.h> // a eliminer

 /*!< XML tags definition */
#define LABEL_XML_SYSTEM                         "system"
#define LABEL_XML_NBCAMS                         "NbCams"
#define LABEL_XML_NOCAM                           "NoCam"
#define LABEL_XML_POSE                             "pose"
#define LABEL_XML_TX                                 "tX"
#define LABEL_XML_TY                                 "tY"
#define LABEL_XML_TZ                                 "tZ"
#define LABEL_XML_THETAUX                          "thuX"
#define LABEL_XML_THETAUY                          "thuY"
#define LABEL_XML_THETAUZ                          "thuZ"

/*!
 \class PER_EXPORT prStereoModelXML prStereoModelXML.h <per/prStereoModelXML.h>
 \brief Class for writing stereo models to XML files and reading them from XML files.
 */
class PER_EXPORT prStereoModelXML : public prSensorModelXML
{
protected:
    /*!
     * \fn static xmlNodePtr putSystem(prCameraModel & rig)
     * \brief Insert the prStereoModel sensor rig parameters in the XML tree
     * \param rig the prStereoModel sensor rig that holds the parameters to save
     * \return the smart pointer to the created XML node
     */
    xmlNodePtr putSystem(prStereoModel & rig);
    
    /*!
     * \fn static int getModel(xmlNodePtr node, prStereoModel & rig)
     * \brief Gets parameters of the prStereoModel sensor rig from the XML tree node
     * \param node the XML node at which to read the prStereoModel parameters
     * \param rig the prStereoModel that holds the parameters to save
     * \return 0 if runned fine
     */
    int getSystem(xmlNodePtr node, prStereoModel & rig);
    
    /*!
     * \fn int readModel(prStereoModel & rig)
     * \brief Reads a stereo rig from the XML tree and saves parameters in a prStereoModel camera model
     * \param rig the prStereoModel sensor rig of the XML tree
     * \return  0 if runned fine
     *         -1 if there is not any stereo rig in the XML file
     *         -2 if the number of sensors is not defined in the XML file
     */
    int readSystem(prStereoModel & rig);
    
    /*!
     * \fn int readNbCams(int & nbSensors)
     * \brief Reads the number of cameras in the stereo rig described in the XML file
     * \param nbSensors output variable that will contain the number of sensors in the rig
     * \return  0 if runned fine
     *         -1 if there is not any stereo rig in the XML file
     *         -2 if the number of sensors is not defined in the XML file
     */
    int readNbCams(int & nbSensors);
    
    /*!
     * \fn int readCameraTypes(std::vector<CameraModelType> & v_sensorTypes)
     * \brief Reads the types of sensors in the stereo rig described in the XML file
     * \param v_sensorTypes output variable that will contain the types of sensors in the rig
     * \return  0 if runned fine
     *         -1 if there is not any stereo rig in the XML file
     *         -2 if the number of sensors is not defined in the XML file
     */
    int readCameraTypes(std::vector<CameraModelType> & v_sensorTypes);
    
    /*!
     * \fn void writeSystem(prStereoModel & rig)
     * \brief Writes a prStereoModel sensors rig to the XML tree
     * \param rig the prStereoModel sensors rig to save
     * \return Nothing
     */
    void writeSystem(prStereoModel & rig);
    
public:
    /*!
     * \fn prStereoModelXML(std::string fic)
     * \brief Constructor of a prStereoModelXML object
     * \param fic the XML filename string
     */
    prStereoModelXML(std::string fic);

    virtual ~prStereoModelXML() override = default;

    /*!
     * \fn void operator<<(prStereoModel & rig)
     * \brief << operator override for writing a prStereoModel sensors rig to a XML file
     *
     * \param rig the prStereoModel sensors rig to write
     * \return Nothing.
     */
    void operator<<(prStereoModel & rig);
    
    /*!
     * \fn void operator>>(std::vector<CameraModelType> & v_sensors)
     * \brief >> operator override for reading the types of sensors constituting the rig described in the XML file
     *
     * \param v_sensors the std::vector of sensor types
     * \return Nothing.
     */
    void operator>>(std::vector<CameraModelType> & v_sensors);
    
    /*!
     * \fn void operator>>(int & nbSensors)
     * \brief >> operator override for reading the number of sensors in the rig described in the XML file
     *
     * \param nbSensors the number of sensors
     * \return Nothing.
     */
    void operator>>(int & nbSensors);
    
    /*!
     * \fn void operator>>(prStereoModel & rig)
     * \brief >> operator override for reading a prStereoModel sensors rig from a XML file
     *
     * \param rig the prStereoModel sensors rig that will hold the loaded parameters
     * \return Nothing.
     */
    void operator>>(prStereoModel & rig);
};

#endif  //_PRSTEREOMODELXML_H
