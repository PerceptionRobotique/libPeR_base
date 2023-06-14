/*!
 \file prCameraModelXML.h
 \brief Header file for the prCameraModelXML class
 \author Guillaume CARON
 \version 0.1
 \date january 2017
 */

#if !defined(_PRCAMERAMODELXML_H)
#define _PRCAMERAMODELXML_H

#include <per/prcommon.h>
#include <per/prSensorModelXML.h>
#include <per/prCameraModel.h>

 /*!< XML tags definition*/
#define LABEL_XML_CAMERA                             "camera"
#define LABEL_XML_CAMERA_NAME                        "name"
#define LABEL_XML_WIDTH                              "image_width"
#define LABEL_XML_HEIGHT                             "image_height"
#define LABEL_XML_MODEL                              "model"
#define LABEL_XML_CAMERA_TYPE                         "type"
#define LABEL_XML_CAMERA_OMNI                         "omni"
#define LABEL_XML_CAMERA_PERSP                        "persp"
#define LABEL_XML_CAMERA_FISHEYE                      "fisheye"
#define LABEL_XML_CAMERA_PARABOLOID                   "paraboloid"
#define LABEL_XML_CAMERA_EQUIRECT                   "equirectangular"
#define LABEL_XML_CAMERA_ORTHO                   "orthographic"
#define LABEL_XML_U0                                 "u0"
#define LABEL_XML_V0                                 "v0"
#define LABEL_XML_AU                                 "au"
#define LABEL_XML_AV                                 "av"
#define LABEL_XML_K0                                "k0"
#define LABEL_XML_K1                                "k1"
#define LABEL_XML_K2                                "k2"
#define LABEL_XML_K3                                "k3"
#define LABEL_XML_K4                                "k4"
#define LABEL_XML_IK0                                "ik0"
#define LABEL_XML_IK1                                "ik1"
#define LABEL_XML_IK2                                "ik2"
#define LABEL_XML_IK3                                "ik3"
#define LABEL_XML_IK4                                "ik4"

/*!
 \class PER_EXPORT prCameraModelXML prCameraModelXML.h <per/prCameraModel.XML.h>
 \brief Base class for writing camera models to XML files and reading them from XML files.
 */
class PER_EXPORT prCameraModelXML : public prSensorModelXML
{
protected:
    /*!
     * \fn int readModel(prCameraModel* cam)
     * \brief Reads a camera model from the XML tree and saves parameters in a prX camera model
     * \param cam the pointer to the last prX camera model (must be allocated before calling readModel) of the XML tree
     * \return the number of cameras in the XML tree
     */
    int readModel(prCameraModel* cam);
    /*!
     * \fn void writeModel(prCameraModel* cam)
     * \brief Writes a prX camera model to the XML tree
     * \param cam the pointer to the prX camera model to save
     * \return Nothing
     */
    void writeModel(prCameraModel* cam);
    
public:
    /*!
     * \fn static xmlNodePtr putModel(prCameraModel* cam)
     * \brief Insert the prX camera model parameters in the XML tree
     * \param cam the pointer to the prX camera model that holds the parameters to save
     * \return the smart pointer to the created XML node
     */
    static xmlNodePtr putModel(prCameraModel*);
    
    /*!
     * \fn static int getModel(xmlDocPtr &doc, xmlNodePtr node,prCameraModel* cam)
     * \brief Gets parameters of the prX camera model cam from the XML tree node
     * \param doc the XML document smart pointer
     * \param node the XML node at which to read the prX camera model parameters
     * \param cam the prX camera model that holds the parameters to save
     * \return 0 if runned fine
     */
    static int getModel(xmlDocPtr &doc, xmlNodePtr node,prCameraModel* cam);

    /*!
     * \fn static int getModelType(xmlDocPtr &doc, xmlNodePtr node_model)
     * \brief Gets the camera model type to know which prX camera model to allocate
     * \param doc the XML document smart pointer
     * \param node the XML node that holds the parameters to read
     * \return -1 if runned bad
     *         the camera type enum int otherwise
     */
    static int getModelType(xmlDocPtr& doc, xmlNodePtr node_model);
    
    /*!
     * \fn prCameraModelXML(std::string fic)
     * \brief Constructor of a prCameraModelXML object
     * \param ficIn the XML filename string
     */
    prCameraModelXML(std::string fic);

    virtual ~prCameraModelXML() override = default;
    
    
    /*!
     * \fn void operator<<(prCameraModel* cam)
     * \brief << operator override for writing a prX camera model to a XML file
     *
     * \param cam the prX camera model to write
     * \return Nothing.
     */
    void operator<<(prCameraModel* cam);
    
    /*!
     * \fn void operator>>(prCameraModel* cam)
     * \brief >> operator override for reading a prX camera model from a XML file
     *
     * \param cam the prX camera model that will hold the loaded parameters
     * \return Nothing.
     */
    void operator>>(prCameraModel* cam);
};

#endif  //_PRCAMERAMODELXML_H
