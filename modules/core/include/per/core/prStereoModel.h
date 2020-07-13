/*!
 \file prStereoModel.h
 \brief Header file for the prStereoModel class
 \author Guillaume CARON
 \version 0.1
 \date january 2017
 */

#if !defined(_PRSTEREOMODEL_H)
#define _PRSTEREOMODEL_H

#include <per/prSensorModel.h>

#include <visp/vpHomogeneousMatrix.h>

/*!
 \class prStereoModel prStereoModel.h <per/prStereoModel.h>
 \brief Class for dealing with rigs of sensor.
 */
class prStereoModel : public prSensorModel
{
public:
    prSensorModel **sen; /*!< Array of pointers to prSensorModel*/
    vpHomogeneousMatrix *sjMr; /*!< Array of frame change matrices from each sensor to the rig origin*/
    
    unsigned int nbsens; /*!< Number of sensors*/
    
    /*!
     * \fn prStereoModel(unsigned int _nbsens)
     * \brief Contructor of the prStereoModel initializing the number of sensors of the rig
     * \param _nbsens the number of sensors
     */
    prStereoModel(unsigned int _nbsens);
    
    /*!
     * \fn ~prStereoModel()
     * \brief Destructor of a prStereoModel object
     */
    ~prStereoModel();
    
    
    /*!
     * \fn void init(unsigned int _nbsens)
     * \brief Initializes the number of sensors of the rig
     * \param _nbsens the number of sensors
     * \return Nothing
     */
    void init(unsigned int _nbsens);
    
    /*!
     * \fn void setSensor(unsigned int j, prSensorModel* _sen)
     * \brief Sets the sensor model of sensor j of the rig
     * \param j sensor index
     * \param _sen pointer to the sensor model
     * \return Nothing
     */
    void setSensor(unsigned int j, prSensorModel* _sen);

    /*!
     * \fn void setsjMr(unsigned int j, vpHomogeneousMatrix & M)
     * \brief Sets the sensor pose in the rig, i.e. the rig to the j-th sensor frame change
     * \param j sensor index
     * \param M the frame change matrix
     * \return Nothing
     */
    void setsjMr(unsigned int j, vpHomogeneousMatrix & M);
    
    /*!
     * \fn int get_nbsens()
     * \brief Accessor to the number of sensors of the rig
     * \return the number of sensors of the rig
     */
    int get_nbsens() {return nbsens;}
};

#endif  //_PRSTEREOMODEL_H
