/*!
 \file prCameraModel.h
 \brief Definition of the base class of camera models
 \author Guillaume CARON
 \version 0.1
 \date january 2017
 */

#if !defined(_PRCAMERAMODEL_H)
#define _PRCAMERAMODEL_H

#include <per/prSensorModel.h>
#include <per/prPointFeature.h>

#include <string>
#include <visp/vpMatrix.h>

/*!
 \enum Enum merging all the considered camera models
 */
typedef enum
{
    Omni, /*!< Unified central camera projection */
    Persp, /*!< classical perspective camera projection */
    Fisheye, /*!< Unified central camera projection */
    Paraboloid, /*!< Adhoc paracatadioptric omnidirectional camera projection */
    Equirectangular, /*!< Equirectangular spherical camera projection */
    Ortho /*!< Orthographic camera projection */
} CameraModelType;


/*!
 \class prCameraModel prCameraModel.h <per/prCameraModel.h>
 \brief Class specifying the prSensorModel one by defining the common attributes and actions between every camera models.
 */
class prCameraModel : public prSensorModel {
public:
    double au,av,u0,v0,inv_av,inv_au; /*!< Minimal set of intrinsic parameters (the ones of K_3x3) */
    
    CameraModelType type; /*!< Type of the camera model among Omni, Persp, Fisheye, ... */
    std::string name; /*!< Name of the camera model */
    bool distorsions; /*!< Are distorsions set? */

    double k[5]; /*!< Distorsion parameters :     
                  * k[0-2] : radial distorsions
                  * k[3-4] : tangential distorsions
                  */
    bool activek[5]; /*!< Which distorsion parameter is considered by the camera model? */

    double ik[5]; /*!< Undistorsion parameters :
                   * k[0-2] : radial undistorsions
                   * k[3-4] : tangential undistorsions
                   */

    int nbActiveParametersBase; /*!< Number of active intrinsic parameters excluding distorsions */
    int nbActiveParameters; /*!< Total number of active intrinsic parameters (including distorsions) */


    /*!
     * \fn prCameraModel(double au=0,double av=0,double u0=0,double v0=0, double k1=0,double k2=0,double k3=0,double k4=0,double k5=0)
     * \brief Contructor of the generic prCameraModel initializing intrinsic parameters common to every camera models deriving prCameraModel
     * \param au the horizontal scale factor from the normalized image plane to the digital image (\alpha_u)
     * \param av the vertical scale factor from the normalized image plane to the digital image (\alpha_v)
     * \param u0 the horizontal coordinate of the principal point (u_0)
     * \param v0 the vertical coordinate of the principal point (v_0)
     * \param k1 the second order radial distorsion parameter (k_1)
     * \param k2 the fourth order radial distorsion parameter (k_2)
     * \param k3 the sixth order radial distorsion parameter (k_3)
     * \param k4 the first tangential distorsion parameter (k_4)
     * \param k5 the second tangential distorsion parameter (k_5)
     */
    prCameraModel(double au=0,double av=0,double u0=0,double v0=0, double k1=0,double k2=0,double k3=0,double k4=0,double k5=0);
    
    /*!
     * \fn ~prCameraModel()
     * \brief Destructor of a prCameraModel object
     */
    ~prCameraModel();
    
    
    /*!
     * \fn void init(double au=0,double av=0,double u0=0,double v0=0, double k1=0,double k2=0,double k3=0,double k4=0,double k5=0)
     * \brief Initialize the intrinsic parameters
     * \param au the horizontal scale factor from the normalized image plane to the digital image (\alpha_u)
     * \param av the vertical scale factor from the normalized image plane to the digital image (\alpha_v)
     * \param u0 the horizontal coordinate of the principal point (u_0)
     * \param v0 the vertical coordinate of the principal point (v_0)
     * \param k1 the second order radial distorsion parameter (k_1)
     * \param k2 the fourth order radial distorsion parameter (k_2)
     * \param k3 the sixth order radial distorsion parameter (k_3)
     * \param k4 the first tangential distorsion parameter (k_4)
     * \param k5 the second tangential distorsion parameter (k_5)
     * \return Nothing
     */
    void init(double au=0,double av=0,double u0=0,double v0=0, double k1=0,double k2=0,double k3=0,double k4=0,double k5=0);

    /*!
     * \fn void init(const prCameraModel& c)
     * \brief Initialize the intrinsic parameters from another prCameraModel object
     * \param c the input prCameraModel object from which intrinsic parameters are going to be copied
     * \return Nothing
     */
    void init(const prCameraModel& c);

    
    /*!
     * \fn void meterPixelConversion(prPointFeature & P)
     * \brief Convert a 2D point from the normalized image plane to the digital image plane
     * \param P the point to convert
     * \return Nothing
     */
    void meterPixelConversion(prPointFeature & P);
    
    /*!
     * \fn void pixelMeterConversion(prPointFeature & P)
     * \brief Convert a 2D point from the digital image plane to the normalized image plane
     * \param P the point to convert
     * \return Nothing
     */
    void pixelMeterConversion(prPointFeature & P);
    
    /*!
     * \fn virtual void project3DImage(prPoint & P)=0
     * \brief Project a 3D point expressed in the sensor frame to the normalized image plane (must be implemented in a specified sensor model)
     * \param P the point to project
     * \return Nothing
     */
    virtual void project3DImage(prPointFeature & P)=0;
    
    /*!
     * \fn CameraModelType getType()
     * \brief Accessor to the type of the current sensor model
     * \return The sensor model type
     */
    /*inline*/ CameraModelType getType(){return type;}

    /*!
     * \fn std::string getName() const
     * \brief Accessor to the name of the current sensor model
     * \return The sensor type name
     */
    /*inline*/ std::string getName() const { return name; }

    /*!
     * \fn double getau() const
     * \brief Accessor to the au (\alpha_u) parameter of the current sensor model
     * \return au
     */
    /*inline*/ double getau() const { return au; }
    
    /*!
     * \fn double getav() const
     * \brief Accessor to the av (\alpha_v) parameter of the current sensor model
     * \return av
     */
    /*inline*/ double getav() const { return av; }
    
    /*!
     * \fn double getu0() const
     * \brief Accessor to the u0 (u_0) parameter of the current sensor model
     * \return u0
     */
    /*inline*/ double getu0() const { return u0; }
    
    /*!
     * \fn double getv0() const
     * \brief Accessor to the v0 (v_0) parameter of the current sensor model
     * \return v0
     */
    /*inline*/ double getv0() const { return v0; }
    
    
    /*!
     * \fn double getk1() const
     * \brief Accessor to the k1 (k_1) radial distorsion parameter of the current sensor model
     * \return k1
     */
    /*inline*/ double getk1() const { return k[0]; }
    
    /*!
     * \fn double getk2() const
     * \brief Accessor to the k2 (k_2) radial distorsion parameter of the current sensor model
     * \return k2
     */
    /*inline*/ double getk2() const { return k[1]; }
    
    /*!
     * \fn double getk3() const
     * \brief Accessor to the k3 (k_3) radial distorsion parameter of the current sensor model
     * \return k3
     */
    /*inline*/ double getk3() const { return k[2]; }
    
    /*!
     * \fn double getk4() const
     * \brief Accessor to the k4 (k_4) tangential distorsion parameter of the current sensor model
     * \return k4
     */
    /*inline*/ double getk4() const { return k[3]; }
    
    /*!
     * \fn double getk5() const
     * \brief Accessor to the k5 (k_5) tangential distorsion parameter of the current sensor model
     * \return k5
     */
    /*inline*/ double getk5() const { return k[4]; }
    
    
    /*!
     * \fn double getik1() const
     * \brief Accessor to the ik1 (k_1') radial undistorsion parameter of the current sensor model
     * \return ik1
     */
    /*inline*/ double getik1() const { return ik[0]; }
    
    /*!
     * \fn double getik2() const
     * \brief Accessor to the ik2 (k_2') radial undistorsion parameter of the current sensor model
     * \return ik2
     */
    /*inline*/ double getik2() const { return ik[1]; }
    
    /*!
     * \fn double getik3() const
     * \brief Accessor to the ik3 (k_3') radial undistorsion parameter of the current sensor model
     * \return ik3
     */
    /*inline*/ double getik3() const { return ik[2]; }
    
    /*!
     * \fn double getik4() const
     * \brief Accessor to the ik4 (k_4') tangential undistorsion parameter of the current sensor model
     * \return ik4
     */
    /*inline*/ double getik4() const { return ik[3]; }
    
    /*!
     * \fn double getik5() const
     * \brief Accessor to the ik5 (k_5') tangential undistorsion parameter of the current sensor model
     * \return ik5
     */
    /*inline*/ double getik5() const { return ik[4]; }
    
    
    /*!
     * \fn int getNbActiveParameters() const
     * \brief Gets the total number of considered instrinsic parameters (including distorsion parameters)
     * \return The number of active parameters
     */
    /*inline*/ int getNbActiveParameters() const { return nbActiveParameters; }

    /*!
     * \fn int getNbActiveDistorsionParameters() const
     * \brief Gets the number of considered distorsion parameters
     * \return The number of active distorsion parameters
     */
    /*inline*/ int getNbActiveDistorsionParameters() const { return nbActiveParameters-nbActiveParametersBase; }

    
    /*!
     * \fn vpMatrix getK() const
     * \brief Gets the K_3x3 affine transformation matrix containing intrinsic parameters au, av, u0 and v0
     * \return K_3x3
     */
    vpMatrix getK() const;
    
    
    /*!
     * \fn void setDistorsions(bool d)
     * \brief Activates or cancels the consideration of distorsion parameters
     * \param d true: distorsions are considered
     *          false : distorsions are not considered
     * \return Nothing
     */
    /*inline*/ void setDistorsions(bool d){ distorsions=d;};
    
    /*!
     * \fn void setType(CameraModelType m)
     * \brief Sets the sensor model type
     * \param m sensor model type
     * \return Nothing
     */
    /*inline*/ void setType(CameraModelType m){ type=m;};
    
    /*!
     * \fn void setName(std::string n)
     * \brief Sets the sensor model type name
     * \param n sensor model type name
     * \return Nothing
     */
    /*inline*/ void setName(std::string n) { name=n; };

    
    /*!
     * \fn void setPixelRatio(double au, double av)
     * \brief Sets the scaling factors au and av (from the normalized image plane to the digital one)
     * \param au the horizontal scale factor from the normalized image plane to the digital image (\alpha_u)
     * \param av the vertical scale factor from the normalized image plane to the digital image (\alpha_v)
     * \return Nothing
     */
    void setPixelRatio(double au, double av);
    
    /*!
     * \fn void setPrincipalPoint(double u0, double v0)
     * \brief Sets the principal point coordinates u0 and v0
     * \param u0 the horizontal coordinate of the principal point (u_0)
     * \param v0 the vertical coordinate of the principal point (v_0)
     * \return Nothing
     */
    void setPrincipalPoint(double u0, double v0);
    
    /*!
     * \fn void setDistorsionParameters(double k1=0,double k2=0,double k3=0,double k4=0,double k5=0)
     * \brief Sets the five distorsion parameters (from the normalized image plane to the distorted normalized image plane)
     * \param k1 the second order radial distorsion parameter (k_1)
     * \param k2 the fourth order radial distorsion parameter (k_2)
     * \param k3 the sixth order radial distorsion parameter (k_3)
     * \param k4 the first tangential distorsion parameter (k_4)
     * \param k5 the second tangential distorsion parameter (k_5)
     * \return Nothing
     */
    void setDistorsionParameters(double k1=0,double k2=0,double k3=0,double k4=0,double k5=0);
    
    /*!
     * \fn void setUndistorsionParameters(double ik1=0,double ik2=0,double ik3=0,double ik4=0,double ik5=0)
     * \brief Sets the five undistorsion parameters (from the distorted normalized image plane to the normalized image plane)
     * \param ik1 the second order radial undistorsion parameter (k_1')
     * \param ik2 the fourth order radial undistorsion parameter (k_2')
     * \param ik3 the sixth order radial undistorsion parameter (k_3')
     * \param ik4 the first tangential undistorsion parameter (k_4')
     * \param ik5 the second tangential undistorsion parameter (k_5')
     * \return Nothing
     */
    void setUndistorsionParameters(double ik1=0,double ik2=0,double ik3=0,double ik4=0,double ik5=0);
    
    /*!
     * \fn void setActiveDistorsionParameters(bool k1=true,bool k2=true,bool k3=true,bool k4=true,bool k5=true);
     * \brief Sets which distorsion/undistorsion parameters are considered in the sensor model
     * \param k1 true: k1 and ik1 are considered
     *           false: k1 and ik1 are not considered
     * \param k2 true: k2 and ik2 are considered
     *           false: k2 and ik2 are not considered
     * \param k3 true: k3 and ik3 are considered
     *           false: k3 and ik3 are not considered
     * \param k4 true: k4 and ik4 are considered
     *           false: k4 and ik4 are not considered
     * \param k5 true: k5 and ik5 are considered
     *           false: k5 and ik5 are not considered
     * \return Nothing
     */
    void setActiveDistorsionParameters(bool k1=true,bool k2=true,bool k3=true,bool k4=true,bool k5=true);
    
    
    /*!
     * \fn virtual prCameraModel& operator=(const prCameraModel& cam)
     * \brief = operator overload for copying a prCameraModel
     *
     * \param cam the source prCameraModel
     * \return the affected prCameraModel.
     */
    virtual prCameraModel& operator=(const prCameraModel& cam);
    
    /*!
     * \fn virtual std::ostream& operator << (std::ostream & os)
     * \brief << operator overload for writing the prCameraModel sensor model parameters to a stream
     *
     * \param os the stream in which to write
     * \return the updated stream.
     */
    virtual std::ostream& operator << (std::ostream & os);
    
};

#endif  //_PRCAMERAMODEL_H
