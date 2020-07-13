/*!
 \file prPointFeature.h
 \brief Header file for the prPointFeature class
 \author Guillaume CARON
 \version 0.1
 \date january 2017
 */

#if !defined(_PRPOINTFEATURE_H)
#define _PRPOINTFEATURE_H

#include <per/prFeature.h>

#include <per/prCartesian3DPointVec.h>
#include <per/prCartesian2DPointVec.h>

#include <visp/vpHomogeneousMatrix.h>
#include <visp/vpColVector.h>

/*!
 \class prPointFeature prPointFeature.h <per/prPointFeature.h>
 \brief Class defining a feature point existing in 2D pixels, normalized 2D, intermediate 3D, camera frame 3D, world frame 3D
 */
class prPointFeature : public prFeature {
public:
    prCartesian3DPointVec wX, sX, ix; /*!< world frame 3D, sensor frame 3D, intermediate projection 3D/2D points*/
    prCartesian2DPointVec x, u; /*!< Metric/Normalized sensor surface, acquisition sample 2D points*/
    
    /*!
     * \fn prPointFeature()
     * \brief Default constructor of a prPointFeature object
     */
    prPointFeature();
    
    /*!
     * \fn ~prPointFeature()
     * \brief Default destructor of a prPointFeature object
     */
    ~prPointFeature();
    
    /*!
     * \fn double get_u()
     * \brief Accessor to the horizontal coordinate u of the point in the digital image plane
     * \return u
     */
    double get_u() const {return u.get_x();}
    
    /*!
     * \fn double get_v()
     * \brief Accessor to the vertical coordinate v of the point in the digital image plane
     * \return v
     */
    double get_v() const {return u.get_y();}
    
    
    /*!
     * \fn void set_u()
     * \brief Sets the horizontal coordinate u of the point in the digital image plane
     * \return Nothing
     */
    /*inline*/ void set_u(const double _u) { u.set_x(_u);}
    
    /*!
     * \fn void set_v()
     * \brief Sets the vertical coordinate v of the point in the digital image plane
     * \return Nothing
     */
    /*inline*/ void set_v(const double _v) { u.set_y(_v);}
    
    
    /*!
     * \fn void set_X()
     * \brief Sets the X coordinate of the 3D point expressed in the sensor frame
     * \return Nothing
     */
    /*inline*/ void set_X(const double X_) { sX.set_X(X_) ; }
    
    /*!
     * \fn void set_Y()
     * \brief Sets the Y coordinate of the 3D point expressed in the sensor frame
     * \return Nothing
     */
    /*inline*/ void set_Y(const double Y_) { sX.set_Y(Y_) ; }
    
    /*!
     * \fn void set_Z()
     * \brief Sets the Z coordinate of the 3D point expressed in the sensor frame
     * \return Nothing
     */
    /*inline*/ void set_Z(const double Z_) { sX.set_Z(Z_) ; }
    
    /*!
     * \fn void set_W()
     * \brief Sets the W coordinate of the 3D point expressed in the sensor frame (generally set to 1, its default value for a Cartesian point)
     * \return Nothing
     */
    /*inline*/ void set_W(const double W_) { sX.set_W(W_) ; }
    

    /*!
     * \fn void set_oX()
     * \brief Sets the X coordinate of the 3D point expressed in the world frame
     * \return Nothing
     */
    /*inline*/ void set_oX(const double X_) { wX.set_X(X_) ; }
    
    /*!
     * \fn void set_oY()
     * \brief Sets the Y coordinate of the 3D point expressed in the world frame
     * \return Nothing
     */
    /*inline*/ void set_oY(const double Y_) { wX.set_Y(Y_) ; }
    
    /*!
     * \fn void set_oZ()
     * \brief Sets the Z coordinate of the 3D point expressed in the world frame
     * \return Nothing
     */
    /*inline*/ void set_oZ(const double Z_) { wX.set_Z(Z_) ; }
    
    /*!
     * \fn void set_oW()
     * \brief Sets the W coordinate of the 3D point expressed in the world frame (generally set to 1, its default value for a Cartesian point)
     * \return Nothing
     */
    /*inline*/ void set_oW(const double W_) { wX.set_W(W_) ; }
    
    
    /*!
     * \fn double get_X()
     * \brief Accessor to the X coordinate of the 3D point expressed in the sensor frame
     * \return X
     */
    double get_X()  const { return sX.get_X() ; }
    
    /*!
     * \fn double get_Y()
     * \brief Accessor to the Y coordinate of the 3D point expressed in the sensor frame
     * \return Y
     */
    double get_Y()  const { return sX.get_Y() ; }
    
    /*!
     * \fn double get_Z()
     * \brief Accessor to the Z coordinate of the 3D point expressed in the sensor frame
     * \return Z
     */
    double get_Z() const  { return sX.get_Z() ; }
    
    /*!
     * \fn double get_W()
     * \brief Accessor to the W coordinate of the 3D point expressed in the sensor frame (generally set to 1, its default value for a Cartesian point)
     * \return W
     */
    double get_W()  const { return sX.get_W() ; }
    

    /*!
     * \fn double get_oX()
     * \brief Accessor to the X coordinate of the 3D point expressed in the world frame
     * \return oX
     */
    double get_oX() const { return wX.get_X() ; }
    
    /*!
     * \fn double get_oY()
     * \brief Accessor to the Y coordinate of the 3D point expressed in the world frame
     * \return oY
     */
    double get_oY() const { return wX.get_Y() ; }
    
    /*!
     * \fn double get_oZ()
     * \brief Accessor to the Z coordinate of the 3D point expressed in the world frame
     * \return oZ
     */
    double get_oZ() const { return wX.get_Z() ; }
    
    /*!
     * \fn double get_oW()
     * \brief Accessor to the W coordinate of the 3D point expressed in the sensor frame (generally set to 1, its default value for a Cartesian point)
     * \return oW
     */
    double get_oW() const { return wX.get_W() ; }
    
    
    /*!
     * \fn void set_x(const double x)
     * \brief Sets the x coordinate of the point expressed in the normalized image plane
     * \param x the horizontal coordinate in the normalized image plane
     * \return Nothing
     */
    /*inline*/ void set_x(const double x_) {  x.set_x(x_) ; }
    
    /*!
     * \fn void set_y(const double y)
     * \brief Sets the y coordinate of the point expressed in the normalized image plane
     * \param y the vertical coordinate in the normalized image plane
     * \return Nothing
     */
    /*inline*/ void set_y(const double y_) {  x.set_y(y_) ; }

    /*!
     * \fn void set_w(const double w)
     * \brief Sets the w coordinate of the point expressed in the normalized image plane
     * \param w the homogeneous coordinate in the normalized image plane (generally set to 1, its default value for a Cartesian point)
     * \return Nothing
     */
    /*inline*/ void set_w(const double w_) {  x.set_w(w_) ; }
    

    /*!
     * \fn double get_x()
     * \brief Accessor to the horizontal coordinate x in the normalized image plane
     * \return x
     */
    double get_x()  const { return x.get_x() ; }
    
    /*!
     * \fn double get_y()
     * \brief Accessor to the vertical coordinate y in the normalized image plane
     * \return y
     */
    double get_y()  const { return x.get_y() ; }
    
    /*!
     * \fn double get_w()
     * \brief Accessor to the homogeneous coordinate w in the normalized image planee (generally set to 1, its default value for a Cartesian point)
     * \return x
     */
    double get_w()  const { return x.get_w() ; }
    

    /*!
     * \fn void setWorldCoordinates(const double ox, const double oy, const double oz)
     * \brief Sets the world coordinates of the 3D point
     * \param oX the X coordinate of the 3D point expressed in the world frame
     * \param oY the Y coordinate of the 3D point expressed in the world frame
     * \param oZ the Z coordinate of the 3D point expressed in the world frame
     * \return Nothing
     */
    void setWorldCoordinates(const double oX,
                             const double oY,
                             const double oZ) ;

    /*!
     * \fn int setWorldCoordinates(const prCartesian3DPointVec &_wX)
     * \brief Sets the world coordinates of the 3D point (the 3D point is automatically set so that W = 1)
     * \param _wX the 4-vector containing the X, Y, Z and W coordinates of the 3D point expressed in the world frame
     * \return  0 if _wX is a homogeneous point
     *         -1 if not (i.e. if it is a vector)
     */
    int setWorldCoordinates(const prCartesian3DPointVec &_wX) ;

    /*!
     * \fn void getWorldCoordinates(double & ox, double & oy, double & oz)
     * \brief Gets the world coordinates of the 3D point
     * \param oX the X coordinate of the 3D point expressed in the world frame (output)
     * \param oY the Y coordinate of the 3D point expressed in the world frame (output)
     * \param oZ the Z coordinate of the 3D point expressed in the world frame (output)
     * \return Nothing
     */
    void getWorldCoordinates(double& oX,
                             double& oY,
                             double& oZ) ;

    /*!
     * \fn void getWorldCoordinates(const prCartesian3DPointVec &_wX)
     * \brief Gets the world coordinates of the 3D point
     * \param wX the 4-vector containing the X, Y, Z and W coordinates of the 3D point expressed in the world frame (output)
     * \return Nothing
     */
    void getWorldCoordinates(prCartesian3DPointVec &_xW) ;
    
    /*!
     * \fn prCartesian3DPointVec getWorldCoordinates(void)
     * \brief Gets the world coordinates of the 3D point
     * \return a 4-vector containing the X, Y, Z and W coordinates of the 3D point expressed in the world frame
     */
    prCartesian3DPointVec getWorldCoordinates(void) ;
    
    /*!
     * \fn void changeFrame(const vpHomogeneousMatrix &sMw, prCartesian3DPointVec &_sX)
     * \brief Computes and outputs the coordinates of the 3D point in the camera frame sX from the internal 3D point expressed in the world frame thanks to the cMo frame change
     * \param sMw frame change from the world to the sensor frame
     * \param _sX the output 4-vector containting the 3D point coordinates in the camera frame
     * \return Nothing
     */
    void changeFrame(const vpHomogeneousMatrix &sMw, prCartesian3DPointVec &_sX) ;
    
    /*!
     * \fn void changeFrame(const vpHomogeneousMatrix &sMw)
     * \brief Computes and update the coordinates of the internal 3D point expressed in the camera frame from the internal 3D point expressed in the world frame thanks to the sMw frame change
     * \param sMw frame change from the world to the sensor frame
     * \return Nothing
     */
    /*inline*/ void changeFrame(const vpHomogeneousMatrix &sMw) {
        
        sX = wX.changeFrame(sMw);
        
    }
    
    
    /*!
     * \fn void getImageMetric(double &x, double &y, double &w)const
     * \brief Gets the normalized image plane point coordinates
     * \param x the horizontal coordinate of the 2D point expressed in the normalized image plane (output)
     * \param y the vertical coordinate of the 2D point expressed in the normalized image plane (output)
     * \param w the homogeneous coordinate of the 2D point expressed in the normalized image plane, generally equal to 1 for a Cartesian point (output)
     * \return Nothing
     */
    void getImageMetric(double &x, double &y, double &w)const;
    
    /*!
     * \fn void setImageMetric(const double x,const double y, const double w)
     * \brief Sets the normalized image plane point coordinates
     * \param x the horizontal coordinate of the 2D point expressed in the normalized image plane
     * \param y the vertical coordinate of the 2D point expressed in the normalized image plane
     * \param w the homogeneous coordinate of the 2D point expressed in the normalized image plane, generally equal to 1 (Cartesian point)
     * \return Nothing
     */
    void setImageMetric(const double x,const double y, const double w);
    
    /*!
     * \fn void getPixUV(double &u, double &v)
     * \brief Gets the digital image plane point coordinates
     * \param u the horizontal coordinate of the 2D point expressed in the digital image plane (output)
     * \param v the vertical coordinate of the 2D point expressed in the digital image plane (output)
     * \return Nothing
     */
    void getPixUV(double &u, double &v)const;
    
    /*!
     * \fn void setPixUV(const double u, const double v)
     * \brief Sets the digital image plane point coordinates
     * \param u the horizontal coordinate of the 2D point expressed in the digital image plane
     * \param v the vertical coordinate of the 2D point expressed in the digital image plane
     * \return Nothing
     */
    void setPixUV(const double u, const double v);
    
    /*!
     * \fn void setObjectPixUV(const double oX, const double oY, const double oZ, const double u, const double v)
     * \brief Sets the 3D point coordinates expressed in the object frame and the digital image plane 2D point coordinates
     * \param oX the X coordinate of the 3D point expressed in the world frame
     * \param oY the Y coordinate of the 3D point expressed in the world frame
     * \param oZ the Z coordinate of the 3D point expressed in the world frame
     * \param u the horizontal coordinate of the 2D point expressed in the digital image plane
     * \param v the vertical coordinate of the 2D point expressed in the digital image plane
     * \return Nothing
     */
    void setObjectPixUV(const double oX, const double oY, const double oZ, const double u, const double v);
    
    /*!
     * \fn prPointFeature& operator=(const prPointFeature&)
     * \brief = operator overload for copying a prPointFeature
     *
     * \param Point the source prPointFeature
     * \return the affected prPointFeature.
     */
    prPointFeature& operator=(const prPointFeature& Point);
    
    double *toDouble(unsigned int place = 0);
};

#endif  //_PRPOINTFEATURE_H
