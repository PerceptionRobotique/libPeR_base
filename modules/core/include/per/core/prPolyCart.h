/*!
 \file prPolyCart.h
 \brief Definition of an omnidirectional camera model (Polynomial Cartesian Model)
 \author Guillaume CARON
 \version 0.1
 \date June 2022
 */

#if !defined(_PRPOLYCART_H)
#define _PRPOLYCART_H

#include <per/prCameraModel.h>

/*!
 \class prPolyCart prPolyCart.h <per/prPolyCart.h>
 \brief Class for dealing with the projection functions of a camera modeled with the polynomial Cartesian model of Scaramuzza.
 */
class PER_EXPORT prPolyCart : public prCameraModel {
public:
    double a[5]; /*!< Polynomial Cartesian model parameters: coefficients of the Taylor expansion of the distorsions */

    /*!
     * \fn prPolyCart(double au=0,double u0=0,double v0=0,double *a=NULL, double k1=0,double k2=0,double k3=0,double k4=0,double k5=0)
     * \brief Contructor of the prPolyCart camera model initializing intrinsic parameters
     * \param au the scale factor for both horizontal and vertical image directions from the normalized image plane to the digital image (\alpha_v = \alpha_u)
     * \param u0 the horizontal coordinate of the principal point (u_0)
     * \param v0 the vertical coordinate of the principal point (v_0)
     * \param a coefficients of the Taylor expansion of the distorsions (a = [a0, a1, a2, a3, a4]) // a1 = 0 always to allow derivation (cf. MCPTAM article)
     * \param k1 the second order radial distorsion parameter (k_1)
     * \param k2 the fourth order radial distorsion parameter (k_2)
     * \param k3 the sixth order radial distorsion parameter (k_3)
     * \param k4 the first tangential distorsion parameter (k_4)
     * \param k5 the second tangential distorsion parameter (k_5)
     */
    prPolyCart(double au=0,double u0=0,double v0=0,double *a=NULL, double k1=0,double k2=0,double k3=0,double k4=0,double k5=0);
    
    /*!
     * \fn ~prPolyCart()
     * \brief Destructor of a prPolyCart object
     */
    virtual ~prPolyCart() override = default;

    
    /*!
     * \fn void init(double au=0,double u0=0,double v0=0,double *a=NULL, double k1=0,double k2=0,double k3=0,double k4=0,double k5=0)
     * \brief Initializing intrinsic parameters
     * \param au the scale factor for both horizontal and vertical image directions from the normalized image plane to the digital image (\alpha_v = \alpha_u)
     
     * \param u0 the horizontal coordinate of the principal point (u_0)
     * \param v0 the vertical coordinate of the principal point (v_0)
     * \param a coefficients of the Taylor expansion of the distorsions (a = [a0, a1, a2, a3, a4])
     * \param k1 the second order radial distorsion parameter (k_1)
     * \param k2 the fourth order radial distorsion parameter (k_2)
     * \param k3 the sixth order radial distorsion parameter (k_3)
     * \param k4 the first tangential distorsion parameter (k_4)
     * \param k5 the second tangential distorsion parameter (k_5)
     * \return Nothing
     */
    void init(double au=0,double u0=0,double v0=0,double *a=NULL, double k1=0,double k2=0,double k3=0,double k4=0,double k5=0);
    
    /*!
     * \fn double geta0()
     * \brief Accessor to the a0 (a_0) parameter of the prPolyCart sensor model
     * \return a0
     */
    double geta0(){return a[0];};
    
    /*!
     * \fn double geta1()
     * \brief Accessor to the a1 (a_1) parameter of the prPolyCart sensor model
     * \return a1
     */
    double geta1(){return a[1];};
    
    /*!
     * \fn double geta2()
     * \brief Accessor to the a2 (a_2) parameter of the prPolyCart sensor model
     * \return a2
     */
    double geta2(){return a[2];};
    
    /*!
     * \fn double geta3()
     * \brief Accessor to the a3 (a_3) parameter of the prPolyCart sensor model
     * \return a3
     */
    double geta3(){return a[3];};
    
    /*!
     * \fn double geta4()
     * \brief Accessor to the a4 (a_4) parameter of the prPolyCart sensor model
     * \return a4
     */
    double geta4(){return a[4];};
    
    /*!
     * \fn void project3DImage(prPointFeature & P)
     * \brief Project a 3D point expressed in the camera frame to the normalized image plane
     * \param P the point to project
     * \return Nothing
     */
    void project3DImage(prPointFeature & P);

    /*!
     * \fn virtual void unProject(prPointFeature & P, double & Depth)=0
     * \brief Compute the 3D point expressed in the sensor frame from the normalized image plane and the depth (must be implemented in a specified sensor model)
     * \param P the point to unproject
     * \param Depth the Depth used to unproject the P (/rho here...)
     * \return true if unprojection could be done
     */
    bool unProject(prPointFeature & P, double & Depth);
    
    /*!
     * \fn void project3DSphere(prPointFeature & P, double & Xs, double & Ys, double & Zs)
     * \brief Project a 3D point expressed in the camera frame to the unit sphere centered at the camera optical center
     * \param P the point to project
     * \param Xs the Cartesian X coordinate of the projected point on the unit sphere
     * \param Ys the Cartesian Y coordinate of the projected point on the unit sphere
     * \param Zs the Cartesian Z coordinate of the projected point on the unit sphere
     * \return Nothing
     */
    void project3DSphere(prPointFeature & P, double & Xs, double & Ys, double & Zs);
    
    /*!
     * \fn virtual void projectImageSphere(prPointFeature & P, double & Xs, double & Ys, double & Zs)
     * \brief Inverse projection from the normalized image plane to the unit sphere centered at the camera optical center
     * \param P the point to backproject
     * \param Xs the Cartesian X coordinate of the projected point on the unit sphere
     * \param Ys the Cartesian Y coordinate of the projected point on the unit sphere
     * \param Zs the Cartesian Z coordinate of the projected point on the unit sphere
     * \return true
     */
    virtual bool projectImageSphere(prPointFeature & P, double & Xs, double & Ys, double & Zs);
    
    /*!
     * \fn template <T> static virtual void projectImageSphere(T & x, T & y, T & *a, T & Xs, T & Ys, T & Zs)
     * \brief Inverse projection from the normalized image plane to the unit sphere centered at the camera optical center
     * \param x the Cartesian x coordinate of the source point in the normalized image plane
     * \param y the Cartesian y coordinate of the source point in the normalized image plane
     * \param a coefficients of the Taylor expansion of the distorsions (a = [a0, a1, a2, a3, a4])
     * \param Xs the Cartesian X coordinate of the projected point on the unit sphere
     * \param Ys the Cartesian Y coordinate of the projected point on the unit sphere
     * \param Zs the Cartesian Z coordinate of the projected point on the unit sphere
     * \return Nothing
     */
    template <typename T>
    static void projectImageSphere(T & u, T & v, T *intrinsics, T & Xs, T & Ys, T & Zs)
    {
        T up = u - intrinsics[1];
        T vp = v - intrinsics[2];
        
        T rho_up = sqrt(up*up + vp*vp);
        T rho_up2 = rho_up*rho_up;
    		T rho_up3 = rho_up2*rho_up;
    		T rho_up4 = rho_up3*rho_up;
        
        T r_rho_up = intrinsics[3] + intrinsics[4]*rho_up + intrinsics[5]*rho_up2 + intrinsics[6]*rho_up3 + intrinsics[7]*rho_up4;
        
        Xs = up / intrinsics[0];
        Ys = vp / intrinsics[0];
        Zs = r_rho_up / intrinsics[0];
    }
    
    /*!
     * \fn prPolyCart& operator=(const prPolyCart& cam)
     * \brief = operator overload for copying a prPolyCart camera model
     *
     * \param cam the source prPolyCart
     * \return the affected prPolyCart.
     */
    prPolyCart& operator=(const prPolyCart& cam);

    /*!
     * \fn void computeSensorJacobian(prPointFeature & P, vpMatrix & LuX)
     * \brief computes the Jacobian of the sensor digital image point coordinates with respect to the 3D coordinates in the camera frame
     * \param P the point for which computing the Jacobian
     * \param LuX the output Jacobian Matrix
     * \return Nothing
     */
    void computeSensorJacobian(prPointFeature & P, vpMatrix & LuX);
};

#endif  //_PRPOLYCART_H
