/*!
 \file prCartesian3DPointVec.h
 \brief Header file for the prCartesian3DPointVec class
 \author Guillaume CARON
 \version 0.1
 \date february 2017
 */

#if !defined(_PRCARTESIAN3DPOINTVEC_H)
#define _PRCARTESIAN3DPOINTVEC_H

#include <per/pr3DPointVec.h>

/*!
 \class PER_EXPORT prCartesian3DPointVec prCartesian3DPointVec.h <per/prCartesian3DPointVec.h>
 \brief Class defining a 3D in homogeneous Cartesian coordinates
 */
class PER_EXPORT prCartesian3DPointVec : public pr3DPointVec
{
public:
    /*!
     * \fn prCartesian3DPointVec(const double X = 0, const double Y = 0, const double Z = 0, const double W = 1)
     * \brief Constructor setting the 3D coordinates of a point (W = 1) or a vector (W = 0)
     * \param X horizontal coordinate
     * \param Y vertical coordinate
     * \param Z depth coordinate
     * \param W homogeneous coordinate
     */
    prCartesian3DPointVec(const double X = 0, const double Y = 0, const double Z = 0, const double W = 1);

    virtual ~prCartesian3DPointVec() override = default;

    /*!
     * \fn void set_X()
     * \brief Sets the X coordinate of the 3D point expressed in the sensor frame
     * \return Nothing
     */
    /*inline*/ void set_X(const double X) { p[0] = X ; }
    
    /*!
     * \fn void set_Y()
     * \brief Sets the Y coordinate of the 3D point expressed in the sensor frame
     * \return Nothing
     */
    /*inline*/ void set_Y(const double Y) { p[1] = Y ; }
    
    /*!
     * \fn void set_Z()
     * \brief Sets the Z coordinate of the 3D point expressed in the sensor frame
     * \return Nothing
     */
    /*inline*/ void set_Z(const double Z) { p[2] = Z ; }
    
    /*!
     * \fn void set_W()
     * \brief Sets the W coordinate of the 3D point expressed in the sensor frame (generally set to 1, its default value for a Cartesian point)
     * \return Nothing
     */
    /*inline*/ void set_W(const double W) { p[3] = W ; }
    
    /*!
     * \fn double get_X()
     * \brief Accessor to the X coordinate of the 3D point expressed in the sensor frame
     * \return X
     */
    double get_X()  const { return p[0] ; }
    
    /*!
     * \fn double get_Y()
     * \brief Accessor to the Y coordinate of the 3D point expressed in the sensor frame
     * \return Y
     */
    double get_Y()  const { return p[1] ; }
    
    /*!
     * \fn double get_Z()
     * \brief Accessor to the Z coordinate of the 3D point expressed in the sensor frame
     * \return Z
     */
    double get_Z() const  { return p[2] ; }
    
    /*!
     * \fn double get_W()
     * \brief Accessor to the W coordinate of the 3D point expressed in the sensor frame (generally set to 1, its default value for a Cartesian point)
     * \return W
     */
    double get_W()  const { return p[3] ; }
    
    /*!
     * \fn prCartesian3DPointVec & changeFrame(const vpHomogeneousMatrix &oMi)
     * \brief Applies the frame change from input frame Fi to output frame Fo
     * \param oMi frame change from the world to the sensor frame
     * \return Modifies the output prCartesian3DPointVec variable
     */
    /*inline*/ prCartesian3DPointVec changeFrame(const vpHomogeneousMatrix & oMi);

    /*!
     * \fn void setPoint(const double & X, const double & Y, const double & Z)
     * \brief Sets the 3D point coordinates, forcing W to 1
     * \param X the point horizontal coordinate
     * \param Y the point vertical coordinate
     * \param Z the point depth coordinate
     * \return Nothing
     */
    void setPoint(const double & X, const double & Y, const double & Z);
    
    //    void setPoint(double X, double Y, double Z);
    
    /*!
     * \fn void setVector(const double & X, const double & Y, const double & Z)
     * \brief Sets the 3D vector coordinates, forcing W to 0
     * \param X the vector horizontal coordinate
     * \param Y the vector vertical coordinate
     * \param Z the vector depth coordinate
     * \return Nothing
     */
    void setVector(const double & X, const double & Y, const double & Z);
    
    /*!
     * \fn void set(const double & X, const double & Y, const double & Z, const double & W)
     * \brief Sets the 3D vector homogeneous coordinates
     * \param X the vector horizontal coordinate
     * \param Y the vector vertical coordinate
     * \param Z the vector depth coordinate
     * \param W the vecotr homogeneous coordinate
     * \return Nothing
     */
    void set(const double & X, const double & Y, const double & Z, const double & W);
    
    /*!
     * \fn int toEuclidean()
     * \brief Transforms the homogeneous coordinates of the point so that its homogeneous coordinate is equal to 1 (thus, only for a point)
     * \return  0 if not a vector
     *         -1 otherwise
     */
    int toEuclidean();
};

#endif  //_PRCARTESIAN3DPOINTVEC_H
