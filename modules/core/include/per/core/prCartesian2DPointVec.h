/*!
 \file prCartesian2DPointVec.h
 \brief Header file for the prCartesian2DPointVec class
 \author Guillaume CARON
 \version 0.1
 \date february 2017
 */

#if !defined(_PRCARTESIAN2DPOINTVEC_H)
#define _PRCARTESIAN2DPOINTVEC_H

#include <per/pr2DPointVec.h>

/*!
 \class PER_EXPORT prCartesian2DPointVec prCartesian2DPointVec.h <per/prCartesian2DPointVec.h>
 \brief Class defining a 2D point in homogeneous Cartesian coordinates
 */
class PER_EXPORT prCartesian2DPointVec : public pr2DPointVec
{
    
public:
    prCartesian2DPointVec();

    virtual ~prCartesian2DPointVec() override = default;
    
    /*!
     * \fn void set_x(const double x)
     * \brief Sets the x coordinate of the point expressed in the normalized image plane
     * \param x the horizontal coordinate in the normalized image plane
     * \return Nothing
     */
    /*inline*/ void set_x(const double x) {  p[0] = x ; }
    
    /*!
     * \fn void set_y(const double y)
     * \brief Sets the y coordinate of the point expressed in the normalized image plane
     * \param y the vertical coordinate in the normalized image plane
     * \return Nothing
     */
    /*inline*/ void set_y(const double y) {  p[1] = y ; }
    
    /*!
     * \fn void set_w(const double w)
     * \brief Sets the w coordinate of the point expressed in the normalized image plane
     * \param w the homogeneous coordinate in the normalized image plane (generally set to 1, its default value for a Cartesian point)
     * \return Nothing
     */
    /*inline*/ void set_w(const double w) {  p[2] = w ; }
    
    
    /*!
     * \fn double get_x()
     * \brief Accessor to the horizontal coordinate x in the normalized image plane
     * \return x
     */
    double get_x()  const { return p[0] ; }
    
    /*!
     * \fn double get_y()
     * \brief Accessor to the vertical coordinate y in the normalized image plane
     * \return y
     */
    double get_y()  const { return p[1] ; }
    
    /*!
     * \fn double get_w()
     * \brief Accessor to the homogeneous coordinate w in the normalized image planee (generally set to 1, its default value for a Cartesian point)
     * \return x
     */
    double get_w()  const { return p[2] ; }

    /*!
     * \fn unsigned int size()
     * \brief Get the number of coordinates of the 2D point
     * \return the number of coordinates
     */
    unsigned int size()
    {
      return 2;
    }
    
    /*!
     * \fn prCartesian2DPointVec & changeFrame(const vpHomogeneousMatrix &oMi)
     * \brief Applies the frame change from input frame Fi to output frame Fo
     * \param oMi frame change from the world to the sensor frame
     * \return Modifies the output prCartesian2DPointVec variable
     */
    /*inline*/ prCartesian2DPointVec changeFrame(const vpHomogeneousMatrix & oMi);

    /*!
     * \fn void setPoint(const double & x, const double & y)
     * \brief Sets the 2D point coordinates, forcing w to 1
     * \param x the point horizontal coordinate
     * \param y the point vertical coordinate
     * \return Nothing
     */
    void setPoint(const double & x, const double & y);
    
    /*!
     * \fn void setVector(const double & x, const double & y)
     * \brief Sets the 2D vector coordinates, forcing w to 0
     * \param x the vector horizontal coordinate
     * \param y the vector vertical coordinate
     * \return Nothing
     */
    void setVector(const double & x, const double & y);
    
    /*!
     * \fn void set(const double & x, const double & y, const double & w)
     * \brief Sets the 2D vector homogeneous coordinates
     * \param x the vector horizontal coordinate
     * \param y the vecotr vertical coordinate
     * \param w the vecotr homogeneous coordinate
     * \return Nothing
     */
    void set(const double & x, const double & y, const double & w);
    
    /*!
     * \fn int toEuclidean()
     * \brief Transforms the homogeneous coordinates of the point so that its homogeneous coordinate is equal to 1 (thus, only for a point)
     * \return  0 if not a vector
     *         -1 otherwise
     */
    int toEuclidean();

    /*!
     * \fn vpColVector& operator-(const prCartesian2DPointVec &f)
     * \brief - operator override to substract two 2D points or vectors
     * *this and f are supposed Euclidean
     * valid for points only
     * \param f the 2D coordinates to substract from the current
     * \return the prPhotometricnnGMS<Ts> object resulting from the difference computation
     */
    vpColVector operator-(const prCartesian2DPointVec &f)
    {
        vpColVector res(2);
        
        res[0] = p[0] - f.p[0];
        res[1] = p[1] - f.p[1];

        return res;
    }
    
};

#endif  //_PRCARTESIAN3DPOINTVEC_H
