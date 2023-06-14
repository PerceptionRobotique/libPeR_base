/*!
 \file prcommon.h
 \brief Header file common to the entire libPR
 \author Guillaume CARON
 \version 0.1
 \date january 2017
 */

#if !defined(_PRCOMMON_H)
#define _PRCOMMON_H

#include <cstring>

#ifdef WIN32
#ifdef PER_SHARED_EXPORT
#define PER_EXPORT __declspec(dllexport)
#else
#define PER_EXPORT __declspec(dllimport)
#endif
#elif __linux__
#ifdef PER_SHARED_EXPORT
#define PER_EXPORT __attribute__((visibility("default")))
#else
#define PER_EXPORT
#endif
#endif

typedef enum{
    IMAGEPLANE_NEARESTNEIGH,
    IMAGEPLANE_BILINEAR
} prInterpType;

/*!
 \class PoseMatrix prcommon.h <per/prcommon.h>
 \brief Class defining a basic 4x4 3D rigid transformation matrix to ensure pr_core being independent
 */
class PER_EXPORT PoseMatrix
{
public:
    /*!
     * \fn PoseMatrix()
     * \brief Default constructor of a PoseMatrix object
     */
    PoseMatrix()
    {
        double *pdata = data;
        for(int i = 0 ; i < 4 ; i++)
        {
            rowPtrs[i] = pdata;
            for(int j = 0 ; j < 4 ; j++, pdata++)
            {
                if(i == j)
                    *pdata = 1.0;
                else
                    *pdata = 0.0;
            }
        }
    }
    
    /*!
     * \fn ~PoseMatrix()
     * \brief Default destructor of a PoseMatrix object
     */
    ~PoseMatrix() {}
    
    /*!
     * \fn double *operator[](int n)
     * \brief [] operator override for writing in a cell of the PoseMatrix
     *
     * \param n row index
     * \return a pointer to the row in which a cell will be modified.
     */
    /*inline*/ double *operator[](int n){return rowPtrs[n];}
    /*!
     * \fn double *operator[](int n)
     * \brief [] operator override for reading a cell of the PoseMatrix
     *
     * \param n row index
     * \return a pointer to the row of which a cell will be read.
     */
    /*inline*/ double *operator[](int n) const {return rowPtrs[n];}
    
    /*!
     * \fn PoseMatrix &operator=(const PoseMatrix &M)
     * \brief = operator override for copying a PoseMatrix
     *
     * \param M the source PoseMatrix
     * \return the affected PoseMatrix.
     */
    PoseMatrix &operator=(const PoseMatrix &M)
    {
        for(int i = 0 ; i < 4 ; i++)
            for(int j = 0 ; j < 4 ; j++)
                rowPtrs[i][j] = M[i][j];
        
        return *this;
    }
    
    
    /*!
     * \fn void inverse()
     * \brief invert the homogeneous matrix
     *
     * [R T]^-1 = [R^T  -R^T T]
     *
     * \return   [R T]^-1
     */
    void inverse()
    {
        int i, j;
        double T[3] = {0.0, 0.0, 0.0};
        
        //transposed rotation bloc (Rt)
        for(i = 0 ; i < 3 ; i++)
            for(j = 0 ; j < 3 ; j++)
                rowPtrs[i][j] = rowPtrs[j][i];
        
        //T : -Rt * translation
        for (i = 0 ; i < 3; i++)
            for(j = 0 ; j < 3 ; j++)
                T[i] -= rowPtrs[i][j]*rowPtrs[j][3];
        
        //Affectating T to the current matrix
        for (i = 0 ; i < 3; i++)
            rowPtrs[i][3] = T[i];
    }
    
    /*!
     * \fn PoseMatrix t()
     * \brief transpose the PoseMatrix
     *
     * \return   M^T
     */
    PoseMatrix t()
    {
        PoseMatrix Mt ;
        int i,j;
        for (i=0;i<4;i++)
            for (j=0;j<4;j++)
                Mt[j][i] = (*this)[i][j];
        
        return Mt;
    }
    
private:
    double *rowPtrs[4]; /*!< Pointers to the PoseMatrix rows*/
public:
    double data[16]; /*!< The 16 cell values of the 4x4 PoseMatrix*/
};

/*!
 \class PoseVector prcommon.h "Common"
 \brief Class defining a basic 6-vector representing a pose (t_X, t_Y, t_Z, \theta w_X, \theta w_Y, \theta w_Z) to ensure pr_core being independent
 */
class PER_EXPORT PoseVector
{
public:
    PoseVector()
    {
        std::memset(data, 0, 6*sizeof(double));
    }
    
    //For writing
    /*inline*/ double & operator[](int n){return data[n];}
    //For reading
    /*inline*/ double operator[](int n) const {return data[n];}
    
    PoseVector &operator=(const PoseVector &v)
    {
        for(int i = 0 ; i < 6 ; i++)
            data[i] = v[i];
        
        return *this;
    }
    
private:
    double data[6];
};

/*!
 \struct s_SImagePoint prcommon.h "Common"
 \brief Basic 2-elements structure representing 2D point coordinates in the digital image plane. Typedefined as SImagePoint.
 */
typedef struct PER_EXPORT s_SImagePoint
{
    float u, v;
}SImagePoint;

#endif  //_PRCOMMON_H
