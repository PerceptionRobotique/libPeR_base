
#if !defined(_PRNEIGHBORHOOD_H)
#define _PRNEIGHBORHOOD_H

#include <per/prRegularlySampledCSImage.h>
// #include <Eigen/Dense>
#include <iostream>
#include <vector>

/*!
 \class prNeighborhood prNeighborhood.h <per/prNeighborhood.h>
 \brief Class for the storage of images neighbors lying on the sphere

 This class aims to create a unique table that contains all the neighbors of an image once projected on the sphere. To do this, a pointer table ```*****Neigh``` is created.
 */
class prNeighborhood
{
public:
    /*!
     * \fn prNeighborhood()
     * \brief Contructor of the prNeighborhood neighbors
     */
    prNeighborhood();

    /*!
     * \fn buildNeighborsCartSphere(int imgWidth, int imgHeight, prCameraModel *camera)
     * \brief Build the neighbors according to the cartesian spherical representation (indexes in the image) according to their proximity on the sphere, given a given projection through a camera type and parameters
     * \param imgWidth width of the considered image
     * \param imgHeight height of the considered image
     * \param camera used sensor to project the image points to the sphere and back, used for computing the neighbors positions
     */
    void buildNeighborsCartSphere(int imgWidth, int imgHeight, prCameraModel *camera);

    /*!
     * \fn deleteNeighbors()
     * \brief Memory delete the neighbors (not implemented yet)
     */
    void deleteNeighbors();

    /*!
     * \fn computeIndexNeighbor(int winSize, double inputIndex)
     * \brief Internal function used to compute the index of the neighbor once back propagated from the sphere to the image plane.
     */
    int computeIndexNeighbor(int winSize, double inputIndex);

    /*!
     * \fn prNeighborhood &operator=(const prNeighborhood &neigh)
     * \brief = operator override to copy neighbor object
     * \TODO: make a clean memory copy of the neighbors to avoid computing them again
     * \param neigh the prNeighborhood object to copy
     * \return the prNeighborhood self reference object
     */
    prNeighborhood &operator=(const prNeighborhood &neigh);

    int *****Neigh;         /*!< Size of the neighborhood table : 3 (for x, y, z directions) x M (number of pixel along first direction) x N (number of pixels along second direction) x 2 (number considered to index a nighbor) x k (number of neighbors, 6 here) */
    bool isNeighSet;        /*!< bool to check if the neighborhood is already set or not */
    double deltaCSSampling; /*!< Sampling size to navigate along the sphere */
    int imgH, imgW, dimS;   /*!< image size */
    int neighSize;          /*!< Total size of the neighborhood table, not used */
    int di, dj;             /*!< frame offset, half the size of considered neighbors (i.e. 3) */
    int nba;
    int nbNeigh;
    double *****pixNeigh;
    std::vector<int> ***quadNeigh;
};

#endif