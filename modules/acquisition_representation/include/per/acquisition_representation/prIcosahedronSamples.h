/*!
 \file prSphericalIcoSamples.h
 \brief Header file for the prSphericalIcoSamples class
 \author Antoine ANDRÉ
 \version 0.1
 \date february 2023
 */

#if !defined(_PRICOSAHEDRONSAMPLES_H)
#define _PRICOSAHEDRONSAMPLES_H

#include <map>
#include <per/prCartesian3DPointVec.h>
#include <per/prXML.h>
#include <vector>

/*!< XML tag definition */
#define LABEL_NB_SAMPLES "nbVertices"
#define LABEL_FACES "faces"
#define LABEL_FACE "face"
#define LABEL_COORDS "coords"
#define LABEL_POINT "pr3DPoint"
#define LABEL_DELAU_POINT "delaunayVertex3DPoint"
#define LABEL_X "x"
#define LABEL_Y "y"
#define LABEL_Z "z"

/*!
 \class prIcosahedronSamples prIcosahedronSamples.h <per/prIcosahedronSamples.h>
 \brief Class for reading and loading samples coordinates coming from a
 subdivided icosahedron
 */
class PER_EXPORT prIcosahedronSamples : public prXML {
public:
  /*!
   * \fn prIcosahedronSamples()
   *
   * \brief Constructor of the icosahedron samples class. Set the path to `.xml`
   * files containing the coordinates of the icosahedron samples
   */
  prIcosahedronSamples(std::string coordsPath = PER_SI_DIR);//"/home/guillaume/Develop/PeR-ws/libPeR/data/subdividedIcosahedron/");

  /*!
   * \fn loadDelaunaySamples(unsigned int subdivLevels)
   *
   * \brief Allocates and fill the array of points coordinates of the Delaunay
   * graph of the Icosahedron as well as the array of "pixels" (not set in that
   * function since needs data)
   *
   * \return  0 if the spherical points coordinates are well loaded
   */
  int loadDelaunaySamples(unsigned int subdivLevels);

  /*!
   * \fn loadVoronoiSamples(unsigned int subdivLevels)
   *
   * \brief Allocate and fill in the array containing the 3D coordinates of the
   * Voronoi cells of the Icosahedron
   *
   * \return  0 if the spherical points coordinates are well loaded
   */
  int loadVoronoiSamples(unsigned int subdivLevels);

  void createNeighborMap(double limitDist);

  prIcosahedronSamples &operator=(const prIcosahedronSamples &icoS) {
    nbVertices = icoS.nbVertices;
    subdivLevels = icoS.subdivLevels;
    coordsPath = icoS.coordsPath;
    neighMap = icoS.neighMap;
    dist = icoS.dist;

    return *this;
  }

  double dist;
  int nbVertices;
  unsigned int subdivLevels;
  std::string coordsPath;
  prCartesian3DPointVec *samplesCoords;
  unsigned int maxSubdivLevels;
  std::map<unsigned int, std::vector<unsigned int>> neighMap;
  std::map<unsigned int, std::vector<prCartesian3DPointVec>> samplesVoronoiMap;
};

#endif //_PRSPHERICALICOSAMPLES_H
