/*!
 \file pointPerspectiveProjection.cpp
 \brief Example program of the PeR core module
 \author Guillaume CARON
 \version 0.1
 \date january 2017
 */

#include <iostream>

#include <per/prPerspective.h>
#include <per/prPointFeature.h>

#include <visp/vpHomogeneousMatrix.h>

/*!
 * \fn main()
 * \brief Main function of the sample program of the PeR core module
 *
 * Perspective projection of a 3D point in a digital image plane
 *
 */
int main()
{
    prPointFeature point; /*!< A feature point */
    prSensorModel *perspcam = NULL; /*!< Generic sensor model */
    vpHomogeneousMatrix cMo; /*!< Sensor pose */
    
    // 3D point coordinates definition in the object frame
    point.setWorldCoordinates(0.1, 0.1, 0.1);
    std::cout << "The 3D point of coordinates oX = " << point.get_oX() << " ; oY = " << point.get_oY() << " ; oZ = " << point.get_oZ() << " in the object frame..." << std::endl;
    
    // Perspective camera definition
    int width = 640, height = 480;
    perspcam = new prPerspective(500, 500, width / 2, height / 2);
    std::cout << "...is observed by a " << ((prCameraModel *)perspcam)->getName() << " camera of intrinsic parameters alpha_u = " << ((prCameraModel *)perspcam)->getau() << " ; alpha_v = " << ((prCameraModel *)perspcam)->getav() << " ; u_0 = " << ((prCameraModel *)perspcam)->getu0() << " ; v_0 = " << ((prCameraModel *)perspcam)->getv0() << "..." << std::endl;
    
    // Object to camera frames transformation definition
    cMo.buildFrom(0.25, 0.0, 0.33, M_PI*0.2, 0.0, 0.0);
    std::cout << "...at camera pose cMo = " << std::endl << cMo << std::endl;

    // Frame change application to get the coordinates in the camera frame
    point.changeFrame(cMo);
    std::cout << "...with respect to what the 3D point coordinates become cX = " << point.get_X() << " ; cY = " << point.get_Y() << " ; cZ = " << point.get_Z() << " in the camera frame..." << std::endl;
    
    // Perspective projection in the normalized image plane
    ((prCameraModel *)perspcam)->project3DImage(point);
    std::cout << "...that projects as x = " << point.get_x() << " ; y = " << point.get_y() << " in the nomalized image plane..." << std::endl;
    
    // Affine transformation to the digital image plane
    ((prCameraModel *)perspcam)->meterPixelConversion(point);
    std::cout << "...corresponding to coordinates u = " << point.get_u() << " ; v = " << point.get_v() << " in the digital image plane" << std::endl;
    
    // Memory free
    delete perspcam;
    
	return 0;
}
