/*!
 \file FisheyeEquidistant2OmniCameraModelConversion.cpp
 \brief Example program to convert a PeR camera model to another
 \author Guillaume CARON
 \version 0.1
 \date march 2022
 */

#include <iostream>

#include <per/prFisheyeEquidistant.h>
#include <per/prOmni.h>
#include <per/prCameraModelConvert.h>

#include <visp/vpMatrix.h>

/*!
 * \fn main()
 * \brief Main function of the sample program to convert a PeR camera model to another
 *
 * Output the omni camera model computed from an equidistant fisheye
 * TODO: load/save xml files
 *
 */
int main()
{
    //prSensorModel *omnicam = NULL, *fisheyeEquidistantcam = NULL; /*!< Generic sensor model */
    prFisheyeEquidistant fisheyeEquidistantcam;
    prOmni omnicam;

    // Equidistant fisheye camera definition
    /*
    // Example of a C-Mount Fujinon fisheye FE185C086HA-1 lens ...
    double f = 2.7e-3; // focal length 2.7mm (Fujinon specs)
    double FoV = 185; // degrees
    // ...mounted on  a Prophesee Gen 3.1 event camera
    double k = 15e-6; // square photodiodes: ku = kv
    double width = 640; // pixels
    double height = 480; // pixels
    */
    /*
    // Example of an M-12 Lensagon fisheye BF10M14522S118 lens ...
    double f = 1.45e-3; // focal length 1.45mm (Lensation specs)
    double FoV = 190; // degrees
    // ...mounted on a eCon e-CAM50 CUNX/NANO 5MP camera (1/2.5" AR0521 CMOS Image sensor from onsemi®)
    double k = 2.2e-6; // square photodiodes: ku = kv
    double width = 2592; // pixels
    double height = 1944; // pixels
	  */
	  /*
    // Example of a M12toCS-Mount Entaniya 280 lens ...
    double f = 1.07e-3; // focal length 1.07mm (Entaniya specs)
    double FoV = 280; // degrees
    // ...mounted on a Flir FL3-U3-13S2C  camera
    double k = 3.63e-6; // square photodiodes: ku = kv
    double width = 1328; // pixels
    double height = 1048; // pixels
    */
    
    /*
    // Example of a C-Mount Fujinon fisheye FE185C086HA-1 lens ...
    double f = 2.7e-3; // focal length 2.7mm (Fujinon specs)
    double FoV = 185; // degrees
    // ...mounted on an IDS UI-3370CP camera
    double k = 5.5e-6; // square photodiodes: ku = kv
    double width = 2048; // pixels
    double height = 2048; // pixels
    */    
    
    // // Example of a C-Mount Fujinon fisheye FE185C086HA-1 lens ...
    // double f = 2.7e-3; // focal length 2.7mm (Fujinon specs)
    // double FoV = 185; // degrees
    // // ...mounted on a Flir FL3-U3-13S2C  camera
    // double k = 5.3e-6; // square photodiodes: ku = kv
    // double width = 1280; // pixels
    // double height = 1024; // pixels

    /*
    // Example of a M12 mount Entaniya Pi Lens 185 fisheye lens ...
    double f = 1.39e-3; // focal length 1.39mm (Entaniya specs)
    double FoV = 185; // degrees
    */
    /*
    // Example of a M12 mount Entaniya 280 fisheye lens ...
    double f = 1.07e-3; // focal length 1.39mm (Entaniya specs)
    double FoV = 280; // degrees
    // ...mounted on a PointGrey FL3-U3-13S2C-CS  camera
    double k = 3.63e-6; // square photodiodes: ku = kv
    double width = 1328; // pixels
    double height = 1048; // pixels
    */

    // Example of a C-Mount Fujinon fisheye FE185C086HA-1 lens ...
    double f = 2.7e-3; // focal length 2.7mm (Fujinon specs)
    double FoV = 185; // degrees
    // ...mounted on a Flir GS3-U3-41C6C-C camera
    double k = 2*5.5e-6; //2* because 5.5 is for non binned 
    double width = 1024; // pixels
    double height = 1024; // pixels
    double u_0 = 506; // center of the image circle
    double v_0 = 490;
    
    fisheyeEquidistantcam.init(f/k, f/k, u_0, v_0);
    std::cout << "The input camera is a " << fisheyeEquidistantcam.getName() << " camera of intrinsic parameters alpha_u = " << fisheyeEquidistantcam.getau() << " ; alpha_v = " << fisheyeEquidistantcam.getav() << " ; u_0 = " << fisheyeEquidistantcam.getu0() << " ; v_0 = " << fisheyeEquidistantcam.getv0() << std::endl;
    
    // Conversion of the equidistant fisheye to an omni camera and export the geometric residuals
    vpMatrix abs_err;
    double residual = prCameraModelConvert::convert(fisheyeEquidistantcam, omnicam, FoV, &abs_err);
    vpMatrix::save("residuals.txt", abs_err); 
    
    std::cout << "The ouput camera is a " << omnicam.getName() << " camera of intrinsic parameters alpha_u = " << omnicam.getau() << " ; alpha_v = " << omnicam.getav() << " ; u_0 = " << omnicam.getu0() << " ; v_0 = " << omnicam.getv0() << " ; xi = " << omnicam.getXi() << std::endl;
    
    std::cout << "An average error of " << residual << " pixel(s) has been made during the conversion" << std::endl;

	//Project a point of the border of the field-of-view to check its coordinates in the image plane
	prPointFeature point; /*!< A feature point */
	vpHomogeneousMatrix cMo; /*!< Sensor pose */
	
  // 3D point coordinates definition in the object frame
  point.setWorldCoordinates(sin(FoV*0.5*M_PI/180.), 0.0, cos(FoV*0.5*M_PI/180.));
  
  // Frame change application to get the coordinates in the camera frame
  point.changeFrame(cMo);
  std::cout << "...with respect to what the 3D point coordinates become cX = " << point.get_X() << " ; cY = " << point.get_Y() << " ; cZ = " << point.get_Z() << " in the camera frame..." << std::endl;
    
  // Fisheye equidistant projection in the normalized image plane
  fisheyeEquidistantcam.project3DImage(point);
  std::cout << "...that projects as x = " << point.get_x() << " ; y = " << point.get_y() << " in the nomalized image plane with the equidistant model..." << std::endl;
    
  // Affine transformation to the digital image plane
  fisheyeEquidistantcam.meterPixelConversion(point);
  std::cout << "...corresponding to coordinates u = " << point.get_u() << " ; v = " << point.get_v() << " in the digital image plane" << std::endl;
  
  // Omni projection in the normalized image plane
  omnicam.project3DImage(point);
  std::cout << "...that projects as x = " << point.get_x() << " ; y = " << point.get_y() << " in the nomalized image plane with the omnidirectional model..." << std::endl;
    
  // Affine transformation to the digital image plane
  omnicam.meterPixelConversion(point);
  std::cout << "...corresponding to coordinates u = " << point.get_u() << " ; v = " << point.get_v() << " in the digital image plane" << std::endl;
    
	return 0;
}
