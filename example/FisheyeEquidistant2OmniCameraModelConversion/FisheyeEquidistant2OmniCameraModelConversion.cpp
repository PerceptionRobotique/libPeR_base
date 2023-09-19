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
    // Example of an M-12 Lensagon fisheye BF10M14522S118 lens ...
    double f = 1.45e-3; // focal length 1.45mm (Lensation specs)
    double FoV = 190; // degrees
    // ...mounted on a eCon e-CAM50 CUNX/NANO 5MP camera (1/2.5" AR0521 CMOS Image sensor from onsemi®)
    double k = 2.2e-6; // square photodiodes: ku = kv
    double width = 2592; // pixels
    double height = 1944; // pixels
    
    fisheyeEquidistantcam.init(f/k, f/k, width*0.5, height*0.5);
    std::cout << "The input camera is a " << fisheyeEquidistantcam.getName() << " camera of intrinsic parameters alpha_u = " << fisheyeEquidistantcam.getau() << " ; alpha_v = " << fisheyeEquidistantcam.getav() << " ; u_0 = " << fisheyeEquidistantcam.getu0() << " ; v_0 = " << fisheyeEquidistantcam.getv0() << std::endl;
    
    // Conversion of the equidistant fisheye to an omni camera
    double residual = prCameraModelConvert::convert(fisheyeEquidistantcam, omnicam, FoV);
    
    std::cout << "The ouput camera is a " << omnicam.getName() << " camera of intrinsic parameters alpha_u = " << omnicam.getau() << " ; alpha_v = " << omnicam.getav() << " ; u_0 = " << omnicam.getu0() << " ; v_0 = " << omnicam.getv0() << " ; xi = " << omnicam.getXi() << std::endl;
    
    std::cout << "An average error of " << residual << " pixel(s) has been made during the conversion" << std::endl;

    
	return 0;
}
