/*!
 \file Omni2PolyCartCameraModelConversion.cpp
 \brief Example program to convert a PeR camera model to another
 \author Guillaume CARON
 \version 0.1
 \date June 2022
 */

#include <iostream>

#include <per/prOmni.h>
#include <per/prCameraModelConvert.h>
#include <per/prPolyCart.h>

#include <visp/vpMatrix.h>

/*!
 * \fn main()
 * \brief Main function of the sample program to convert a PeR camera model to another
 *
 * Output the prPolyCart camera model computed from an omni
 * TODO: load/save xml files
 *
 */
int main()
{
    //prSensorModel *omnicam = NULL, *fisheyeEquidistantcam = NULL; /*!< Generic sensor model */    
    prOmni omnicam;
		prPolyCart polycartcam;

    // Omni camera definition
    
    // // Example of a XXX catadioptric optics...
    // double au = 195; // 
    // double av = 195; // 
    // // ...mounted on  
    // double xi = 0.8; // unitless
    // double width = 430; // pixels
    // double height = 430; // pixels
    // double u0 = width*0.5; // pixels
    // double v0 = height*0.5; // pixels
    
    // Example of the PanoraMIS dataset, parameters are given with the dataset (Sequence 6)
    double au = 231.4616771; // 
    double av = 232.4224827; // 
    // ...mounted on  
    double xi = 0.958288499; // unitless
    double u0 = 319.7040894 ; // pixels
    double v0 = 310.9443395; // pixels
    double theoreticalFOV = 210;

    // // Example of the OCamCalib toolbox dataset (VMRImage*.jpg), calibrated with the Barreto's model using MIXEDVISION
    // double au = 259.8888419; // 
    // double av = 259.3345743; // 
    // // ...mounted on  
    // double xi = 0.9751280249; // unitless
    // double u0 = 514.1675494; // pixels
    // double v0 = 382.7969728; // pixels
    // double theoreticalFOV = 180;

    omnicam.init(au, av, u0, v0, xi);
    std::cout << "The input camera is a " << omnicam.getName() << " camera of intrinsic parameters alpha_u = " << omnicam.getau() << " ; alpha_v = " << omnicam.getav() << " ; u_0 = " << omnicam.getu0() << " ; v_0 = " << omnicam.getv0() << " ; xi = " << omnicam.getXi() << std::endl;
    
    // Conversion of the equidistant fisheye to an omni camera and export the algebraic residuals
    vpMatrix abs_err;
    double residual = prCameraModelConvert::convert(omnicam, polycartcam, theoreticalFOV, &abs_err);
    vpMatrix::save("residuals.txt", abs_err); 

    std::cout << "The ouput camera is a " << polycartcam.getName() << " camera of intrinsic parameters alpha_u = alpha_v = " << polycartcam.getau() << " ; u_0 = " << polycartcam.getu0() << " ; v_0 = " << polycartcam.getv0() << " ; a_0 = " << polycartcam.geta0() << " ; a_1 = " << polycartcam.geta1() << " ; a_2 = " << polycartcam.geta2() << " ; a_3 = " << polycartcam.geta3() << " ; a_4 = " << polycartcam.geta4() << std::endl;
    
    std::cout << "An average error of " << residual << " pixel(s) has been made during the conversion" << std::endl;

    
	return 0;
}
