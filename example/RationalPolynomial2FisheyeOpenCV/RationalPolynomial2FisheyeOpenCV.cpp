/*!
 \file RationalPolynomial2FisheyeOpenCV.cpp
 \brief Example program to convert a PeR camera model to another
 \author Eva Goichon
 \version 0.1
 \date July 2024
 */

#include <iostream>

#include <per/prPerspective.h>
#include <per/prCameraModelConvert.h>

#include <visp/vpMatrix.h>
// #include <prCameraModelConvert.h>

/*!
 * \fn main()
 * \brief Example program's main function for converting a set of distortion parameters from a polynomial model to OpenCV's Fisheye model
 *
 * Output the Fisheye OpenCV camera model
 * TODO: load/save xml files
 *
 */
int main()
{
    // Perspective camera definition with rational polynomial model of distortions
    // Example of an Azure Kinect Depth camera 
    //double au = 251.93849182128906; //fx 
    //double av = 252.07272338867188; //fy
    //double u0 = 254.28904724121094; //cx
    //double v0 = 255.1665802001953;  //cy

    // Example of an Azure Kinect IR camera 
    double au = 503.8769836425781; //fx 
    double av = 504.14544677734375; //fy
    double u0 = 509.0780944824219; //cx
    double v0 = 510.8331604003906; //cy
    
    double k1 = 0.4452361762523651;
    double k2 = -0.027260301634669304;
    double k3 = -0.0019093812443315983;
    double k4 = 0.00011894194904016331; // tangential p1
    double k5 = 2.8838716389145702e-05; // tangential p2
    double k6 = 0.7864969968795776;
    double k7 = 0.04874652251601219;
    double k8 = -0.011641541495919228;
    
    double FoV = 120; // degrees; horizontal
    
    prPerspective inperspcam(au, av, u0, v0, k1, k2, k3, k4, k5, k6, k7, k8);
    prPerspective outfisheye(au, av, u0, v0, k1, k2, k3, k4);


    std::cout << "The input camera is a " << inperspcam.getName() << " camera of intrinsic parameters " << std::endl;
    inperspcam.operator<<(std::cout);

    vpMatrix abs_err;
    double residual = prCameraModelConvert::conversion(&inperspcam, &outfisheye, FoV, &abs_err);
    vpMatrix::save("residuals.txt", abs_err);

    // std::cout << "The ouput camera is a Fisheye camera of intrinsic parameters " << std::endl;
    // outfisheye.operator<<(std::cout);

    std::cout << "The ouput camera is  OpencvFisheye camera of intrinsic parameters " << std::endl;
    std::cout << "Camera parameters for OpencvFisheye projection without distortion:" << std::endl;
    std::cout << "  au = " << outfisheye.getau() <<"\t av = "<< outfisheye.getav() << std::endl ;
    std::cout << "  u0 = " << outfisheye.getu0() <<"\t v0 = "<< outfisheye.getv0() << std::endl ;
    std::cout <<  "with distorsions parameters for undistorted to distorted transformation :" << std::endl;
    std::cout << "k1 = " << outfisheye.getk1() <<"\t k2 = " << outfisheye.getk2() <<"\t k3 = " << outfisheye.getk3() <<"\t k4 = " << outfisheye.getk4()<< std::endl;


    
    std::cout << "An average error of " << residual << " pixel(s) has been made during the conversion" << std::endl;

    
	return 0;
}
