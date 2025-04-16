/*!
 \file CartesianPolynomial2FisheyeOpenCV.cpp
 \brief Example program to convert a PeR camera model to another
 \author Eva GOICHON    
 \version 0.1
 \date August 2024
 */

#include <iostream>

#include <per/prPerspective.h>
#include <per/prCameraModelConvert.h>
#include <per/prPolyCart.h>

#include <visp/vpMatrix.h>

/*!
 * \fn main()
 * \brief Main function of the sample program to convert a PeR camera model to another
 *
 * Output the Fisheye OpenCV camera model computed from an prPolyCart
 * TODO: load/save xml files
 *
 */
int main()
{
    prPolyCart polycartcam;
    
    double au = 1;
    double av = 1; 
    double a[5];
    
    // //py_ocamcalib Sequence6 Panoramis
    // double u0 = 316.11746197386447; 
    // double v0 = 310.8244883081968;
    // a[0]= 123.14453786699245;
    // a[1] = 0;
    // a[2] = -0.003153164438933958;
    // a[3] = 7.209622221439549e-06;
    // a[4] = -1.2886257196977669e-08;
    // double FoV = 210;
    // double image_size = 620;
    
    //Calibrated with ocamcalib (4th order): Sequence 6 Panoramis
    double u0 = 321.502861; 
    double v0 = 311.665234;
    a[0]= 121.1861;
    a[1] = 0;
    a[2] = -2.791683e-03 ;
    a[3] = 4.565693e-06;
    a[4] = -7.412085e-09;
    double FoV = 210;
    double image_size = 620;

    // //Parameters of Zurich Cata dataset https://rpg.ifi.uzh.ch/fov.html
    // double u0 = 320.0;
    // double v0 = 240.0;
    // a[0] = 70.25584688010954;
    // a[1] = 0;
    // a[2] = -0.0005433216164761568;
    // a[3] = -2.102936744599082e-05;
    // a[4] = 8.54310806364036e-09;
    // double FoV = 210;
    // double image_size = 640;

    // //Parameters of Zurich fisheye dataset https://rpg.ifi.uzh.ch/fov.html
    // double u0 = 320.0;
    // double v0 = 240.0;
    // a[0] = 179.471829787234;
    // a[1] = 0;
    // a[2] = -0.002316743975;
    // a[3] = 3.635968439375e-06;
    // a[4] = -2.0546506810625e-08;
    // double FoV = 180;
    // double image_size = 640;

    // //Scaramuzza parameters for his dataset
    // double u0 =519.4290787359728;//516.4379;//311.875204; 
    // double v0 = 377.0391374256972;//383.0140;//321.364040;
    // a[0]= 130.72918883362064;// 131.0074;//1.211861e+02;
    // a[1] = 0;
    // a[2] = -0.002093358194106016;//-0.0018;//-2.791683e-03;
    // a[3] = 1.5028352840062567e-06;//4.565693e-06; HEAD
    // a[4] = -2.44444779088339e-09;//-7.412085e-09;
    
    double k1 = 0;
    double k2 = 0; 
    double k3 = 0;
    double k4 = 0;

    polycartcam.init(au, u0, v0, a);
     
    prPerspective outfisheye(au, av, u0, v0, k1, k2, k3, k4);

    std::cout << "The input camera is a "<< polycartcam.getName() << " camera of intrinsic parameters alpha_u = alpha_v = " 
            << polycartcam.getau() << " ; u_0 = " << polycartcam.getu0() << " ; v_0 = " << polycartcam.getv0() << " ; a_0 = " 
            << polycartcam.geta0() << " ; a_1 = " << polycartcam.geta1() << " ; a_2 = " << polycartcam.geta2() << " ; a_3 = " 
            << polycartcam.geta3() << " ; a_4 = " << polycartcam.geta4() << std::endl;
    
    // Conversion of the cartesian polynomial model to OpenCV's Fisheye model and export the algebraic residuals
    vpMatrix abs_err;
    
    //conversion with reverse approach
    double residual = prCameraModelConvert::FCVreverseConvert(&polycartcam, &outfisheye, FoV, &abs_err,image_size);

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

