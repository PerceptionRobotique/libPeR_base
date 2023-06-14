/*!
 \file saveAndLoadStereocameraModelXML.cpp
 \brief Example program of the PeR io module
 \author Guillaume CARON
 \version 0.1
 \date january 2017
 */

#include <iostream>

#include <per/prPerspective.h>
#include <per/prOmni.h>
#include <per/prStereoModel.h>

#include <per/prStereoModelXML.h>

/*!
 * \fn main()
 * \brief Main function of the sample program of the PeR io module
 *
 * Saving a stereovision system made of an omni cam and a perspective one to an XML file and load
 *
 */
int main()
{
    prStereoModel stereoCam(2);

    prSensorModel *omnicam = NULL, *perspcam = NULL; /*!< Generic sensor model */
    vpHomogeneousMatrix cpMco; /*!< omni cam to perspective one relative pose */

    // Omnidirectional camera definition
    omnicam = new prOmni(350, 350, 320, 240, 0.92);
    std::cout << "The stereo rig is made of a " << ((prCameraModel *)omnicam)->getName() << " camera of intrinsic parameters alpha_u = " << ((prCameraModel *)omnicam)->getau() << " ; alpha_v = " << ((prCameraModel *)omnicam)->getav() << " ; u_0 = " << ((prCameraModel *)omnicam)->getu0() << " ; v_0 = " << ((prCameraModel *)omnicam)->getv0() << " ; xi = " << ((prOmni *)omnicam)->getXi() << "..." << std::endl;
    
    // Perspective camera definition
    perspcam = new prPerspective(500, 500, 320, 240);
    std::cout << "... and a " << ((prCameraModel *)perspcam)->getName() << " camera of intrinsic parameters alpha_u = " << ((prCameraModel *)perspcam)->getau() << " ; alpha_v = " << ((prCameraModel *)perspcam)->getav() << " ; u_0 = " << ((prCameraModel *)perspcam)->getu0() << " ; v_0 = " << ((prCameraModel *)perspcam)->getv0() << "..." << std::endl;
    
    // Setting the relative pose between cameras
    cpMco.buildFrom(-0.12, 0.0, 0.37, -M_PI*0.5, 0.0, 0.0);
    std::cout << "...at camera pose cpMco = " << std::endl << cpMco << std::endl << " relatively to the first camera..." << std::endl;

    // Making the stereo rig
    stereoCam.setSensor(0, omnicam);
    stereoCam.setSensor(1, perspcam);
    stereoCam.setsjMr(1, cpMco);

    // Save the stereo rig to XML file
    {
        prStereoModelXML toFile("myStereoRig.xml");
        
        toFile << stereoCam;
    }
    std::cout << "... and after saving the rig to an XML file..." << std::endl;
    
    //Forget the rig and its elements
    stereoCam.init(0);
    // Memory free
    delete omnicam;
    delete perspcam;
    
    //Recreate an empty rig
    stereoCam.init(2);
    
    // Load the stereo rig parameters from the XML file
    {
        prStereoModelXML fromFile("myStereoRig.xml");
        
        fromFile >> stereoCam;
    }
    std::cout << "... and loading the XML file to an empty rig..." << std::endl;
    
    // If a sensor is loaded, print its parameters
    if(stereoCam.get_nbsens() >= 1)
    {
        std::cout << "the stereo rig is made of a " << ((prCameraModel *)stereoCam.sen[0])->getName() << " camera of intrinsic parameters alpha_u = " << ((prCameraModel *)stereoCam.sen[0])->getau() << " ; alpha_v = " << ((prCameraModel *)stereoCam.sen[0])->getav() << " ; u_0 = " << ((prCameraModel *)stereoCam.sen[0])->getu0() << " ; v_0 = " << ((prCameraModel *)stereoCam.sen[0])->getv0();
        if(((prCameraModel *)stereoCam.sen[0])->getType() == Omni)
            std::cout << " ; xi = " << ((prOmni *)stereoCam.sen[0])->getXi();
        std::cout << "..." << std::endl;
    }

    // If a second sensor is loaded, print its parameters and its pose relatively to the first camera
    if(stereoCam.get_nbsens() >= 2)
    {
        std::cout << "... and a " << ((prCameraModel *)stereoCam.sen[1])->getName() << " camera of intrinsic parameters alpha_u = " << ((prCameraModel *)stereoCam.sen[1])->getau() << " ; alpha_v = " << ((prCameraModel *)stereoCam.sen[1])->getav() << " ; u_0 = " << ((prCameraModel *)stereoCam.sen[1])->getu0() << " ; v_0 = " << ((prCameraModel *)stereoCam.sen[1])->getv0();
        if(((prCameraModel *)stereoCam.sen[1])->getType() == Omni)
            std::cout << " ; xi = " << ((prOmni *)stereoCam.sen[1])->getXi();
        std::cout << "..." << std::endl;
   
        std::cout << "...at camera pose cpMco = " << std::endl << stereoCam.sjMr[1] << std::endl << " relative to the first camera." << std::endl;
    }
    
	return 0;
}
