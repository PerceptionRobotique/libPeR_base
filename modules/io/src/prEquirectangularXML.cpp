#include <per/prEquirectangularXML.h>

prEquirectangularXML::prEquirectangularXML(std::string fic):prCameraModelXML(fic)
{
    
}

void prEquirectangularXML::operator<<(prEquirectangular &cam)
{
    prCameraModelXML::operator<<((prCameraModel*)&cam);
}

void prEquirectangularXML::operator>>(prEquirectangular &cam)
{
    prCameraModelXML::operator>>((prCameraModel*)&cam);
}
