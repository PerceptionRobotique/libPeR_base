#include <per/prPerspectiveXML.h>

prPerspectiveXML::prPerspectiveXML(std::string fic):prCameraModelXML(fic)
{
    
}

void prPerspectiveXML::operator<<(prPerspective &cam)
{
    prCameraModelXML::operator<<((prCameraModel*)&cam);
}

void prPerspectiveXML::operator>>(prPerspective &cam)
{
    prCameraModelXML::operator>>((prCameraModel*)&cam);
}
