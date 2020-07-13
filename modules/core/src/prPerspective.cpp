#include <per/prPerspective.h>

prPerspective::prPerspective(double au,double av,double u0,double v0, double k1,double k2,double k3,double k4,double k5):prCameraModel(au,av,u0,v0, k1,k2,k3,k4,k5)
{
    name = "Perspective";
    type = Persp;
}

prPerspective::~prPerspective()
{
    
}

void prPerspective::project3DImage(prPointFeature &pt)
{
    pt.set_x(pt.get_X()/pt.get_Z());
    pt.set_y(pt.get_Y()/pt.get_Z());
    pt.set_w(1.0);
}

std::ostream& prPerspective::operator<<(std::ostream & os)
{
    os << "Camera parameters for perspective projection without distortion:" << std::endl;
    prCameraModel::operator<<(os);
    
    return os;
}

