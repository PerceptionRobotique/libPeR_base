#include <per/prPerspective.h>

prPerspective::prPerspective(double au,double av,double u0,double v0, double k1,double k2,double k3,double k4,double k5,double k6,double k7,double k8):prCameraModel(au,av,u0,v0, k1,k2,k3,k4,k5,k6,k7,k8)
{
    //std::cout << "camera perspective created" << std::endl;
    name = "Perspective";
    type = Persp;
}

void prPerspective::project3DImage(prPointFeature &pt)
{
    pt.set_x(pt.get_X()/pt.get_Z());
    pt.set_y(pt.get_Y()/pt.get_Z());
    pt.set_w(1.0);
}

bool prPerspective::unProject(prPointFeature &pt, double &Depth)
{
    //std::cout << "prPerspective::unProject" << std::endl;
    pt.set_X(pt.get_x()*Depth);
    pt.set_Y(pt.get_y()*Depth);
    pt.set_Z(Depth);

    return true;
}

void prPerspective::project3DSphere(prPointFeature& P, double& Xs, double& Ys, double& Zs)
{
    double inv_norme;

    Xs = P.get_X();
    Ys = P.get_Y();
    Zs = P.get_Z();

    inv_norme = 1.0 / sqrt(Xs * Xs + Ys * Ys + Zs * Zs);

    Xs *= inv_norme;
    Ys *= inv_norme;
    Zs *= inv_norme;
}

bool prPerspective::projectImageSphere(prPointFeature &P, double &Xs, double &Ys, double &Zs) {}

std::ostream& prPerspective::operator<<(std::ostream & os)
{
    os << "Camera parameters for perspective projection without distortion:" << std::endl;
    prCameraModel::operator<<(os);
    
    return os;
}

// du / dX
void prPerspective::computeSensorJacobian(prPointFeature & P, vpMatrix & LuX)
{
    double x = P.get_x(), y = P.get_y(), iZ = 1.0 / P.get_Z();
    //get du / dx from mother class
    vpMatrix Lux(2,2);
    prCameraModel::computeSensorJacobian(P, Lux);

    vpMatrix LxX(2,3);
    LxX[0][0] = iZ; LxX[0][1] = 0;  LxX[0][2] = -x*iZ;
    LxX[1][0] = 0;  LxX[1][1] = iZ; LxX[1][2] = -y*iZ;

    LuX = Lux*LxX;
}
