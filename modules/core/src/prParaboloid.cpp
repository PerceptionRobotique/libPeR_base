#include <per/prParaboloid.h>

prParaboloid::prParaboloid(double au,double av,double u0,double v0,double hh, double k1, double k2, double k3, double k4, double k5):prCameraModel(au,av,u0,v0, k1, k2, k3, k4, k5)
{
    type = Paraboloid;
    name = "Paraboloid";
    h = hh;
}

void prParaboloid::init(double au,double av,double u0,double v0,double hh, double k1, double k2, double k3, double k4, double k5)
{
    prCameraModel::init(au,av,u0,v0, k1, k2, k3, k4, k5);
    type = Paraboloid;
    name = "Paraboloid";
    h=hh;
}

void prParaboloid::project3DImage(prPointFeature & P)
{
    double X = P.get_X(), Y = P.get_Y(), Z = P.get_Z(), fact;
    fact = h / (sqrt(X*X+Y*Y+Z*Z)+Z);
    
    P.set_x(fact * X);
    P.set_y(fact * Y);
    //P.set_z(fact * Z); //sur le miroir...
}

bool prParaboloid::unProject(prPointFeature & P, double & Depth)
{
    return false;
    //TODO
    //return true;
}

void prParaboloid::projectImageMiroir(prPointFeature & P, double & Xm, double & Ym, double & Zm)
{
    double x = P.get_x(), y = P.get_y();
    
    Xm = x;
    Ym = y;
    Zm = (h*h-x*x-y*y)/(2*h);//((-(h*h)+(x*x+y*y))/(2*h));
    //Zm = (x*x+y*y)/(4*h*h);
}

void prParaboloid::project3DMiroir(prPointFeature & P, double & Xm, double & Ym, double & Zm)
{
    double X = P.get_X(), Y = P.get_Y(), Z = P.get_Z(), fact;
    fact = h / (sqrt(X*X+Y*Y+Z*Z)+Z);
    
    Xm = fact * X;
    Ym = fact * Y;
    Zm = fact * Z;
}

prParaboloid& prParaboloid::operator=(const prParaboloid& cam)
{
    prCameraModel::operator=(cam);
    this->h = cam.h;
    return *this;
}

// du / dX
void prParaboloid::computeSensorJacobian(prPointFeature & P, vpMatrix & LuX)
{
    //TODO
}
