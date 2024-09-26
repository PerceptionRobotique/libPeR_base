#include <per/prFisheyeEquidistant.h>

prFisheyeEquidistant::prFisheyeEquidistant(double au,double av,double u0,double v0, double k1,double k2,double k3,double k4,double k5):prCameraModel(au,av,u0,v0, k1, k2, k3, k4, k5)
{
    type = FisheyeEquidistant;
    name = "FisheyeEquidistant";
}

void prFisheyeEquidistant::init(double _au,double _av,double _u0,double _v0, double _k1,double _k2,double _k3,double _k4,double _k5)
{
    prCameraModel::init(_au,_av,_u0,_v0,_k1,_k2,_k3, _k4, _k5);
    
    type = FisheyeEquidistant;
    name = "FisheyeEquidistant";
    
}

void prFisheyeEquidistant::project3DImage(prPointFeature & P)
{
		double X = P.get_X(), Y = P.get_Y(), Z = P.get_Z();
		double phi = acos(Z/sqrt(X*X+Y*Y+Z*Z)); // angle between Z-axis and [X,Y,Z]
		double theta = atan2(Y, X); //
		 
    P.set_x(phi*cos(theta));
    P.set_y(phi*sin(theta));
}

bool prFisheyeEquidistant::unProject(prPointFeature & P, double & Depth)
{
    //to do
    return false;
}

//idem prOmni
void prFisheyeEquidistant::project3DSphere(prPointFeature & P, double & Xs, double & Ys, double & Zs)
{
    double inv_norme;
    
    Xs = P.get_X();
    Ys = P.get_Y();
    Zs = P.get_Z();
    
    inv_norme = 1.0 / sqrt(Xs*Xs+Ys*Ys+Zs*Zs);
    
    Xs *= inv_norme;
    Ys *= inv_norme;
    Zs *= inv_norme;
}

bool prFisheyeEquidistant::projectImageSphere(prPointFeature & P, double & Xs, double & Ys, double & Zs)
{
	std::cout << "prFisheyeEquidistant::projectImageSphere not implemented" << std::endl;

    return -1;
		/*
    double x = P.get_x(), y = P.get_y(), fact;
    fact = (xi + sqrt(1.0+(1.0-xi*xi)*(x*x + y*y))) / (x*x + y*y + 1.0);
    
    Xs = fact * x;
    Ys = fact * y;
    Zs = fact - xi;
    */
}

prFisheyeEquidistant& prFisheyeEquidistant::operator=(const prFisheyeEquidistant& cam)
{
    prCameraModel::operator=(cam);
    return *this;
}

// du / dX
void prFisheyeEquidistant::computeSensorJacobian(prPointFeature & P, vpMatrix & LuX)
{
    //TODO
}
