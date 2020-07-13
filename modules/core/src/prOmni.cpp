#include <per/prOmni.h>

prOmni::prOmni(double au,double av,double u0,double v0,double xii, double k1,double k2,double k3,double k4,double k5):prCameraModel(au,av,u0,v0, k1, k2, k3, k4, k5)
{
    type = Omni;
    name = "Omni";
    xi = xii;
    
    //ATTENTION : mise en commentaire juste pour fixer xi a 1
    nbActiveParameters++;
    nbActiveParametersBase++;
}

prOmni::~prOmni()
{
}

void prOmni::init(double _au,double _av,double _u0,double _v0,double _xi, double _k1,double _k2,double _k3,double _k4,double _k5)
{
    prCameraModel::init(_au,_av,_u0,_v0,_k1,_k2,_k3, _k4, _k5);
    /*	au = _au;
     av = _av;
     u0 = _u0;
     v0 = _v0;
     inv_au = 1.0/au;
     inv_av = 1.0/av;*/
    
    type = Omni;
    name = "Omni";
    xi=_xi;
    
}

void prOmni::project3DImage(prPointFeature & P)
{
    double X = P.get_X(), Y = P.get_Y(), Z = P.get_Z(), den;
    den = Z + xi * sqrt(X*X+Y*Y+Z*Z);
    
    P.set_x(X/den);
    P.set_y(Y/den);
}

void prOmni::project3DSphere(prPointFeature & P, double & Xs, double & Ys, double & Zs)
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

int prOmni::projectImageSphere(prPointFeature & P, double & Xs, double & Ys, double & Zs)
{
    double x = P.get_x(), y = P.get_y(), fact;

    fact = 1.0+(1.0-xi*xi)*(x*x + y*y);
    if(fact >= 0.0)
    {
      fact = (xi + sqrt(fact)) / (x*x + y*y + 1.0);
    
      Xs = fact * x;
      Ys = fact * y;
      Zs = fact - xi;

      return 0;
    }
    return -1;
}

prOmni& prOmni::operator=(const prOmni& cam)
{
    prCameraModel::operator=(cam);
    this->xi = cam.xi;
    return *this;
}
