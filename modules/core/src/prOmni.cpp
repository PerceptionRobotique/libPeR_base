#include <per/prOmni.h>

prOmni::prOmni(double au,double av,double u0,double v0,double xii, double k1,double k2,double k3,double k4,double k5):prCameraModel(au,av,u0,v0, k1, k2, k3, k4, k5)
{
    //std::cout << "camera omni created" << std::endl;

    type = Omni;
    name = "Omni";
    xi = xii;
    
    //ATTENTION : mise en commentaire juste pour fixer xi a 1
    nbActiveParameters++;
    nbActiveParametersBase++;
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

bool prOmni::unProject(prPointFeature & P, double & Depth)
{
    double Xs, Ys, Zs;
    //to do
    if(!projectImageSphere(P, Xs, Ys, Zs))
        return false;

    P.set_X(Xs*Depth);
    P.set_Y(Ys*Depth);
    P.set_Z(Zs*Depth);
    
    return true;
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

bool prOmni::projectImageSphere(prPointFeature & P, double & Xs, double & Ys, double & Zs)
{
    double x = P.get_x(), y = P.get_y();
    double ra = 1.0+(1.0-xi*xi)*(x*x + y*y);
    if(ra < 0)
    {
        Xs = 0;
        Ys = 0;
        Zs = 1;

        return false;        
    }

    double fact = (xi + sqrt(ra)) / (x*x + y*y + 1.0);
    
    Xs = fact * x;
    Ys = fact * y;
    Zs = fact - xi;

    return true;
}

prOmni& prOmni::operator=(const prOmni& cam)
{
    prCameraModel::operator=(cam);
    this->xi = cam.xi;
    return *this;
}

// du / dX
void prOmni::computeSensorJacobian(prPointFeature & P, vpMatrix & LuX)
{
    double x = P.get_x(), y = P.get_y();
    double x2=x*x, y2=y*y;
    double xi2=xi*xi;
    double irho = 1./sqrt(pow(P.get_X(),2)+pow(P.get_Y(),2)+pow(P.get_Z(),2));
    double ra = 1+(1.-xi2)*(x2+y2);
    if(ra < 0)
    {
        LuX.resize(2,3,true);
    }
    else
    {
        double gamma = sqrt(ra);
        //get du / dx from mother class
        vpMatrix Lux(2,2);
        prCameraModel::computeSensorJacobian(P, Lux);

        vpMatrix LxX(2,3);
        LxX[0][0] = (gamma+xi*y2-xi2*x2*gamma)*irho/(1+xi*gamma); LxX[0][1] = -x*y*xi*irho;                                 LxX[0][2] = -x*gamma*irho;
        LxX[1][0] = LxX[0][1];                                    LxX[1][1] = (gamma+xi*x2-xi2*y2*gamma)*irho/(1+xi*gamma); LxX[1][2] = -y*gamma*irho;

        LuX = Lux*LxX;
    }
}
