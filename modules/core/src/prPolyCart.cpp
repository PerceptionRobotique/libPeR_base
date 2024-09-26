#include <per/prPolyCart.h>

prPolyCart::prPolyCart(double au,double u0,double v0,double *aa, double k1,double k2,double k3,double k4,double k5):prCameraModel(au,au,u0,v0, k1, k2, k3, k4, k5)
{
    type = PolyCart;
    name = "PolyCart";
    
    if(aa != NULL)
    {
    	for(unsigned int i = 0 ; i < 5 ; i++)
    		a[i] = aa[i];
    }
    else
    {
    	for(unsigned int i = 0 ; i < 5 ; i++)
    		a[i] = 0;
    }
    
    //ATTENTION : mise en commentaire juste pour fixer xi a 1
    nbActiveParameters+=4; //a_1 is not optimized since assumed = 0 always
    nbActiveParametersBase+=4;
}

void prPolyCart::init(double _au,double _u0,double _v0,double *_a, double _k1,double _k2,double _k3,double _k4,double _k5)
{
    prCameraModel::init(_au,_au,_u0,_v0,_k1,_k2,_k3, _k4, _k5);
    
    type = PolyCart;
    name = "PolyCart";
    
    if(_a != NULL)
    {
    	for(unsigned int i = 0 ; i < 5 ; i++)
    		a[i] = _a[i];
    }
    else
    {
    	for(unsigned int i = 0 ; i < 5 ; i++)
    		a[i] = 0;
    }    
}

void prPolyCart::project3DSphere(prPointFeature & P, double & Xs, double & Ys, double & Zs)
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

void prPolyCart::project3DImage(prPointFeature & P)
{
	std::cout << "prPolyCart::project3DImage(prPointFeature & P) not implemented" << std::endl;
	/*
    double X = P.get_X(), Y = P.get_Y(), Z = P.get_Z(), den;
    den = Z + xi * sqrt(X*X+Y*Y+Z*Z);
    
    P.set_x(X/den);
    P.set_y(Y/den);
   */
}

bool prPolyCart::unProject(prPointFeature & P, double & Depth)
{
    return false;
    //TODO
    //return true;
}

bool prPolyCart::projectImageSphere(prPointFeature & P, double & Xs, double & Ys, double & Zs)
{
    double up = P.get_u()-u0, vp = P.get_v()-v0;
    
    double rho_up = sqrt(up*up + vp*vp);
    double rho_up2 = rho_up*rho_up;
    double rho_up3 = rho_up2*rho_up;
    double rho_up4 = rho_up3*rho_up;
        
    double r_rho_up = a[0] + a[1]*rho_up + a[2]*rho_up2 + a[3]*rho_up3 + a[4]*rho_up4;
        
    Xs = up / au;
    Ys = vp / au;
    Zs = r_rho_up / au;
    
    return true;
}

prPolyCart& prPolyCart::operator=(const prPolyCart& cam)
{
    prCameraModel::operator=(cam);
    for(unsigned int i = 0 ; i < 5 ; i++)
    	this->a[i] = cam.a[i];   
    	
    return *this;
}

// du / dX
void prPolyCart::computeSensorJacobian(prPointFeature & P, vpMatrix & LuX)
{
    //TODO
}
