#include <per/prPointFeature.h>


prPointFeature::prPointFeature()
{

}

prPointFeature::~prPointFeature()
{
    
}


//! set the point world coordinates
void
prPointFeature::setWorldCoordinates(const double wX_,
                            const double wY_,
                            const double wZ_)
{
    wX.setPoint(wX_, wY_, wZ_);
}


int
prPointFeature::setWorldCoordinates(const prCartesian3DPointVec & _wX)
{
    if(_wX.get_W() ==  0)
        return -1;
    
    wX.set(_wX.get_X(), _wX.get_Y(), _wX.get_Z(), _wX.get_W());
    wX.toEuclidean();
    
    return 0;
}

void
prPointFeature::getWorldCoordinates(double& wX_,
                            double& wY_,
                            double& wZ_)
{
    wX_ = wX.get_X();
    wY_ = wX.get_Y();
    wZ_ = wX.get_Z();
}


void
prPointFeature::getWorldCoordinates(prCartesian3DPointVec & _wX)
{
    _wX.setPoint(wX.get_X(), wX.get_Y(), wX.get_Z());
}


prCartesian3DPointVec
prPointFeature::getWorldCoordinates(void)
{
    return this->wX;
}

//! Compute the new 3D coordinates of the point in the new camera frame.
void
prPointFeature::changeFrame(const vpHomogeneousMatrix &sMw, prCartesian3DPointVec &_sP)
{
    _sP = wX.changeFrame(sMw);
}

void prPointFeature::getImageMetric(double &x, double &y, double &w)const
{
    x = get_x();
    y = get_y();
    w = get_w();
}

void prPointFeature::setImageMetric(const double x,const double y, const double w)
{
    set_x(x);
    set_y(y);
    set_w(w);
}

void prPointFeature::getPixUV(double &ui, double &vi)const
{
    ui = get_u();
    vi = get_v();
}

void prPointFeature::setPixUV(const double u, const double v)
{
    set_u(u);
    set_v(v);
}

void prPointFeature::setObjectPixUV(const double oX, const double oY, const double oZ, const double u, const double v)
{
    set_oX(oX);
    set_oY(oY);
    set_oZ(oZ);
    setPixUV(u,v);
}

prPointFeature& prPointFeature::operator=(const prPointFeature& pt)
{
    set_u(pt.get_u());
    set_v(pt.get_v());
    set_oX(pt.get_oX());
    set_oY(pt.get_oY());
    set_oZ(pt.get_oZ());
    set_X(pt.get_X());
    set_Y(pt.get_Y());
    set_Z(pt.get_Z());
    set_x(pt.get_x());
    set_y(pt.get_y());
    
    return *this;
}

double *prPointFeature::toDouble(unsigned int place)
{
    return NULL;
}