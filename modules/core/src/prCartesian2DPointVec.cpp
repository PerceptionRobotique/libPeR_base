#include <per/prCartesian2DPointVec.h>

prCartesian2DPointVec::prCartesian2DPointVec()
{
    p.resize(3);
}

prCartesian2DPointVec prCartesian2DPointVec::changeFrame(const vpHomogeneousMatrix & oMi)
{
    prCartesian2DPointVec ox;
    double x = get_x(), y = get_y(), w = get_w();
    
    ox.set(oMi[0][0]*x+ oMi[0][1]*y+ oMi[0][3]*w,
           oMi[1][0]*x+ oMi[1][1]*y+ oMi[1][3]*w,
           oMi[3][0]*x+ oMi[3][1]*y+ oMi[3][3]*w);

    ox.toEuclidean();
    
    return ox;
}

void prCartesian2DPointVec::setPoint(const double & x, const double & y)
{
    p[0] = x;
    p[1] = y;
    p[2] = 1.0;
}

void prCartesian2DPointVec::setVector(const double & x, const double & y)
{
    p[0] = x;
    p[1] = y;
    p[2] = 0.0;
}

void prCartesian2DPointVec::set(const double & x, const double & y, const double & w)
{
    p[0] = x;
    p[1] = y;
    p[2] = w;
}

int prCartesian2DPointVec::toEuclidean()
{
    double d = get_w();
   if(d == 0)
       return -1;
    
    d = 1.0/d;
    set(get_x()*d, get_y()*d, 1.0);
    
    return 0;
}
