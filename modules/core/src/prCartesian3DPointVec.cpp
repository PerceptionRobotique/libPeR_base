#include <per/prCartesian3DPointVec.h>

prCartesian3DPointVec::prCartesian3DPointVec(const double X, const double Y, const double Z, const double W)
{
    p.resize(4);
    
    set(X, Y, Z, W);
    toEuclidean();
}

prCartesian3DPointVec prCartesian3DPointVec::changeFrame(const vpHomogeneousMatrix & oMi)
{
    prCartesian3DPointVec oX;
    double X = get_X(), Y = get_Y(), Z = get_Z(), W = get_W();
    
    oX.set(oMi[0][0]*X+ oMi[0][1]*Y+ oMi[0][2]*Z+ oMi[0][3]*W,
           oMi[1][0]*X+ oMi[1][1]*Y+ oMi[1][2]*Z+ oMi[1][3]*W,
           oMi[2][0]*X+ oMi[2][1]*Y+ oMi[2][2]*Z+ oMi[2][3]*W,
           oMi[3][0]*X+ oMi[3][1]*Y+ oMi[3][2]*Z+ oMi[3][3]*W);

    oX.toEuclidean();
    
    return oX;
}

void prCartesian3DPointVec::setPoint(const double & X, const double & Y, const double & Z)
{
    p[0] = X;
    p[1] = Y;
    p[2] = Z;
    p[3] = 1;
}


/*void prCartesian3DPointVec::setPoint(double X, double Y, double Z)
{
    p[0] = X;
    p[1] = Y;
    p[2] = Z;
    p[3] = 1;
}*/

void prCartesian3DPointVec::setVector(const double & X, const double & Y, const double & Z)
{
    p[0] = X;
    p[1] = Y;
    p[2] = Z;
    p[3] = 0;
}

void prCartesian3DPointVec::set(const double & X, const double & Y, const double & Z, const double & W)
{
    p[0] = X;
    p[1] = Y;
    p[2] = Z;
    p[3] = W;
}

int prCartesian3DPointVec::toEuclidean()
{
    double d = get_W();
    if(d == 0)
        return -1;
    
    d = 1.0/d;
    set(get_X()*d, get_Y()*d, get_Z()*d, 1.0);
    
    return 0;
}
