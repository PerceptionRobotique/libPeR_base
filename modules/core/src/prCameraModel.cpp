#include <per/prCameraModel.h>

prCameraModel::prCameraModel(double m_au,double m_av,double m_u0,double m_v0, double m_k1, double m_k2, double m_k3, double m_k4, double m_k5)
{
    ik[0] = ik[1] = ik[2] = ik[3] = ik[4] = 0;
    
    init(m_au, m_av, m_u0, m_v0, m_k1, m_k2, m_k3, m_k4, m_k5);
    
    nbActiveParameters = nbActiveParametersBase = 4;
    
    if((m_k1!=0.0) || (m_k2!=0.0) || (m_k3!=0.0) || (m_k4!=0.0) || (m_k5!=0.0) )
    {
        distorsions = true;
        
        if(m_k1!=0.0)
        {
            k[0] = m_k1;
            activek[0] = true;
            nbActiveParameters++;
        }
        if(m_k2!=0.0)
        {
            k[1] = m_k2;
            activek[1] = true;
            nbActiveParameters++;
        }
        if(m_k3!=0.0)
        {
            k[2] = m_k3;
            activek[2] = true;
            nbActiveParameters++;
        }
        if(m_k4!=0.0)
        {
            k[3] = m_k4;
            activek[3] = true;
            nbActiveParameters++;
        }
        if(m_k5!=0.0)
        {
            k[4] = m_k5;
            activek[4] = true;
            nbActiveParameters++;
        }
    }
    else {
        distorsions = false;
        
        for(int i = 0 ; i < 5 ; i++)
        {
            k[i] = 0.0 ;
            activek[i] = false;
        }
    }
}

void prCameraModel::init(double m_au,double m_av,double m_u0,double m_v0, double m_k1, double m_k2, double m_k3, double m_k4, double m_k5)
{
    name = "Perspective";
    au = m_au;
    av = m_av;
    u0 = m_u0;
    v0 = m_v0;
    inv_au = 1.0/au;
    inv_av = 1.0/av;
    
    k[0] = m_k1;
    k[1] = m_k2;
    k[2] = m_k3;
    k[3] = m_k4;
    k[4] = m_k5;
    
}

void prCameraModel::init(const prCameraModel &c)
{
    name = "Perspective";
    *this = c;
}

void prCameraModel::meterPixelConversion(prPointFeature & P)
{
    if(distorsions)
    {
        double x = P.get_x(), y = P.get_y(), r2, distR, xy, xd, yd;
        r2 = x*x+y*y;
        distR = 1;
        xy = x*y;
        if(activek[0])
            distR += k[0]*r2;
        if(activek[1])
            distR += k[1]*r2*r2; //r4
        if(activek[2])
            distR += k[2]*r2*r2*r2; //r6
        
        xd = x*distR;
        yd = y*distR;
        
        if(activek[3])
        {
            xd += 2*k[3]*xy;
            yd += k[3]*(r2 + 2*y*y);
        }
        if(activek[4])
        {
            xd += k[4]*(r2 + 2*x*x);
            yd += 2*k[4]*xy;
        }
        
        P.set_u(xd * au + u0);
        P.set_v(yd * av + v0);
    }
    else
    {
        P.set_u(P.get_x() * au + u0);
        P.set_v(P.get_y() * av + v0);
    }
}

void prCameraModel::pixelMeterConversion(prPointFeature & P)
{
    if(distorsions)
    {
        /*
         //!!!!!! Attention : temporaire !!!!!!
         // il faudra estimer les paramètres de
         // distorsion distorted->undistorted
         // pour pouvoir faire cette opération
         //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         P.set_x((P.get_u() - u0) * inv_au);
         P.set_y((P.get_v() - v0) * inv_av);*/
        
        //!!!!!! Attention : temporaire !!!!!!
        // il faudra estimer les paramètres de
        // distorsion distorted->undistorted
        // pour pouvoir faire cette opération
        //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        double rd2, idistR, xdyd, xd, yd, x, y;
        xd = (P.get_u() - u0) * inv_au;
        yd = (P.get_v() - v0) * inv_av;
        rd2 = xd*xd+yd*yd;
        idistR = 1;
        xdyd = xd*yd;
        if(activek[0])
            idistR += ik[0]*rd2;
        if(activek[1])
            idistR += ik[1]*rd2*rd2; //r4
        if(activek[2])
            idistR += ik[2]*rd2*rd2*rd2; //r6
        
        x = xd*idistR;
        y = yd*idistR;
        
        if(activek[3])
        {
            x += 2*ik[3]*xdyd;
            y += ik[3]*(rd2 + 2*yd*yd);
        }
        if(activek[4])
        {
            x += ik[4]*(rd2 + 2*xd*xd);
            y += 2*ik[4]*xdyd;
        }
        
        P.set_x(x);
        P.set_y(y);
    }
    else
    {
        P.set_x((P.get_u() - u0) * inv_au);
        P.set_y((P.get_v() - v0) * inv_av);
    }
    
}

void prCameraModel::setPrincipalPoint(double u0, double v0)
{
    this->u0    = u0 ;
    this->v0    = v0 ;
}

void prCameraModel::setPixelRatio(double au, double av)
{
    this->au    = au ;
    this->av    = av ;
    
    this->inv_au = 1./au;
    this->inv_av = 1./av;
}

void prCameraModel::setDistorsionParameters(double k1,double k2,double k3,double k4,double k5)
{
    this->k[0] = k1;
    this->k[1] = k2;
    this->k[2] = k3;
    this->k[3] = k4;
    this->k[4] = k5;
}

void prCameraModel::setUndistorsionParameters(double ik1,double ik2,double ik3,double ik4,double ik5)
{
    this->ik[0] = ik1;
    this->ik[1] = ik2;
    this->ik[2] = ik3;
    this->ik[3] = ik4;
    this->ik[4] = ik5;
}

void prCameraModel::setActiveDistorsionParameters(bool k1,bool k2,bool k3,bool k4,bool k5)
{
    nbActiveParameters = nbActiveParametersBase;
    distorsions = false;
    
    activek[0] = k1;
    activek[1] = k2;
    activek[2] = k3;
    activek[3] = k4;
    activek[4] = k5;
    
    for(int i = 0 ; i < 5 ; i++)
        if(activek[i])
        {
            nbActiveParameters++;
            distorsions = true;
        }
}

prCameraModel& prCameraModel::operator=(const prCameraModel& cam)
{
    au = cam.au ;
    av = cam.av ;
    u0 = cam.u0 ;
    v0 = cam.v0 ;
    
    inv_au = cam.inv_au;
    inv_av = cam.inv_av;
    
    distorsions = cam.distorsions;
    nbActiveParameters = cam.nbActiveParameters;
    nbActiveParametersBase = cam.nbActiveParametersBase;
    for(int i = 0 ; i < 5 ; i++)
    {
        k[i] = cam.k[i] ;
        ik[i] = cam.ik[i] ;
        activek[i] = cam.activek[i];
    }
    
    return *this;
}

std::ostream& prCameraModel::operator<<(std::ostream & os)
{
    //os << "Camera parameters for generic model :" << std::endl ;
    os << "  au = " << getau() <<"\t av = "<< getav() << std::endl ;
    os << "  u0 = " << getu0() <<"\t v0 = "<< getv0() << std::endl ;
    
    if(distorsions)
    {
        os << "with distorsions parameters for undistorted to distorted transformation :" << std::endl ;
        if(activek[0])
            os << "  k1 = " << getk1();
        if(activek[1])
            os <<"\t k2 = "<< getk2();
        if(activek[2])
            os << "\t  k3 = " << getk3();
        if(activek[3])
            os <<"\t k4 = "<< getk4();
        if(activek[4])
            os << "\t  k5 = " << getk5() << std::endl ;
        os << "with undistorsions parameters for distorted to undistorted transformation :" << std::endl ;
        if(activek[0])
            os << "  ik1 = " << getik1();
        if(activek[1])
            os <<"\t ik2 = "<< getik2();
        if(activek[2])
            os << "\t  ik3 = " << getik3();
        if(activek[3])
            os <<"\t ik4 = "<< getik4();
        if(activek[4])
            os << "\t  ik5 = " << getik5() << std::endl ;
    }
    
    return os;
}

vpMatrix prCameraModel::getK() const
{
    vpMatrix K;
    
    K.resize(3,3) ;
    K = 0.0 ;
    K[0][0] = au ;
    K[1][1] = av ;
    K[0][2] = u0 ;
    K[1][2] = v0 ;
    K[2][2] = 1.0 ;
    
    return K; 
}
