#include <per/prStereoModel.h>

prStereoModel::prStereoModel(unsigned int _nbsens)
{
    init(_nbsens);
}

void prStereoModel::init(unsigned int _nbsens)
{
    nbsens = _nbsens;
    if(nbsens)
    {
        sen = new prSensorModel *[nbsens];
        sjMr = new vpHomogeneousMatrix[nbsens]; //attention sjMr[0] est inutilise
        sjMr[0].eye(); // au cas où...
    }
}

prStereoModel::~prStereoModel()
{
    if(nbsens)
    {
        // marche pô		delete [] sen; // ou delete sen;
        //		delete [] sjMr;
    }
}

void prStereoModel::setSensor(unsigned int i, prSensorModel* _sen)
{
    sen[i] = _sen;
}

void prStereoModel::setsjMr(unsigned int i, vpHomogeneousMatrix & M)
{
    if( (nbsens > 1) && (i > 0) && (i < nbsens) )
        sjMr[i] = M;
}

prStereoModel &prStereoModel::operator=(const prStereoModel &cam)
{
    this->nbsens = cam.nbsens;
    for (int i = 0; i < cam.nbsens; i++)
    {
        prSensorModel::operator=(*cam.sen[i]);
        this->sjMr[i] = cam.sjMr[i];
    }

    return *this;
}
