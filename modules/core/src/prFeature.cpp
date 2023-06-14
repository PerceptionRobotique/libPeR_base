#include <per/prFeature.h>

prFeature::prFeature() : doubleFeature(0)
{
    
}

prFeature::~prFeature()
{
    if(doubleFeature != 0)
        delete [] doubleFeature;
}

