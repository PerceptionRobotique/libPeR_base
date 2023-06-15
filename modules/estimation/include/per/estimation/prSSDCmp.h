/*!
 \file prSSDCmp.h
 \brief Definition of the SSD class feature sets comparator
 \author Guillaume CARON
 \version 0.1
 \date august 2017
 */
#if !defined(_PRSSDCMP_H)
#define _PRSSDCMP_H

#include <per/prFeaturesComparator.h>

#include <vector>

#include <visp/vpColVector.h>
#include <visp/vpRobust.h>

/*!
 \class PER_EXPORT prSSDCmp<T> prSSDCmp.h <per/prSSDCmp.h>
 \brief Class for comparing sets of features using the Sum of Squared Differences (SSD)
 */
template<typename Ts, class T>
class PER_EXPORT prSSDCmp : public prFeaturesComparator<T>
{
public:
    /*!
     * \fn prSSDCmp(prFeaturesSet<T> & fSet1, prFeaturesSet<T> & fSet2, bool _robustEstimation = false)
     * \brief Constructor of a prSSDCmp object that directly compare both provided sets of features (both must have the same number of features)
     * \param fSet1 the first set of features considered in the comparison
     * \param fSet2 the second set of features considered in the comparision
     * \param _robustEstimation boolean to set if a robust estimator must be considered to compute confidence weights simultaneously to the error computation
     */
    template< template<typename Taet> class Ta>
    prSSDCmp(prFeaturesSet<Ts, T, Ta> & fSet1, prFeaturesSet<Ts, T, Ta> & fSet2, bool _robustEstimation = false) : prFeaturesComparator<T>(), nbFeatures(0), robustEstimation(_robustEstimation)
    {
        compareFeatures(fSet1, fSet2);
    }
    
    /*!
     * \fn ~prSSDCmp()
     * \brief Default Destructor of a prSSDCmp object
     */
    virtual ~prSSDCmp() override
    {
        
    }
    
    /*!
     * \fn int compareFeatures(prFeaturesSet<T> & fSet1, prFeaturesSet<T> & fSet2)
     * \brief Implements the compareFeatures method: computes the difference between features of the first set and the ones of the second, feature per feature, at corresponding indexes. If robust estimation has been set, computes the confidence weights as well.
     * \param fSet1 the first set of features considered in the comparison
     * \param fSet2 the second set of features considered in the comparision
     * \return 0 if the feature sets might be compared
     */
    template< template<typename Taet> class Ta>
    int compareFeatures(prFeaturesSet<Ts, T, Ta> & fSet1, prFeaturesSet<Ts, T, Ta> & fSet2)
    {
        typename std::vector<T>::iterator it_fSet1 = fSet1.set.begin();
        typename std::vector<T>::iterator it_fSet2 = fSet2.set.begin();
        
        featuresComparison.clear();
        for(; it_fSet1 != fSet1.set.end() ; it_fSet1++, it_fSet2++)
        {
//            std::cout << it_fSet1->toDouble()[0] << std::endl;
//            std::cout << it_fSet2->toDouble()[0] << std::endl;
            featuresComparison.push_back(*it_fSet1-*it_fSet2);

//            std::cout << "nbFeat " << featuresComparison.size() << " / " << fSet1.set.size() << " | " << fSet1.set.size() << std::endl;
        }
        
        nbFeatures = featuresComparison.size();
        
        if(robustEstimation)
        {
            vpRobust robust(nbFeatures);
            
            featuresCmp2VpColVector();
            weights.resize(nbFeatures);
            weights=1; //imperatif pour TUKEY et HUBER !!!!
            robust.MEstimator(vpRobust::CAUCHY, residues, weights);
        }
        
        return 0;
    }
    
    /*!
     * \fn void featuresCmp2VpColVector()
     * \brief Data reorganization of the feature sets differences to a column vector.
     * \return Nothing
     */
    void featuresCmp2VpColVector()
    {
        if(nbFeatures != 0)
        {
            residues.resize(nbFeatures);
            unsigned long i = 0;
            typename std::vector<T>::iterator it_featuresComparison = featuresComparison.begin();
            double *pt_data = residues.data;
            for( ; i < nbFeatures ; i++, it_featuresComparison++, pt_data++)
                *pt_data = (it_featuresComparison->toDouble())[0]; //only GMS for this feature
        }
    }
    
    /*!
     * \fn T getCost()
     * \brief Accessor to get the SSD cost value: computes the SSD itself, as well, from the sets differences.
     * \return the SSD cost between the considered feature sets, in the considered feature type container
     */
    T getCost()
    {
        T cost;
        unsigned long i = 0;
        for(typename std::vector<T>::iterator it_fC = featuresComparison.begin();
            i < nbFeatures ; it_fC++, i++)
        {
            cost += (*it_fC)*(*it_fC);
        }
        
        return cost;
    }
    
    /*!
     * \fn T getRobustCost()
     * \brief Accessor to get the robust SSD cost value: computes the weighted SSD itself, as well, from the sets differences and weights got from the robust estimation, if set.
     * \return the robust SSD cost between the considered feature sets, in the considered feature type container
     */
    T getRobustCost()
    {
        if(robustEstimation)
        {
            T cost;
            unsigned long i = 0;
            double *pt_data = weights.data;
            for(typename std::vector<T>::iterator it_fC = featuresComparison.begin();
                i < nbFeatures ; it_fC++, i++, pt_data++)
            {
                cost += ((*it_fC)*(*pt_data))*((*it_fC)*(*pt_data));
            }
            
            return cost;
        }
        else
            return getCost();
    }
    
    /*!
     * \fn vpColVector getRobustWeigths()
     * \brief Accessor to the robust weights column vector.
     * \return the column vector containing weights if their computation were set
     *         a column vector of one element equal to -1 if the computation of weights were not set
     */
    vpColVector getRobustWeigths()
    {
        if(robustEstimation)
            return weights;
        else
        {
            vpColVector noRobust(1);
            noRobust = -1;
            return noRobust;
        }
    }
    
    /*!
     * \fn vpColVector getFeaturesCmpVpColVector()
     * \brief Accessor to the feature sets differences column vector.
     * \return the column vector containing the features sets differences
     */
    vpColVector getFeaturesCmpVpColVector()
    {
        if(!robustEstimation)
            featuresCmp2VpColVector();
        
        return residues;
    }
    
    
    //a implementer :
    //toErrorVector ?
//private:
    std::vector<T> featuresComparison;  /*!< the vector to store the feature sets differences */
    vpColVector residues;  /*!< the column vector storing the feature sets differences */
    vpColVector weights;  /*!< the column vector storing the feature sets differences weights */
    unsigned long nbFeatures;   /*!< number of features in the comparator */
    bool robustEstimation;   /*!< boolean to know is the computation of confidence weights is required (robustEstimation == true) or not (robustEstimation == false) */
};

#endif  //_PRSSDCMP_H
