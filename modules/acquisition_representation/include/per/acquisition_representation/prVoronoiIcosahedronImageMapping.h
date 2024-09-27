/*!
 \file prVoronoiIcosahedronImageMapping.h
 \brief Header file for the prVoronoiIcosahedronImageMapping class
 \author Antoine ANDRE
 \version 0.1
 \date may 2023
 */

#if !defined(_PRVORONOIICOSAHEDRONIMAGEMAPPING_H)
#define _PRVORONOIICOSAHEDRONIMAGEMAPPING_H

#include <per/prAcquisitionModel.h>
#include <per/prCartesian3DPointVec.h>

#include <per/prEquirectangular.h>
#include <per/prOmni.h>
#include <per/prPointFeature.h>
#include <per/prStereoModel.h>

// #include <per/prcommon.h>

#include <per/prIcosahedronSamples.h>
#include <per/prRegularlySampledCSImage.h>
// #include </usr/include/eigen3/Eigen/Dense>
#include <Eigen/Dense>

#include <opencv2/core/core.hpp>
#include <opencv2/core/types.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>

#include <visp/vpDisplayX.h>
#include <visp/vpImage.h>
#include <visp/vpImageIo.h>
#include <visp/vpImageTools.h>
#include <visp/vpTime.h>

template <typename T>
class PER_EXPORT prVoronoiIcosahedronImageMapping : public prRegularlySampledCSImage<T> {
public:
  prVoronoiIcosahedronImageMapping(unsigned int _subdivLevels = 3) {
    this->subdivLevels = _subdivLevels;
    this->isCellMapInit = false;
    loadSphere();
  }

  int loadSphere() {
    prIcosahedronSamples icoSampler;

    if (this->subdivLevels > icoSampler.maxSubdivLevels)
      return -1;

    icoSampler.loadVoronoiSamples(this->subdivLevels);
    this->nbSamples = icoSampler.nbVertices;
    this->bitmap = new T[this->nbSamples];

    this->ge = new prCartesian3DPointVec[this->nbSamples];

    // copy samples coordinates
    prCartesian3DPointVec *pt_XS = (prCartesian3DPointVec *)this->ge;
    const prCartesian3DPointVec *pt_XS_src = icoSampler.samplesCoords;
    voroSamples = icoSampler.samplesVoronoiMap;
    float inorme = 1.0f;
    for (unsigned long i = 0; i < this->nbSamples; i++, pt_XS++, pt_XS_src++) {
      // renormalize coordinates (might be necessary due to the text exportation
      inorme = 1.0f / sqrt(pt_XS_src->get_X() * pt_XS_src->get_X() +
                           pt_XS_src->get_Y() * pt_XS_src->get_Y() +
                           pt_XS_src->get_Z() * pt_XS_src->get_Z());
      pt_XS->set(pt_XS_src->get_X() * inorme, pt_XS_src->get_Y() * inorme,
                 pt_XS_src->get_Z() * inorme, pt_XS_src->get_W());
    }

    return 0;
  }

  ~prVoronoiIcosahedronImageMapping() {}

  int createImageCellsMap(prStereoModel &stereoCam) {
    vpHomogeneousMatrix c2Rc1 = stereoCam.sjMr[1];
    prCartesian3DPointVec XSs;
    prCartesian3DPointVec *pt_XS_s = (prCartesian3DPointVec *)this->ge;
    prPointFeature P, Ps;

    double epsilon = 1e-6;
    unsigned int icam = 0, i, j;
    double u, v;

    centerCoords = Eigen::ArrayXXd(this->nbSamples, 2);

    for (unsigned long ns = 0; ns < this->nbSamples; ns++, pt_XS_s++) {
      Ps.sX = *pt_XS_s;

      if (Ps.get_Z() > 0) {
        icam = 0;
      } else {
        icam = 1;
        Ps.sX = Ps.sX.changeFrame(c2Rc1);
      }

      ((prOmni *)(stereoCam.sen[icam]))->project3DImage(Ps);

      ((prOmni *)(stereoCam.sen[icam]))->meterPixelConversion(Ps);

      int uCenter = Ps.get_u();
      int vCenter = Ps.get_v();
      centerCoords(ns, 0) = uCenter;
      centerCoords(ns, 1) = vCenter;

      int nbSommets = voroSamples[ns].size();

      bool icam0 = true;
      bool icam1 = true;
      int icamCenter = icam;

      Eigen::MatrixXd voroContour(0, 2);

      for (j = 0; j < nbSommets; j++) {
        prCartesian3DPointVec ptXs = voroSamples[ns][j];
        P.sX = ptXs;

        if (P.get_Z() > 0) {
          icam = 0;
        } else {
          icam = 1;
          P.sX = P.sX.changeFrame(c2Rc1);
        }

        ((prOmni *)(stereoCam.sen[icam]))->project3DImage(P);

        ((prOmni *)(stereoCam.sen[icam]))->meterPixelConversion(P);

        u = P.get_u();
        v = P.get_v();

        if (icamCenter == icam) {
          voroContour.conservativeResize(voroContour.rows() + 1, 2);
          voroContour(voroContour.rows() - 1, 0) = u;
          voroContour(voroContour.rows() - 1, 1) = v;
        }
      }

      int bboxMinX = (int)floor(voroContour.col(0).minCoeff());
      int bboxMinY = (int)floor(voroContour.col(1).minCoeff());
      int bboxMaxX = (int)ceil(voroContour.col(0).maxCoeff());
      int bboxMaxY = (int)ceil(voroContour.col(1).maxCoeff());

      double area = 0;
      double correctArea = 0;
      Eigen::Matrix3d trglTemp;

      if (voroContour.rows() >= 5) {
        for (int c = 0; c < voroContour.rows() - 1; c++) {
          trglTemp << voroContour(c, 0), voroContour(c, 1), 1,
              voroContour(c + 1, 0), voroContour(c + 1, 1), 1, uCenter, vCenter,
              1;
          correctArea += 0.5 * abs(trglTemp.determinant());
        }
        trglTemp << voroContour(0, 0), voroContour(0, 1), 1,
            voroContour(voroContour.rows() - 1, 0),
            voroContour(voroContour.rows() - 1, 1), 1, uCenter, vCenter, 1;
        correctArea += 0.5 * abs(trglTemp.determinant());
      } else if (voroContour.rows() == 3) {
        trglTemp << voroContour(0, 0), voroContour(0, 1), 1, voroContour(1, 0),
            voroContour(1, 1), 1, voroContour(2, 0), voroContour(2, 1), 1;
        correctArea += 0.5 * abs(trglTemp.determinant());
      } else if (voroContour.rows() == 4) {
        trglTemp << voroContour(0, 0), voroContour(0, 1), 1, voroContour(1, 0),
            voroContour(1, 1), 1, voroContour(2, 0), voroContour(2, 1), 1;
        correctArea += 0.5 * abs(trglTemp.determinant());

        trglTemp << voroContour(0, 0), voroContour(0, 1), 1, voroContour(2, 0),
            voroContour(2, 1), 1, voroContour(3, 0), voroContour(3, 1), 1;
        correctArea += 0.5 * abs(trglTemp.determinant());
      }

      Eigen::ArrayXXd pointsInVoro(0, 2);
      for (int row = bboxMinY; row <= bboxMaxY; row++) {
        for (int col = bboxMinX; col <= bboxMaxX; col++) {
          area = 0;
          for (int c = 0; c < voroContour.rows() - 1; c++) {
            trglTemp << voroContour(c, 0), voroContour(c, 1), 1,
                voroContour(c + 1, 0), voroContour(c + 1, 1), 1, col, row, 1;
            area += 0.5 * abs(trglTemp.determinant());
          }

          trglTemp << voroContour(0, 0), voroContour(0, 1), 1,
              voroContour(voroContour.rows() - 1, 0),
              voroContour(voroContour.rows() - 1, 1), 1, col, row, 1;
          area += 0.5 * abs(trglTemp.determinant());
          if (fabs(area - correctArea) < epsilon) {
            pointsInVoro.conservativeResize(pointsInVoro.rows() + 1, 2);
            pointsInVoro(pointsInVoro.rows() - 1, 0) = row;
            pointsInVoro(pointsInVoro.rows() - 1, 1) = col;
          }
        }
      }

      verticeHexaVec.push_back(pointsInVoro);
      vecVoroContours.push_back(voroContour);
    }

    isCellMapInit = true;
    return 0;
  }

  int createImageCellsMap(prEquirectangular &equiCam) {

    prCartesian3DPointVec XSs;
    prCartesian3DPointVec *pt_XS_s = (prCartesian3DPointVec *)this->ge;
    unsigned int icam = 0, i, j;
    prPointFeature P, Ps;
    double epsilon = 1e-6;
    double u, v, du, dv, unmdu, unmdv;

    // cv::Mat fullBlackImg(imHeight, imWidth, CV_8UC3);
    // fullBlackImg.setTo(0);
    centerCoords = Eigen::ArrayXXd(0, 2);

    // #pragma omp parallel for
    for (unsigned long ns = 0; ns < this->nbSamples; ns++) {
      std::vector<Eigen::MatrixXd> vecTempVoroContours;
      Ps.sX = *pt_XS_s;

      equiCam.project3DImage(Ps);

      equiCam.meterPixelConversion(Ps);

      int uCenter = Ps.get_u();
      int vCenter = Ps.get_v();
      centerCoords.conservativeResize(centerCoords.rows() + 1, 2);
      centerCoords(centerCoords.rows() - 1, 0) = uCenter;
      centerCoords(centerCoords.rows() - 1, 1) = vCenter;

      int nbSommets = voroSamples[ns].size();

      bool icam0 = true;
      bool icam1 = true;
      int icamCenter = icam;

      Eigen::MatrixXd voroContour(0, 2);

      bool pushVoro = true;

      for (j = 0; j < nbSommets; j++) {
        prCartesian3DPointVec ptXs = voroSamples[ns][j];
        P.sX = ptXs;

        equiCam.project3DImage(P);

        equiCam.meterPixelConversion(P);

        u = P.get_u();
        v = P.get_v();

        voroContour.conservativeResize(voroContour.rows() + 1, 2);
        voroContour(voroContour.rows() - 1, 0) = u;
        voroContour(voroContour.rows() - 1, 1) = v;
      }
      for (j = 0; j < voroContour.rows() - 1; j++) {
        if (voroContour(j, 0) < imWidth / 4 &&
            voroContour(j + 1, 0) > 3 * imWidth / 4) {
          pushVoro = false;
        }
        if (voroContour(j, 0) > 3 * imWidth / 4 &&
            voroContour(j + 1, 0) < imWidth / 4) {
          pushVoro = false;
        }
      }
      if (voroContour(0, 0) > 3 * imWidth / 4 &&
          voroContour(voroContour.rows() - 1, 0) < imWidth / 4) {
        pushVoro = false;
      }

      if (pushVoro) {
        vecTempVoroContours.push_back(voroContour);
        vecVoroContours.push_back(voroContour);
      } else {

        Eigen::MatrixXd voroContourLeft(0, 2);
        for (j = 0; j < nbSommets; j++) {
          if (voroContour(j, 0) <= imWidth / 2) {
            voroContourLeft.conservativeResize(voroContourLeft.rows() + 1, 2);
            voroContourLeft(voroContourLeft.rows() - 1, 0) =
                voroContour(voroContourLeft.rows() - 1, 0);
            voroContourLeft(voroContourLeft.rows() - 1, 1) =
                voroContour(voroContourLeft.rows() - 1, 1);
          } else {
            voroContourLeft.conservativeResize(voroContourLeft.rows() + 1, 2);
            voroContourLeft(voroContourLeft.rows() - 1, 0) = 0;
            voroContourLeft(voroContourLeft.rows() - 1, 1) =
                voroContour(voroContourLeft.rows() - 1, 1);
          }
        }

        Eigen::MatrixXd voroContourRight(0, 2);
        for (j = 0; j < nbSommets; j++) {
          if (voroContour(j, 0) >= imWidth / 2) {
            voroContourRight.conservativeResize(voroContourRight.rows() + 1, 2);
            voroContourRight(voroContourRight.rows() - 1, 0) =
                voroContour(voroContourRight.rows() - 1, 0);
            voroContourRight(voroContourRight.rows() - 1, 1) =
                voroContour(voroContourRight.rows() - 1, 1);
          } else {
            voroContourRight.conservativeResize(voroContourRight.rows() + 1, 2);
            voroContourRight(voroContourRight.rows() - 1, 0) = imWidth - 1;
            voroContourRight(voroContourRight.rows() - 1, 1) =
                voroContour(voroContourRight.rows() - 1, 1);
          }
        }

        // vecVoroContours.push_back(voroContourLeft);
        vecVoroContours.push_back(voroContourRight);
        vecTempVoroContours.push_back(voroContourLeft);
        vecTempVoroContours.push_back(voroContourRight);
      }
      vecTotVoroContours.push_back(vecTempVoroContours);
      pt_XS_s++;
    }

    // #pragma omp parallel for
    for (int i = 0; i < vecTotVoroContours.size(); i++) {
      Eigen::ArrayXXd pointsInVoro(0, 2);
      for (int j = 0; j < vecTotVoroContours[i].size(); j++) {

        Eigen::MatrixXd voroContour = vecTotVoroContours[i][j];
        int uCenter, vCenter;
        if (vecTotVoroContours[i].size() == 1) {
          uCenter = centerCoords(i, 0);
          vCenter = centerCoords(i, 1);
        } else {
          uCenter = voroContour.col(0).mean();
          vCenter = voroContour.col(1).mean();
        }

        // std::cout << voroContour << std::endl;

        int bboxMinX = (int)floor(voroContour.col(0).minCoeff());
        int bboxMinY = (int)floor(voroContour.col(1).minCoeff());
        int bboxMaxX = (int)ceil(voroContour.col(0).maxCoeff());
        int bboxMaxY = (int)ceil(voroContour.col(1).maxCoeff());

        double area = 0;
        double correctArea = 0;
        Eigen::Matrix3d trglTemp;

        if (voroContour.rows() >= 5) {
          for (int c = 0; c < voroContour.rows() - 1; c++) {
            trglTemp << voroContour(c, 0), voroContour(c, 1), 1,
                voroContour(c + 1, 0), voroContour(c + 1, 1), 1, uCenter,
                vCenter, 1;
            correctArea += 0.5 * abs(trglTemp.determinant());
          }
          trglTemp << voroContour(0, 0), voroContour(0, 1), 1,
              voroContour(voroContour.rows() - 1, 0),
              voroContour(voroContour.rows() - 1, 1), 1, uCenter, vCenter, 1;
          correctArea += 0.5 * abs(trglTemp.determinant());
        } else if (voroContour.rows() == 3) {
          trglTemp << voroContour(0, 0), voroContour(0, 1), 1,
              voroContour(1, 0), voroContour(1, 1), 1, voroContour(2, 0),
              voroContour(2, 1), 1;
          correctArea += 0.5 * abs(trglTemp.determinant());
        } else if (voroContour.rows() == 4) {
          trglTemp << voroContour(0, 0), voroContour(0, 1), 1,
              voroContour(1, 0), voroContour(1, 1), 1, voroContour(2, 0),
              voroContour(2, 1), 1;
          correctArea += 0.5 * abs(trglTemp.determinant());

          trglTemp << voroContour(0, 0), voroContour(0, 1), 1,
              voroContour(2, 0), voroContour(2, 1), 1, voroContour(3, 0),
              voroContour(3, 1), 1;
          correctArea += 0.5 * abs(trglTemp.determinant());
        }

        // std::cout << "---" << std::endl;

        for (int row = bboxMinY; row <= bboxMaxY; row++) {
          for (int col = bboxMinX; col <= bboxMaxX; col++) {
            if (row >= 0 && col >= 0 && row < imHeight - 1 &&
                col < imWidth - 1) {
              area = 0;
              for (int c = 0; c < voroContour.rows() - 1; c++) {
                trglTemp << voroContour(c, 0), voroContour(c, 1), 1,
                    voroContour(c + 1, 0), voroContour(c + 1, 1), 1, col, row,
                    1;
                area += 0.5 * abs(trglTemp.determinant());
              }
              trglTemp << voroContour(0, 0), voroContour(0, 1), 1,
                  voroContour(voroContour.rows() - 1, 0),
                  voroContour(voroContour.rows() - 1, 1), 1, col, row, 1;
              area += 0.5 * abs(trglTemp.determinant());
              if (fabs(area - correctArea) < epsilon) {
                pointsInVoro.conservativeResize(pointsInVoro.rows() + 1, 2);
                pointsInVoro(pointsInVoro.rows() - 1, 0) = row;
                pointsInVoro(pointsInVoro.rows() - 1, 1) = col;
              }
            }
          }
        }
      }
      verticeHexaVec.push_back(pointsInVoro);

      // std::cout << i << std::endl;
    }

    isCellMapInit = true;
    return 0;
  }

  int buildFromTwinOmni(vpImage<T> &I, prStereoModel &stereoCam,
                        vpImage<unsigned char> *Mask = NULL) {
    this->I_req = I;

    if (this->nbSamples == 0)
      return -1;

    if (this->bitmapf != NULL) {
      delete[] this->bitmapf;
      this->bitmapf = NULL;
    }

    this->bitmapf = new float[this->nbSamples];
    // this->bitmapfUpdated = new float[this->nbSamples];
    this->isBitmapfSet = true;

    T *pt_bitmap = this->bitmap;
    float *pt_bitmapf = this->bitmapf;
    // float *pt_bitmapf_updated = this->bitmapfUpdated;
    double u, v, du, dv, unmdu, unmdv;
    unsigned int i, j;

    imWidth = I.getWidth();
    imHeight = I.getHeight();

    if (!isCellMapInit || verticeHexaVec.empty()) {
      createImageCellsMap(stereoCam);
    }

    sumTotColors = 0;
    double cumulColor = 0;
    // #pragma omp parallel
    for (unsigned long ns = 0; ns < this->nbSamples; ns++) {
      *pt_bitmap = 0;
      *pt_bitmapf = 0.;
      cumulColor = 0;
      if (verticeHexaVec[ns].rows() > 1) {
        // #pragma omp parallel for reduction(+ : cumulColor)
        for (j = 0; j < verticeHexaVec[ns].rows(); j++) {
          int uCumul = verticeHexaVec[ns](j, 0);
          int vCumul = verticeHexaVec[ns](j, 1);
          cumulColor += I[uCumul][vCumul];
        }

        cumulColor = (double)cumulColor / (double)verticeHexaVec[ns].rows();
      }

      *pt_bitmap = (int)cumulColor;
      *pt_bitmapf = (float)cumulColor;
      // *pt_bitmapf_updated = (float)cumulColor;
      sumTotColors += (float)cumulColor;
      pt_bitmap++;
      pt_bitmapf++;
      // pt_bitmapf_updated++;
    }

    return 0;
  }

  int buildFromTwinOmni(cv::Mat &I, prStereoModel &stereoCam,
                        vpImage<unsigned char> *Mask = NULL) {
    // this->I_req = I;

    if (this->nbSamples == 0)
      return -1;

    if (this->bitmapf != NULL) {
      delete[] this->bitmapf;
      this->bitmapf = NULL;
    }

    this->bitmapf = new float[this->nbSamples];
    this->bitmapfUpdated = new float[this->nbSamples];
    this->isBitmapfSet = true;

    T *pt_bitmap = this->bitmap;
    float *pt_bitmapf = this->bitmapf;
    float *pt_bitmapf_updated = this->bitmapfUpdated;
    double u, v, du, dv, unmdu, unmdv;
    unsigned int i, j;

    imWidth = I.cols;
    imHeight = I.rows;

    if (!isCellMapInit || verticeHexaVec.empty()) {
      createImageCellsMap(stereoCam);
    }

    sumTotColors = 0;
    double cumulColor = 0;
    uint8_t pixelVal;
    uint8_t *pixelPtrIn = (uint8_t *)I.data;

    // #pragma omp parallel
#pragma omp parallel for
    for (unsigned long ns = 0; ns < this->nbSamples; ns++) {
      *pt_bitmap = 0;
      *pt_bitmapf = 0.;
      cumulColor = 0;
      if (verticeHexaVec[ns].rows() > 1) {
        // #pragma omp parallel for reduction(+ : cumulColor)
        for (j = 0; j < verticeHexaVec[ns].rows(); j++) {
          int uCumul = verticeHexaVec[ns](j, 0);
          int vCumul = verticeHexaVec[ns](j, 1);

          cumulColor += pixelPtrIn[uCumul * imWidth + vCumul];

          // cumulColor += I.at<unsigned char>(uCumul, vCumul);
        }

        cumulColor = (double)cumulColor / (double)verticeHexaVec[ns].rows();
      }

      *pt_bitmap = (int)cumulColor;
      *pt_bitmapf = (float)cumulColor;
      *pt_bitmapf_updated = (float)cumulColor;
      sumTotColors += (float)cumulColor;
      pt_bitmap++;
      pt_bitmapf++;
      pt_bitmapf_updated++;
    }

    return 0;
  }

  int toTwinOmni(vpImage<T> &I_r, vpPoseVector &r, prStereoModel &stereoCam,
                 vpImage<unsigned char> *Mask = NULL) {
    if (this->nbSamples == 0)
      return -1;

    T *pt_bitmap = this->bitmap;
    float *pt_bitmapf = this->bitmapf;
    vpHomogeneousMatrix c2Rc1 = stereoCam.sjMr[1];
    prCartesian3DPointVec XSs;
    prCartesian3DPointVec *pt_XS = (prCartesian3DPointVec *)this->ge;
    int j;
    prPointFeature P;

    vpHomogeneousMatrix dMc;
    dMc.buildFrom(r);
    dMc = dMc.inverse(); // a justifier clairement dans la doc
    Eigen::Vector3d vecAB, vecAC, vecBC, vecAP, vecBP, vecCP;
    double x_a, y_a, z_a, x_b, y_b, z_b, x_c, y_c, z_c, x_p, y_p, z_p;
    double areaAPC, areaAPB, areaBPC, areaTot;
    double wa, wb, wc;
    double weightedInt;
    cv::Mat I_out_CV(I_r.getHeight(), I_r.getWidth(), CV_8UC1);
    I_out_CV.setTo(0);

    for (unsigned long ns = 0; ns < this->nbSamples; ns++, pt_XS++) {
      P.sX = *pt_XS;
      P.sX = P.sX.changeFrame(dMc);
      // int localIndex = getMinDistance3D(P.sX);
      std::vector<unsigned int> closestPoint = this->get3ClosestPoints(P.sX);

      x_a = ((prCartesian3DPointVec *)this->ge)[closestPoint[0]].get_X();
      y_a = ((prCartesian3DPointVec *)this->ge)[closestPoint[0]].get_Y();
      z_a = ((prCartesian3DPointVec *)this->ge)[closestPoint[0]].get_Z();

      x_b = ((prCartesian3DPointVec *)this->ge)[closestPoint[1]].get_X();
      y_b = ((prCartesian3DPointVec *)this->ge)[closestPoint[1]].get_Y();
      z_b = ((prCartesian3DPointVec *)this->ge)[closestPoint[1]].get_Z();

      x_c = ((prCartesian3DPointVec *)this->ge)[closestPoint[2]].get_X();
      y_c = ((prCartesian3DPointVec *)this->ge)[closestPoint[2]].get_Y();
      z_c = ((prCartesian3DPointVec *)this->ge)[closestPoint[2]].get_Z();

      x_p = P.sX.get_X();
      y_p = P.sX.get_Y();
      z_p = P.sX.get_Z();

      vecAB << x_b - x_a, y_b - y_a, z_b - z_a;
      vecAC << x_c - x_a, y_c - y_a, z_c - z_a;
      vecBC << x_c - x_b, y_c - y_b, z_c - z_b;

      vecAP << x_p - x_a, y_p - y_a, z_p - z_a;
      vecBP << x_p - x_b, y_p - y_b, z_p - z_b;
      vecCP << x_p - x_c, y_p - y_c, z_p - z_c;

      // double areaTot = 0.5 * vecAB.cross(vecAC).squaredNorm();
      areaAPC = 0.5 * vecAC.cross(vecAP).squaredNorm();
      areaAPB = 0.5 * vecAB.cross(vecAP).squaredNorm();
      areaBPC = 0.5 * vecBC.cross(vecBP).squaredNorm();
      areaTot = areaAPC + areaAPB + areaBPC;

      wa = areaBPC / areaTot;
      wb = areaAPC / areaTot;
      wc = areaAPB / areaTot;

      weightedInt = wa * this->bitmap[closestPoint[0]] +
                    wb * this->bitmap[closestPoint[1]] +
                    (1.0 - wa - wb) * this->bitmap[closestPoint[2]];

      // int localIndex = ns;
      for (j = 0; j < verticeHexaVec[ns].rows(); j++) {
        int uCumul = verticeHexaVec[ns](j, 0);
        int vCumul = verticeHexaVec[ns](j, 1);
        I_r[uCumul][vCumul] = (T)weightedInt;
      }

      // vpImageConvert::convert(I_r, I_out_CV);
      // cv::imshow("i out temp", I_out_CV);
      // cv::waitKey(1);
    }

    return 0;
  }

  int toImage(vpImage<T> &I_r, vpImage<unsigned char> *Mask = NULL) {
    if (this->nbSamples == 0)
      return -1;

    T *pt_bitmap = this->bitmap;
    float *pt_bitmapf = this->bitmapf;
    prCartesian3DPointVec *pt_XS = (prCartesian3DPointVec *)this->ge;
    int j;

    for (unsigned long ns = 0; ns < this->nbSamples; ns++, pt_bitmap++) {
      // #pragma omp parallel for
      for (j = 0; j < verticeHexaVec[ns].rows(); j++) {
        int uCumul = verticeHexaVec[ns](j, 0);
        int vCumul = verticeHexaVec[ns](j, 1);
        I_r[uCumul][vCumul] = *pt_bitmap;
      }
    }

    return 0;
  }

  int toImage(cv::Mat &I_r, vpImage<unsigned char> *Mask = NULL) {
    if (this->nbSamples == 0)
      return -1;

    T *pt_bitmap = this->bitmap;
    // float *pt_bitmapf = this->bitmapf;
    // prCartesian3DPointVec *pt_XS = (prCartesian3DPointVec *)this->ge;
    int j;
    uint8_t *pixelPtrIn = (uint8_t *)I_r.data;

#pragma omp parallel for
    for (unsigned long ns = 0; ns < this->nbSamples; ns++) {
      // #pragma omp parallel for
      for (j = 0; j < verticeHexaVec[ns].rows(); j++) {
        int uCumul = verticeHexaVec[ns](j, 0);
        int vCumul = verticeHexaVec[ns](j, 1);

        pixelPtrIn[uCumul * I_r.cols + vCumul] = (uint8_t)*pt_bitmap;
        // I_r.at<T>(uCumul, vCumul) = *pt_bitmap;
      }
      pt_bitmap++;
    }

    return 0;
  }

  int toImageF(cv::Mat &I_r, vpImage<unsigned char> *Mask = NULL) {
    if (this->nbSamples == 0)
      return -1;

    T *pt_bitmap = this->bitmap;
    float *pt_bitmapf = this->bitmapf;
    prCartesian3DPointVec *pt_XS = (prCartesian3DPointVec *)this->ge;
    int j;

    for (unsigned long ns = 0; ns < this->nbSamples;
         ns++, pt_bitmap++, pt_bitmapf++) {
      // #pragma omp parallel for
      for (j = 0; j < verticeHexaVec[ns].rows(); j++) {
        int uCumul = verticeHexaVec[ns](j, 0);
        int vCumul = verticeHexaVec[ns](j, 1);
        I_r.at<float>(uCumul, vCumul) = *pt_bitmapf;
      }
    }

    return 0;
  }
  int toImageRotate(cv::Mat &I_r, vpPoseVector &r) {

    if (this->nbSamples == 0)
      return -1;

    T *pt_bitmap = this->bitmap;
    float *pt_bitmapf = this->bitmapf;
    prCartesian3DPointVec XSs;
    prCartesian3DPointVec *pt_XS = (prCartesian3DPointVec *)this->ge;
    int j;
    prPointFeature P;

    vpHomogeneousMatrix dMc;
    dMc.buildFrom(r);
    dMc = dMc.inverse(); // a justifier clairement dans la doc
    Eigen::Vector3d vecAB, vecAC, vecBC, vecAP, vecBP, vecCP;
    double x_a, y_a, z_a, x_b, y_b, z_b, x_c, y_c, z_c, x_p, y_p, z_p;
    double areaAPC, areaAPB, areaBPC, areaTot;
    double wa, wb, wc;
    double weightedInt;

    for (unsigned long ns = 0; ns < this->nbSamples; ns++, pt_XS++) {
      P.sX = *pt_XS;
      P.sX = P.sX.changeFrame(dMc);
      // int localIndex = getMinDistance3D(P.sX);
      std::vector<unsigned int> closestPoint = this->get3ClosestPoints(P.sX);

      x_a = ((prCartesian3DPointVec *)this->ge)[closestPoint[0]].get_X();
      y_a = ((prCartesian3DPointVec *)this->ge)[closestPoint[0]].get_Y();
      z_a = ((prCartesian3DPointVec *)this->ge)[closestPoint[0]].get_Z();

      x_b = ((prCartesian3DPointVec *)this->ge)[closestPoint[1]].get_X();
      y_b = ((prCartesian3DPointVec *)this->ge)[closestPoint[1]].get_Y();
      z_b = ((prCartesian3DPointVec *)this->ge)[closestPoint[1]].get_Z();

      x_c = ((prCartesian3DPointVec *)this->ge)[closestPoint[2]].get_X();
      y_c = ((prCartesian3DPointVec *)this->ge)[closestPoint[2]].get_Y();
      z_c = ((prCartesian3DPointVec *)this->ge)[closestPoint[2]].get_Z();

      x_p = P.sX.get_X();
      y_p = P.sX.get_Y();
      z_p = P.sX.get_Z();

      vecAB << x_b - x_a, y_b - y_a, z_b - z_a;
      vecAC << x_c - x_a, y_c - y_a, z_c - z_a;
      vecBC << x_c - x_b, y_c - y_b, z_c - z_b;

      vecAP << x_p - x_a, y_p - y_a, z_p - z_a;
      vecBP << x_p - x_b, y_p - y_b, z_p - z_b;
      vecCP << x_p - x_c, y_p - y_c, z_p - z_c;

      // double areaTot = 0.5 * vecAB.cross(vecAC).squaredNorm();
      areaAPC = 0.5 * vecAC.cross(vecAP).squaredNorm();
      areaAPB = 0.5 * vecAB.cross(vecAP).squaredNorm();
      areaBPC = 0.5 * vecBC.cross(vecBP).squaredNorm();
      areaTot = areaAPC + areaAPB + areaBPC;

      wa = areaBPC / areaTot;
      wb = areaAPC / areaTot;
      wc = areaAPB / areaTot;

      weightedInt = wa * this->bitmapf[closestPoint[0]] +
                    wb * this->bitmapf[closestPoint[1]] +
                    (1.0 - wa - wb) * this->bitmapf[closestPoint[2]];

      // int localIndex = ns;
      for (j = 0; j < verticeHexaVec[ns].rows(); j++) {
        int uCumul = verticeHexaVec[ns](j, 0);
        int vCumul = verticeHexaVec[ns](j, 1);
        I_r.at<T>(uCumul, vCumul) = (T)weightedInt;
      }
    }

    return 0;
  }

  int buildFromEquiRect(vpImage<T> &I, prEquirectangular &equiCam,
                        vpImage<unsigned char> *Mask = NULL) {
    this->I_req = I;

    if (this->nbSamples == 0)
      return -1;

    if (this->bitmapf != NULL) {
      delete[] this->bitmapf;
      this->bitmapf = NULL;
    }

    this->bitmapf = new float[this->nbSamples];
    this->bitmapfUpdated = new float[this->nbSamples];
    this->isBitmapfSet = true;

    T *pt_bitmap = this->bitmap;
    float *pt_bitmapf = this->bitmapf;
    float *pt_bitmapf_updated = this->bitmapfUpdated;
    double u, v, du, dv, unmdu, unmdv;
    unsigned int i, j;

    imWidth = I.getWidth();
    imHeight = I.getHeight();

    if (!isCellMapInit || verticeHexaVec.empty()) {
      createImageCellsMap(equiCam);
    }

    for (unsigned long ns = 0; ns < this->nbSamples;
         ns++, pt_bitmap++, pt_bitmapf++, pt_bitmapf_updated++) {
      *pt_bitmap = 0;
      *pt_bitmapf = 0.;
      double cumulColor = 0;

      if (verticeHexaVec[ns].rows() > 1) {
        for (j = 0; j < verticeHexaVec[ns].rows(); j++) {
          int uCumul = verticeHexaVec[ns](j, 0);
          int vCumul = verticeHexaVec[ns](j, 1);
          cumulColor += I[uCumul][vCumul];
        }

        cumulColor = (double)cumulColor / (double)verticeHexaVec[ns].rows();
      }

      *pt_bitmap = (int)cumulColor;
      *pt_bitmapf = (float)cumulColor;
      *pt_bitmapf_updated = (float)cumulColor;
    }

    return 0;
  }

  int buildFromEquiRect(cv::Mat I, prEquirectangular &equiCam,
                        vpImage<unsigned char> *Mask = NULL) {
    // this->I_req = I;

    if (this->nbSamples == 0)
      return -1;

    if (this->bitmapf != NULL) {
      delete[] this->bitmapf;
      this->bitmapf = NULL;
    }

    this->bitmapf = new float[this->nbSamples];
    this->bitmapfUpdated = new float[this->nbSamples];
    this->isBitmapfSet = true;

    T *pt_bitmap = this->bitmap;
    float *pt_bitmapf = this->bitmapf;
    float *pt_bitmapf_updated = this->bitmapfUpdated;
    double u, v, du, dv, unmdu, unmdv;
    unsigned int i, j;

    imWidth = I.cols;
    imHeight = I.rows;

    if (!isCellMapInit || verticeHexaVec.empty()) {
      createImageCellsMap(equiCam);
    }

    for (unsigned long ns = 0; ns < this->nbSamples;
         ns++, pt_bitmap++, pt_bitmapf++, pt_bitmapf_updated++) {
      *pt_bitmap = 0;
      *pt_bitmapf = 0.;
      double cumulColor = 0;

      if (verticeHexaVec[ns].rows() > 1) {
        for (j = 0; j < verticeHexaVec[ns].rows(); j++) {
          int uCumul = verticeHexaVec[ns](j, 0);
          int vCumul = verticeHexaVec[ns](j, 1);
          cumulColor += I.at<unsigned char>(uCumul, vCumul);
        }

        cumulColor = (double)cumulColor / (double)verticeHexaVec[ns].rows();
      }

      *pt_bitmap = (int)cumulColor;
      *pt_bitmapf = (float)cumulColor;
      *pt_bitmapf_updated = (float)cumulColor;
    }

    return 0;
  }

  std::map<unsigned int, std::vector<prCartesian3DPointVec>> voroSamples;
  std::vector<Eigen::ArrayXXd> verticeHexaVec;
  Eigen::ArrayXXd centerCoords;
  std::vector<Eigen::MatrixXd> vecVoroContours;
  std::vector<std::vector<Eigen::MatrixXd>> vecTotVoroContours;
  double sumTotColors;
  bool isCellMapInit;
  int imWidth, imHeight;
};

#endif
