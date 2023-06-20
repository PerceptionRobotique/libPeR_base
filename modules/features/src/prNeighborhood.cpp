#include <per/prNeighborhood.h>

prNeighborhood::prNeighborhood()
{
    isNeighSet = false;

    nba = 2;
    nbNeigh = 6;
    di = 3;
    dj = 3;
}

void prNeighborhood::buildNeighborsCartSphere(int imgWidth, int imgHeight, prCameraModel *camera)
{

    if (isNeighSet)
    {
        return;
    }
    isNeighSet = true;
    // std::cout << "u0 " << camera->u0 << " v0 " << camera->v0 << std::endl;

    // Voisinage 3D
    Neigh = new int ****[3];

    int r, c, p, h, l, n; // Better initialise them in memory here, less work than defining it in every loop

    int nbri = imgHeight - 2 * di;
    int nbrj = imgWidth - 2 * dj;

    double x, y, Xs, Ys, Zs, u, v, vXs, vYs, vZs;
    double Xs0, Ys0, Zs0, Xs1, Ys1, Zs1;
    prPointFeature P;

    // std::cout << "--" << std::endl;
    // Computing the actual DeltaCSSampling with camera projection parameters
    P.set_u(camera->u0);
    P.set_v(camera->v0);
    camera->projectImageSphere(P, Xs0, Ys0, Zs0);
    P.set_u(camera->u0 + 1);
    P.set_v(camera->v0 + 1);
    camera->pixelMeterConversion(P);

    camera->projectImageSphere(P, Xs, Ys, Zs);

    // std::cout << "----" << std::endl;

    int u0 = camera->u0;
    int v0 = camera->v0;
    deltaCSSampling = sqrt(vpMath::sqr(Xs - Xs0) + vpMath::sqr(Ys - Ys0) + vpMath::sqr(Zs - Zs0));

    imgH = imgHeight;
    imgW = imgWidth;
    dimS = 2;
    neighSize = 3 * imgHeight * imgWidth * nba * nbNeigh;

    // First, initialise the arrays
    Neigh[0] = new int ***[imgHeight];

    for (h = 0; h < imgHeight; h++)
    {
        Neigh[0][h] = new int **[imgWidth];
        for (l = 0; l < imgWidth; l++)
        {
            Neigh[0][h][l] = new int *[nba];
            for (n = 0; n < nba; n++)
            {
                Neigh[0][h][l][n] = new int[nbNeigh];
            }
        }
    }

    Neigh[1] = new int ***[imgHeight];
    for (h = 0; h < imgHeight; h++)
    {
        Neigh[1][h] = new int **[imgWidth];
        for (l = 0; l < imgWidth; l++)
        {
            Neigh[1][h][l] = new int *[nba];
            for (n = 0; n < nba; n++)
            {
                Neigh[1][h][l][n] = new int[nbNeigh];
            }
        }
    }

    Neigh[2] = new int ***[imgHeight];
    for (h = 0; h < imgHeight; h++)
    {
        Neigh[2][h] = new int **[imgWidth];
        for (l = 0; l < imgWidth; l++)
        {
            Neigh[2][h][l] = new int *[nba];
            for (n = 0; n < nba; n++)
            {
                Neigh[2][h][l][n] = new int[nbNeigh];
            }
        }
    }

    // mise à zéro des coordonnées des voisins qui ne seront pas utilisées (bords des images)
    // bord supérieur
    for (r = 0; r < di; r++)
        for (c = 0; c < imgWidth; c++)
            for (p = 0; p < nbNeigh; p++)
                Neigh[0][r][c][0][p] = Neigh[0][r][c][1][p] = Neigh[1][r][c][0][p] = Neigh[1][r][c][1][p] = Neigh[2][r][c][0][p] = Neigh[2][r][c][1][p] = 0;

    // bord gauche
    for (r = di; r < imgHeight; r++)
        for (c = 0; c < dj; c++)
            for (p = 0; p < nbNeigh; p++)
                Neigh[0][r][c][0][p] = Neigh[0][r][c][1][p] = Neigh[1][r][c][0][p] = Neigh[1][r][c][1][p] = Neigh[2][r][c][0][p] = Neigh[2][r][c][1][p] = 0;

    // bord droit
    for (r = di; r < imgHeight; r++)
        for (c = dj + nbrj; c < imgWidth; c++)
            for (p = 0; p < nbNeigh; p++)
                Neigh[0][r][c][0][p] = Neigh[0][r][c][1][p] = Neigh[1][r][c][0][p] = Neigh[1][r][c][1][p] = Neigh[2][r][c][0][p] = Neigh[2][r][c][1][p] = 0;

    // bord inférieur
    for (r = di + nbri; r < imgHeight; r++)
        for (c = dj; c < imgWidth - nbrj; c++)
            for (p = 0; p < nbNeigh; p++)
                Neigh[0][r][c][0][p] = Neigh[0][r][c][1][p] = Neigh[1][r][c][0][p] = Neigh[1][r][c][1][p] = Neigh[2][r][c][0][p] = Neigh[2][r][c][1][p] = 0;

    // Once the arrays are defined, they need to be filled with the actual pixel intensities from the image
    // Now it's time to cross the whole image (ie. along the variables r and c)
    int i, j;
    prPointFeature P_temp;
    bool isProjPoint;

    for (r = di; r < (di + nbri); r++)
    {
        for (c = dj; c < (dj + nbrj); c++)
        {
            P_temp.set_u(c);
            P_temp.set_v(r);
            camera->pixelMeterConversion(P_temp);
            isProjPoint = camera->projectImageSphere(P_temp, Xs, Ys, Zs);

            if (isProjPoint)
            {

                // Voisinage "linÃ©aire" en Xs
                vXs = Xs - deltaCSSampling * nbNeigh / 2.0;
                vYs = Ys;
                vZs = Zs;

                i = 0;
                for (int dtheta = -nbNeigh / 2; dtheta < 0; dtheta++, vXs += deltaCSSampling, i++)
                {
                    P_temp.set_X(vXs);
                    P_temp.set_Y(vYs);
                    P_temp.set_Z(vZs);

                    camera->project3DImage(P_temp);
                    camera->meterPixelConversion(P_temp);
                    u = P_temp.get_u();
                    v = P_temp.get_v();

                    Neigh[0][r][c][1][i] = computeIndexNeighbor(imgWidth, u);
                    Neigh[0][r][c][0][i] = computeIndexNeighbor(imgHeight, v);

                    if (!isProjPoint)
                    {
                        Neigh[0][r][c][1][i] = 0;
                        Neigh[0][r][c][0][i] = 0;
                    }
                }

                vXs = Xs + deltaCSSampling;
                i = nbNeigh / 2;
                for (int dtheta = 1; dtheta <= nbNeigh / 2; dtheta++, vXs += deltaCSSampling, i++)
                {
                    P_temp.set_X(vXs);
                    P_temp.set_Y(vYs);
                    P_temp.set_Z(vZs);
                    camera->project3DImage(P_temp);
                    camera->meterPixelConversion(P_temp);
                    u = P_temp.get_u();
                    v = P_temp.get_v();

                    Neigh[0][r][c][1][i] = computeIndexNeighbor(imgWidth, u);
                    Neigh[0][r][c][0][i] = computeIndexNeighbor(imgHeight, v);

                    if (!isProjPoint)
                    {
                        Neigh[0][r][c][1][i] = 0;
                        Neigh[0][r][c][0][i] = 0;
                    }
                }

                // Voisinage "linÃ©aire" en Ys
                vXs = Xs;
                vYs = Ys - deltaCSSampling * nbNeigh / 2.0;
                vZs = Zs;
                i = 0;
                for (int dtheta = -nbNeigh / 2; dtheta < 0; dtheta++, vYs += deltaCSSampling, i++)
                {
                    P_temp.set_X(vXs);
                    P_temp.set_Y(vYs);
                    P_temp.set_Z(vZs);

                    camera->project3DImage(P_temp);
                    camera->meterPixelConversion(P_temp);
                    u = P_temp.get_u();
                    v = P_temp.get_v();

                    Neigh[1][r][c][1][i] = computeIndexNeighbor(imgWidth, u);
                    Neigh[1][r][c][0][i] = computeIndexNeighbor(imgHeight, v);

                    if (!isProjPoint)
                    {
                        Neigh[1][r][c][1][i] = 0;
                        Neigh[1][r][c][0][i] = 0;
                    }
                }
                vYs = Ys + deltaCSSampling;
                i = nbNeigh / 2;
                for (int dtheta = 1; dtheta <= nbNeigh / 2; dtheta++, vYs += deltaCSSampling, i++)
                {
                    P_temp.set_X(vXs);
                    P_temp.set_Y(vYs);
                    P_temp.set_Z(vZs);

                    camera->project3DImage(P_temp);
                    camera->meterPixelConversion(P_temp);
                    u = P_temp.get_u();
                    v = P_temp.get_v();

                    Neigh[1][r][c][1][i] = computeIndexNeighbor(imgWidth, u);
                    Neigh[1][r][c][0][i] = computeIndexNeighbor(imgHeight, v);

                    if (!isProjPoint)
                    {
                        Neigh[1][r][c][1][i] = 0;
                        Neigh[1][r][c][0][i] = 0;
                    }
                }

                // Voisinage "linÃ©aire" en Zs
                vXs = Xs;
                vYs = Ys;
                vZs = Zs - deltaCSSampling * nbNeigh / 2.0;
                i = 0;
                for (int dtheta = -nbNeigh / 2; dtheta < 0; dtheta++, vZs += deltaCSSampling, i++)
                {
                    P_temp.set_X(vXs);
                    P_temp.set_Y(vYs);
                    P_temp.set_Z(vZs);

                    camera->project3DImage(P_temp);
                    camera->meterPixelConversion(P_temp);
                    u = P_temp.get_u();
                    v = P_temp.get_v();

                    Neigh[2][r][c][1][i] = computeIndexNeighbor(imgWidth, u);
                    Neigh[2][r][c][0][i] = computeIndexNeighbor(imgHeight, v);

                    if (!isProjPoint)
                    {
                        Neigh[2][r][c][1][i] = 0;
                        Neigh[2][r][c][0][i] = 0;
                    }
                }
                vZs = Zs + deltaCSSampling;
                i = nbNeigh / 2;
                for (int dtheta = 1; dtheta <= nbNeigh / 2; dtheta++, vZs += deltaCSSampling, i++)
                {
                    P_temp.set_X(vXs);
                    P_temp.set_Y(vYs);
                    P_temp.set_Z(vZs);

                    camera->project3DImage(P_temp);
                    camera->meterPixelConversion(P_temp);
                    u = P_temp.get_u();
                    v = P_temp.get_v();

                    Neigh[2][r][c][1][i] = computeIndexNeighbor(imgWidth, u);
                    Neigh[2][r][c][0][i] = computeIndexNeighbor(imgHeight, v);

                    if (!isProjPoint)
                    {
                        Neigh[2][r][c][1][i] = 0;
                        Neigh[2][r][c][0][i] = 0;
                    }
                }
            }
        }
    }
}

int prNeighborhood::computeIndexNeighbor(int winSize, double inputIndex)
{
    int index = 0;

    if (inputIndex > winSize - 4 || inputIndex < 4)
    {
        index = 0;
    }
    else
    {
        index = (int)round(inputIndex);
    }

    return index;
}

void prNeighborhood::deleteNeighbors()
{
    if (!isNeighSet)
        return;

    isNeighSet = false;

    for (int numDim = 0; numDim < dimS; numDim++)
    {
        for (int h = 0; h < imgH; h++)
        {
            for (int l = 0; l < imgW; l++)
            {
                for (int n = 0; n < nba; n++)
                {
                    delete[] Neigh[numDim][h][l][n];
                }
                delete[] Neigh[numDim][h][l];
            }
            delete[] Neigh[numDim][h];
        }
        delete[] Neigh[numDim];
    }
    delete[] Neigh;
}

prNeighborhood &prNeighborhood::operator=(const prNeighborhood &neigh)
{
    deltaCSSampling = neigh.deltaCSSampling;
    imgH = neigh.imgH;
    imgW = neigh.imgW;
    dimS = neigh.dimS;
    neighSize = neigh.neighSize;

    return *this;
}
