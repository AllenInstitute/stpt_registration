/*=========================================================================

  idpProjectionLocalAlignmentDispatcher.cxx

  Copyright (c) Allen Institute for Brain Science. All rights reserved.

=========================================================================*/
#ifndef __idpProjectionLocalAlignmentDispatcher_txx
#define __idpProjectionLocalAlignmentDispatcher_txx

#include "idpProjectionLocalAlignmentDispatcher.h"
#include "tinyxml.h"
#include "idpRegistrationUtilities.h"
#include "itkSliceImageConstIterator.h"
#include "idpXMLUtilities.h"
#include "idpImageSeries.h"
#include "idpAffineCorrelationVolumeRegistration.h"
#include "idpMutualInfoVersorRigidVolumeRegistration.h"
#include "itkIdentityTransform.h"
#include "idpMutualInfoCenteredAffineRegularGradientDescentRegistration.h"
#include <sstream>

#include <iostream>
#include <fstream>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <vector>
#include <algorithm>
#include <functional>
#include <unistd.h>
#include <string.h>
#include <stdarg.h>
#include <unistd.h>
#include "vnl/vnl_matrix.txx"
#include "vnl/vnl_c_vector.txx"

namespace itk {
    namespace idp {

        /**
         * Constructor
         */
        template<class TPixel>
        ProjectionLocalAlignmentDispatcher<TPixel>::
        ProjectionLocalAlignmentDispatcher() {

            this->SetSampleNum(64);

            this->SetImageSize0(256);
            this->SetImageSize1(256);
            this->SetImageSize2(256);

            this->SetMetricValue0(0.0);
            this->SetMetricValue1(0.0);

            this->SetNumBin(32);
            this->SetMaxFixed(255);
            this->SetMaxMoving(255);
        }

        /**
         * Destructor
         */
        template<class TPixel>
        ProjectionLocalAlignmentDispatcher<TPixel>::
        ~ProjectionLocalAlignmentDispatcher() {

        }

        /**
         * PrintSelf
         */
        template<class TPixel>
        void
        ProjectionLocalAlignmentDispatcher<TPixel>::
        PrintSelf(std::ostream& os, Indent indent) const {
            Superclass::PrintSelf(os, indent);

        }

        template<class TPixel>
        void
        ProjectionLocalAlignmentDispatcher<TPixel>::
        minimumInd(float *numbers, float &value, int &index, int length) {
            value = numbers[0];
            int i;
            index = 0;
            for (i = 1; i < length; i++) {
                if (numbers[i] < value) {
                    index = i;
                    value = numbers[i];
                }
            }
        }

        template<class TPixel>
        void
        ProjectionLocalAlignmentDispatcher<TPixel>::
        writeOutput(float* data, char* name, int length) {

            std::ofstream ofs1(name, std::ofstream::binary);
            ofs1.write((char *) data, length * sizeof (float));
            ofs1.close();
        }

        template<class TPixel>
        void
        ProjectionLocalAlignmentDispatcher<TPixel>::
        DfmfldModulation(float* dfmfld, int* ordered, int* parents, int stepsize) {
            int xsize = this->GetImageSize0();
            int ysize = this->GetImageSize1();
            int zsize = this->GetImageSize2();
            int xgridsize = xsize / stepsize;
            int ygridsize = ysize / stepsize;
            int zgridsize = zsize / stepsize;

            int num_vertices = xgridsize * ygridsize*zgridsize;

            int num_neighbours = 6;
            float* edgecost = new float[num_vertices * num_neighbours];
            int* index_neighbours = new int[num_vertices * num_neighbours];
            for (int i = 0; i < num_vertices * num_neighbours; i++) {
                edgecost[i] = 0.0;
                index_neighbours[i] = -1;
            }

            int xdelta[6] = {-1, 1, 0, 0, 0, 0};
            int ydelta[6] = {0, 0, -1, 1, 0, 0};
            int zdelta[6] = {0, 0, 0, 0, -1, 1};
            int x, y, z, xnew, ynew, znew;

            //calculate weights 
            for (int k = 0; k < zgridsize; k++) {
                for (int j = 0; j < ygridsize; j++) {
                    for (int i = 0; i < xgridsize; i++) {
                        for (int nb = 0; nb < num_neighbours; nb++) {
                            if ((i + ydelta[nb]) >= 0 & (i + ydelta[nb]) < xgridsize & (j + xdelta[nb]) >= 0 & (j + xdelta[nb]) < ygridsize & (k + zdelta[nb]) >= 0 & (k + zdelta[nb]) < zgridsize) {
                                index_neighbours[i + j * xgridsize + k * xgridsize * ygridsize + nb * num_vertices] = i + ydelta[nb]+(j + xdelta[nb]) * xgridsize + (k + zdelta[nb]) * xgridsize*ygridsize;
                                for (int k1 = 0; k1 < stepsize; k1++) {
                                    for (int j1 = 0; j1 < stepsize; j1++) {
                                        for (int i1 = 0; i1 < stepsize; i1++) {
                                            x = j * stepsize + j1;
                                            y = i * stepsize + i1;
                                            z = k * stepsize + k1;
                                            xnew = (j + xdelta[nb]) * stepsize + j1;
                                            ynew = (i + ydelta[nb]) * stepsize + i1;
                                            znew = (k + zdelta[nb]) * stepsize + k1;
                                            edgecost[i + j * xgridsize + k * xgridsize * ygridsize + nb * num_vertices] += fabs(dfmfld[y + x * xsize + z * xsize * ysize] - dfmfld[ynew + xnew * xsize + znew * xsize * ysize]);
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }

            int root = xgridsize / 2 + ygridsize / 2 * xgridsize + zgridsize / 2 * xgridsize*ygridsize;

            std::vector<Relation> sorted;
            bool* vertices = new bool[num_vertices];
            int* level = new int[num_vertices];
            for (int i = 0; i < num_vertices; i++) {
                vertices[i] = false;
                parents[i] = -1;
            }

            level[root] = 0;
            int last = root;
            vertices[root] = true;
            Relation edgeout = Relation(0.0, -1, -1);
            Relation minedge = Relation(0.0, -1, -1);
            float cost = 0.0;
            
            for (int i = 0; i < num_vertices - 1; i++) { // have all grids added
                for (int j = 0; j < num_neighbours; j++) {
                    int n = index_neighbours[last + j * num_vertices];
                    if (n >= 0) {
                        sorted.push_back(Relation(edgecost[last + j * num_vertices], last, n));
                        push_heap(sorted.begin(), sorted.end());
                    }
                }
                last = -1;
                //find lowest
                while (last == -1) {
                    minedge = sorted.front();
                    pop_heap(sorted.begin(), sorted.end());
                    sorted.pop_back();
                    last = newRelation(minedge, edgeout, vertices); //return next valid 
                }
                cost += edgeout.weight;
                vertices[last] = true;
                level[edgeout.vert2] = level[edgeout.vert1] + 1;
                parents[edgeout.vert2] = edgeout.vert1;
            }

            //re-caculate
            int maxlevel = 0;
            for (int i = 0; i < num_vertices; i++) {
                if (level[i] > maxlevel)
                    maxlevel = level[i];
            }
            maxlevel++;
            int* leveloffset = new int[maxlevel];
            int* levelcount = new int[maxlevel];
            for (int i = 0; i < maxlevel; i++) {
                leveloffset[i] = 0;
                levelcount[i] = 0;
            }
            for (int i = 0; i < num_vertices; i++) {
                if (level[i] < maxlevel - 1)
                    leveloffset[level[i] + 1]++; //counting number 
            }
            for (int i = 1; i < maxlevel; i++) {
                leveloffset[i] += leveloffset[i - 1];
            }
            for (int i = 0; i < num_vertices; i++) {
                int num = leveloffset[level[i]] + levelcount[level[i]];
                levelcount[level[i]]++;
                ordered[num] = i;
            }

            sorted.clear();

            delete edgecost;
            delete index_neighbours;
            delete levelcount;
            delete leveloffset;
            delete vertices;
            delete level;

        }

        template<class TPixel>
        void
        ProjectionLocalAlignmentDispatcher<TPixel>::
        energy(float *val, int* index, int length, float offset, int k, int* v, float* z, float* f, int* indexinv) {
            float INF = 1e10;

            int j = 0;
            z[0] = -INF;
            z[1] = INF;
            v[0] = 0;
            for (int q = 1; q < length; q++) {
                float s = ((val[q * k] + powf((float) q + offset, 2.0))-(val[v[j] * k] + powf((float) v[j] + offset, 2.0))) / (2.0 * (float) (q - v[j]));
                while (s <= z[j]) {
                    j--;
                    s = ((val[q * k] + powf((float) q + offset, 2.0))-(val[v[j] * k] + powf((float) v[j] + offset, 2.0))) / (2.0 * (float) (q - v[j]));
                }
                j++;
                v[j] = q;
                z[j] = s;
                z[j + 1] = INF;
            }
            j = 0;
            for (int q = 0; q < length; q++) {
                f[q] = val[q * k];
                indexinv[q] = index[q * k];
            }
            for (int q = 0; q < length; q++) {
                while (z[j + 1] < q) {
                    j++;
                }
                index[q * k] = indexinv[v[j]];
                val[q * k] = powf((float) q - ((float) v[j] + offset), 2.0) + f[v[j]];
            }

        }
        

        template<class TPixel>
        void
        ProjectionLocalAlignmentDispatcher<TPixel>::
        regulate(float* r, int* indr, int rl, float dx, float dy, float dz) {

            for (int i = 0; i < rl * rl * rl; i++) {
                indr[i] = i;
            }
            int* v = new int[rl];
            float* z = new float[rl + 1];
            float* f = new float[rl];
            int* i1 = new int[rl];

            for (int k = 0; k < rl; k++) {
                for (int i = 0; i < rl; i++) {
                    energy(r + i + k * rl*rl, indr + i + k * rl*rl, rl, -dx, rl, v, z, f, i1);
                }
            }
            for (int k = 0; k < rl; k++) {
                for (int j = 0; j < rl; j++) {
                    energy(r + j * rl + k * rl*rl, indr + j * rl + k * rl*rl, rl, -dy, 1, v, z, f, i1);
                }
            }

            for (int j = 0; j < rl; j++) {
                for (int i = 0; i < rl; i++) {
                    energy(r + i + j*rl, indr + i + j*rl, rl, -dz, rl*rl, v, z, f, i1);
                }
            }
            delete []i1;
            delete []f;

            delete []v;
            delete []z;


        }

        template<class TPixel>
        void
        ProjectionLocalAlignmentDispatcher<TPixel>::
        interpolate3d(float* interp, float* input, float* x1, float* y1, float* z1, int xsize, int ysize, int zsize, int xsizenew, int ysizenew, int zsizenew, bool flag) {
            using namespace std;

            for (int k = 0; k < zsize; k++) {
                for (int j = 0; j < ysize; j++) {
                    for (int i = 0; i < xsize; i++) {
                        int x = (int) floor(x1[i + j * xsize + k * xsize * ysize]);
                        int y = (int) floor(y1[i + j * xsize + k * xsize * ysize]);
                        int z = (int) floor(z1[i + j * xsize + k * xsize * ysize]);
                        float deltax = x1[i + j * xsize + k * xsize * ysize] - x;
                        float deltay = y1[i + j * xsize + k * xsize * ysize] - y;
                        float deltaz = z1[i + j * xsize + k * xsize * ysize] - z;

                        if (flag) {
                            x += j;
                            y += i;
                            z += k;
                        }
                        interp[i + j * xsize + k * xsize * ysize] = (1.0 - deltax)*(1.0 - deltay)*(1.0 - deltaz) * input[min(max(y, 0), xsizenew - 1) + min(max(x, 0), ysizenew - 1) * xsizenew + min(max(z, 0), zsizenew - 1) * xsizenew * ysizenew]+
                                (1.0 - deltax) * deltay * (1.0 - deltaz) * input[min(max(y + 1, 0), xsizenew - 1) + min(max(x, 0), ysizenew - 1) * xsizenew + min(max(z, 0), zsizenew - 1) * xsizenew * ysizenew] +
                                deltax * (1.0 - deltay)*(1.0 - deltaz) * input[min(max(y, 0), xsizenew - 1) + min(max(x + 1, 0), ysizenew - 1) * xsizenew + min(max(z, 0), zsizenew - 1) * xsizenew * ysizenew]+
                                (1.0 - deltax)*(1.0 - deltay) * deltaz * input[min(max(y, 0), xsizenew - 1) + min(max(x, 0), ysizenew - 1) * xsizenew + min(max(z + 1, 0), zsizenew - 1) * xsizenew * ysizenew] +
                                deltax * deltay * (1.0 - deltaz) * input[min(max(y + 1, 0), xsizenew - 1) + min(max(x + 1, 0), ysizenew - 1) * xsizenew + min(max(z, 0), zsizenew - 1) * xsizenew * ysizenew]+
                                (1.0 - deltax) * deltay * deltaz * input[min(max(y + 1, 0), xsizenew - 1) + min(max(x, 0), ysizenew - 1) * xsizenew + min(max(z + 1, 0), zsizenew - 1) * xsizenew * ysizenew] +
                                deltax * (1.0 - deltay) * deltaz * input[min(max(y, 0), xsizenew - 1) + min(max(x + 1, 0), ysizenew - 1) * xsizenew + min(max(z + 1, 0), zsizenew - 1) * xsizenew * ysizenew] +
                                deltax * deltay * deltaz * input[min(max(y + 1, 0), xsizenew - 1) + min(max(x + 1, 0), ysizenew - 1) * xsizenew + min(max(z + 1, 0), zsizenew - 1) * xsizenew * ysizenew];
                    }
                }
            }
        }


        template<class TPixel>
        void
        ProjectionLocalAlignmentDispatcher<TPixel>::
        InterpDeformationFld(float* fldx_in, float* fldy_in, float* fldz_in, float* fldx, float* fldy, float* fldz, int xsize_in, int ysize_in, int zsize_in, int xsize, int ysize, int zsize) {


            float scale_m = (float) xsize_in / (float) xsize;
            float scale_n = (float) ysize_in / (float) ysize;
            float scale_o = (float) zsize_in / (float) zsize;

            float* x1 = new float[xsize_in * ysize_in * zsize_in];
            float* y1 = new float[xsize_in * ysize_in * zsize_in];
            float* z1 = new float[xsize_in * ysize_in * zsize_in];
            for (int k = 0; k < zsize_in; k++) {
                for (int j = 0; j < ysize_in; j++) {
                    for (int i = 0; i < xsize_in; i++) {
                        x1[i + j * xsize_in + k * xsize_in * ysize_in] = j / scale_n;
                        y1[i + j * xsize_in + k * xsize_in * ysize_in] = i / scale_m;
                        z1[i + j * xsize_in + k * xsize_in * ysize_in] = k / scale_o;
                    }
                }
            }

            interpolate3d(fldx_in, fldx, x1, y1, z1, xsize_in, ysize_in, zsize_in, xsize, ysize, zsize, false);
            interpolate3d(fldy_in, fldy, x1, y1, z1, xsize_in, ysize_in, zsize_in, xsize, ysize, zsize, false);
            interpolate3d(fldz_in, fldz, x1, y1, z1, xsize_in, ysize_in, zsize_in, xsize, ysize, zsize, false);

            delete []x1;
            delete []y1;
            delete []z1;

        }

        template<class TPixel>
        void
        ProjectionLocalAlignmentDispatcher<TPixel>::
        fldinv(float* fldx_in, float* fldy_in, float* fldz_in, float* fldx, float* fldy, float* fldz, int m, int n, int o) {
            float* un = new float[m * n * o];
            float* vn = new float[m * n * o];
            float* wn = new float[m * n * o];
            float* uin = new float[m * n * o];
            float* vin = new float[m * n * o];
            float* win = new float[m * n * o];

            for (int i = 0; i < m * n * o; i++) {
                uin[i] = -fldx[i];
                vin[i] = -fldy[i];
                win[i] = -fldz[i];
                fldx_in[i] = 0.0;
                fldy_in[i] = 0.0;
                fldz_in[i] = 0.0;
            }
            for (int it = 0; it < 10; it++) {
                interpolate3d(un, uin, fldx_in, fldy_in, fldz_in, m, n, o, m, n, o, true);
                interpolate3d(vn, vin, fldx_in, fldy_in, fldz_in, m, n, o, m, n, o, true);
                interpolate3d(wn, win, fldx_in, fldy_in, fldz_in, m, n, o, m, n, o, true);
                for (int i = 0; i < m * n * o; i++) {
                    fldx_in[i] = un[i];
                    fldy_in[i] = vn[i];
                    fldz_in[i] = wn[i];

                }
            }
            delete []uin;
            delete []vin;
            delete []win;
            delete []un;
            delete []vn;
            delete []wn;
        }

        template<class TPixel>
        void
        ProjectionLocalAlignmentDispatcher<TPixel>::
        combineDeformation(float* fldx, float* fldy, float* fldz, float* fldx1, float* fldy1, float* fldz1, float* fldx2, float* fldy2, float* fldz2, int m, int n, int o) {
            float* uc = new float[m * n * o];
            float* vc = new float[m * n * o];
            float* wc = new float[m * n * o];

            interpolate3d(uc, fldx1, fldx2, fldy2, fldz2, m, n, o, m, n, o, true);
            interpolate3d(vc, fldy1, fldx2, fldy2, fldz2, m, n, o, m, n, o, true);
            interpolate3d(wc, fldz1, fldx2, fldy2, fldz2, m, n, o, m, n, o, true);

            for (int i = 0; i < m * n * o; i++) {
                fldx[i] = uc[i] + fldx2[i];
                fldy[i] = vc[i] + fldy2[i];
                fldz[i] = wc[i] + fldz2[i];

            }
            delete []uc;
            delete []vc;
            delete []wc;

        }

        template<class TPixel>
        void
        ProjectionLocalAlignmentDispatcher<TPixel>::
        smoothDfm(float* fldx, float* fldy, float *fldz, int m, int n, int o, int expsteps, int factor) {
            float* fldx1 = new float[m * n * o];
            float* fldy1 = new float[m * n * o];
            float* fldz1 = new float[m * n * o];

            float* fldx2 = new float[m * n * o];
            float* fldy2 = new float[m * n * o];
            float* fldz2 = new float[m * n * o];

            float scale = 1.0 / (float) factor;
            float coeff = 1.0 / (float) powf(2.0, expsteps);
            for (int i = 0; i < m * n * o; i++) {
                fldx1[i] = coeff * fldx[i] * scale;
                fldy1[i] = coeff * fldy[i] * scale;
                fldz1[i] = coeff * fldz[i] * scale;

            }
            for (int it = 0; it < expsteps; it++) {
                combineDeformation(fldx2, fldy2, fldz2, fldx1, fldy1, fldz1, fldx1, fldy1, fldz1, m, n, o);
                for (int i = 0; i < m * n * o; i++) {
                    fldx1[i] = fldx2[i];
                    fldy1[i] = fldy2[i];
                    fldz1[i] = fldz2[i];
                }
            }
            for (int i = 0; i < m * n * o; i++) {
                fldx[i] = fldx2[i]*(float) factor;
                fldy[i] = fldy2[i]*(float) factor;
                fldz[i] = fldz2[i]*(float) factor;

            }
            delete []fldx2;
            delete []fldy2;
            delete []fldz2;
            delete []fldx1;
            delete []fldy1;
            delete []fldz1;

        }

        template<class TPixel>
        void
        ProjectionLocalAlignmentDispatcher<TPixel>::
        symmetricMapping(float* xforward, float* yforward, float* zforward, float* xbackward, float* ybackward, float* zbackward, int m, int n, int o, int factor) {
            float factor1 = 1.0 / (float) factor;
            float* usym = new float[m* n * o];
            float* vsym = new float[m* n * o];
            float* wsym = new float[m* n * o];
            float* usym2 = new float[m* n * o];
            float* vsym2 = new float[m* n * o];
            float* wsym2 = new float[m* n * o];

            for (int i = 0; i < m * n * o; i++) {
                usym[i] = xforward[i]*(factor1 * 0.5);
                vsym[i] = yforward[i]*(factor1 * 0.5);
                wsym[i] = zforward[i]*(factor1 * 0.5);
                usym2[i] = xbackward[i]*(factor1 * 0.5);
                vsym2[i] = ybackward[i]*(factor1 * 0.5);
                wsym2[i] = zbackward[i]*(factor1 * 0.5);

            }
            float* ui = new float[m * n * o];
            float* vi = new float[m * n * o];
            float* wi = new float[m * n * o];
            float* u2i = new float[m * n * o];
            float* v2i = new float[m * n * o];
            float* w2i = new float[m * n * o];
            fldinv(ui, vi, wi, usym, vsym, wsym, m, n, o);
            fldinv(u2i, v2i, w2i, usym2, vsym2, wsym2, m, n, o);

            combineDeformation(xforward, yforward, zforward, u2i, v2i, w2i, usym, vsym, wsym, m, n, o);
            combineDeformation(xbackward, ybackward, zbackward, ui, vi, wi, usym2, vsym2, wsym2, m, n, o);

            for (int i = 0; i < m * n * o; i++) {
                xforward[i] *= (float) factor;
                yforward[i] *= (float) factor;
                zforward[i] *= (float) factor;
                xbackward[i] *= (float) factor;
                ybackward[i] *= (float) factor;
                zbackward[i] *= (float) factor;
            }

            delete []ui;
            delete []vi;
            delete []wi;
            delete []u2i;
            delete []v2i;
            delete []w2i;
            delete []usym;
            delete []vsym;
            delete []wsym;
            delete []usym2;
            delete []vsym2;
            delete []wsym2;

        }

        template<class TPixel>
        void
        ProjectionLocalAlignmentDispatcher<TPixel>::
        readRaw(char str2[], float*& pixels, int SZ) {
            FILE * pFile;
            int i;
            size_t result;

            pFile = fopen(str2, "rb");
            if (pFile == NULL) {
                fputs("File error", stderr);
                exit(1);
            }

            fseek(pFile, 0, 0);

            size_t lSize;
            float * buffer;
            lSize = SZ;
            buffer = (float*) malloc(sizeof (float)*lSize);
            if (buffer == NULL) {
                fputs("File error", stderr);
                exit(2);
            }

            // copy the file into the buffer:
            result = fread(buffer, sizeof (float), lSize, pFile);
            if (result != lSize) {
                fputs("Reading error", stderr);
                exit(3);
            }

            for (i = 0; i < (SZ); i++)
                pixels[i] = buffer[i];
            free(buffer);
            fclose(pFile);
        }

        template<class TPixel>
        void
        ProjectionLocalAlignmentDispatcher<TPixel>::
        genDfm(float* fldx, float* fldy, float* fldz, float* fldx0, float* fldy0, float* fldz0, float* costMatrix, float lambda, int hw, int stepsize, float gridmin, int* sorted, int* relations) {

            using namespace std;

            int sizex = this->GetImageSize0();
            int sizey = this->GetImageSize1();
            int sizez = this->GetImageSize2();

            int gridsizex = sizex / stepsize;
            int gridsizey = sizey / stepsize;
            int gridsizez = sizez / stepsize;


            int gridsizevol = gridsizex * gridsizey*gridsizez;

            // deformation field
            int searchsize = hw * 2 + 1;
            float* xs = new float[searchsize * searchsize * searchsize];
            float* ys = new float[searchsize * searchsize * searchsize];
            float* zs = new float[searchsize * searchsize * searchsize];

            for (int i = 0; i < searchsize; i++) {
                for (int j = 0; j < searchsize; j++) {
                    for (int k = 0; k < searchsize; k++) {
                        xs[i + j * searchsize + k * searchsize * searchsize] = (j - hw) * gridmin;
                        ys[i + j * searchsize + k * searchsize * searchsize] = (i - hw) * gridmin;
                        zs[i + j * searchsize + k * searchsize * searchsize] = (k - hw) * gridmin;
                    }
                }
            }
            int searchvol = searchsize * searchsize*searchsize;

            int *labels = new int[gridsizevol];
            short *indexall = new short[gridsizevol * searchvol];
            float *energy = new float[searchvol];
            float *values = new float[searchvol];
            int *index = new int[searchvol];


            for (int i = 0; i < searchvol; i++) {
                energy[i] = 0;
            }

            int frac = (int) (gridsizevol / 25);

            //calculate cost
            for (int i = (gridsizevol - 1); i > 0; i--) { //do for each control point
                if ((i % frac) == 0) {
                    cout << "x" << flush;
                }
                int next = sorted[i];
                int last = relations[sorted[i]];
                int znext = next / (gridsizex * gridsizey);
                int xnext = (next - znext * gridsizex * gridsizey) / gridsizex;
                int ynext = next - znext * gridsizex * gridsizey - xnext*gridsizex;
                
                int zlast = last / (gridsizex * gridsizey);
                int xlast = (last - zlast * gridsizex * gridsizey) / gridsizex;
                int ylast = last - zlast * gridsizex * gridsizey - xlast*gridsizex;

                for (int l = 0; l < searchvol; l++) {
                    energy[l] = costMatrix[next + l * gridsizevol];
                }

                float deltax = (fldx0[ylast + xlast * gridsizex + zlast * gridsizex * gridsizey] - fldx0[ynext + xnext * gridsizex + znext * gridsizex * gridsizey]) / (float) gridmin;
                float deltay = (fldy0[ylast + xlast * gridsizex + zlast * gridsizex * gridsizey] - fldy0[ynext + xnext * gridsizex + znext * gridsizex * gridsizey]) / (float) gridmin;
                float deltaz = (fldz0[ylast + xlast * gridsizex + zlast * gridsizex * gridsizey] - fldz0[ynext + xnext * gridsizex + znext * gridsizex * gridsizey]) / (float) gridmin;

                regulate(energy, index, searchsize, deltax, deltay, deltaz);

                for (int l = 0; l < searchvol; l++) {
                    indexall[ynext + xnext * gridsizex + znext * gridsizex * gridsizey + l * gridsizevol] = index[l];
                }

                for (int l = 0; l < searchvol; l++) {
                    costMatrix[last + l * gridsizevol] += energy[l];
                }

            }

            //compuate label
            int i = 0;
            int oroot = sorted[i];
            int z1 = oroot / (gridsizex * gridsizey);
            int x1 = (oroot - z1 * gridsizex * gridsizey) / gridsizex;
            int y1 = oroot - z1 * gridsizex * gridsizey - x1*gridsizex;
            for (int l = 0; l < searchvol; l++) {
                energy[l] = costMatrix[oroot + l * gridsizevol];
            }
            float value;
            int ind;
            minimumInd(energy, value, ind, searchvol);
            for (int l = 0; l < searchvol; l++) {
                indexall[y1 + x1 * gridsizex + z1 * gridsizex * gridsizey + l * gridsizevol] = l;
            }
            labels[y1 + x1 * gridsizex + z1 * gridsizex * gridsizey] = ind;
            fldx[y1 + x1 * gridsizex + z1 * gridsizex * gridsizey] = xs[ind] + fldx0[y1 + x1 * gridsizex + z1 * gridsizex * gridsizey];
            fldy[y1 + x1 * gridsizex + z1 * gridsizex * gridsizey] = ys[ind] + fldy0[y1 + x1 * gridsizex + z1 * gridsizex * gridsizey];
            fldz[y1 + x1 * gridsizex + z1 * gridsizex * gridsizey] = zs[ind] + fldz0[y1 + x1 * gridsizex + z1 * gridsizex * gridsizey];


            //compute displacements 
            for (int i = 1; i < gridsizevol; i++) {
                int next = sorted[i];
                int last = relations[sorted[i]];
                int z1 = next / (gridsizex * gridsizey);
                int x1 = (next - z1 * gridsizex * gridsizey) / gridsizex;
                int y1 = next - z1 * gridsizex * gridsizey - x1*gridsizex;
                int z2 = last / (gridsizex * gridsizey);
                int x2 = (last - z2 * gridsizex * gridsizey) / gridsizex;
                int y2 = last - z2 * gridsizex * gridsizey - x2*gridsizex;

                ind = indexall[y1 + x1 * gridsizex + z1 * gridsizex * gridsizey + labels[y2 + x2 * gridsizex + z2 * gridsizex * gridsizey] * gridsizex * gridsizey * gridsizez];
                labels[y1 + x1 * gridsizex + z1 * gridsizex * gridsizey] = ind;
                fldx[y1 + x1 * gridsizex + z1 * gridsizex * gridsizey] = xs[ind] + fldx0[y1 + x1 * gridsizex + z1 * gridsizex * gridsizey];
                fldy[y1 + x1 * gridsizex + z1 * gridsizex * gridsizey] = ys[ind] + fldy0[y1 + x1 * gridsizex + z1 * gridsizex * gridsizey];
                fldz[y1 + x1 * gridsizex + z1 * gridsizex * gridsizey] = zs[ind] + fldz0[y1 + x1 * gridsizex + z1 * gridsizex * gridsizey];

            }

            delete energy;
            delete values;
            delete index;
            delete indexall;
            delete labels;

            //return 0;
        }

        template<class TPixel>
        void
        ProjectionLocalAlignmentDispatcher<TPixel>::
        optimization(float* imgFix, float* imgMov, float* energy, float lambda, int hw, float gridstep, float mingrid) {
            using namespace std;

            float* fixed;
            float* moving = imgMov;

            bool subpixel = false;
            if (mingrid == 0.5)
                subpixel = true;

            float gamma = (float) gridstep / (lambda * (float) mingrid);

            int m = this->GetImageSize0();
            int n = this->GetImageSize1();
            int o = this->GetImageSize2();
            int sz = m*n*o;

            int gridx = m / gridstep;
            int gridy = n / gridstep;
            int gridz = o / gridstep;
            int sz1 = gridx*gridy*gridz;

            int len = hw * 2 + 1;
            int lengthfull = (int) pow(len, 3.0);
            float* xs = new float[lengthfull];
            float* ys = new float[lengthfull];
            float* zs = new float[lengthfull];

            for (int i = 0; i < len; i++) {
                for (int j = 0; j < len; j++) {
                    for (int k = 0; k < len; k++) {
                        xs[i + j * len + k * len * len] = (float) ((j - hw) * mingrid);
                        ys[i + j * len + k * len * len] = (float) ((i - hw) * mingrid);
                        zs[i + j * len + k * len * len] = (float) ((k - hw) * mingrid);

                    }
                }
            }

            int hw_local;
            if (subpixel)
                hw_local = hw;
            else
                hw_local = hw * (int) mingrid;

            float* movingi;
            int* movingiInd; //array of the index of bin for the histgram computing
            int* fixedInd;

            const float binSizeFixed = this->GetMaxFixed() / this->GetNumBin();
            const float binSizeMoving = this->GetMaxMoving() / this->GetNumBin();

            int mi = m;
            int ni = n;
            int oi = o;

            if (subpixel) {
                //interpolation 
                mi *= 2;
                ni *= 2;
                oi *= 2;
                int szi = mi*ni*oi;
                float* x1 = new float[szi];
                float* y1 = new float[szi];
                float* z1 = new float[szi];
                movingi = new float[szi];
                movingiInd = new int[szi];
                for (int k = 0; k < oi; k++) {
                    for (int j = 0; j < ni; j++) {
                        for (int i = 0; i < mi; i++) {
                            x1[i + j * mi + k * mi * ni] = 0.5 * (float) j;
                            y1[i + j * mi + k * mi * ni] = 0.5 * (float) i;
                            z1[i + j * mi + k * mi * ni] = 0.5 * (float) k;
                        }
                    }
                }
                interpolate3d(movingi, moving, x1, y1, z1, mi, ni, oi, m, n, o, false);
                delete []x1;
                delete []y1;
                delete []z1;
                for (int i = 0; i < lengthfull; i++) {
                    xs[i] *= 2.0;
                    ys[i] *= 2.0;
                    zs[i] *= 2.0;
                }


            } else {
                movingi = new float[sz];
                movingiInd = new int[sz];
                for (int i = 0; i < sz; i++) {
                    movingi[i] = moving[i];
                }


            }

            fixed = new float[m * n * o];
            fixedInd = new int[m * n * o];
            for (int i = 0; i < m * n * o; i++) {

                fixed[i] = imgFix[i];
                fixedInd[i] = (int) floor(fixed[i] / binSizeFixed);
            }

            int samples = this->GetSampleNum();
            bool randommode = samples < pow(gridstep, 3);
            int maxsamp;
            if (randommode) {
                maxsamp = samples;
            } else {
                maxsamp = (int) pow(gridstep, 3);
            }
            float* cost1 = new float[lengthfull];
            float* costcount = new float[lengthfull];

            int frac = (int) (sz1 / 25);

            float epsilon = gamma / (float) maxsamp;
            int xind, yind, zind;

            {
                float* normalizerA = new float[lengthfull];
                float* normalizerB = new float[lengthfull];

                for (int i = 0; i < sz1; i++) {
                    if ((i % frac) == 0) {
                        cout << "x" << flush;
                    }
                    int z1 = i / (gridx * gridy);
                    int x1 = (i - z1 * gridx * gridy) / gridx;
                    int y1 = i - z1 * gridx * gridy - x1*gridx;

                    z1 *= gridstep;
                    x1 *= gridstep;
                    y1 *= gridstep;

                    bool boundaries = true; //check boundaries 
                    if (subpixel) {
                        if (x1 * 2 + (gridstep - 1)*2 + hw_local >= ni | y1 * 2 + (gridstep - 1)*2 + hw_local >= mi | z1 * 2 + (gridstep - 1)*2 + hw_local >= oi)
                            boundaries = false;
                        if (x1 * 2 - hw_local < 0 | y1 * 2 - hw_local < 0 | z1 * 2 - hw_local < 0)
                            boundaries = false;
                    }
                    else {
                        if (x1 + (gridstep - 1) + hw_local >= ni | y1 + (gridstep - 1) + hw_local >= mi | z1 + (gridstep - 1) + hw_local >= oi)
                            boundaries = false;
                        if (x1 - hw_local < 0 | y1 - hw_local < 0 | z1 - hw_local < 0)
                            boundaries = false;
                    }


                    for (int l = 0; l < lengthfull; l++) {
                        cost1[l] = 0.0;
                        normalizerA[l] = 0.0;
                        normalizerB[l] = 0.0;
                    }


                    for (int j1 = 0; j1 < maxsamp; j1++) {
                        int i1;
                        if (randommode)
                            i1 = (int) (rand() * pow(gridstep, 3) / float(RAND_MAX));
                        else
                            i1 = j1;
                        int zz = i1 / (gridstep * gridstep);
                        int xx = (i1 - zz * gridstep * gridstep) / gridstep;
                        int yy = i1 - zz * gridstep * gridstep - xx*gridstep;

                        xx += x1;
                        yy += y1;
                        zz += z1;

                        for (int l = 0; l < lengthfull; l++) {
                            if (not(boundaries)) {
                                if (subpixel) {
                                    xind = max(min(xx * 2 + (int) xs[l], ni - 1), 0);
                                    yind = max(min(yy * 2 + (int) ys[l], mi - 1), 0);
                                    zind = max(min(zz * 2 + (int) zs[l], oi - 1), 0);
                                } else {
                                    xind = max(min(xx + (int) (xs[l]), ni - 1), 0);
                                    yind = max(min(yy + (int) (ys[l]), mi - 1), 0);
                                    zind = max(min(zz + (int) (zs[l]), oi - 1), 0);
                                }
                            } else {
                                if (subpixel) {
                                    xind = xx * 2 + (int) xs[l];
                                    yind = yy * 2 + (int) ys[l];
                                    zind = zz * 2 + (int) zs[l];
                                } else {
                                    xind = xx + (int) xs[l];
                                    yind = yy + (int) ys[l];
                                    zind = zz + (int) zs[l];
                                }
                            }
                            cost1[l] += fixed[yy + xx * m + zz * m * n] * movingi[yind + xind * mi + zind * mi * ni];
                            normalizerA[l] += powf(fixed[yy + xx * m + zz * m * n], 2.0);
                            normalizerB[l] += powf(movingi[yind + xind * mi + zind * mi * ni], 2.0);
                        }
                    }

                    for (int l = 0; l < lengthfull; l++) {
                        energy[i + l * sz1] = epsilon * (1 - pow(double(cost1[l] * cost1[l] / normalizerA[l] / normalizerB[l]), 2));
                    }
                }

                delete []normalizerA;
                delete []normalizerB;
            }

            delete []movingi;
            delete []fixed;
            delete []movingiInd;
            delete []fixedInd;

            delete []cost1;
            delete []costcount;
            delete []xs;
            delete []ys;
            delete []zs;

            //return 0;
        }

        template<class TPixel>
        void
        ProjectionLocalAlignmentDispatcher<TPixel>::
        updateImage(float* warped, float* imgFix, float* imgMov, float* u1, float* v1, float* w1) {
            int m = this->GetImageSize0();
            int n = this->GetImageSize1();
            int o = this->GetImageSize2();


            float metric = 0;
            float metric_start = 0;


            interpolate3d(warped, imgFix, u1, v1, w1, m, n, o, m, n, o, true);

            for (int i = 0; i < m; i++) {
                for (int j = 0; j < n; j++) {
                    for (int k = 0; k < o; k++) {
                        metric += pow(imgMov[i + j * m + k * m * n] - warped[i + j * m + k * m * n], 2);
                        metric_start += pow(imgMov[i + j * m + k * m * n] - imgFix[i + j * m + k * m * n], 2);
                    }
                }
            }

            metric /= m * n*o;
            metric_start /= m * n*o;

            this->SetMetricValue0(metric_start);
            this->SetMetricValue1(metric);
        }

        template<class TPixel>
        void
        ProjectionLocalAlignmentDispatcher<TPixel>::
        deformableReg(std::string modelDirectory, std::string outputDirectory, int debugLevel) {


            // prepare 
            this->PopulateVolume();
            typename Affine3DTransformType::Pointer vtoa = this->GetAtlasToVolumeTransform();
            vtoa->Print(std::cout);

            typename RegistrationImageType::Pointer sagittalAtlas;
            std::string fname;
            fname = modelDirectory + "/P56/atlasVolume/average_template_25.nrrd";
            itk::idp::ReadImage<RegistrationImageType>(fname.c_str(), sagittalAtlas);

            // resample 	
            typedef itk::ResampleImageFilter< VolumeType, VolumeType > ResampleFilterType;
            typename ResampleFilterType::Pointer resampler = ResampleFilterType::New();
            resampler->SetInput(this->GetVolume());

            resampler->SetTransform(vtoa.GetPointer());
            resampler->SetSize(sagittalAtlas->GetLargestPossibleRegion().GetSize());
            resampler->SetOutputOrigin(sagittalAtlas->GetOrigin());
            resampler->SetOutputSpacing(sagittalAtlas->GetSpacing());
            resampler->SetOutputDirection(sagittalAtlas->GetDirection());
            resampler->SetDefaultPixelValue(0);
            resampler->Update();

            typedef itk::CastImageFilter< VolumeType, RegistrationImageType > CastFilterType;
            typename CastFilterType::Pointer castFilter = CastFilterType::New();
            castFilter->SetInput(resampler->GetOutput());
            castFilter->Update();


            // output image
            if (debugLevel > 2) {
                fname = outputDirectory;
                fname += "/alignmentInput.nii.gz";
                itk::idp::WriteImage< RegistrationImageType >(fname.c_str(), castFilter->GetOutput());
            }


            using namespace std;

            std::string outputflow1, outputdef;

            outputflow1 = outputDirectory + "/grids.raw";
            outputdef = outputDirectory + "/deformed.nii.gz";

            char* outputflow = new char[200];
            strcpy(outputflow, outputflow1.c_str());

            float lambda = 1e-5;
            int samplenum = 1024;
            int numLevels = 4;
            int* grid_size = new int[numLevels];
            grid_size[0] = 30;

            for (int i = 1; i < numLevels; i++) {
                grid_size[i] = ceil(grid_size[i - 1] / 2);
            }


            this->SetSampleNum(samplenum);

            float* imgMov;
            float* imgFix;
            int xsize, ysize, zsize;

            // get image size
            const RegistrationImageType::SizeType sizeOfImage = sagittalAtlas->GetLargestPossibleRegion().GetSize();
            xsize = sizeOfImage[0];
            ysize = sizeOfImage[1];
            zsize = sizeOfImage[2];
            std::cout << "Input image size: " << xsize << "*" << ysize << "*" << zsize << std::endl;

            // get image
            imgMov = new float [xsize * ysize * zsize];
            imgFix = new float [xsize * ysize * zsize];
            RegistrationImageType::IndexType Index;
            for (long i = 0; i < zsize; i++) {
                for (long j = 0; j < ysize; j++) {
                    for (long k = 0; k < xsize; k++) {
                        Index[0] = k;
                        Index[1] = j;
                        Index[2] = i;
                        imgFix[k + j * xsize + i * xsize * ysize] = sagittalAtlas->GetPixel(Index);
                        imgMov[k + j * xsize + i * xsize * ysize] = castFilter->GetOutput()->GetPixel(Index);
                    }
                }
            }


            this->SetImageSize0(xsize);
            this->SetImageSize1(ysize);
            this->SetImageSize2(zsize);

            int xsize_start = xsize;
            int ysize_start = ysize;
            int zsize_start = zsize;
            int volsize = xsize_start*ysize_start*zsize_start;

            float *warpedImgFix = new float[xsize_start * ysize_start * zsize_start];
            float *warpedImgMov = new float[xsize_start * ysize_start * zsize_start];


            //==========================================================================================
            int* grid_hw = new int[numLevels];
            for (int i = 0; i < numLevels; i++) {
                grid_hw[i] = ceil(grid_size[i] / 2);
            }

            float* label_size_min = new float[numLevels];
            for (int i = 0; i < numLevels; i++) {
                label_size_min[i] = 1.0;
            }
            label_size_min[numLevels - 1] = 0.5;

            int stepsize;
            int hwidth;
            float minsize;

            // initialize the dfmfld
            float* xfld = new float[volsize];
            float* yfld = new float[volsize];
            float* zfld = new float[volsize];
            float* xfldinv = new float[volsize];
            float* yfldinv = new float[volsize];
            float* zfldinv = new float[volsize];
            for (int i = 0; i < volsize; i++) {
                xfld[i] = 0.0;
                yfld[i] = 0.0;
                zfld[i] = 0.0;
                xfldinv[i] = 0.0;
                yfldinv[i] = 0.0;
                zfldinv[i] = 0.0;
            }
            
            int xgridsize, ygridsize, zgridsize, gridvolsize;
            int xsizethis, ysizethis, zsizethis, volsizethis;
            xgridsize = xsize_start / grid_size[0];
            ygridsize = ysize_start / grid_size[0];
            zgridsize = zsize_start / grid_size[0];
            gridvolsize = xgridsize*ygridsize*zgridsize;
            float* xlabel = new float[gridvolsize];
            float* ylabel = new float[gridvolsize];
            float* zlabel = new float[gridvolsize];
            float* xlabelinv = new float[gridvolsize];
            float* ylabelinv = new float[gridvolsize];
            float* zlabelinv = new float[gridvolsize];
            for (int i = 0; i < gridvolsize; i++) {
                xlabel[i] = 0.0;
                ylabel[i] = 0.0;
                zlabel[i] = 0.0;
                xlabelinv[i] = 0.0;
                ylabelinv[i] = 0.0;
                zlabelinv[i] = 0.0;
            }

            //==========================================================================================	
            for (int level = 0; level < numLevels; level++) {


                updateImage(warpedImgFix, imgMov, imgFix, xfld, yfld, zfld);
                updateImage(warpedImgMov, imgFix, imgMov, xfldinv, yfldinv, zfldinv);

                minsize = label_size_min[level];
                stepsize = grid_size[level];
                hwidth = grid_hw[level];

                int gridvolsize = (int) pow(hwidth * 2 + 1, 3.0);
                xsizethis = xsize_start / stepsize;
                ysizethis = ysize_start / stepsize;
                zsizethis = zsize_start / stepsize;
                volsizethis = xsizethis * ysizethis*zsizethis;

                // spacing
                float* xfldstart = new float[volsizethis];
                float* yfldstart = new float[volsizethis];
                float* zfldstart = new float[volsizethis];
                float* xfldstartinv = new float[volsizethis];
                float* yfldstartinv = new float[volsizethis];
                float* zfldstartinv = new float[volsizethis];
                InterpDeformationFld(xfldstart, yfldstart, zfldstart, xlabel, ylabel, zlabel, xsizethis, ysizethis, zsizethis, xgridsize, ygridsize, zgridsize);
                InterpDeformationFld(xfldstartinv, yfldstartinv, zfldstartinv, xlabelinv, ylabelinv, zlabelinv, xsizethis, ysizethis, zsizethis, xgridsize, ygridsize, zgridsize);

                cout << "***********************************************************************************\n";
                cout << "Level #" << level << "   Grid size: " << stepsize << "\n";

                // release some memory
                delete xlabel;
                delete ylabel;
                delete zlabel;
                delete xlabelinv;
                delete ylabelinv;
                delete zlabelinv;

                xlabel = new float[volsizethis];
                ylabel = new float[volsizethis];
                zlabel = new float[volsizethis];
                xlabelinv = new float[volsizethis];
                ylabelinv = new float[volsizethis];
                zlabelinv = new float[volsizethis];

                // data structure 
                int* optimizedForward = new int[volsizethis];
                int* relationForward = new int[volsizethis];
                DfmfldModulation(imgFix, optimizedForward, relationForward, stepsize);
                int* optimizedBackward = new int[volsizethis];
                int* relationBackward = new int[volsizethis];
                DfmfldModulation(imgMov, optimizedBackward, relationBackward, stepsize);

                cout << "Begin computation ...\n";
                //cout<<"***********************************************************************************\n";		

                // optimization 
                float* costforward = new float[volsizethis * gridvolsize];
                float* costbackward = new float[volsizethis * gridvolsize];

                optimization(imgFix, warpedImgFix, costforward, lambda, hwidth, stepsize, minsize);
                optimization(imgMov, warpedImgMov, costbackward, lambda, hwidth, stepsize, minsize);

                genDfm(xlabel, ylabel, zlabel, xfldstart, yfldstart, zfldstart, costforward, lambda, hwidth, stepsize, minsize, optimizedForward, relationForward);
                genDfm(xlabelinv, ylabelinv, zlabelinv, xfldstartinv, yfldstartinv, zfldstartinv, costbackward, lambda, hwidth, stepsize, minsize, optimizedBackward, relationBackward);

                // smooth transformations 
                smoothDfm(xlabel, ylabel, zlabel, xsizethis, ysizethis, zsizethis, 4, stepsize);
                smoothDfm(xlabelinv, ylabelinv, zlabelinv, xsizethis, ysizethis, zsizethis, 4, stepsize);
                symmetricMapping(xlabel, ylabel, zlabel, xlabelinv, ylabelinv, zlabelinv, xsizethis, ysizethis, zsizethis, stepsize);


                //interp deformations from grid-resolution to high-resolution
                InterpDeformationFld(xfld, yfld, zfld, xlabel, ylabel, zlabel, xsize_start, ysize_start, zsize_start, xsizethis, ysizethis, zsizethis);
                InterpDeformationFld(xfldinv, yfldinv, zfldinv, xlabelinv, ylabelinv, zlabelinv, xsize_start, ysize_start, zsize_start, xsizethis, ysizethis, zsizethis);

                updateImage(warpedImgFix, imgMov, imgFix, xfld, yfld, zfld);
                cout << "\n" << "before this level: " << this->GetMetricValue0() << "\n" << "after this level: " << this->GetMetricValue1() << "\n";
                xgridsize = xsizethis;
                ygridsize = ysizethis;
                zgridsize = zsizethis;
                cout << "\n";

                // change the sampling rate per level
                this->SetSampleNum((int) this->GetSampleNum() / 4);
                if (this->GetSampleNum() < 64) this->SetSampleNum(64);

                delete []xfldstart;
                delete []yfldstart;
                delete []zfldstart;
                delete []xfldstartinv;
                delete []yfldstartinv;
                delete []zfldstartinv;
                delete []costforward;
                delete []costbackward;
                delete []relationForward;
                delete []optimizedForward;
                delete []relationBackward;
                delete []optimizedBackward;

            }

            // write out the results
            float* flow = new float[volsizethis * 3];
            for (int i = 0; i < volsizethis; i++) {
                flow[i] = xlabel[i];
                flow[i + volsizethis] = ylabel[i];
                flow[i + volsizethis * 2] = zlabel[i];
            }
            writeOutput(flow, outputflow, volsizethis * 3);

            // write deformed image 
            RegistrationImageType::Pointer warpedImg = RegistrationImageType::New();
            warpedImg->SetRegions(sagittalAtlas->GetLargestPossibleRegion());
            warpedImg->CopyInformation(sagittalAtlas);
            warpedImg->Allocate();
            for (long i = 0; i < zsize; i++) {
                for (long j = 0; j < ysize; j++) {
                    for (long k = 0; k < xsize; k++) {
                        Index[0] = k;
                        Index[1] = j;
                        Index[2] = i;
                        warpedImg->SetPixel(Index, warpedImgFix[k + j * xsize + i * xsize * ysize]);
                    }
                }
            }

            // cast to io image type
            typedef itk::CastImageFilter< RegistrationImageType, IOImageType > InvCastFilterType;
            typename InvCastFilterType::Pointer castFilterInv = InvCastFilterType::New();
            castFilterInv->SetInput(warpedImg);

            typedef itk::ImageFileWriter< IOImageType > WriterType;
            WriterType::Pointer writer = WriterType::New();

            writer->SetFileName(outputdef.c_str());
            writer->SetInput(castFilterInv->GetOutput());
            writer->Update();

            cout << "Registration finished!\n";

            delete []grid_size;
            delete []grid_hw;
            delete []label_size_min;

            // end of registration
        }

        template<class TPixel>
        void
        ProjectionLocalAlignmentDispatcher<TPixel>::
        ResampleDeformationFld(char* controlGrid, std::string modelDirectory, const RegistrationImageType::Pointer fldX, const RegistrationImageType::Pointer fldY, const RegistrationImageType::Pointer fldZ) {

            // read atlas 
            typename RegistrationImageType::Pointer sagittalAtlas;
            std::string fname;
            fname = modelDirectory + "/P56/atlasVolume/average_template_25.nrrd";
            itk::idp::ReadImage<RegistrationImageType>(fname.c_str(), sagittalAtlas);

            using namespace std;

            FILE * pFile;
            long sizefile;

            // open grid file and get the size
            pFile = fopen(controlGrid, "rb");
            if (pFile == NULL) perror("Error opening file");
            else {
                fseek(pFile, 0, SEEK_END);
                sizefile = ftell(pFile);
                fclose(pFile);
            }
            sizefile /= 4;
            float* grids = new float[sizefile];
            readRaw(controlGrid, grids, sizefile);
            int channelSZ = sizefile / 3;
            float* u1 = new float[channelSZ];
            float* v1 = new float[channelSZ];
            float* w1 = new float[channelSZ];
            for (int i = 0; i < channelSZ; i++) {
                u1[i] = grids[i];
                v1[i] = grids[i + channelSZ];
                w1[i] = grids[i + channelSZ * 2];
            }

            RegistrationImageType::SizeType size = sagittalAtlas->GetLargestPossibleRegion().GetSize();
            RegistrationImageType::SpacingType spacing = sagittalAtlas->GetSpacing();

            int m = size[0];
            int n = size[1];
            int o = size[2];
            int sz = m * n*o;
            int step1 = round(powf((float) sz / (float) channelSZ, 0.3333333));
            cout << "grid-size: " << step1 << "\n";
            cout << "tempalte size: " << size[0] << "x" << size[1] << "x" << size[2] << "\n";
            int m1 = m / step1;
            int n1 = n / step1;
            int o1 = o / step1;

            // dfm-fields
            float* ux = new float[sz];
            float* vx = new float[sz];
            float* wx = new float[sz];
            InterpDeformationFld(ux, vx, wx, u1, v1, w1, m, n, o, m1, n1, o1);

            // add the spacing information 
            RegistrationImageType::IndexType Index;
            fldX->SetRegions(sagittalAtlas->GetLargestPossibleRegion());
            fldY->SetRegions(sagittalAtlas->GetLargestPossibleRegion());
            fldZ->SetRegions(sagittalAtlas->GetLargestPossibleRegion());
            fldX->CopyInformation(sagittalAtlas);
            fldY->CopyInformation(sagittalAtlas);
            fldZ->CopyInformation(sagittalAtlas);
            fldX->Allocate();
            fldY->Allocate();
            fldZ->Allocate();

            for (long i = 0; i < o; i++) {
                for (long j = 0; j < n; j++) {
                    for (long k = 0; k < m; k++) {
                        Index[0] = k;
                        Index[1] = j;
                        Index[2] = i;
                        fldX->SetPixel(Index, vx[k + j * m + i * m * n] * spacing[0]);
                        fldY->SetPixel(Index, ux[k + j * m + i * m * n] * spacing[1]);
                        fldZ->SetPixel(Index, wx[k + j * m + i * m * n] * spacing[2]);
                    }
                }
            }
        }

    } // end namespace idp
} //end namespace itk

#endif
