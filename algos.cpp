#include <bits/stdc++.h>
#include <omp.h>

#include "algos.h"

using namespace std;

const double TwoPi = 4 * acosf(0);

// последовательно
void my_algos::Seq(double *AVal, double *FTvl, int len) {

    int i, j, n, m, Mmax, Istp;
    double Tmpr, Tmpi, Wtmp, Theta;
    double Wpr, Wpi, Wr, Wi;

    n = len * 2;

    for (i = 0; i < n; i += 2) {
        FTvl[i] = 0;
        FTvl[i + 1] = AVal[i / 2];
    }

    i = 1;
    j = 1;

    while (i < n) {
        if (j > i) {
            Tmpr = FTvl[i];
            FTvl[i] = FTvl[j];
            FTvl[j] = Tmpr;

            Tmpr = FTvl[i + 1];
            FTvl[i + 1] = FTvl[j + 1];
            FTvl[j + 1] = Tmpr;
        }
        i += 2;
        m = len;

        while ((m >= 2) && (j > m)) {
            j -= m;
            m >>= 1;
        }
        j += m;
    }
    Mmax = 2;

    while (Mmax < n) {

        Theta = -TwoPi / Mmax;
        Wpi = std::sin(Theta);
        Wtmp = sin(Theta / 2);
        Wpr = Wtmp * Wtmp * 2;
        Istp = Mmax * 2;
        Wr = 1;
        Wi = 0;
        m = 1;

        while (m < Mmax) {
            i = m;
            m = m + 2;
            Tmpr = Wr;
            Tmpi = Wi;
            Wr = Wr - Tmpr * Wpr - Tmpi * Wpi;
            Wi = Wi + Tmpr * Wpi - Tmpi * Wpr;

            while (i < n) {
                j = i + Mmax;
                Tmpr = Wr * FTvl[j] - Wi * FTvl[j - 1];
                Tmpi = Wi * FTvl[j] + Wr * FTvl[j - 1];

                FTvl[j] = FTvl[i] - Tmpr;
                FTvl[j - 1] = FTvl[i - 1] - Tmpi;
                FTvl[i] = FTvl[i] + Tmpr;
                FTvl[i - 1] = FTvl[i - 1] + Tmpi;
                i = i + Istp;
            }
        }
        Mmax = Istp;
    }
}


void my_algos::Vec(double *AVal, double *FTvl, int len) {

    int i, j, n, m, Istp;
    double Tmpr, Tmpi, Wtmp, Theta;
    double Wpr, Wpi, Wr, Wi;

    n = len * 2;

    for (i = 0; i < n; i += 2) {
        FTvl[i] = 0;
        FTvl[i + 1] = AVal[i / 2];
    }

    i = 1;
    j = 1;

    while (i < n) {
        if (j > i) {
            Tmpr = FTvl[i];
            FTvl[i] = FTvl[j];
            FTvl[j] = Tmpr;

            Tmpr = FTvl[i + 1];
            FTvl[i + 1] = FTvl[j + 1];
            FTvl[j + 1] = Tmpr;
        }
        i += 2;
        m = len;

        while ((m >= 2) && (j > m)) {
            j -= m;
            m >>= 1;
        }
        j += m;
    }

    int to_mmax[25], Z;
    for (i = 0; i < 25; i++) {
        to_mmax[i] = (1 << i);
        if (to_mmax[i] >= n) {
            Z = i;
            break;
        }
    }

#pragma omp parallel for
    for (int y = 1; y < Z; y++) {

        int Mmax = to_mmax[y];

        Theta = -TwoPi / Mmax;
        Wpi = sin(Theta);
        Wtmp = sin(Theta / 2);
        Wpr = Wtmp * Wtmp * 2;
        Istp = Mmax * 2;
        Wr = 1;
        Wi = 0;

        for (m = 1; m < Mmax; m += 2) {

            Tmpr = Wr;
            Tmpi = Wi;
            Wi += Tmpr * Wpi - Tmpi * Wpr;
            Wr -= Tmpr * Wpr + Tmpi * Wpi;

            for (i = m; i < n; i += Istp) {
                j = i + Mmax;
                Tmpr = Wr * FTvl[j] - Wi * FTvl[j - 1];
                Tmpi = Wi * FTvl[j] + Wr * FTvl[j - 1];

                FTvl[j] = FTvl[i] - Tmpr;
                FTvl[j - 1] = FTvl[i - 1] - Tmpi;
                FTvl[i] = FTvl[i] + Tmpr;
                FTvl[i - 1] = FTvl[i - 1] + Tmpi;
            }

        }
    }
}