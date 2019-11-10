#include <bits/stdc++.h>
#include <omp.h>

#include "algos.h"

using namespace std;

const double TwoPi = 4 * acosf(0);

// последовательно
void my_algos::Seq(double *AVal, double *FTvl, int len) {

    int i, j, n, m, Mmax, Istp;
    double Wtmp, Theta, Tmpr;
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
        Wpi = sin(Theta);
        Wtmp = sin(Theta / 2);
        Wpr = Wtmp * Wtmp * 2;
        Istp = Mmax * 2;
        Wr = 1;
        Wi = 0;
        m = 1;

        while (m < Mmax) {
            i = m;
            m = m + 2;
            double Tmpr = Wr;
            double Tmpi = Wi;
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

void as_bin(int u, int Z) {
    /* бинарное представление числа u в контексте Z бит */
    for (int i = Z - 1; i >= 0; i--) {
        if ((u & (1 << i)) != 0)
            cout << 1;
        else
            cout << 0;
    }
    cout << endl;
}

int compute(int kp1, int Z) {
    /* вспомогательная функция для вычисления второго итератора */
    int res = 0;
    for (int i = 0; i < Z; i++) {
        if ((1 << (Z - i - 1)) & kp1)
            res |= (1 << i);
    }
    return res + 1;
}


void my_algos::Vec(double *AVal, double *FTvl, int len) {

    int n = len * 2;
    int to_mmax[30], Z;

    for (int i = 0; i < 30; i++) {
        to_mmax[i] = (1 << i);
        if (to_mmax[i] >= n) {
            Z = i;
            break;
        }
    }

    for (int i = 0; i < n; i += 2) {
        FTvl[i] = 0;
        FTvl[i + 1] = AVal[i / 2];
    }

#pragma omp parallel for
    for (int k = 0; k < len; k++) {
        int i = 2 * k + 3;
        int j = compute(k + 1, Z);

        if (j > i) {
            double Tmpr = FTvl[i];
            FTvl[i] = FTvl[j];
            FTvl[j] = Tmpr;

            Tmpr = FTvl[i + 1];
            FTvl[i + 1] = FTvl[j + 1];
            FTvl[j + 1] = Tmpr;
        }
    }

    omp_set_nested(0);

    for (int y = 1; y < Z; y++) {

        int Mmax = to_mmax[y];

        double Theta = -TwoPi / Mmax;
        double Wpi = sin(Theta);
        double Wtmp = sin(Theta / 2);
        double Wpr = Wtmp * Wtmp * 2;
        int Istp = Mmax * 2;
        double Wr = 1;
        double Wi = 0;

        for (int m = 1; m < Mmax; m += 2) {
            double Tmpr = Wr;
            double Tmpi = Wi;
            Wr -= Tmpr * Wpr + Tmpi * Wpi;
            Wi += Tmpr * Wpi - Tmpi * Wpr;

            #pragma omp parallel for private(Tmpr) private(Tmpi)
            for (int i = m; i < n; i += Istp) {
                int j = i + Mmax;
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