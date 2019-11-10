/*
 * Напишите программу для вычисления быстрого преобразование Фурье.
 * Реализуйте алгоритм последовательно, последовательно с использованием
 * векторизации, с использованием параллельных вычислений на центральном
 * процессоре, с использованием вычислений на графическом процессоре.
 * Для каждого способа проверьте корректность преобразование и замерьте
 * время выполнения. Измерьте среднее время выполнения (от ста запусков),
 * проанализируйте полученные результаты. Сформулируйте зависимость времени
 * выполнения каждого способа от размера входной последовательности.
 * */

#include <bits/stdc++.h>
#include "algos.h"

using namespace std;
using namespace chrono;

void profiler(void (*func)(double *, double *, int),
              double *in, double *out, int len,
              char *alias, int times = 1) {

    auto start(high_resolution_clock::now());

    while (times--)
        (*func)(in, out, len);

    auto end(high_resolution_clock::now());
    auto duration(duration_cast<milliseconds>(end - start));

    printf("%s (%ld ms):\n", alias, duration.count());
}

void presenter(int len, double *output) {
    int i, j;
    for (i = j = 0; i < len; i++, j += 2) {
        printf("%.2f %c %.4fi%s",
               output[j + 1],
               output[j] >= 0 ? '+' : '-',
               abs(output[j]),
               i + 1 == len ? "." : ", "
        );
    }
    printf("\n\n\n");
}

int main() {

    freopen("input.txt", "rt", stdin);
    freopen("output.txt", "wt", stdout);

    int len, i;
    cin >> len;
    auto *input = new double[len];
    auto *output = new double[2 * len];

    for (i = 0; i < len; ++i) {
        cin >> input[i];
    }

    profiler(my_algos::Seq, input, output, len, (char *) "Seq");
    presenter(len, output);

    profiler(my_algos::Vec, input, output, len, (char *) "Vec");
    presenter(len, output);

    return 0;
}
