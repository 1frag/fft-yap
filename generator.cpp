//
// Created by ifrag on 10.11.2019.
//

#include <bits/stdc++.h>
#include "generator.h"

using namespace std;

int main(int argc, char *argv[]) {
    int count = 10;

    if (argc > 1) {
        count = atoi(argv[1]);
    }

    freopen("input.txt", "wt", stdout);

    count = pow(2, count);
    cout << count << endl;

    for (int i = 0; i < count; i++) {
        cout << rand() % 2 << " ";
    }

}
