#include <cstdio>
#include <iostream>
#include <cmath>
using namespace std;

double spectral_distance(double C[], double sp[], int n) {
    double d = 0;
    for (int i = 0; i < n; i++) {
        d += (C[i] - sp[i]) * (C[i] - sp[i]);
    }
    return sqrt(d);
}
