#include <cstdio>
#include <cstdlib>
#include <iostream>

using namespace std;

extern int bco_srg(int, int, int);

int main() {
    for (int s = 1; s <= 100; s++) {
        bco_srg(s, 6, 30);
    }
    return 0;
}