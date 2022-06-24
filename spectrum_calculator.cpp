#include <cstdio>
#include <cstdlib>
#include <iostream>

using namespace std;

extern void jacobi_eigenvalue(int, double[], int, double[], double[], int *, int *);

int main() {
    double *v = (double *)malloc(sizeof(double) * (400));
    if (v == NULL) {
        printf("Malloc");
        return 1;
    }
    double *sp = (double *)malloc(sizeof(double) * 20);
    if (sp == NULL) {
        printf("malloc");
        return 1;
    }
    int it_num;
    int rot_num;
    double A[400];
    for (int i = 0; i < 400; i++)
    {
        A[i] = 0;
    }
    
    for (int i = 1; i < 20; i++)
    {
        A[i*20] = 1;
        A[i] = 1;
    }
       

    jacobi_eigenvalue(20, A, 100000, v, sp, &it_num, &rot_num);

    for (int i = 0; i < 20; i++) {
        cout << sp[i] << ' ';
    }
    cout << endl;
    /*for (int i = 0; i < 20; i++) {
        for (int j = 0; j < 20; j++) {
            cout << A[i * 20 + j] << ' ';
        }
        cout << endl;
    }
    */
    cout << endl;

    return 0;
}