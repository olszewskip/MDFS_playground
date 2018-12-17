#include <iostream>
using namespace std;

int*** matrix3D(size_t dim0, size_t dim1, size_t dim2) {
// use size_t or int?

    int*** M = NULL;
    M = new int**[dim0];

    for(int i=0; i<dim0; i++) {
        M[i] = new int*[dim1];
            for(int j=0; j<dim1; j++) {
                M[i][j] = new int[dim2];
            }
    }

    return M;
}

void del3D(int***M, size_t dim0, size_t dim1, size_t dim2) {

    for (int i=0; i<dim0; i++) {
        for (int j=0; j<dim1; j++) {
            delete[] M[i][j];
        }
        delete[] M[i];
    }
    delete[] M;
}

// why doesn't it work?
void spawn3D(int*** M, int dim0, int dim1, int dim2) {

    M = new int**[dim0];
    
    for(int i=0; i<dim0; i++) {
        M[i] = new int*[dim1];
            for(int j=0; j<dim1; j++) {
                M[i][j] = new int[dim2];
            }
    }
}


int main(int argc, char **argv) {

    cout << "1)\n";
    int*** A = NULL;

    cout << "2)\n";
    A = matrix3D(2,3,1);

    cout << "3)\n";
    A[1][2][0] = 123;

    cout << "4)\n";
    printf("%d\n", A[1][2][0]);

    cout << "5)\n";
    del3D(A, 2,3,1);

    //cout << "6)\n";
    //printf("%d\n", A[1][2][0]); //error

    cout << "7)\n";
    int*** B = NULL;

    cout << "8)\n";
    spawn3D(B, 2,3,1);
    
    cout << "9)\n";
    B[1][2][0] = 123; //why error?

    cout << "10)\n";
    printf("%d\n", B[1][2][0]);

    return 0;
} 