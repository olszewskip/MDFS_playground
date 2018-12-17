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
            for (int j=0; j<dim1; j++) {
                M[i][j] = new int[dim2];
            }
    }
}


int main(int argc, char **argv) {

    size_t dim0, dim1, dim2;
    dim0 = 2;
    dim1 = 3;
    dim2 = 1;

    int ***A = matrix3D(dim0, dim1, dim2);

    A[1][2][0] = 123;
    printf("%d\n", A[1][2][0]);

    del3D(A, 2,3,1);

    int*** B = NULL;
    spawn3D(B, 2,3,1);
    B[1][2][0] = 123; // error

    return 0;
} 