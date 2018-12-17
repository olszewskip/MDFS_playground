#include <stdlib.h>
#include <stdio.h>

int*** matrix3D(size_t dim0, size_t dim1, size_t dim2) {

    int ***M = NULL;
    M = malloc(dim0 * sizeof(*M));

    for (int i=0; i<dim0; i++){
        M[i] = malloc(dim1 * sizeof(**M));
        for (int j=0; j<dim1; j++) {
            M[i][j] = malloc(dim1 * sizeof(***M));
        }
    }

    return M; 
}

void del3D(int***M, size_t dim0, size_t dim1, size_t dim2) {

    for (int i=0; i<dim0; i++) {
        for (int j=0; j<dim1; j++) {
            free(M[i][j]);
        }
        free(M[i]);
    }
    free(M);
}

// why doesn't it work?
void spawn3D(int***M, size_t dim0, size_t dim1, size_t dim2) {

    M = malloc(dim0 * sizeof(*M));

    for (int i=0; i<dim0; i++){
        M[i] = malloc(dim1 * sizeof(**M));
        for (int j=0; j<dim1; j++) {
            M[i][j] = malloc(dim1 * sizeof(***M));
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

    del3D(A, dim0, dim1, dim2);

    int ***B = NULL;
    spawn3D(B, dim0, dim1, dim2);
    B[1][2][0] = 123; // error

    return 0;
}