#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int isNumber(char const *str) {

    int len = strlen(str);
    for (int i = 0; i < len; i++) {
        if(!(str[i] >= '0' && str[i] <= '9'))
            return 0;
    }
    
    return 1;
}

double** createMatrix(int r, int c) {
    double** matrix = (double**) malloc(sizeof(double*) * r);
    for (int i = 0; i < r; i++) {
        matrix[i] = (double*) malloc(c * sizeof(double));
    }
    
    for (int i = 0; i < r; i++) {
        for (int j = 0; j < c; j++) {
            matrix[i][j] = rand() % 6;
        }
    }
    
    return matrix;
}

double* createVector(int col) {
    double* vector = (double*) malloc(col * sizeof(double));
    for (int i = 0; i < col; i++) {
        vector[i] = rand() % 6;
    }
    
    return vector;
}

void printMatrix(double** matrix, int row, int col) {
    for (int r = 0; r < row; r++) {
        for (int c = 0; c < col; c++) {
            printf("%f ", matrix[r][c]);
        }
        printf("\n");
    }
    printf("\n");
}

void printVector(double* vector, int col) {
    for (int i = 0; i < col; i++) {
        printf("%f\n", vector[i]);
    }
    printf("\n");
}

double* multiply(double** matrix, double* vector, int row, int col) {
    double* result = (double*) malloc(row * sizeof(double));
    
    for (int i = 0; i < row; i++) {
        double a = 0;
        for (int j = 0; j < col; j++) {
            a += matrix[i][j] * vector[j];
        }
        result[i] = a;
    }
    
    return result;
}

void printFile(char const* fileName, double* result, int row) {
    FILE *out = fopen(fileName, "w");
    for (int i = 0; i < row; i++) {
        fprintf(out, "%f\n", result[i]);
    }
    fclose(out);
}

int main(int argc, char const *argv[]) {
    
    if ( argc != 4 ) {
        printf("Unexpected number of arguments!\n");
        exit(1);
    }
    
    if ( !isNumber(argv[1]) || !isNumber(argv[2]) ) {
        printf("Please enter first two arguments as a number!\n");
        exit(1);
    }
    
    int row = atoi(argv[1]);
    int col = atoi(argv[2]);
    
    srand(111032);
    
    double** matrix = createMatrix(row, col);
    double* vector = createVector(col);
    
    printMatrix(matrix, row, col);
    printVector(vector, col);
    
    double* result = multiply(matrix, vector, row, col);
    
    printVector(result, row);
    printFile(argv[3], result, row);
    
    free(matrix);
    free(vector);
    free(result);
    
    return 0;
}