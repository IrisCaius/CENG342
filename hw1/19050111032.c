/*
@author     İbrahim Gülmez 19050111032
@details    In this program matrix - vector multiplication is implemented with parallel programming.
            In the implementation OpenMPI's mpi.h library is used.
            Communicatios between processes is done with MPI_Send and MPI_Recv functions.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>
#include <sys/time.h>

double getTime() {
    struct timeval t;
    gettimeofday(&t, NULL);
    return (double) t.tv_sec + (double) t.tv_usec / 1000000.0;
}

/// @brief      Checks whether the given string is numeric
/// @param str  The string will be controlled
/// @return     Returns 1 if the is numeric else returns 0
int isNumeric(char const *str) {
    int len = strlen(str);
    if( str[0] == '0' ) return 0;
    for (int i = 0; i < len; i++) {
        if(!(str[i] >= '0' && str[i] <= '9'))
            return 0;
    }
    return 1;
}

void asssign_random_numbers_to_matrix(double *matrix, int n, int r, int c) {
    /// @authors https://cboard.cprogramming.com/c-programming/17939-random-double-numbers-post114817.html#post114817
    for (int i = 0; i < r; i++) {
        for (int j = 0; j < c; j++) {
            matrix[i*n + j] = ( (double)rand() * ( 10 -(-10) ) ) / (double)RAND_MAX + (-10);
        }
    }
}

/// @brief          Creates vector that its values random double values between -10 and 10.
/// @param length   Length of the vector
/// @return         1D array with random double values.
void asssign_random_numbers_to_vector(double *vector, int length) {
    /// @authors https://cboard.cprogramming.com/c-programming/17939-random-double-numbers-post114817.html#post114817
    for (int i = 0; i < length; i++) {
        vector[i] = ( (double)rand() * ( 10 -(-10) ) ) / (double)RAND_MAX + (-10);
    }
}

void printMatrix(double** matrix, int row, int col) {
    for (int r = 0; r < row; r++) {
        for (int c = 0; c < col; c++) {
            printf("%f\t", matrix[r][c]);
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

void matrix_vector_multiplication(double *matrix, double* vector, double* result, int n, int row, int col, int rank, int size) {
    int low = (int) (rank * (double)row / (double)size);  
    int high = (int) ((rank + 1) * (double)row / (double)size);
    for (int i = low; i < high; i++) {
        double a = 0;
        for (int j = 0; j < col; j++) {
            a += matrix[i*n + j] * vector[j];
        }
        result[i] = a;
    }
}

void printFile(char const* fileName, double* result1, double *result2, int n1, int n2) {
    FILE *out = fopen(fileName, "w");
    for (int i = 0; i < n1; i++) {
        fprintf(out, "%f\n", result1[i]);
    }
    fprintf(out, "\n------------------\n");
    for (int i = 0; i < n2; i++) {
        fprintf(out, "%f\n", result2[i]);
    }
    fclose(out);
}

void parallel_matrix_vector_multiplication(double *matrix, double *vector, double *result, int n, int rank, int size) {
    if ( rank == 0 ) {
        asssign_random_numbers_to_matrix(matrix, n, n, n);
        asssign_random_numbers_to_vector(vector, n);
        
        // Send generated matrixes and vectors to the other processes
        for (int i = 1; i < size; i++) {
            MPI_Send(matrix, n*n, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
            MPI_Send(vector, n, MPI_DOUBLE, i, 1, MPI_COMM_WORLD);
        }
        
        // Calculate matrix vector multiplication for process 0's area
        matrix_vector_multiplication(matrix, vector, result, n, n, n, rank, size);
        
        for (int i = 1; i < size; i++) {
            double local_result[n];
            
            // Receiving other process' results
            MPI_Recv(local_result, n, MPI_DOUBLE, i, i * 100 + 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            
            // Assigning other process' results to process 0's result
            int low = (int) (i * (double)n / (double)size);  
            int high = (int) ((i + 1) * (double)n / (double)size);
            for (int i = low; i < high; i++) {
                result[i] = local_result[i];
            }
        }
        
        // printVector(result, n);
    } 
    else {
        // Receive generated matrixes and vectors from process 0  to the other processes
        MPI_Recv(matrix, n*n, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(vector, n, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        
        // Calculate matrix vector multiplication for other processes' area
        matrix_vector_multiplication(matrix, vector, result, n, n, n, rank, size);
        
        // Send calculated matrix vector multiplication results to process 0
        MPI_Send(result, n, MPI_DOUBLE, 0, rank * 100 + 1, MPI_COMM_WORLD);
    }
}

int main(int argc, char const *argv[]) {
    
    if ( argc != 4 ) {
        printf("Unexpected number of arguments!\n");
        exit(1);
    }
    
    // ./19050111032 1003 1003 output.txt
    if ( !isNumeric(argv[1]) || !isNumeric(argv[2]) ) {
        printf("Please enter third, fifth, sixth arguments as a number!\n");
        exit(1);
    }
    
    int rank, size;
    MPI_Init(NULL, NULL);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    int n1 = atoi(argv[1]);
    int n2 = atoi(argv[2]);
    
    // Since I got type error (passed array -double (*)[n1]) (function parameter - double**),
    // I used 1D array as 2D array. This technique has gotten from Pacheo's examples
    double matrix1[n1*n1];
    double matrix2[n2*n2];
    double vector1[n1], vector2[n2];
    double result1[n1], result2[n2];
    
    MPI_Barrier( MPI_COMM_WORLD );
    double start1 = getTime();
    parallel_matrix_vector_multiplication(matrix1, vector1, result1, n1, rank, size);
    double finish1 = getTime();
    
    MPI_Barrier( MPI_COMM_WORLD );
    double start2 = getTime();
    parallel_matrix_vector_multiplication(matrix2, vector2, result2, n2, rank, size);
    double finish2 = getTime();
    
    if ( rank == 0 ) {
        printf("Elapsed time is %.3f seconds for parallel mx1v1 with %d processes\n", finish1 - start1, size);
        printf("Elapsed time is %.3f seconds for parallel mx2v2 with %d processes\n", finish2 - start2, size);
        
        printFile(argv[3], result1, result2, n1, n2);
    }
    
    MPI_Finalize();
    return 0;
}