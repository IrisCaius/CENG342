/*
@author     İbrahim Gülmez 19050111032
@details    In this program matrix - vector multiplication is implemented with parallel programming.
            In the implementation OpenMPI's mpi.h library is used.
            Communicatios between processes is done with MPI_Send and MPI_Recv functions.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include </usr/lib/x86_64-linux-gnu/openmpi/include/mpi.h>
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

void asssign_random_numbers_to_matrix(double *matrix, int r, int c) {
    /// @authors https://cboard.cprogramming.com/c-programming/17939-random-double-numbers-post114817.html#post114817
    for (int i = 0; i < r; i++) {
        for (int j = 0; j < c; j++) {
            matrix[i*r + j] = ( (double)rand() * ( 10 -(-10) ) ) / (double)RAND_MAX + (-10);
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

void matrix_vector_multiplication(double *matrix, double *vector, double* result, int row, int col, int start_row) {
    // printf("Start Row: %d, Row: %d\n", start_row, row);
    for (int i = 0; i < row; i++) {
        double a = 0;
        for (int j = 0; j < col; j++) {
            a += matrix[i*row + j] * vector[j];
        }
        result[start_row + i] = a;
    }
}

void printFile(char const* fileName, double* result, int n1, int n2) {
    FILE *out = fopen(fileName, "w");
    for (int i = 0; i < n1; i++) {
        fprintf(out, "%f\n", result[i]);
    }
    fclose(out);
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
    
    double result[n2];
    
    int start_row = (int) (rank * (double)n1 / (double)size);  
    int end_row = (int) ((rank + 1) * (double)n1 / (double)size);
    int local_num_rows = end_row - start_row + 1;
    
    // printf("Rank: %d, [%d, %d)\n", rank, start_row, end_row);
    
    MPI_Barrier( MPI_COMM_WORLD );
    double start;
    if ( rank == 0 ) {
        // Since I got type error (passed array -double (*)[n1]) (function parameter - double**),
        // I used 1D array as 2D array. This technique has gotten from Pacheo's examples
        double *matrix = (double*) malloc(n1*n2*sizeof(double));
        double *vector = (double*) malloc(n2*sizeof(double));
        
        asssign_random_numbers_to_matrix(matrix, n1, n2);
        asssign_random_numbers_to_vector(vector, n2);
        
        start = getTime();
        
        // Send generated matrixes and vectors to the other processes
        for (int i = 1; i < size; i++) {
            int low = (int) (i * (double)n1 / (double)size);  
            int high = (int) ((i + 1) * (double)n1 / (double)size);
            int num_of_elements = (high - low + 1) * n2;
            // printf("Lİne 133 %d\n", num_of_elements);
            MPI_Send(&matrix[low], num_of_elements, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
            MPI_Send(vector, n2, MPI_DOUBLE, i, 1, MPI_COMM_WORLD);
        }
;
        // Calculate matrix vector multiplication for process 0's area
        matrix_vector_multiplication(matrix, vector, result, local_num_rows, n2, start_row);
        
        for (int i = 1; i < size; i++) {
            double local_result[n2];
            int low = (int) (i * (double)n1 / (double)size);  
            int high = (int) ((i + 1) * (double)n1 / (double)size);
            
            // Receiving other process' results
            MPI_Recv(local_result, n2, MPI_DOUBLE, i, i * 100 + 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            // Assigning other process' results to process 0's result
            for (int i = low; i < high; i++) {
                result[i] = local_result[i];
            }
        }
        
        // printVector(result, n1);
        free(matrix);
        free(vector);
    } 
    else {
        double matrix[local_num_rows * n2], vector[n2];
        // Receive generated matrixes and vectors from process 0  to the other processes
        MPI_Recv(matrix, local_num_rows * n2, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(vector, n2, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        
        // Calculate matrix vector multiplication for other processes' area
        matrix_vector_multiplication(matrix, vector, result, local_num_rows, n2, start_row);
        
        // Send calculated matrix vector multiplication results to process 0
        MPI_Send(result, n2, MPI_DOUBLE, 0, rank * 100 + 1, MPI_COMM_WORLD);
    }
    
    
    double finish = getTime();
    
    
    if ( rank == 0 ) {
        printf("Elapsed time is %.6f seconds for parallel mxv with %d processes\n", finish - start, size);
        printFile(argv[3], result, n1, n2);
    }
    
    MPI_Finalize();
    return 0;
}