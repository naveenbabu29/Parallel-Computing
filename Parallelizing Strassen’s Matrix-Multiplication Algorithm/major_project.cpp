#define __is_convertible(F, T) is_convertible<F, T>::value

#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <time.h>
#include <ctime>
#include <vector>
#include <math.h>
#include <cmath>
#include <limits.h>
#include <omp.h>

using namespace std;

// Terminal Matrix size
int matrix_limit;

int **allocate_matrix(int row, int col)
{
    int **row_ptrs = new int *[row];  // Allocate memory for rows
    int *mat_data = new int[row *col];   // Allocate memory for each element

    // Assign pointer to each row
    for(int row_idx = 0; row_idx < row; row_idx++)
    {
        row_ptrs[row_idx] = mat_data + (row_idx *col);
    }

    return row_ptrs;
}

void deallocate_memory(int **matrix_arr)
{
    delete[] (*matrix_arr); // Free matrix elements
    delete matrix_arr; // Free row pointers
}

void normal_matrix_mul(int matrix_size, 
					   int **matrix_A,
					   int **matrix_B, 
					   int **matrix_C)
{    
    for(int i = 0; i < matrix_size; i++)
    {   
        for(int j = 0; j < matrix_size; j++)
        {
            int sum = 0;

            for(int k = 0; k < matrix_size; k++)
            {
                sum = sum + matrix_A[i][k] * matrix_B[k][j];
            }

            matrix_C[i][j] = sum;
        }
    }
}


void print_results(int error_count, 
				   int size,
				   int k1,
				   auto omp_get_max_threads,
				   auto elapsed_time)
{
    if(error_count == 0)
    {
		printf("Success: Matrix size = %d, K' = %d, No. of threads = %d, Execution time = %lf seconds, Errors = %d \n", size, k1, omp_get_max_threads, elapsed_time, error_count);
	}
	else
    {
		printf("Fail: Matrix size = %d, K' = %d, No. of threads = %d, Execution time = %lf seconds, Errors = %d \n", size, k1, omp_get_max_threads, elapsed_time, error_count);
	}    
}


void copy_matrix(int **destination, 
				 int dest_rows, 
				 int **source, 
				 int mf, 
				 int nf)
{
    for(int i = 0; i < dest_rows; i++)
    {
        destination[i] = source[mf + i] + nf;
    }
}


void add_matrix(int **result, 
				int rows, 
				int cols, 
				int **matrix_X, 
				int **matrix_Y)
{
    for(int row = 0; row < rows; row++)
    {
        for(int col = 0; col < cols; col++)
        {
            result[row][col] = matrix_X[row][col] + matrix_Y[row][col];
        }
    }
}


void subtract_matrix(int **result, 
					 int rows,
					 int cols, 
					 int **matrix_X, 
					 int **matrix_Y)
{
    for (int row = 0; row < rows; row++)
    {
        for (int col = 0; col < cols; col++)
        {
            *(*(result + row) + col) = *(*(matrix_X + row) + col) - *(*(matrix_Y + row) + col);
        }
    }
}


void strassen_multiplication(int matrix_size, 
							 int **first_matrix, 
							 int **second_matrix, 
							 int **C)
{
    // Standard algorithm for multiplication - Terminal condition
    if( matrix_size <= matrix_limit )
    {
        for(int i = 0; i < matrix_size; i++)
        {
            for(int j = 0; j < matrix_size; j++)
            {
                int sum = 0;
                for( int k = 0; k < matrix_size; k++ )
                {
                    sum = sum + first_matrix[i][k] * second_matrix[k][j];
                }

                C[i][j] = sum;
            }
        }
    }

    // Strassen's Recursive algorithm - regular case
    else
    {
        int x = matrix_size / 2;
        int y = matrix_size / 2;
        int z = matrix_size / 2;
        int sub_matrix_size = matrix_size / 2;

        int **m1 = allocate_matrix(x, y);
        int **weight_am_1 = allocate_matrix(x,z);
        int **weight_bm_6 = allocate_matrix(z, y);
		    
        int **weight_bm_1 = allocate_matrix(z, y);
        int **m2 = allocate_matrix(x, y);
        int **weight_am_7 = allocate_matrix(x, z);
		    
        int **weight_am_2 = allocate_matrix(x, z);
        int **weight_bm_7 = allocate_matrix(z, y);
        int **m3 = allocate_matrix(x, y);
		    
        int **weight_bm_3 = allocate_matrix(z, y);
        int **m4 = allocate_matrix(x, y);
        int **weight_bm_4 = allocate_matrix(z, y);
		    
        int **m5 = allocate_matrix(x, y);
        int **weight_am_5 = allocate_matrix(x, z);
        int **m6 = allocate_matrix(x, y);
		    
        int **m7 = allocate_matrix(x, y);
        int **weight_am_6 = allocate_matrix(x, z);

        int **first_11 = new int*[x];
        int **second_11 = new int*[z];
        int **C11 = new int*[x];
        double aa = 1;
        aa++;

        int **C12 = new int*[x];
        int **second_12 = new int*[z];
        int **first_12 = new int*[x];
		auto ab = 123;

        int **first_21 = new int*[x];
        int **second_21 = new int*[z];
        int **C21 = new int*[x];
		aa++; 
		ab--;

        int **C22 = new int*[x];
        int **second_22 = new int*[z];
        int **first_22 = new int*[x];

        copy_matrix(first_11, x, first_matrix,  0,  0);
        copy_matrix(C11, x, C,  0,  0);

        copy_matrix(second_21, z, second_matrix, z,  0);

        copy_matrix(second_12, z, second_matrix,  0, y);
        copy_matrix(C12, x, C,  0, y);
        
        copy_matrix(first_22, x, first_matrix, x, z);

        copy_matrix(second_11, z, second_matrix,  0,  0);
        copy_matrix(C21, x, C, x,  0);
        
        copy_matrix(first_21, x, first_matrix, x,  0);

        copy_matrix(second_22, z, second_matrix, z, y);
        copy_matrix(C22, x, C, x, y);
        
        copy_matrix(first_12, x, first_matrix,  0, z);

#pragma omp task
        {
            // M1 = (A11 + A22)*(B11 + B22)
            add_matrix(weight_am_1, x, z, first_11, first_22);
            add_matrix(weight_bm_1, z, y, second_11, second_22);
            strassen_multiplication(sub_matrix_size, weight_am_1, weight_bm_1, m1);
        }

#pragma omp task
		{
            // M2 = (A21 + A22) * B11
			add_matrix(weight_am_2, x, z, first_21, first_22);
			strassen_multiplication(sub_matrix_size, weight_am_2, second_11, m2);
		}

#pragma omp task
		{
            // M3 = A11 * (B12 - B22)
			subtract_matrix(weight_bm_3, z, y, second_12, second_22);
			strassen_multiplication(sub_matrix_size, first_11, weight_bm_3, m3);
		}

#pragma omp task
		{
            // M4 = A22* (B21 - B11)
			subtract_matrix(weight_bm_4, z, y, second_21, second_11);
			strassen_multiplication(sub_matrix_size, first_22, weight_bm_4, m4);
		}

#pragma omp task
		{
            // M5 = (A11 + A12) * B22
			add_matrix(weight_am_5, x, z, first_11, first_12);
			strassen_multiplication(sub_matrix_size, weight_am_5, second_22, m5);
		}

#pragma omp task
		{
            // M6 = (A21 - A11)*(B11 + B12)
			subtract_matrix(weight_am_6, x, z, first_21, first_11);
			add_matrix(weight_bm_6, z, y, second_11, second_12);
			strassen_multiplication(sub_matrix_size, weight_am_6, weight_bm_6, m6);
		}

#pragma omp task
		{
            // M7 = (A12 - A22)*(B21 + B22)
			subtract_matrix(weight_am_7, x, z, first_12, first_22);
			add_matrix(weight_bm_7, z, y, second_21, second_22);
			strassen_multiplication(sub_matrix_size, weight_am_7, weight_bm_7, m7);
		}

#pragma omp taskwait
#pragma omp parallel for
        for(int i = 0; i < x; i++)
		{
            for(int j = 0; j < y; j++)
            {
                C11[i][j] = m1[i][j] + m4[i][j] - m5[i][j] + m7[i][j];
                C21[i][j] = m2[i][j] + m4[i][j];
                C12[i][j] = m3[i][j] + m5[i][j];
                C22[i][j] = m1[i][j] - m2[i][j] + m3[i][j] + m6[i][j];
            }
		}

		deallocate_memory(m1);
		deallocate_memory(m2);
		deallocate_memory(m3);
		deallocate_memory(m4);
		deallocate_memory(m5);
		deallocate_memory(m6);
		deallocate_memory(m7);

		deallocate_memory(weight_am_1);
		deallocate_memory(weight_am_2);
		deallocate_memory(weight_am_5);
		deallocate_memory(weight_am_6);
		deallocate_memory(weight_am_7);

		deallocate_memory(weight_bm_1);
		deallocate_memory(weight_bm_3);
		deallocate_memory(weight_bm_4);
		deallocate_memory(weight_bm_6);
		deallocate_memory(weight_bm_7);

		delete[] first_11; 
		delete[] C21;
		delete[] second_11; 
		delete[] C11;
		delete[] first_22; 
		delete[] second_12;
		delete[] first_12; 
		delete[] second_21; 
		delete[] C12;
		delete[] C22; 
		delete[] first_21; 
		delete[] second_22;
    }

}

int count_error(int **A, int **B, int size)
{
    int error_count = 0;

    for(int i = 0; i < size; i++)
    {
        for(int j = 0; j < size; j++)
        {
            if(A[i][j] != B[i][j])
            {
                error_count++;
            }
        }
    }
    return error_count;
}


int main(int argc, char *argv[])
{
    int k = atoi(argv[1]);
    int k1 = atoi(argv[2]);

    if(k - k1 < 0)
    {
        printf("\nk' must be less than or equal to k, setting matrix limit to 1 \n");
		matrix_limit = 1;
    }
    else
    {
        matrix_limit = pow(2, (k-k1));
    }

    int num_threads = atoi(argv[3]);
    num_threads = pow(2, num_threads);
    
    int size = pow(2, k);

    int **mat_A = allocate_matrix(size, size);
    int **mat_B = allocate_matrix(size, size);
    int **mat_result = allocate_matrix(size, size);
    int **mat_test = allocate_matrix(size, size);

    srand((unsigned)time(NULL));

    for(int i = 0; i < size; i++)
	{
        for(int j = 0; j < size; j++)
        {
            mat_A[i][j] = rand()%100;
            mat_B[i][j] = rand()%100;
        }
	}

    // Disable dynamic adjustment of number of threads
    omp_set_dynamic(0);

    // Set number of threads
    omp_set_num_threads(num_threads);

    double start_time = omp_get_wtime();

#pragma omp parallel
    {
#pragma omp single
        {
            strassen_multiplication(size, mat_A, mat_B, mat_result);
        }
    }

    auto time_taken = omp_get_wtime() - start_time;

    normal_matrix_mul(size, mat_A, mat_B, mat_test);

    int err = count_error(mat_result, mat_test, size);
    print_results(err, size, k1, omp_get_max_threads(), time_taken);

    deallocate_memory(mat_A);
    deallocate_memory(mat_B);
    deallocate_memory(mat_result);
    deallocate_memory(mat_test);

    return 0;
}