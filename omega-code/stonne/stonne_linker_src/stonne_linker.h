#include "../include/Config.h"
#include <iostream>

#ifndef __stonne_linker__
#define __stonne_linker__
void simulateDenseConvForward(std::string layer_name, float* input, float* weight, float* output, int R, int S, int C, int K, int G, int N, int X, int Y, int X_, int Y_, int strides, int pad_x, int pad_y, std::string path_to_tile, Config stonne_cfg);

//This function performs the prunning on its own and gets the bitmaps and sparse representation according to that prunning configuration. The prunning is done by prunning the sparsity_level% lowest amount of data in the STA matrix. 
void simulateSparseGemmForward(std::string layer_name, float* KN_matrix_raw, float* MK_matrix_raw, float* output_raw, int N, int G, int gemm_M, int gemm_K, int gemm_N, float sparsity_level, Config stonne_cfg, Dataflow dataflow);

//This function already gets the bitmaps and the matrices in a sparse representaion. 
void simulateSparseGemmWithBitmapsForward(std::string layer_name, float* KN_matrix_raw, float* MK_matrix_raw, float* output_raw, int N, int G, int gemm_M, int gemm_K, int gemm_N, unsigned int* MK_bitmap, unsigned int* KN_bitmap, Config stonne_cfg, Dataflow dataflow);

void simulateDenseGemmForward(std::string layer_name, float* KN_matrix_raw, float* MK_matrix_raw, float* output_raw, int N, int G, int gemm_M, int gemm_K, int gemm_N, std::string path_to_tile, Config stonne_cfg);

//Sparse Dense GEMM with CSR as encoding technique
void simulateSparseDenseGemm(std::string layer_name, float* MK_sparse_matrix, float* KN_dense_matrix, float* output_raw, int M, int K, int N, unsigned int* MK_col_id, unsigned int* MK_row_pointer, int T_N, int T_K, Config stonne_cfg);
#endif
