#include <iostream>
#include "STONNEModel.h"
#include "types.h"
#include <chrono>
#include <assert.h>
#include "testbench.h"
#include <string>
#include <math.h>
#include <utility.h>

using namespace std;


int main(int argc, char *argv[]) {
    float EPSILON=0.05;
    unsigned int MAX_RANDOM=10; //Variable used to generate the random values
    /** Generating the inputs and outputs **/

    //Layer parameters (See MAERI paper to find out the taxonomy meaning)
    std::string layer_name="TestLayer";
    unsigned int K=8;
    unsigned int M=4; 
    unsigned int N=4;
    Config stonne_cfg; //Hardware parameters
    stonne_cfg.m_MSNetworkCfg.ms_size=8;
    stonne_cfg.m_SDMemoryCfg.n_read_ports=8;
    stonne_cfg.m_SDMemoryCfg.n_write_ports=8;

    //Calculating output parameters
    unsigned int O_rows=M;
    unsigned int O_columns = N;

    //Creating arrays to store the ifmap ofmap and weights

    unsigned int MK_size=M*K;
    unsigned int NK_size=N*K;
    unsigned int output_size=M*N;

    float MK_matrix[] = {1.0, 3.0, 5.0, 2.0, 3.0, 1.0, 5.0, 3.0, 2.0, 1.0, 2.0, 5.0, 1.0, 2.0, 3.0, 4.0, 2.0, 3.0, 5.0};
    float KN_matrix[] = {1.0 ,3.0, 1.0, 2.0, 3.0, 1.0, 2.0, 1.0, 2.0, 1.0, 2.0, 3.0, 1.0, 2.0, 4.0, 1.0, 4.0};
    unsigned int metadata_MK[] = {1,1,0,0,1,0,0,0,
	                          0,0,1,1,1,1,1,1,
				  1,1,1,1,1,1,1,0,
				  0,0,0,1,1,1,0,0}; //Actually these are bits

    unsigned int metadata_KN[]={0,1,0,1,
	                        1,1,0,1,
				1,0,0,1,
				1,0,1,0,
                                1,1,1,1,
                                1,0,0,0,
                                0,0,1,1,
                                1,0,0,0};
    float output_matrix[output_size]; //Used to store the CPU computed values to compare with the simulator version
    unsigned int metadata_output[output_size];

    //computing CPU version
    //sequential_layer(R, S, C, K, G, N, X, Y, strides, ifmap, filter, ofmap_cpu); 

    //Computing the CNN Layer with the simulator
    Stonne* stonne_instance = new Stonne(stonne_cfg); //Creating instance of the simulator
    stonne_instance->loadGEMM("GEMM_test", N, K, M, MK_matrix, KN_matrix, metadata_MK, metadata_KN, output_matrix, metadata_output, MK_STA_KN_STR ); //Loading GEMM
    stonne_instance->run(); //Running the simulator 
    //Printing the results
    std::cout << "MK bitmap:" << std::endl;
    for(int i=0; i<M; i++) {
        for(int j=0; j<K; j++) {
		std::cout << metadata_MK[i*K+j] << ", ";
        }

	std::cout << std::endl;
    }

    std::cout << "MK matrix: " << std::endl;
    for(int i=0; i<19; i++) {
        std::cout << MK_matrix[i] << ", ";
    }

    std::cout << std::endl << std::endl;

        std::cout << "KN bitmap:" << std::endl;
    for(int i=0; i<K; i++) {
        for(int j=0; j<N; j++) {
                std::cout << metadata_MK[i*N+j] << ", ";
        }

	std::cout << std::endl;
    }

    std::cout << "KN matrix: " << std::endl;
    for(int i=0; i<17; i++) {
        std::cout << KN_matrix[i] << ", ";
    }

    std::cout << std::endl << std::endl;
    std::cout << "Resulting matrix:" << std::endl;
    for(int i=0; i<M; i++) {
        for(int j=0; j<N; j++) {
            std::cout << output_matrix[i*N+j] << ", ";
        }
	std::cout << std::endl;
    }

}


