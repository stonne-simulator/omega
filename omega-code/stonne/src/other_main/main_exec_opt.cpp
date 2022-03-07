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

void configConvParameters(int argc, char *argv[], Config &stonne_cfg, std::string &layer_name, unsigned int &R, unsigned int &S, unsigned int &C, unsigned int &K, unsigned int &G, unsigned int &N, unsigned int &X, unsigned int &Y, unsigned int &strides,
                      unsigned int &T_R, unsigned int &T_S, unsigned int &T_C, unsigned int &T_K, unsigned int &T_G, unsigned int &T_N, unsigned int &T_X_, unsigned int &T_Y_);

void configSparseGEMMParameters(int argc, char *argv[], Config &stonne_cfg, std::string &layer_name, unsigned int &M, unsigned int &N, unsigned int &K, unsigned int &MK_sparsity, unsigned int &KN_sparsity, Dataflow &dataflow, unsigned int &optimize);

bool runConvCommand(int argc, char *argv[]);
bool runSparseGEMMCommand(int argc, char *argv[]);
bool runHelpCommand();
float* generateMatrixDense(unsigned int rows, unsigned int cols, unsigned int sparsity);

void generateSparseDense(unsigned int rows, unsigned int cols, unsigned int sparsity);

unsigned int* generateBitMapFromDense(float* denseMatrix, unsigned int rows, unsigned int cols, GENERATION_TYPE gen_type);

float* generateMatrixSparseFromDense(float* denseMatrix, unsigned int* bitmap, unsigned int rows, unsigned int cols, GENERATION_TYPE gen_type);

void printDenseMatrix(float* matrix, unsigned int rows, unsigned int cols);
void printBitMap(unsigned int* bitmap, unsigned int rows, unsigned int cols);
void printSparseMatrix(float* sparseMatrix, unsigned int* bitmap, unsigned int rows, unsigned int cols);

int main(int argc, char *argv[]) {
    if(argc > 1) { //IF there is at least one parameter, -h is checked
        string arg = argv[1];
        if(arg=="-h") {
            runHelpCommand();
        }
      
        else if(arg=="-CONV") {
            runConvCommand(argc, argv);
        }

	else if(arg=="-SparseGEMM") {
            runSparseGEMMCommand(argc, argv);
	}

	else {
	    std::cout << "How to use STONNE User Interface: ./" << argv[0] << " -h" << std::endl;
	}
    }

    else {
        std::cout << "How to use STONNE User Interface: ./" << argv[0] << " -h" << std::endl;
    }
}

bool runConvCommand(int argc, char *argv[]) {
    float EPSILON=0.05;
    unsigned int MAX_RANDOM=10; //Variable used to generate the random values
    /** Generating the inputs and outputs **/

    //Layer parameters (See MAERI paper to find out the taxonomy meaning)
    std::string layer_name="TestLayer";
    unsigned int R=1;                                  // R
    unsigned int S=3;                                  // S
    unsigned int C=1;                                  // C
    unsigned int K=1;                                  // K
    unsigned int G=1;                                  // G
    unsigned int N=1;                                  // N
    unsigned int X=1;                                  // X //TODO CHECK X=1 and Y=1
    unsigned int Y=3;                                  // Y
    unsigned int strides=1;                            // Strides
 
    //Tile parameters (See MAERI paper to find out the taxonomy meaning)
    unsigned int T_R=1;                                // T_R
    unsigned int T_S=3;                                // T_S
    unsigned int T_C=1;                                // T_C
    unsigned int T_K=1;                                // T_K
    unsigned int T_G=1;                                // T_G
    unsigned int T_N=1;                                // T_N
    unsigned int T_X_=1;                               // T_X
    unsigned int T_Y_=1;                               // T_Y   
    Config stonne_cfg; //Hardware parameters
//    stonne_cfg.m_MSNetworkCfg.ms_size=128;
    configConvParameters(argc, argv, stonne_cfg, layer_name, R, S, C, K, G, N, X, Y, strides, T_R, T_S, T_C, T_K, T_G, T_N, T_X_, T_Y_); //Modify stonne_cfg and the variables according to user arguments

    //Calculating output parameters
    unsigned int X_= (X - R + strides) / strides;      // X_
    unsigned int Y_=(Y - S + strides) / strides;       // Y_
    std::cout << "Output Size: (X'=" << X_ << ", Y'=" << Y_ << std::endl; 


    //Creating arrays to store the ifmap ofmap and weights
    unsigned int ifmap_size=N*X*Y*C;
    unsigned int filter_size=R*S*(C/G)*K;
    unsigned int ofmap_size=N*X_*Y_*K;
    float* ifmap = new float[ifmap_size];
    float* filter = new float[filter_size];
    float* ofmap = new float[ofmap_size];
    float* ofmap_cpu = new float[ofmap_size]; //Used to store the CPU computed values to compare with the simulator version

    //Filling the arrays with random values
    for(int i=0; i<ifmap_size; i++) {
        ifmap[i]=rand()%MAX_RANDOM;
    }

    for(int i=0;i<filter_size; i++) {
        filter[i]=rand()%MAX_RANDOM;
    }

    //computing CPU version
    sequential_layer(R, S, C, K, G, N, X, Y, strides, ifmap, filter, ofmap_cpu); 

    /** END of generating the inputs and outputs **/
    //
    //
    //
    /** Configuring and running the accelerator  **/
    
    //Computing the CNN Layer with the simulator
    Stonne* stonne_instance = new Stonne(stonne_cfg); //Creating instance of the simulator
    stonne_instance->loadDNNLayer(CONV, layer_name, R, S, C, K, G, N, X, Y, strides, ifmap, filter, ofmap, CNN_DATAFLOW); //Loading the layer
    stonne_instance->loadTile(T_R, T_S, T_C, T_K, T_G, T_N, T_X_, T_Y_); //Loading the tile
    stonne_instance->run(); //Running the simulator 

    /** END of configuring and running the accelerator  **/
    //
    //
    //
    /** CHECKING the results to make sure that the output is correct  **/

    //Comparing the results
    for(int i=0;i<ofmap_size; i++) {
        float difference=fabs(ofmap[i]-ofmap_cpu[i]);
        if(difference > EPSILON) {
            std::cout << "ERROR position " << i <<  ": Value ofmap simulator: " << ofmap[i] << ". Value ofmap CPU: " << ofmap_cpu[i] << std::endl;
            std::cout << "\033[1;31mT test not passed\033[0m" << std::endl;
            delete[] ifmap;
            delete[] filter;
            delete[] ofmap;
            delete[] ofmap_cpu;
            delete stonne_instance;
            assert(false); //Always false
            
        }
    }


    //If the code does not stop then the TEST is correct
    std::cout << "\033[1;32mTest passed correctly \033[0m" << std::endl << std::endl;

    delete[] ifmap;
    delete[] filter;
    delete[] ofmap;
    delete[] ofmap_cpu;
    delete stonne_instance; 
    return true;
}

bool runSparseGEMMCommand(int argc, char *argv[]) {
    //Layer parameters (See SIGMA paper to find out the taxonomy meaning)
    float EPSILON=0.05;
    std::string layer_name="GEMMTestLayer";
    unsigned int M=4;                                  // M
    unsigned int N=4;                                  // N
    unsigned int K=8;                                  // K
    unsigned int optimize = 0; //False

    unsigned int MK_sparsity=20;
    unsigned int KN_sparsity=20;
    Dataflow dataflow = MK_STA_KN_STR;

    Config stonne_cfg;
    stonne_cfg.m_SDMemoryCfg.mem_controller_type=SIGMA_SPARSE_GEMM;

    configSparseGEMMParameters(argc, argv, stonne_cfg, layer_name, M, N, K, MK_sparsity, KN_sparsity, dataflow, optimize);

    //Creating MK matrix
    float* MK_dense_matrix_no_organized = generateMatrixDense(M, K, MK_sparsity);
    float* MK_dense_matrix = new float[M*K];
    
    //KN matrix
    float* KN_dense_matrix_no_organized = generateMatrixDense(K, N, KN_sparsity);
    float* KN_dense_matrix = new float[K*N];

    for(int i=0; i<M*K; i++) {
        MK_dense_matrix[i]=MK_dense_matrix_no_organized[i];
    }

    for(int i=0; i<K*N; i++) {
        KN_dense_matrix[i]=KN_dense_matrix_no_organized[i];

    }


    //See if it is necessary to reorganize
    unsigned int* order_table;
    if(optimize & (dataflow == MK_STA_KN_STR)) {


        order_table = calculateOrdering (MK_dense_matrix_no_organized, M,  K, GEN_BY_ROWS);
        organizeMatrix (MK_dense_matrix, M, K, order_table, GEN_BY_ROWS);

    }

    else if(optimize & (dataflow==MK_STR_KN_STA)) {
        order_table = calculateOrdering (KN_dense_matrix_no_organized, K,  N, GEN_BY_COLS);
        organizeMatrix (KN_dense_matrix, K, N, order_table, GEN_BY_COLS);


    }

    //Creating outputs
    float* cpu_output = new float[M*N];
    float* acc_output = new float[M*N];
    unsigned int* acc_bitmap = new unsigned int[M*N]; //Currently is not generated by the accelerator


    //Generating bitmaps
    unsigned int* MK_bitmap = generateBitMapFromDense(MK_dense_matrix, M, K, GEN_BY_ROWS);
    unsigned int* KN_bitmap = generateBitMapFromDense(KN_dense_matrix, K, N, GEN_BY_COLS);

    //Generating sparse matrix
    float* MK_sparse_matrix = generateMatrixSparseFromDense(MK_dense_matrix, MK_bitmap, M, K, GEN_BY_ROWS);
    float* KN_sparse_matrix = generateMatrixSparseFromDense(KN_dense_matrix, KN_bitmap, K, N, GEN_BY_COLS);
/*
    std::cout << "MK_matrix_no_organized: " << std::endl;
    printDenseMatrix(MK_dense_matrix_no_organized, M, K);
    std::cout << "MK_matrix: " << std::endl;
    printDenseMatrix(MK_dense_matrix, M, K);
     std::cout << "order_table" << std::endl;
    for(int i=0; i<M; i++) {
        std::cout << order_table[i] << std::endl;
    }
*/
    /*
    std::cout << "KN Matrix: " << std::endl;
    printDenseMatrix(KN_dense_matrix, K, N);
    std::cout << "MK_bitmap: " << std::endl;
    printBitMap(MK_bitmap, M, K);
    std::cout << "KN_bitmap: " << std::endl;
    printBitMap(KN_bitmap, K,N);
    std::cout << "MK_sparse_matrix: " << std::endl;
    printSparseMatrix(MK_sparse_matrix, MK_bitmap, M, K);
    std::cout << "KN_sparse_matrix: " << std::endl;
    printSparseMatrix(KN_sparse_matrix, KN_bitmap, K, N);
    */
    //Running STONNE
    Stonne* stonne_instance = new Stonne(stonne_cfg); //Creating instance of the simulator
    stonne_instance->loadGEMM(layer_name, N, K, M, MK_sparse_matrix, KN_sparse_matrix, MK_bitmap, KN_bitmap, acc_output, acc_bitmap, dataflow ); //Loading GEMM
    stonne_instance->run(); //Running the simulator
    if(optimize && (dataflow==MK_STA_KN_STR)) {
        organizeMatrixBack (acc_output, M,  N, order_table, GEN_BY_ROWS);
    }

    else if(optimize && (dataflow==MK_STR_KN_STA)) {
        organizeMatrixBack(acc_output, M, N, order_table, GEN_BY_COLS);
    }

    /** CHECKING the results to make sure that the output is correct  **/
    std::cout << "Running CPU version to compare results" << std::endl;
    //Generating cpu output
    cpu_gemm(MK_dense_matrix_no_organized, KN_dense_matrix_no_organized, cpu_output, M, N, K);
/*
    std::cout << "Output matrix generated by CPU: " << std::endl;
    printDenseMatrix(cpu_output, M, N);
    std::cout << "Output matrix generated by STONNE: " << std::endl;
    printDenseMatrix(acc_output, M, N);
*/

    //Comparing the results
    for(int i=0;i<M; i++) {
        for(int j=0; j<N; j++) {
            float difference=fabs(cpu_output[i*N+j]-acc_output[i*N+j]);
            if(difference > EPSILON) {
                std::cout << "ERROR position (" << i << "," << j <<  "): Value out simulator: " << acc_output[i*N+j] << ". Value out CPU: " << cpu_output[i*N+j] << std::endl;
                std::cout << "\033[1;31mT test not passed\033[0m" << std::endl;
		delete[] MK_dense_matrix;
                delete[] KN_dense_matrix;
		delete[] MK_dense_matrix_no_organized;
		delete[] KN_dense_matrix_no_organized;
                delete[] MK_bitmap;
                delete[] KN_bitmap;
                delete[] MK_sparse_matrix;
                delete[] KN_sparse_matrix;
                delete[] cpu_output;
                delete[] acc_output;
                delete[] acc_bitmap;
		if(optimize)
		    delete[] order_table;
                delete stonne_instance;

                assert(false); //Always false

            }

        }
    }


    //If the code does not stop then the TEST is correct
    std::cout << "\033[1;32mTest passed correctly \033[0m" << std::endl << std::endl;


    delete[] MK_dense_matrix_no_organized;
    delete[] KN_dense_matrix_no_organized;
    delete[] MK_dense_matrix;
    delete[] KN_dense_matrix;
    delete[] MK_bitmap;
    delete[] KN_bitmap;
    delete[] MK_sparse_matrix;
    delete[] KN_sparse_matrix;
    delete[] cpu_output;
    delete[] acc_output;
    delete[] acc_bitmap;
    if(optimize)
        delete[] order_table;
    
    delete stonne_instance;
    return true;    



}

bool runHelpCommand() {
    std::cout << "Help: " << std::endl;
    std::cout << "EXECUTION OF DENSE CNNs AND DENSE GEMMS" << std::endl;
    std::cout << "********************************************************************************" << std::endl;
    std::cout << "********************************************************************************" << std::endl;
    std::cout << "********************************************************************************" << std::endl << std::endl;
    std::cout << "./stonne [-h, -CONV, -SparseGEMM]" << std::endl;
    std::cout << "Hardware parameters" << std::endl;
    std::cout << "---------------------------------------------------------------------------------" << std::endl << std::endl;
    std::cout << "-num_ms=[x]: Number of multiplier switches (must be power of 2)" << std::endl;
    std::cout << "-dn_bw=[x]: Number of read ports in the SDMemory (must be power of 2)" << std::endl;
    std::cout << "-rn_bw=[x]: Number of write ports in the SDMemory (must be power of 2)" << std::endl;
    std::cout << "-rn_type=[0=ASNETWORK, 1=FENETWORK]: type of the ReduceNetwork to be used" << std::endl;
    std::cout << "-print_stats=[0,1]: Flag that enables the printing of the statistics" << std::endl;
    std::cout << "Layer configuration parameters" << std::endl;
    std::cout << "---------------------------------------------------------------------------------" << std::endl << std::endl;
    std::cout << "-layer_name: Name of the layer to run. The output statistic file will be called by this name" << std::endl; 
    std::cout << "-R: Number of flter rows" << std::endl;
    std::cout << "-S: Number of filter columns" << std::endl;
    std::cout << "-C: Number of filter and input channels" << std::endl;
    std::cout << "-K: Number of filters and output channels" << std::endl;
    std::cout << "-G: Number of groups" << std::endl;
    std::cout << "-N: Number of inputs (Only 1 is supported so far)" << std::endl;
    std::cout << "-X: Number of input rows" << std::endl;
    std::cout << "-Y: Number of input columns" << std::endl;
    std::cout << "-strides: Stride value used in the layer" << std::endl;
    std::cout << "Tile configuration parameters" << std::endl;
    std::cout << "---------------------------------------------------------------------------------" << std::endl << std::endl;
    std::cout << "-T_R: Number of flter rows mapped at a time" << std::endl;
    std::cout << "-T_S: Number of filter columns mapped at a time" << std::endl;
    std::cout << "-T_C: Number of filter and input channels per group mapped at a time" << std::endl;
    std::cout << "-T_K: Number of filters and output channels per group mapped at a time" << std::endl;
    std::cout << "-T_G: Number of groups mappd at a time" << std::endl;
    std::cout << "-T_N: Number of inputs mapped at a time (Only 1 is supported so far)" << std::endl;
    std::cout << "-T_X_: Number of input rows mapped at a time" << std::endl;
    std::cout << "-T_Y_: Number of input columns mapped a time" << std::endl;
    std::cout << "** Please take into consideration that: **" << std::endl;
    std::cout << "Virtual Neuron Size (VN_Size) will be T_R*T_S*T_C" << std::endl;
    std::cout << "Number of Virtual Neurons mapped (Num_VNs) will be T_K*T_G*T_N*T_X_*T_T_" << std::endl;
    std::cout << "The minimum number of MSwitches needed will be at least VN_Size*Num_VNs" << std::endl;
    std::cout << "However, if folding (iteration over the same VN) is enabled, 1 extra MSwitch per VN will be needed to manage the psum" << std::endl;
    std::cout << "In this case, the minimum number of MSwitches needed will be at least (VN_Size+1)*Num_VNs" << std::endl;
    std::cout << "Folding (iteration over the same virtual neuron) will be enabled if (R/T_S)*(S/T_S)*(C/T_C) > 1" << std::endl;        

    std::cout << std::endl << std::endl << std::endl;
    std::cout << "EXECUTION OF SPARSE GEMMs" << std::endl;
    std::cout << "********************************************************************************" << std::endl;
    std::cout << "********************************************************************************" << std::endl;
    std::cout << "********************************************************************************" << std::endl << std::endl;
    std::cout << "Hardware parameters" << std::endl;
    std::cout << "---------------------------------------------------------------------------------" << std::endl << std::endl;
    std::cout << "-num_ms=[x]: Number of multiplier switches (must be power of 2)" << std::endl;
    std::cout << "-dn_bw=[x]: Number of read ports in the SDMemory (must be power of 2)" << std::endl;
    std::cout << "-rn_bw=[x]: Number of write ports in the SDMemory (must be power of 2)" << std::endl;
    std::cout << "-print_stats=[0,1]: Flag that enables the printing of the statistics" << std::endl;
    std::cout << "Layer configuration parameters" << std::endl;
    std::cout << "---------------------------------------------------------------------------------" << std::endl << std::endl;
    std::cout << "-layer_name: Name of the layer to run. The output statistic file will be called by this name" << std::endl;
    std::cout << "-M=Number of rows MK matrix" << std::endl;
    std::cout << "-N=Number of columns KN matrix" << std::endl;
    std::cout << "-K=Number of columns MK and rows KN matrix (cluster size)" << std::endl;
    std::cout << "-MK_sparsity=Percentage of sparsity MK matrix (0-100)" << std::endl;
    std::cout << "-KN_sparsity=Percentahe of sparsity KN matrix (0-100)" << std::endl;
    std::cout << "-dataflow=MK_STA_KN_STR or MK_STR_KN_STA" << std::endl;
    std::cout << "-optimize=[0,1]: apply compiler-based optimizations" << std::endl;

    exit(0);
    return true; //Never executed
}

//This function modifies the default values of the parameters according to user arguments.
void configConvParameters(int argc, char *argv[], Config &stonne_cfg, std::string &layer_name, unsigned int &R, unsigned int &S, unsigned int &C, unsigned int &K, unsigned int &G, unsigned int &N, unsigned int &X, unsigned int &Y, unsigned int &strides,
                      unsigned int &T_R, unsigned int &T_S, unsigned int &T_C, unsigned int &T_K, unsigned int &T_G, unsigned int &T_N, unsigned int &T_X_, unsigned int &T_Y_) {

    //Parsing
    for(int i=2; i<argc; i++) { //0 is the name of the program and 1 is the execution command type
        string arg = argv[i];
        //Spliting using = character
        string::size_type pos = arg.find('=');
        if(arg.npos != pos) {
            string value_str=arg.substr(pos+1);
            string name=arg.substr(0, pos);
            unsigned int value;
            if((name != "-layer_name") && (name != "-rn_type")) { //string parameters
                value=stoi(value_str);
            }
            //Checking parameter name
            if(name=="-num_ms") {
                if(!ispowerof2(value)) {   //Checking that the num_ms is power of 2
                    std::cout << "Error: -num_ms must be power of 2" << std::endl;
                    exit(1);
                }
                std::cout << "Changing num_ms to " << value << std::endl; //To debug
                stonne_cfg.m_MSNetworkCfg.ms_size=value;
            }

            else if(name=="-dn_bw") {
                if(!ispowerof2(value)) {
                    std::cout << "Error: -dn_bw must be power of 2" << std::endl;
                    exit(1);
                }
                std::cout << "Changing dn_bw to " << value << std::endl; //To debug
                stonne_cfg.m_SDMemoryCfg.n_read_ports=value;
            }

            else if(name=="-rn_bw") {
                if(!ispowerof2(value)) {
                    std::cout << "Error: -rn_bw must be power of 2" << std::endl; 
                    exit(1);
                }
                std::cout << "Changing rn_bw to " << value << std::endl;
                stonne_cfg.m_SDMemoryCfg.n_write_ports=value;
            }



            else if(name=="-print_stats") {
                if((value != 0) && (value != 1)) {
                    std::cout << "Error: -print_stats only supports 0 or 1" << std::endl;
                    exit(1);
                }
                std::cout << "Changing print_stats to " << value << std::endl;
                stonne_cfg.print_stats_enabled=value; 
            }

            else if(name=="-rn_type") {
                std::cout << "Changing rn_type to " << value_str << std::endl;
                stonne_cfg.m_ASNetworkCfg.reduce_network_type=get_type_reduce_network_type(value_str);
            }
            //Running configuration parameters (layer and tile)
   
           //Layer parameters
           else if(name=="-layer_name") {
               std::cout << "Changing layer_name to " << value_str << std::endl;
               layer_name=value_str; 
           }
        
           else if(name=="-R") {
                std::cout << "Changing R to " << value << std::endl;
                R=value;
           }

           else if(name=="-S") {
                std::cout << "Changing S to " << value << std::endl;
                S=value;
           }

           else if(name=="-C") {
                std::cout << "Changing C to " << value << std::endl;
                C=value;
           } 
        
           else if(name=="-K") {
                std::cout << "Changing K to " << value << std::endl;
                K=value;
           }
  
           else if(name=="-G") {
               std::cout << "Changing G to " << value << std::endl;
               G=value;
           }
  
           else if(name=="-N") {
                std::cout << "Changing N to " << value << std::endl;
                N=value;
           }

           else if(name=="-X") {
                std::cout << "Changing X to " << value << std::endl;
                X=value;
           }

           else if(name=="-Y") {
                std::cout << "Changing Y to " << value << std::endl;
                Y=value;
           }

           else if(name=="-strides") {
               std::cout << "Changing strides to " << value << std::endl;
               strides=value;
           }

           //Tile parameters
           else if(name=="-T_R") {
                std::cout << "Changing T_R to " << value << std::endl;
                T_R=value;
           } 

           else if(name=="-T_S") {
                std::cout << "Changing T_S to " << value << std::endl;
                T_S=value;
           }

           else if(name=="-T_C") {
                std::cout << "Changing T_C to " << value << std::endl;
                T_C=value;
           }

           else if(name=="-T_K") {
                std::cout << "Changing T_K to " << value << std::endl;
                T_K=value;
           }

           else if(name=="-T_G") {
               std::cout << "Changing T_G to " << value << std::endl;
               T_G=value;
           }

           else if(name=="-T_N") {
                std::cout << "Changing T_N to " << value << std::endl;
                T_N=value;
           }

           else if(name=="-T_X_") {
                std::cout << "Changing T_X_ to " << value << std::endl;
                T_X_=value;
           }

           else if(name=="-T_Y_") {
                std::cout << "Changing T_Y_ to " << value << std::endl;
                T_Y_=value;
           }

           //Parameter is not recognized
           else {
                std::cout << "Error: parameter " << name << " does not exist" << std::endl;
                exit(1);
            }

 
    
           

        }
        else {

            std::cout << "Error: parameter " << arg << " does not exist" << std::endl;
            exit(1);

        }
    }
}

void configSparseGEMMParameters(int argc, char *argv[], Config &stonne_cfg, std::string &layer_name, unsigned int &M, unsigned int &N, unsigned int &K, unsigned int &MK_sparsity, unsigned int &KN_sparsity, Dataflow &dataflow, unsigned int &optimize) {
    //Parsing
    for(int i=2; i<argc; i++) { //0 is the name of the program and 1 is the execution command type
        string arg = argv[i];
        //Spliting using = character
        string::size_type pos = arg.find('=');
        if(arg.npos != pos) {
            string value_str=arg.substr(pos+1);
            string name=arg.substr(0, pos);
            unsigned int value;
            if((name != "-layer_name") && (name != "-dataflow")) { //string parameters
                value=stoi(value_str);
            }
            //Checking parameter name
            if(name=="-num_ms") {
                if(!ispowerof2(value)) {   //Checking that the num_ms is power of 2
                    std::cout << "Error: -num_ms must be power of 2" << std::endl;
                    exit(1);
                }
                std::cout << "Changing num_ms to " << value << std::endl; //To debug
                stonne_cfg.m_MSNetworkCfg.ms_size=value;
            }

            else if(name=="-dn_bw") {
                if(!ispowerof2(value)) {
                    std::cout << "Error: -dn_bw must be power of 2" << std::endl;
                    exit(1);
                }
                std::cout << "Changing dn_bw to " << value << std::endl; //To debug
                stonne_cfg.m_SDMemoryCfg.n_read_ports=value;
            }

            else if(name=="-rn_bw") {
                if(!ispowerof2(value)) {
                    std::cout << "Error: -rn_bw must be power of 2" << std::endl; 
                    exit(1);
                }
                std::cout << "Changing rn_bw to " << value << std::endl;
                stonne_cfg.m_SDMemoryCfg.n_write_ports=value;
            }



            else if(name=="-print_stats") {
                if((value != 0) && (value != 1)) {
                    std::cout << "Error: -print_stats only supports 0 or 1" << std::endl;
                    exit(1);
                }
                std::cout << "Changing print_stats to " << value << std::endl;
                stonne_cfg.print_stats_enabled=value; 
            }

            //Running configuration parameters (layer)
   
           //Layer parameters
           else if(name=="-layer_name") {
               std::cout << "Changing layer_name to " << value_str << std::endl;
               layer_name=value_str; 
           }
        
           else if(name=="-M") {
                std::cout << "Changing M to " << value << std::endl;
                M=value;
           }

           else if(name=="-N") {
                std::cout << "Changing N to " << value << std::endl;
                N=value;
           }

           else if(name=="-K") {
                std::cout << "Changing K to " << value << std::endl;
                K=value;
           } 
        
           else if(name=="-MK_sparsity") {
                std::cout << "Changing MK_sparsity to " << value << std::endl;
                MK_sparsity=value;
           }
  
           else if(name=="-KN_sparsity") {
               std::cout << "Changing KN_sparsity to " << value << std::endl;
               KN_sparsity=value;
           }
  
	   else if(name=="-dataflow") {
                std::cout << "Changing dataflow to " << value_str << std::endl;
                dataflow=get_type_dataflow_type(value_str);
            }

	    else if(name=="-optimize") {
               std::cout << "Changing optimize " << value << std::endl;
               optimize=value;
           }



           //Parameter is not recognized
           else {
                std::cout << "Error: parameter " << name << " does not exist" << std::endl;
                exit(1);
            }

 
    
           

        }
        else {

            std::cout << "Error: parameter " << arg << " does not exist" << std::endl;
            exit(1);

        }
    }
}


