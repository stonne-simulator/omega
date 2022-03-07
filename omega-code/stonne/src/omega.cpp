#include <iostream>
#include "STONNEModel.h"
#include "types.h"
#include <chrono>
#include <assert.h>
#include "testbench.h"
#include <string>
#include <math.h>
#include <utility.h>
#include<fstream>
#include<sstream>

void configParameters(int argc, char *argv[], int &V, int &F, int &G, int &E, int &T_Va, int &T_N, int &T_Fa, int &T_Vc, int &T_G, int &T_Fc, int &pe_agg, int &pe_cmb, int &dn_bw_agg, int &dn_bw_cmb, int &rn_bw_agg, int &rn_bw_cmb, std::string &edge_path, std::string &vertex_path)
{
	    //Parsing all the parameters from the command line
    for(int i=1; i<argc; i++) { //0 is the name of the program
        std::string arg = argv[i];
        //Spliting using = character
        std::string::size_type pos = arg.find('=');
        if(arg.npos != pos) {
            std::string value_str=arg.substr(pos+1);
            std::string name=arg.substr(0, pos);
            unsigned int value;
            if((name!="-vertex_path") && (name!="-edge_path")) { //string parameters
                value=stoi(value_str);
            }
            //Checking parameter name
            if(name=="-pe_agg") {
                if(!ispowerof2(value)) {   //Checking that the num_ms is power of 2
                    std::cout << "Error: -num_ms must be power of 2" << std::endl;
                    exit(1);
                }
                std::cout << "Changing pe_agg to " << value << std::endl; //To debug
                pe_agg=value;
            }

            else if(name=="-pe_cmb") {
                if(!ispowerof2(value)) {   //Checking that the num_ms is power of 2
                    std::cout << "Error: -num_ms must be power of 2" << std::endl;
                    exit(1);
                }
                std::cout << "Changing pe_cmb to " << value << std::endl; //To debug
                pe_cmb=value;
            }
 			else if(name=="-dn_bw_agg") {
                if(!ispowerof2(value)) {   //Checking that the num_ms is power of 2
                    std::cout << "Error: -dn_bw must be power of 2" << std::endl;
                    exit(1);
                }
                std::cout << "Changing dn_bw_agg to " << value << std::endl; //To debug
                dn_bw_agg=value;
            }

            else if(name=="-dn_bw_cmb") {
                if(!ispowerof2(value)) {   //Checking that the num_ms is power of 2
                    std::cout << "Error: -dn_bw must be power of 2" << std::endl;
                    exit(1);
                }
                std::cout << "Changing dn_bw_cmb to " << value << std::endl; //To debug
                dn_bw_cmb=value;
            }

            else if(name=="-rn_bw_agg") {
                if(!ispowerof2(value)) {   //Checking that the num_ms is power of 2
                    std::cout << "Error: -rn_bw must be power of 2" << std::endl;
                    exit(1);
                }
                std::cout << "Changing rn_bw_agg to " << value << std::endl; //To debug
                rn_bw_agg=value;
            }

            else if(name=="-rn_bw_cmb") {
                if(!ispowerof2(value)) {   //Checking that the num_ms is power of 2
                    std::cout << "Error: -rn_bw must be power of 2" << std::endl;
                    exit(1);
                }
                std::cout << "Changing rn_bw_cmb to " << value << std::endl; //To debug
                rn_bw_cmb=value;
            }


            //Running configuration parameters (layer)
   
           //Layer parameters
       
           else if(name=="-V") {
                std::cout << "Changing V to " << value << std::endl;
                V=value;
           }

           else if(name=="-F") {
                std::cout << "Changing F to " << value << std::endl;
                F=value;
           }

           else if(name=="-G") {
                std::cout << "Changing G to " << value << std::endl;
                G=value;
           } 
        
           else if(name=="-E") {
                std::cout << "Changing E to " << value << std::endl;
                E=value;
           }
  
           else if(name=="-T_Va") {
               std::cout << "Changing T_Va to " << value << std::endl;
               T_Va=value;
           }

	 else if(name=="-T_N") {
               std::cout << "Changing T_N to " << value << std::endl;
               T_N=value;
           }

	   else if(name=="-T_Fa") {
               std::cout << "Changing T_Fa to " << value << std::endl;
               T_Fa=value;
           }

	   else if(name=="-T_Vc") {
               std::cout << "Changing T_Vc to " << value << std::endl;
               T_Vc=value;
           }

	   else if(name=="-T_G") {
               std::cout << "Changing T_G to " << value << std::endl;
               T_G=value;
           }

	else if(name=="-T_Fc") {
               std::cout << "Changing T_Fc to " << value << std::endl;
               T_Fc=value;
           }
  
	   else if(name=="-edge_path") {
                std::cout << "Changing edge_path to " << value_str << std::endl;
                edge_path=value_str;
            }

          else if(name=="-vertex_path") {
                std::cout << "Changing vertex_path to " << value_str << std::endl;
                vertex_path=value_str;
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

//end of parsing code

enum GRAN{ELEMENT, ROW, COLUMN, SEQUENTIAL};
enum KERNELORDER{AC,CA};
int main(int argc, char *argv[])
{

	//EXTERNAL INPUTS - Loop order, Kernel order, dimensions - V,F,G. Tile sizes - T_Va, T_Fa, T_N, T_Vc, T_Fc, T_G, adjacency matrix
	GRAN granularity;
	//KERNELORDER order;
	KERNELORDER order = AC;
	//char loop_order_agg[4], loop_order_cmb[4];
	int V,F,G,E;
	V=3786;F=29;G=2;E=14456;
	int T_Va, T_N, T_Fa, T_Vc, T_G, T_Fc;
	T_Va=17;T_N=1;T_Fa=15;T_Vc=85;T_G=1;T_Fc=3;
	int T_Vp, T_Fp;
	int pe_agg=256, pe_cmb=256, dn_bw_agg=256, dn_bw_cmb=256, rn_bw_agg=256, rn_bw_cmb=256;
	std::string vertex_path = "sample_graphs/vertex_proteins_batch64.txt";
	std::string edge_path = "sample_graphs/edge_proteins_batch64.txt";
	char loop_order_agg[4]="VFN";
	char loop_order_cmb[4]="VGF";
	configParameters(argc, argv, V, F, G, E, T_Va, T_N, T_Fa, T_Vc, T_G, T_Fc, pe_agg, pe_cmb, dn_bw_agg, dn_bw_cmb, rn_bw_agg, rn_bw_cmb,  edge_path, vertex_path);
	if(order==AC)			//TODO: FUTURE WORK: Support modelling all the loop orders. The loop order depend on that of the current cost model which is STONNE. We aim to generalize omega's model to wider range of cost models in the future
	{
		if(loop_order_agg[0]==loop_order_cmb[0] && loop_order_agg[1]==loop_order_cmb[1] && (loop_order_agg[0]=='V' || loop_order_agg[0]=='F'))
			granularity = ELEMENT;
		else if(loop_order_agg[0]==loop_order_cmb[0] && loop_order_agg[0]=='V')
			granularity = ROW;
		else if(loop_order_agg[0]==loop_order_cmb[0] && loop_order_agg[0]=='F')
			granularity = COLUMN;
		else
			granularity = SEQUENTIAL;

	}
	else if(order==CA)	//TODO: Future Work: Support modelling all phase orders and loop orders. Currently parallelism strategy, Inter-phase dataflows and mappings are tunable
	{
		if(((loop_order_agg[0]=='N' && loop_order_cmb[0]=='V') && (loop_order_agg[1]=='F' && loop_order_cmb[1]=='G'))||((loop_order_agg[0]=='F' && loop_order_cmb[0]=='G') && (loop_order_agg[1]=='N' && loop_order_cmb[1]=='V')))
			granularity = ELEMENT;
		else if(loop_order_agg[0]=='N' && loop_order_cmb[0]=='V')
			granularity = ROW;
		else if(loop_order_agg[0]=='F' && loop_order_cmb[0]=='G')
			granularity = COLUMN;
		else
			granularity = SEQUENTIAL;
	}

	
	
	//std::cout<<"\n"<<loop_order_agg<<"\n";
	std::cout<<"\n"<<loop_order_agg<<"\t"<<loop_order_cmb<<"\t"<<granularity<<"\n";
	
	if(granularity==ELEMENT)				//TODO: Future work: When more loop orders are supported in STONNE.
	{
		int n_cycles=0;
		T_Vp = (T_Va>T_Vc) ? T_Va : T_Vc;
		T_Fp = (T_Fa>T_Fc) ? T_Fa : T_Fc;
		int n_buffers = 2*T_Vp*T_Fp;
		//Run Stonne in loop, Update n_cycles
		//Need to do loop order manipulations
		std::cout<<"Number of cycles = "<<n_cycles<<"\nNumber of buffers = "<<n_buffers;
	}
	
	if(granularity==ROW)				// We currently focus on this for PP. Will add support as more loop orders get supported
	{
		int n_cycles=0;
		T_Vp = (T_Va>T_Vc) ? T_Va : T_Vc;
		int n_buffers = 2*T_Vp*F;
		
		float EPSILON=0.05;
    		std::string layer_name="SparseDenseTestLayer";
		
		/*
		//float* MK_dense_matrix_no_organized = generateMatrixDense(V, V, 99.8);
		float* MK_dense_matrix_no_organized =  new float[V*V];
		std::ifstream MK_mat("eval_set/adj_reddit_batch32.txt");
		std::string line,temp;
		

		for(int i=0;i<V;i++){
		std::getline(MK_mat, line);
		std::stringstream ss(line);
			for(int j=0;j<V;j++){
				std::getline(ss,temp,',');
				MK_dense_matrix_no_organized[i*V+j]= (unsigned int) stoi(temp);
				//std::cout<<"\t"<<temp;
				}
		}
		
		    MK_mat.close();
		
        	float* MK_dense_matrix = new float[V*V];
		
		 for(int i=0; i<V*V; i++) {
       		 MK_dense_matrix[i]=MK_dense_matrix_no_organized[i];

    		}
//	
//		

    		    int nnz=0;
*/

		//Reading the meta data of the graph since the structure of the graph influences the performance

		int* MK_col_id = new int[E]; //generateMinorIDFromDense(MK_dense_matrix, V, V, nnz, GEN_BY_ROWS);
		
		std::ifstream MK_col(edge_path);	//Reading the column-id of adjancency matrix from the edge file.
		std::string line,temp;
		

		std::getline(MK_col, line);
		std::stringstream ss(line);
			for(int j=0;j<E;j++){
				std::getline(ss,temp,',');
				MK_col_id[j]= (int) stoi(temp);
				//std::cout<<"\t"<<temp;
				}
		
		    MK_col.close();
		

    		int* MK_row_pointer = new int[V+1];//generateMajorPointerFromDense(MK_dense_matrix, V , V, GEN_BY_ROWS);

		std::ifstream MK_row(vertex_path);	//Reading the row pointer of the adjacency matrix from the vertex file
		std::string line1,temp1;
		
		std::getline(MK_row, line1);
		std::stringstream ss1(line1);
			for(int j=0;j<V+1;j++){
				std::getline(ss1,temp1,',');
				MK_row_pointer[j]= (int) stoi(temp1);
				//std::cout<<"\t"<<temp;
				}
		
		    MK_row.close();

    		//int* MK_col_id = generateMinorIDFromDense(MK_dense_matrix, V, V, nnz, GEN_BY_ROWS);
    		//int* MK_row_pointer = generateMajorPointerFromDense(MK_dense_matrix, V , V, GEN_BY_ROWS);
		
		
//		std::cout<<"\nMK Dense matrix - \n";
//	
//	for(int i=0; i<V; i++) {
//        for(int j=0; j<V; j++) {
//            std::cout << MK_dense_matrix[i*V+j] << "\t";
//	}
//	std::cout << "\n";
//    }

		std::cout << "\n\n-------MK Row pointer-------------------------------------------------------------------------------------------------------------- - \n"; 
		for(int i=0;i<=V;i++)
    		{
    			std::cout<<MK_row_pointer[i]<<"\t";
    		}
    		std::cout << "\n\n";	



		
		std::cout << "\n\n----------MK Col ID-------------------------------------------------------------------------------------------------------------\n"; 
		for(int i=0;i<E;i++)
    		{
    			std::cout<<MK_col_id[i]<<"\t";
    		}
    		std::cout << "\n\n";	
		
    		//Generating sparse matrix
    		//float* MK_sparse_matrix = generateMatrixSparseFromDenseNoBitmap(MK_dense_matrix, V, V, GEN_BY_ROWS);
    
		float* MK_sparse_matrix = new float[E]; //generateMatrixSparseFromDenseNoBitmap(MK_dense_matrix, V, V, GEN_BY_ROWS);
		
		for(int i=0;i<E;i++)
			MK_sparse_matrix[i]=1;

    		float* KN_dense_matrix_no_organized = generateMatrixDense(V, F, 0);		//Generating dense matrices from the dimensions
    		float* KN_dense_matrix = new float[V*F];

		


   		 for(int i=0; i<V*F; i++) {
       		 KN_dense_matrix[i]=KN_dense_matrix_no_organized[i];
    		}
    		
    		unsigned int* clocked_op = new unsigned int [V*F];				//clocked_op is used to collect timestamps. Please search for the term to code modified in STONNE to collect the timestamps.
    
    		float* cpu_output = new float[V*F];
    		float* acc_output = new float[V*F];
    
    		std::cout<<"\n\n\n\n------Simulating SpMM-----------------------------------------------------------------------------------------------------------\n\n\n";
	///////////////////////////////////////////////////////////////////////////// OMEGA flow Step 1: Runnning SpMM //////////////////////////////////////////////////////////////////////////////////////////	
		//Configure STONNE for the SpMM phase
		Config stonne_cfg1, stonne_cfg2;
		
				stonne_cfg1.m_MSNetworkCfg.ms_size=pe_agg;
                stonne_cfg1.m_SDMemoryCfg.n_read_ports=dn_bw_agg;					
                stonne_cfg1.m_SDMemoryCfg.n_write_ports=rn_bw_agg;  				
                stonne_cfg1.m_ASNetworkCfg.accumulation_buffer_enabled=1;			
                stonne_cfg1.m_SDMemoryCfg.mem_controller_type=MAGMA_SPARSE_DENSE;
                stonne_cfg1.print_stats_enabled=0; 
                if(T_N>1)
               		stonne_cfg1.m_ASNetworkCfg.reduce_network_type=get_type_reduce_network_type("ASNETWORK");
               	else
                	stonne_cfg1.m_ASNetworkCfg.reduce_network_type=get_type_reduce_network_type("TEMPORALRN");
                Stonne* stonne_instance = new Stonne(stonne_cfg1);
                

    	stonne_instance->loadSparseDense(layer_name, F, V, V, MK_sparse_matrix, KN_dense_matrix, (unsigned int*)MK_col_id, (unsigned int*) MK_row_pointer, acc_output, T_Fa, T_N); //Loading Sparse Dense
    stonne_instance->loadClocking(clocked_op);
    stonne_instance->run(); //Running the simulator to get the intra-phase stats and the timestamps
    unsigned int a_cycles = stonne_instance->getCycles(); //Cycles returned by the simulator for the SpMM phase
	std::cout<<"\n\n------Timestamps for simulated SpMM-----------------------------------------------------------------------------------------------------------\n";
    for(int i=F-1;i<V*F;i+=F)
    	std::cout<<clocked_op[i]<<"\t";	//Timestamps are saved in this variable when the simulator is running

    std::cout<<"\n\n------Statistics for simulated SpMM-----------------------------------------------------------------------------------------------------------\n";
    ///////////// Other stats
    MSNetworkStats mulstatssp = stonne_instance->getMultiplierNetworkStats();
    std::cout<<"\n\nRF Stats:\n";
    std::cout<<"RF weight reads = "<<mulstatssp.n_l1_weight_reads<<"\n";
    std::cout<<"RF weight writes = "<<mulstatssp.n_l1_weight_writes<<"\n";
    std::cout<<"RF input reads = "<< mulstatssp.n_l1_input_reads<<"\n";
    std::cout<<"RF input writes = "<< mulstatssp.n_l1_input_writes<<"\n";
    std::cout<<"RF psum reads = "<<mulstatssp.n_l1_psum_reads<<"\n";
    std::cout<<"RF psum writes = "<<mulstatssp.n_l1_psum_writes<<"\n";
    //std::cout<<"Number of multiplications = "<<mulstatssp.n_multiplications<<"\n";
    //std::cout<<"Number of link traversals = "<<mulstatssp.n_local_network_traversals<<"\n";
    
    std::cout<<"\n\nAdder stats: \n";
    if(T_N==1)
	std::cout<<"Reduction network is linear since reduction is temporal\n";
    else
    {
    	ASNetworkStats adderstatssp = stonne_instance->getASNetworkStats();
    	std::cout<<"Number of link traversals = "<<adderstatssp.n_total_traversals<<"\n";
    }    

    DSNetworkStats diststatssp = stonne_instance->getDSNetworkStats();
    std::cout<<"\n\nDistribution stats: \n";
    std::cout<<"Number of link traversals = "<<diststatssp.n_total_traversals<<"\n";
    
    SDMemoryStats memstatssp = stonne_instance->getMemoryStats();
    std::cout<<"\n\nGlobal Buffer stats: \n";
    std::cout << "Global SRAM weight reads = " << memstatssp.n_SRAM_weight_reads <<"\n";
    std::cout << "Global SRAM input reads = " << memstatssp.n_SRAM_input_reads <<"\n";
    std::cout << "Global SRAM psum reads = " << memstatssp.n_SRAM_psum_reads << "\n";
    std::cout << "Global SRAM psum writes = " << memstatssp.n_SRAM_psum_writes << "\n";    


    /** CHECKING the results to make sure that the output is correct  **/
	//std::cout<<¨\n\n¨;
 /*   std::cout << "Running CPU version to compare results" << std::endl;
    //Generating cpu output
    cpu_gemm(MK_dense_matrix_no_organized, KN_dense_matrix_no_organized, cpu_output, V, F, V);

    //Comparing the results
    for(int i=0;i<V; i++) {
        for(int j=0; j<F; j++) {
            float difference=fabs(cpu_output[i*F+j]-acc_output[i*F+j]);
            if(difference > EPSILON) {
                std::cout << "ERROR position (" << i << "," << j <<  "): Value out simulator: " << acc_output[i*F+j] << ". Value out CPU: " << cpu_output[i*F+j] << std::endl;
                std::cout << "\033[1;31mT test not passed\033[0m" << std::endl;
		delete[] MK_dense_matrix;
                delete[] KN_dense_matrix;
		delete[] MK_dense_matrix_no_organized;
		delete[] KN_dense_matrix_no_organized;
                delete[] MK_col_id;
                delete[] MK_row_pointer;
                delete[] MK_sparse_matrix;
                delete[] cpu_output;
                delete[] acc_output;
                delete stonne_instance;	
                
                
                
             
	}
	}}
	
	
    //If the code does not stop then the TEST is correct
    std::cout << "\033[1;32mTest passed correctly \033[0m" << std::endl << std::endl;
*/

	//Validation of simulated matmul
       std::cout << "Running CPU version to compare results" << std::endl;
    //Generating cpu output
    //cpu_gemm(MK_dense_matrix_no_organized, KN_dense_matrix_no_organized, cpu_output, V, F, V);

  for(int v=0;v<V;v++)
	for(int f=0;f<F;f++)
		cpu_output[v*F+f]=0;

   for(int v=0;v<V;v++) {
	int Nb=MK_row_pointer[v+1]-MK_row_pointer[v];    
	    for(int f=0;f<F;f++)  
	      for(int n=0;n<Nb;n++)
	            cpu_output[v*F+f]+=MK_sparse_matrix[MK_row_pointer[v]+n]*KN_dense_matrix[MK_col_id[MK_row_pointer[v]+n]*F+f];
				}

    //Comparing the results
    for(int i=0;i<V; i++) {
        for(int j=0; j<F; j++) {
            float difference=fabs(cpu_output[i*F+j]-acc_output[i*F+j]);
            if(difference > EPSILON) {
                std::cout << "ERROR position (" << i << "," << j <<  "): Value out simulator: " << acc_output[i*F+j] << ". Value out CPU: " << cpu_output[i*F+j] << std::endl;
                std::cout << "\033[1;31mT test not passed\033[0m" << std::endl;
		//delete[] MK_dense_matrix;
                delete[] KN_dense_matrix;
		//delete[] MK_dense_matrix_no_organized;
		delete[] KN_dense_matrix_no_organized;
                delete[] MK_col_id;
                delete[] MK_row_pointer;
                delete[] MK_sparse_matrix;
                delete[] cpu_output;
                delete[] acc_output;
                delete stonne_instance;	
                
                
                
             
	}
	}}
	
    //If the code does not stop then the TEST is correct
    std::cout << "\033[1;32mTest passed correctly \033[0m" << std::endl << std::endl;
    //delete[] MK_dense_matrix_no_organized;
    delete[] KN_dense_matrix_no_organized;
    //delete[] MK_dense_matrix;
    delete[] KN_dense_matrix;
    delete[] MK_col_id;
    delete[] MK_row_pointer;
    delete[] MK_sparse_matrix;
    delete[] cpu_output;
    delete[] acc_output;
    delete stonne_instance;
	
	//////////////////////////////////////////////////////////////////////// Step 2: DenseGEMM phase ///////////////////////////////////////////////////////////////////////////////////////////
	     		std::cout<<"\n\n\n\n------Simulating DenseGEMM-----------------------------------------------------------------------------------------------------------\n\n\n";

	 layer_name="DenseTestLayer";
		
		float* MK_dense_matrix_no_organized2 = generateMatrixDense(V, F, 0);
        	float* MK_dense_matrix2 = new float[V*F];
		
		 for(int i=0; i<V*F; i++) {
       		 MK_dense_matrix2[i]=MK_dense_matrix_no_organized2[i];

    		}
		
    		
    		float* KN_dense_matrix_no_organized2 = generateMatrixDense(F, G, 0);
    		float* KN_dense_matrix2 = new float[F*G];

		


   		 for(int i=0; i<F*G; i++) {
       		 KN_dense_matrix2[i]=KN_dense_matrix_no_organized2[i];
    		}
    		
    		unsigned int* clocked_op2 = new unsigned int [V*G];
    
    		float* output_cpu2 = new float[V*G];
    		float* acc_output2 = new float[V*G];
                
                
                //Configuring STONNE to run DenseGEMM.
                stonne_cfg2.m_MSNetworkCfg.ms_size=pe_cmb;
                stonne_cfg2.m_SDMemoryCfg.n_read_ports=dn_bw_cmb;
                stonne_cfg2.m_SDMemoryCfg.n_write_ports=rn_bw_cmb; 
                stonne_cfg2.m_ASNetworkCfg.accumulation_buffer_enabled=1;
                stonne_cfg2.m_SDMemoryCfg.mem_controller_type=MAERI_DENSE_WORKLOAD;
                stonne_cfg2.print_stats_enabled=0; 
                if(T_Fc>1)
               		stonne_cfg2.m_ASNetworkCfg.reduce_network_type=get_type_reduce_network_type("ASNETWORK");
               	else
                	stonne_cfg2.m_ASNetworkCfg.reduce_network_type=get_type_reduce_network_type("TEMPORALRN");
                        		

		Stonne* stonne_instance2 = new Stonne(stonne_cfg2);
		


    	stonne_instance2->loadDenseGEMM(layer_name, G, F, V, MK_dense_matrix2, KN_dense_matrix2, acc_output2, CNN_DATAFLOW); //Loading Sparse Dense
    stonne_instance2->loadClocking(clocked_op2);
        stonne_instance2->loadGEMMTile(T_G, T_Fc, T_Vc); //Loading the tile
    stonne_instance2->run(); //Running the simulator
    
    unsigned int* clocked_op_cp = new unsigned int[V*G];         
    for(int i=0;i<V*G;i++)
    	  clocked_op_cp[i]=clocked_op2[i];
    
    for(int i=0;i<V;i++)
    	for(int j=0;j<G;j++)
    	    	clocked_op2[i*G+j]=clocked_op_cp[j*V+i];
    
    unsigned int c_cycles = stonne_instance2->getCycles();
   std::cout<<"\n\n------Timestamps for simulated GEMM-----------------------------------------------------------------------------------------------------------\n";
    for(int i=G-1;i<V*G;i+=G)
    	std::cout<<clocked_op2[i]<<"\t";
    /** CHECKING the results to make sure that the output is correct  **/
//std::cout<<¨\n\n¨;   
std::cout << "Running CPU version to compare results" << std::endl;
    //Generating cpu output
        sequential_layer(1, F, 1, V, 1, 1, G, F, 1, KN_dense_matrix2, MK_dense_matrix2, output_cpu2); //Supposes that MK=inputs (M=batch size) and KN=filters (N=number of filters)
    //Comparing the results

		  for(int i=0;i<V*G; i++) {
        float difference=fabs(acc_output2[i]-output_cpu2[i]);
        if(difference > EPSILON) {
            std::cout << "ERROR position " << i <<  ": Value ofmap simulator: " << acc_output2[i] << ". Value ofmap CPU: " << output_cpu2[i] << std::endl;
            std::cout << "\033[1;31mT test not passed\033[0m" << std::endl;
            delete[] MK_dense_matrix2;
            delete[] KN_dense_matrix2;
            delete[] acc_output2;
            delete[] output_cpu2;
            delete stonne_instance;
            assert(false); //Always false
            
        }
    }


    //If the code does not stop then the TEST is correct
    std::cout << "\033[1;32mTest passed correctly \033[0m" << std::endl << std::endl;


   std::cout<<"\n\n------Statistics for simulated GEMM-----------------------------------------------------------------------------------------------------------\n";

   MSNetworkStats mulstats = stonne_instance2->getMultiplierNetworkStats();
    std::cout<<"\n\nRF Stats:\n";
    std::cout<<"RF weight reads = "<<mulstats.n_l1_weight_reads<<"\n";
    std::cout<<"RF weight writes = "<<mulstats.n_l1_weight_writes<<"\n";
    std::cout<<"RF input reads = "<< mulstats.n_l1_input_reads<<"\n";
    std::cout<<"RF input writes = "<< mulstats.n_l1_input_writes<<"\n";
    std::cout<<"RF psum reads = "<<mulstats.n_l1_psum_reads<<"\n";
    std::cout<<"RF psum writes = "<<mulstats.n_l1_psum_writes<<"\n";
    //std::cout<<"Number of multiplications = "<<mulstats.n_multiplications<<"\n";
    //std::cout<<"Number of link traversals = "<<mulstats.n_local_network_traversals<<"\n";
    
    std::cout<<"\n\nAdder stats: \n";
    if(T_Fc==1)
	std::cout<<"Reduction network is linear since reduction is temporal\n";
    else
    {
    	ASNetworkStats adderstats = stonne_instance2->getASNetworkStats();
    	std::cout<<"Number of link traversals = "<<adderstats.n_total_traversals<<"\n";
    }

    DSNetworkStats diststats = stonne_instance2->getDSNetworkStats();
    std::cout<<"\n\nDistribution stats: \n";
    std::cout<<"Number of link traversals = "<<diststats.n_total_traversals<<"\n";
    
    SDMemoryStats memstats = stonne_instance2->getMemoryStats();
    std::cout<<"\n\nGlobal Buffer stats: \n";
    std::cout << "Global SRAM weight reads = " << memstats.n_SRAM_weight_reads <<"\n";
    std::cout << "Global SRAM input reads = " << memstats.n_SRAM_input_reads <<"\n";
    std::cout << "Global SRAM psum reads = " << memstats.n_SRAM_psum_reads << "\n";
    std::cout << "Global SRAM psum writes = " << memstats.n_SRAM_psum_writes << "\n";


    delete[] MK_dense_matrix2;
    delete[] KN_dense_matrix2;
    delete[] acc_output2;
    delete[] output_cpu2;
    delete stonne_instance2;

///////////////////////////////////////////////////////////////////// Step 3: Analytical Modelling //////////////////////////////////////////////////////////////////////////////////////////////////////////////
	if(T_Va==1){ 
                // This calculation is for PP and involves estimation of runtime through timestamps.
                n_cycles = clocked_op[T_Vp*F-1];
    	        n_cycles += (clocked_op[2*T_Vp*F-1]-clocked_op[T_Vp*F-1]>clocked_op2[T_Vp*G-1])?clocked_op[2*T_Vp*F-1]-clocked_op[T_Vp*F-1]:clocked_op2[T_Vp*G-1];
    	        for(int i=2;i<V/T_Vp;i++)
    	        	n_cycles+=(clocked_op[(i*F+F)*T_Vp-1]-clocked_op[((i-1)*F+F)*T_Vp-1]>clocked_op2[((i-1)*G+G)*T_Vp-1]-clocked_op2[((i-2)*G+G)*T_Vp-1])?clocked_op[T_Vp*(i*F+F)-1]-clocked_op[((i-1)*F+F)*T_Vp-1]:clocked_op2[((i-1)*G+G)*T_Vp-1]-clocked_op2[((i-2)*G+G)*T_Vp-1];
		
		//////// Add the zero remaindered tile if it exists    	        
                if(V%T_Vp!=0)
		{
			n_cycles+=(clocked_op[V*F-1]-clocked_op[(V/T_Vp)*T_Vp*F-1] > clocked_op2[(V/T_Vp)*T_Vp*G-1]-clocked_op2[(V/T_Vp-1)*T_Vp*G-1]) ? clocked_op[V*F-1]-clocked_op[(V/T_Vp)*T_Vp*F-1] : clocked_op2[(V/T_Vp)*T_Vp*G-1]-clocked_op2[(V/T_Vp-1)*T_Vp*G-1];
		
		n_cycles+=clocked_op2[V*G-1]-clocked_op2[(V/T_Vp)*T_Vp*G-1];
		}	
                else
			n_cycles+=clocked_op2[V*G-1]-clocked_op2[(V-T_Vp)*G-1];
		std::cout<<"\n\n------------------------------------------------Inter-phase stats--------------------------------------------------\n";
		std::cout<<"\n\n------------------------------------------------PP_AC(VFN, VGF), numPEs=pe_agg+pe_cmb------------------------------\n";
		std::cout<<"\n\n\nNumber of PP cycles = "<<n_cycles<<"\n\n\nNumber of PP buffers utilized = "<<n_buffers<<"\n\n\n";
	}

	
	int* clk2 = new int[V/T_Va+(V%T_Va!=0)];		////////This becomes ceil(V/T_Va)	

	
	        int max;
		int m_cycles = 0;
		for(int i=0;i<V/T_Va;i++)	
		{
			max=0;
			for(int j=0;j<T_Va;j++)
			{
				int v_index=i*T_Va+j;
				if(v_index==0)
					max=clocked_op[F-1];
				else{
				if(clocked_op[v_index*F+F-1]-clocked_op[(v_index-1)*F+F-1]>max)
					max=clocked_op[v_index*F+F-1]-clocked_op[(v_index-1)*F+F-1];
				}
			}
		m_cycles+=max;
		clk2[i]=m_cycles;
		//std::cout<<"\t"<<max;
		}
		
               if(V%T_Va!=0)                                 ////////If V%T_Va!=0, then we also need to consider the max of last few and add it to clk2 
		{
			for(int i=T_Va*(V/T_Va);i<V;i++)
			{
				max=0;
				if(clocked_op[i*F+F-1]-clocked_op[(i-1)*F+F-1]>max)
					max=clocked_op[i*F+F-1]-clocked_op[(i-1)*F+F-1];
			}
		m_cycles+=max;
		clk2[V/T_Va]=m_cycles;
		//std::cout<<"\t rem - "<<max;
		}
	
			if(T_Va!=1){ /////////////////// SpMM always takes T_M=1, we calculate the exact runtime for T_Va>1 using an analytical model ///////////////////////
                
                n_cycles = clk2[T_Vp/T_Va-1];
    	        n_cycles += (clk2[2*T_Vp/T_Va-1]-clk2[T_Vp/T_Va-1]>clocked_op2[T_Vp*G-1])?clk2[2*T_Vp/T_Va-1]-clk2[T_Vp/T_Va-1]:clocked_op2[T_Vp*G-1];
    	        for(int i=2;i<V/T_Vp;i++)
    	        	n_cycles+=(clk2[(i+1)*T_Vp/T_Va-1]-clk2[i*T_Vp/T_Va-1]>clocked_op2[((i-1)*G+G)*T_Vp-1]-clocked_op2[((i-2)*G+G)*T_Vp-1])?clk2[(i+1)*T_Vp/T_Va-1]-clk2[i*T_Vp/T_Va-1]:clocked_op2[((i-1)*G+G)*T_Vp-1]-clocked_op2[((i-2)*G+G)*T_Vp-1];
		
		//////// Add the zero remaindered tile if it exists    	        
                if(V%T_Vp!=0)
		{
			////////
			n_cycles+=(clk2[V/T_Va-1+(V%T_Va!=0)]-clk2[(V/T_Vp)*T_Vp/T_Va-1]>clocked_op2[(V/T_Vp)*T_Vp*G-1]-clocked_op2[(V/T_Vp-1)*G*T_Vp-1])?clk2[V/T_Va-1+(V%T_Va!=0)]-clk2[(V/T_Vp)*T_Vp/T_Va-1]:clocked_op2[(V/T_Vp)*T_Vp*G-1]-clocked_op2[(V/T_Vp-1)*G*T_Vp-1];
			n_cycles+=clocked_op2[V*G-1]-clocked_op2[(V/T_Vp)*T_Vp*G-1];
			
		}			
                else    	        
		n_cycles+=clocked_op2[V*G-1]-clocked_op2[(V-T_Vp)*G-1];
		std::cout<<"\n\n------------------------------------------------Inter-phase stats--------------------------------------------------\n";
		std::cout<<"\n\n------------------------------------------------PP_AC(VFN, VGF), numPEs=pe_agg+pe_cmb------------------------------\n";
		std::cout<<"\n\n\nNumber of PP cycles = "<<n_cycles<<"\n\nNumber of PP buffers utilized = "<<n_buffers<<"\n\n";
	}

		//////// Modifying SP. SP-Optimized requires (VFN,VFG). We also estimate the psum overhead and add it to the stats computed by (VFN,VGF)
		std::cout<<"\n\n------------------------------------------------Seq_AC(VFN, VGF), numPEs=pe_agg=pe_cmb------------------------------\n";
		std::cout<<"\n\nSimulated Aggregation Cycles (T_Va=1) = "<<a_cycles<<"\n\n";
		std::cout<<"\n\nCombination Cycles = "<<c_cycles<<"\n\n";
		std::cout<<"\n\nAggregation Cycles = "<<m_cycles<<"\n\n";
		int s_cycles=m_cycles+c_cycles;
		int s_buffers=V*F;
		std::cout<<"\n\nSeqential cycles = "<<s_cycles<<"\t Intermediate buffering for sequential = "<<s_buffers<<"\n\n";
		int l_cycles=(V/T_Vc+(V%T_Vc!=0))*(F/T_Fc+(F%T_Fc!=0))*pe_cmb/dn_bw_cmb;

		std::cout<<"\n\n--------------------------------SP_AC(VFN, VFG)(Psum overhead of VFG considered), numPEs=pe_agg=pe_cmb------------------------------\n";
		std::cout<<"\n\nRedistribution Cycles saved= "<<l_cycles<<"\n\n"; //////////////// Redistribution cycles saved by SP /////////////////
		int sp_buffers = G*T_Vc*(F/T_Fc+(F%T_Fc!=0));
		int p_cycles = ((sp_buffers/pe_cmb)+(sp_buffers%pe_cmb!=0))*(V/T_Vc+(V%T_Vc!=0));	// Adding the partial sum overhead of VFN,VFG
		if(F==T_Fc)
		p_cycles=0;
		std::cout<<"\n\nPsum overhead Cycles = "<<p_cycles<<"\n\n";
		int sp_cycles = s_cycles-l_cycles+p_cycles;
		std::cout<<"\n\nSP dataflow cycles = "<<sp_cycles<<"\t Intermediate buffering for interleaving = "<<sp_buffers<<"\n\n";

		////////////Printing buffer stats and analytical model here
		
		std::cout<<"------------------------------------Global Buffer stats-----------------------------------------------------------------------------------\n";
		int sp_psum_accesses=2*G*V*(F/T_Fc+(F%T_Fc!=0));
	        if(F==T_Fc)
			sp_psum_accesses=0;
		std::cout<<"\n\nSP Psum accesses = \n\n"<<sp_psum_accesses;
		int pp_global_accesses=memstats.n_SRAM_weight_reads+memstats.n_SRAM_input_reads+memstats.n_SRAM_psum_reads+memstats.n_SRAM_psum_writes+
			memstatssp.n_SRAM_weight_reads+memstatssp.n_SRAM_input_reads+memstatssp.n_SRAM_psum_reads+memstatssp.n_SRAM_psum_writes;
		std::cout<<"\n\nPP and SEQ Global Buffer Accesses = \n\n"<<pp_global_accesses;
		int sp_global_accesses=memstats.n_SRAM_input_reads+memstats.n_SRAM_psum_reads+memstats.n_SRAM_psum_writes+
			memstatssp.n_SRAM_weight_reads+memstatssp.n_SRAM_input_reads+sp_psum_accesses;
		
		std::cout<<"\n\nSP Global Buffer Accesses = \n\n"<<sp_global_accesses;	
		std::cout<<"\n\n\n\n";	

		

	}
	
	
	
	
	
	
	if(granularity==COLUMN)				//TODO: Future work
	{
		int n_cycles=0;
		T_Fp = (T_Fa>T_Fc) ? T_Fa : T_Fc;
		int n_buffers = 2*V*T_Fp;
		//Run Stonne in loop, Update n_cycles
		//Need to do loop order manipulations
		std::cout<<"Number of cycles = "<<n_cycles<<"\nNumber of buffers = "<<n_buffers;
		
	}
	

	return 0;
	
}
