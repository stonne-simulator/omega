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

enum GRAN{ELEMENT, ROW, COLUMN, SEQUENTIAL};
enum KERNELORDER{AC,CA};
int main()
{

	//EXTERNAL INPUTS - Loop order, Kernel order, dimensions - V,F,G. Tile sizes - T_Va, T_Fa, T_N, T_Vc, T_Fc, T_G, adjacency matrix
	GRAN granularity;
	//KERNELORDER order;
	KERNELORDER order = AC;
	//char loop_order_agg[4], loop_order_cmb[4];
	int V,F,G,D;
	V=1313;F=136;G=2,D=12;
	int T_Va, T_N, T_Fa, T_Vc, T_G, T_Fc;
	T_Va=1;T_N=1;T_Fa=136;T_Vc=1;T_G=1;T_Fc=136;
	int T_Vp, T_Fp;
	int pe_agg=512, pe_cmb=512;
	char loop_order_agg[4]="VFN";
	char loop_order_cmb[4]="VGF";
	if(order==AC)
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
	else if(order==CA)
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
	
	if(granularity==ELEMENT)
	{
		int n_cycles=0;
		T_Vp = (T_Va>T_Vc) ? T_Va : T_Vc;
		T_Fp = (T_Fa>T_Fc) ? T_Fa : T_Fc;
		int n_buffers = 2*T_Vp*T_Fp;
		//Run Stonne in loop, Update n_cycles
		//Need to do loop order manipulations
		std::cout<<"Number of cycles = "<<n_cycles<<"\nNumber of buffers = "<<n_buffers;
	}
	
	if(granularity==ROW)
	{
		int n_cycles=0;
		T_Vp = (T_Va>T_Vc) ? T_Va : T_Vc;
		int n_buffers = 2*T_Vp*F;
		
		float EPSILON=0.05;
    		std::string layer_name="SparseDenseTestLayer";
		
		
		float* MK_dense_matrix_no_organized = generateMatrixDenseSampled(V, V, D);
		/*float* MK_dense_matrix_no_organized =  new float[V*V];
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

*/
		
        	float* MK_dense_matrix = new float[V*V];
		
		 for(int i=0; i<V*V; i++) {
       		 MK_dense_matrix[i]=MK_dense_matrix_no_organized[i];

    		}
//	
//		

    		    int nnz=0;
    		int* MK_col_id = generateMinorIDFromDense(MK_dense_matrix, V, V, nnz, GEN_BY_ROWS);
    		int* MK_row_pointer = generateMajorPointerFromDense(MK_dense_matrix, V , V, GEN_BY_ROWS);
		
		
//		std::cout<<"\nMK Dense matrix - \n";
//	
//	for(int i=0; i<V; i++) {
//        for(int j=0; j<V; j++) {
//            std::cout << MK_dense_matrix[i*V+j] << "\t";
//	}
//	std::cout << "\n";
//    }

		std::cout << "\n\nMK Row pointer - \n"; 
		for(int i=0;i<=V;i++)
    		{
    			std::cout<<MK_row_pointer[i]<<"\t";
    		}
    		std::cout << "\n\n";		
		
    		//Generating sparse matrix
    		float* MK_sparse_matrix = generateMatrixSparseFromDenseNoBitmap(MK_dense_matrix, V, V, GEN_BY_ROWS);
    
    		float* KN_dense_matrix_no_organized = generateMatrixDense(V, F, 0);
    		float* KN_dense_matrix = new float[V*F];

		


   		 for(int i=0; i<V*F; i++) {
       		 KN_dense_matrix[i]=KN_dense_matrix_no_organized[i];
    		}
    		
    		unsigned int* clocked_op = new unsigned int [V*F];
    
    		float* cpu_output = new float[V*F];
    		float* acc_output = new float[V*F];
    
    		
		
		//Configure STONNE
		Config stonne_cfg1, stonne_cfg2;
		
		stonne_cfg1.m_MSNetworkCfg.ms_size=pe_agg;
                stonne_cfg1.m_SDMemoryCfg.n_read_ports=pe_agg;
                stonne_cfg1.m_SDMemoryCfg.n_write_ports=pe_agg; 
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
    stonne_instance->run(); //Running the simulator
    unsigned int a_cycles = stonne_instance->getCycles();
    for(int i=F-1;i<V*F;i+=F)
    	std::cout<<clocked_op[i]<<"\t";
    /** CHECKING the results to make sure that the output is correct  **/
	//std::cout<<¨\n\n¨;
   /* std::cout << "Running CPU version to compare results" << std::endl;
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

    delete[] MK_dense_matrix_no_organized;
    delete[] KN_dense_matrix_no_organized;
    delete[] MK_dense_matrix;
    delete[] KN_dense_matrix;
    delete[] MK_col_id;
    delete[] MK_row_pointer;
    delete[] MK_sparse_matrix;
    delete[] cpu_output;
    delete[] acc_output;
    delete stonne_instance;
	
	
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
                
                
                
                stonne_cfg2.m_MSNetworkCfg.ms_size=pe_cmb;
                stonne_cfg2.m_SDMemoryCfg.n_read_ports=pe_cmb;
                stonne_cfg2.m_SDMemoryCfg.n_write_ports=pe_cmb; 
                stonne_cfg2.m_ASNetworkCfg.accumulation_buffer_enabled=1;
                stonne_cfg2.m_SDMemoryCfg.mem_controller_type=MAERI_DENSE_WORKLOAD;
                stonne_cfg2.print_stats_enabled=0; 
                if(T_Fc>1)
               		stonne_cfg2.m_ASNetworkCfg.reduce_network_type=get_type_reduce_network_type("ASNETWORK");
               	else
                	stonne_cfg2.m_ASNetworkCfg.reduce_network_type=get_type_reduce_network_type("TEMPORALRN");
                        		

		Stonne* stonne_instance2 = new Stonne(stonne_cfg2);
		
                
//TO BE CONTINUED
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
    for(int i=T_Vp*G-1;i<V*G;i+=T_Vp*G)
    	std::cout<<clocked_op2[i]<<"\t";
    /** CHECKING the results to make sure that the output is correct  **/
//std::cout<<¨\n\n¨;   
/*std::cout << "Running CPU version to compare results" << std::endl;
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
*/
    delete[] MK_dense_matrix2;
    delete[] KN_dense_matrix2;
    delete[] acc_output2;
    delete[] output_cpu2;
    delete stonne_instance2;


	if(T_Va==1){ 
                
                n_cycles = clocked_op[T_Vp*F-1];
    	        n_cycles += (clocked_op[2*T_Vp*F-1]-clocked_op[T_Vp*F-1]>clocked_op2[T_Vp*G-1])?clocked_op[2*T_Vp*F-1]-clocked_op[T_Vp*F-1]:clocked_op2[T_Vp*G-1];
    	        for(int i=2;i<V/T_Vp;i++)
    	        	n_cycles+=(clocked_op[(i*F+F)*T_Vp-1]-clocked_op[((i-1)*F+F)*T_Vp-1]>clocked_op2[((i-1)*G+G)*T_Vp-1]-clocked_op2[((i-2)*G+G)*T_Vp-1])?clocked_op[T_Vp*(i*F+F)-1]-clocked_op[((i-1)*F+F)*T_Vp-1]:clocked_op2[((i-1)*G+G)*T_Vp-1]-clocked_op2[((i-2)*G+G)*T_Vp-1];
    	        	n_cycles+=clocked_op2[V*G-1]-clocked_op2[(V-T_Vp)*G-1];
		std::cout<<"\n\n\nNumber of PP cycles = "<<n_cycles<<"\n\n\nNumber of PP buffers utilized = "<<n_buffers<<"\n\n\n";
	}

	
	int* clk2 = new int[V/T_Va];

	
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
		clk2[i]=max;
		std::cout<<"\t"<<max;
		}

			if(T_Va!=1){ 
                
                n_cycles = clk2[T_Vp/T_Va-1];
    	        n_cycles += (clk2[2*T_Vp/T_Va-1]>clocked_op2[T_Vp*G-1])?clk2[2*T_Vp/T_Va-1]:clocked_op2[T_Vp*G-1];
    	        for(int i=2;i<V/T_Vp;i++)
    	        	n_cycles+=(clk2[(i+1)*T_Vp/T_Va-1]>clocked_op2[((i-1)*G+G)*T_Vp-1]-clocked_op2[((i-2)*G+G)*T_Vp-1])?clk2[(i+1)*T_Vp/T_Va-1]:clocked_op2[((i-1)*G+G)*T_Vp-1]-clocked_op2[((i-2)*G+G)*T_Vp-1];
    	        	n_cycles+=clocked_op2[V*G-1]-clocked_op2[(V-T_Vp)*G-1];
		std::cout<<"\n\n\nNumber of PP cycles = "<<n_cycles<<"\n\n\nNumber of PP buffers utilized = "<<n_buffers<<"\n\n\n";
	}


		std::cout<<"\n\nA Cycles = "<<a_cycles<<"\n\n";
		std::cout<<"\n\nC Cycles = "<<c_cycles<<"\n\n";
		std::cout<<"\n\nAggregation Cycles = "<<m_cycles<<"\n\n";
		int s_cycles=m_cycles+c_cycles;
		int s_buffers=V*F;
		std::cout<<"\n\nSeqential cycles = "<<s_cycles<<"\t Intermediate buffering for sequential = "<<s_buffers<<"\n\n";
		int l_cycles=V/T_Vc*F/T_Fc;
		std::cout<<"\n\nL Cycles = "<<l_cycles<<"\n\n";
		int sp_buffers = G*T_Vc*F/T_Fc;
		int p_cycles = ((sp_buffers/pe_cmb)+(sp_buffers%pe_cmb!=0))*V/T_Vc;
		if(F==T_Fc)
		p_cycles=0;
		std::cout<<"\n\nPsum Cycles = "<<p_cycles<<"\n\n";
		int sp_cycles = s_cycles-l_cycles+p_cycles;
		std::cout<<"\n\nInterleaving cycles = "<<sp_cycles<<"\t Intermediate buffering for interleaving = "<<sp_buffers<<"\n\n";
	}
	
	
	
	
	
	
	if(granularity==COLUMN)
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
