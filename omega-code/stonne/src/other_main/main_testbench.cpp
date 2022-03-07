//Created by Francisco Munoz-Martinez on 17/06/2019                             

#include <iostream>
#include "MAERIModel.h"
#include "types.h"
#include <chrono>
#include <assert.h>
#include "testbench.h"

using namespace std;
int main(int argc, char** argv) {
//    hand_tests();                                                             
//    run_simple_tests();                                                       
      unsigned int num_ms=32;
      run_maeri_architecture_tests(LATE_SYNTHETIC, 32);
}

