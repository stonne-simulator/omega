//Created by Francisco Munoz-Martinez on 17/06/2019                             

#include <iostream>
#include "MAERIModel.h"
#include "types.h"
#include <chrono>
#include <assert.h>
#include "testbench.h"
#include "Tile.h"

using namespace std;
int main(int argc, char** argv) {

    string name_file="tile_configuration.txt";
    Tile tile(name_file);
    //tile is supposed to have the values
    std::cout << "T_R: " << tile.get_T_R() << std::endl;
    std::cout << "T_S: " << tile.get_T_S() << std::endl;
    std::cout << "T_C: " << tile.get_T_C() << std::endl;
    std::cout << "T_K: " << tile.get_T_K() << std::endl;
    std::cout << "T_N: " << tile.get_T_N() << std::endl;
    std::cout << "T_X': " << tile.get_T_X_() << std::endl;
    std::cout << "T_Y_': " << tile.get_T_Y_() << std::endl;
   


}

