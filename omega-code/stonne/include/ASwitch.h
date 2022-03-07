//Created 19/06/2019

#ifndef __ASwitch__h
#define __ASwitch__h

#include "types.h"
#include "DataPackage.h"
#include "Connection.h"
#include <vector>
#include "Fifo.h"
#include "Unit.h"
#include "Config.h"
#include "Stats.h"
/*
*/

class ASwitch : public Unit {
private:
    unsigned int level;                              // Level where the Adder is set in the tree
    unsigned int num_in_level;                       // Number of the Adder in the level
    unsigned int num_ms;                             // These three parameters are for routing. In hardware it is not neccesary since it is used a bit vector
    unsigned int input_ports;                        // Input port per branch (left, right).
    unsigned int output_ports;                       // output port 
    unsigned int forwarding_ports;                  //Intermedium links between some nodes of the same level
    unsigned int buffers_capacity;                   //Buffers size in bytes
    unsigned int port_width;                         //Bit width of each port

    unsigned int busID;                              //CollectionBus connected to this ASwitch
    unsigned int inputID;                            //Number of input of the Collection Bus busID connected to this AS.
    unsigned int accumulationBufferID;               //Number of the accumulation buffer if used
    //Inputs fifos. These Fifos are not in the picture of figure 3. However, we build it for simplicity and extensibility
    Fifo* input_psum_left_fifo;                            // Array of packages that are received 
    Fifo* input_psum_right_fifo;                           // Array of packages that are received
    Fifo* input_fw_fifo;                          // Array of packages that are received though the fw link
   
    // Output Fifos
    Fifo* output_psum_fifo;                      // Output fifo to the parent
    Fifo* output_fw_fifo;                        // Ouptut fifo to the fw link

    std::vector<DataPackage*> psums_created;         // Array of packages generated by this specific AS. Used to delete the packages at the end of exec
     
    unsigned int current_capacity;                   // the capacity must not exceed the input bw of the connection
    Connection* inputLeftConnection;                 // This is the input left connection of the Adder
    Connection* inputRightConnection;                // This is the input right connection of the adder
    Connection* outputConnection;                    // This is the output connection of the adder
    Connection* memoryConnection;                    // This is the connection of the adder with the memory collect bus to write the psum if neccesary

    // Related intermedimum links (forwardings). This links are NULL if the Adder has no connection with neighbour. (Look up the paper of MAERI for more details)
    Connection* forwardingConnection;                // connection from the neighbour Adder. In hardware we use two. One for inputs and another for outputs and the adder slect the correct one
                                                     // using the configuration. However, in software is enough with one link and a flag that indicates if it is going to receive or send data

    cycles_t latency;                                  // Number of cycles to compute a sum. This is configurable since can vary depending on the implementation technology and frequency

    //Configuration signals. Note the direction of the fl is also a configuration signal. 
    fl_t fl_direction;                               // RECEIVE or SEND. Each adder will be configured differently
    bool fw_enabled; //Set if the fw is enabled
    adderconfig_t  config_mode;                       // Configuration mode of the adder. ADD_2_1, ADD_3_1, ADD_1_1_PLUS_FW_1_1 or  FW_2_2
    adderoperation_t operation_mode;                  // Operation to perform by the AS (ADDER, COMPARATOR). DEFAULT=ADDER
    bool left_child_enabled;                          // Indicates if the as receives data from the left child
    bool right_child_enabled;                         // Indicates if the as receives data from the right child
    bool forward_to_memory;                           // To optimize. If the psum is already completed send directly to the collection bus

    //Operation functions. This functions can be changed in order to perform different types of length operations
    DataPackage* perform_operation_2_operands(DataPackage* pck_left, DataPackage* pck_right);    //Perform 2:1 sum
    DataPackage* perform_operation_3_operands(DataPackage* pck_left, DataPackage* pck_right, DataPackage* pck_forward); //Perform 3:1 sum
    //Routing functions. Based on the configuration of the switch, the routing is done in a different way.
    void route_2_1_config(); //2:1 routing. This perform a sum with left and right inputs and send the output to the parent link or fw link.
    void route_3_1_config(); //3;2 routing. This perform a 3:1 sum with left, right and fw link inputs and send the output to the parent
    void route_1_1_plus_fw_1_1_config(); //1:1 ADD plus 1:1 fw. Send one branch to the parent and the other to the fw link
    void route_fw_2_2_config(); //2:2 fw. This perform a forwarding of the 2 left and right inputs sending it through the parent link.

    cycles_t local_cycle;
    ASwitchStats aswitchStats; //To track the behaviour of the ASwitch


public:
    ASwitch(id_t id, std::string name, unsigned int level, unsigned int num_in_level, Config stonne_cfg); 
    ASwitch(id_t id, std::string name, unsigned int level, unsigned int num_in_level, Config stonne_cfg, Connection* inputLeftConnection, 
                   Connection* inputRightConnection, Connection* forwardingConnection, Connection* outputConnection, Connection* memoryConnection);
    ~ASwitch();

    //Connection setters
    void setInputLeftConnection(Connection* inputLeftConnection);                       // Set the input left connection of the Adder
    void setInputRightConnection(Connection* inputRightConnection);                     // Set the input right connection of the Adder
    void setOutputConnection(Connection* outputConnection);                             // Set the output connection of the Adder
    void setForwardingConnection(Connection* forwardingConnection);                     // Set the forwarding connection of the Adder 
    void setMemoryConnection(Connection* memoryConnection, unsigned int busID, unsigned int inputID); // Set the memory connection of the adder. Indicates the line connected to be shown in the stats.
    void resetSignals(); //Reset all configuration signals

    // Getters
    const unsigned int getLevel()      const {return this->level;}
    const unsigned int getNumInLevel() const {return this->num_in_level;}
    const unsigned int getInputPorts() const {return this->input_ports;}                            // Get the input bw per branch
    const unsigned int getOutputPorts() const {return this->output_ports;}                          // Get the output ports
    const unsigned int getForwardingPorts() const {return this->forwarding_ports;}                  // Get the forwarding ports
    const cycles_t getLatency()    const {return this->latency;}                        // Get the number of cycles to compute
    const fl_t getFlDirection()    const {return this->fl_direction;}                   // Get the direction of the fl for the current configuration. RECEIVE or SEND
    bool isFwEnabled()             const {return this->fw_enabled;}
    unsigned int getAccumulationBufferID()    const {return this->accumulationBufferID;}
    const adderconfig_t getConfig()           const {return this->config_mode;} // Get the current configuration of the adder. ADD_2_1, ADD_3_1, ADD_1_1_PLUS_FW_1_1 or  FW_2_2
    const adderoperation_t getOperationMode() const {return this->operation_mode;}             // Get the current operation configured (ADD or POOL) TODO POOL NOT SUPPORTED YET

    //Configuration settings (control signals)
    void setForwardingLinkDirection(fl_t fl_direction); // Set if the configuration of the FL for a certain AS is RECEIVE or SEND. This is a control signal
    void setConfigurationMode(adderconfig_t config_mode);    // Set the configuration mode of the adder. This depends on the compiler and it is a control signal. 
    void setChildsEnabled(bool left_child_enabled, bool right_child_enabled); //Indicates if the as receives data from left and right
    void setForwardingToMemoryEnabled(bool forwarding_to_memory);
    void setOperationMode(adderoperation_t operation_mode);


    // Functionality
    void send(); //Packages of data to be sent depending on routing. 
    void receive_childs(); //take data from inputLeftConnection and inputRightConnection and save it in its corresponding buffers 
    void receive_fwlink();

                     //Buffers seem not to appear in the Adder switch implemented in hardware (just there is a buffer in the output). However, they are just a logical unit used in the software model. 
                     //This is done like this to make easier possible future implementations.

    void cycle();   //Computing a cycle. Based on routing the AS decides where the data goes. 

    void printConfiguration(std::ofstream& out, unsigned int indent);  //This function prints the configuration of ASwitch such as  the operation mode, augmented link enabled, etc
    void printStats(std::ofstream& out, unsigned int indent);
    void printEnergy(std::ofstream& out, unsigned int indent);
    ASwitchStats getStats() {return this->aswitchStats;}
	

};

#endif

