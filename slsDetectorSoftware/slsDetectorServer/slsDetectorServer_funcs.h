#ifndef SERVER_FUNCS_H
#define SERVER_FUNCS_H

#include "sls_receiver_defs.h"
#include <stdlib.h>

// initialization functions
void basictests();
void init_detector(int);
int decode_function(int);
const char* getFunctionName(enum detFuncs func);
void function_table();
int M_nofunc(int);

// functions called by client
int exec_command(int);
int get_error(int);
int get_detector_type(int);
int set_number_of_modules(int);
int get_max_number_of_modules(int);
int set_external_signal_flag(int);
int set_external_communication_mode(int);
int get_id(int);
int digital_test(int);
int analog_test(int);
int enable_analog_out(int);
int calibration_pulse(int);
int set_dac(int);
int get_adc(int);
int write_register(int);
int read_register(int);
int write_memory(int);
int read_memory(int);
int set_channel(int);
int get_channel(int);
int set_all_channels(int);
int set_chip(int);
int get_chip(int);
int set_all_chips(int);
int set_module(int);
int get_module(int);
int set_settings(int);
int get_threshold_energy(int);
int set_threshold_energy(int);
int start_acquisition(int);
int stop_acquisition(int);
int start_readout(int);
int get_run_status(int);
int start_and_read_all(int);
int read_frame(int);
int read_all(int);
int set_timer(int);
int get_time_left(int);
int set_dynamic_range(int);
int set_readout_flags(int);
int set_roi(int);
int set_speed(int);
int execute_trimming(int);
int exit_server(int);
int lock_server(int);
int get_last_client_ip(int);
int set_port(int);
int update_client(int);
int send_update(int);
int configure_mac(int);
int load_image(int);
int set_master(int);
int set_synchronization(int);
int read_counter_block(int);
int reset_counter_block(int);
int calibrate_pedestal(int);
int enable_ten_giga(int);
int set_all_trimbits(int);
int set_ctb_pattern(int);
int write_adc_register(int);
int set_counter_bit(int);
int pulse_pixel(int);
int pulse_pixel_and_move(int);
int pulse_chip(int);
int set_rate_correct(int);
int get_rate_correct(int);
int set_network_parameter(int);
int program_fpga(int);
int reset_fpga(int);
int power_chip(int);
int set_activate(int);
int prepare_acquisition(int);


#endif
