
/**
 * @author Ian Johnson
 * @version 1.0
 */


//daq register definitions
#define DAQ_REG_CTRL                  1
#define DAQ_REG_CHIP_CMDS             2
#define DAQ_REG_STATIC_BITS           3
#define DAQ_REG_CLK_ROW_CLK_NTIMES    3
#define DAQ_REG_SHIFT_IN_32           3
#define DAQ_REG_READOUT_NROWS         3
#define DAQ_REG_SEND_N_TESTPULSES     3

#define DAQ_REG_NEXPOSURES            3
#define DAQ_REG_EXPOSURE_TIMER        4 // == (31 downto 3) * 10^(2 downto 0) 
#define DAQ_REG_EXPOSURE_REPEAT_TIMER 5 // == (31 downto 3) * 10^(2 downto 0) 
#define DAQ_REG_SUBFRAME_EXPOSURES    6
#define DAQ_REG_STATUS                7 //also pg and fifo status register

#define DAQ_CTRL_RESET              0x80000000
#define DAQ_CTRL_START              0x40000000
#define ACQ_CTRL_START              0x50000000 //this is 0x10000000 (acq) | 0x40000000 (daq)
#define DAQ_CTRL_STOP               0x00000000

//direct chip commands to the DAQ_REG_CHIP_CMDS register
#define DAQ_SET_STATIC_BIT          0x00000001
#define DAQ_RESET_COMPLETELY        0x0000000E
#define DAQ_RESET_PERIPHERY         0x00000002
#define DAQ_RESET_PIXEL_COUNTERS    0x00000004
#define DAQ_RESET_COLUMN_SELECT     0x00000008

#define DAQ_STORE_IMAGE             0x00000010
#define DAQ_RELEASE_IMAGE_STORE     0x00000020

#define DAQ_SEND_A_TOKEN_IN         0x00000040
#define DAQ_CLK_ROW_CLK_NTIMES      0x00000080
#define DAQ_SERIALIN_SHIFT_IN_32    0x00000100
#define DAQ_LOAD_16ROWS_OF_TRIMBITS 0x00000200

#define DAQ_IGNORE_INITIAL_CRAP     0x00000400 //crap before readout
#define DAQ_READOUT_NROWS           0x00000800
#define DAQ_CLKOUT_LAST_4_BITS_AND_RETURN_TO_START 0x00001000 //last 4 bit of data in the last frame

#define DAQ_RELEASE_IMAGE_STORE_AFTER_READOUT  0x00002000
#define DAQ_RESET_PIXEL_COUNTERS_AFTER_READOUT 0x00004000

#define DAQ_CLK_ROW_CLK_TO_SELECT_NEXT_ROW     0x00008000
#define DAQ_CLK_MAIN_CLK_TO_SELECT_NEXT_PIXEL  0x00010000
#define DAQ_SEND_N_TEST_PULSES                 0x00020000

#define DAQ_CHIP_CONTROLLER_HALF_SPEED         0x00040000 //everything at 100 MHz (50MHz ddr readout)
#define DAQ_CHIP_CONTROLLER_QUARTER_SPEED      0x00080000 //everything at  50 MHz (25MHz ddr readout)
#define DAQ_CHIP_CONTROLLER_SUPER_SLOW_SPEED   0x000c0000 //everything at  ~200 kHz (200 kHz MHz ddr readout)

#define DAQ_FIFO_ENABLE                        0x00100000

//direct chip commands to the DAQ_REG_CHIP_CMDS register
#define DAQ_NEXPOSURERS_SAFEST_MODE_ROW_CLK_BEFORE_MODE 0x00200000 //row clk is before main clk readout sequence
#define DAQ_NEXPOSURERS_NORMAL_NONPARALLEL_MODE         0x00400000 //expose ->readout ->expose -> ..., with store is always closed
#define DAQ_NEXPOSURERS_PARALLEL_MODE                   0x00600000 //parallel acquire/read mode

//DAQ_NEXPOSURERS_READOUT_COMPLETE_IMAGES is old now hard-wired in the firmware that every image comes with a header
//#define DAQ_NEXPOSURERS_READOUT_COMPLETE_IMAGES    0x00800000 //DAQ_IGNORE_INITIAL_CRAP and DAQ_CLKOUT_LAST_4_BITS_AND_RETURN_TO_START

#define DAQ_NEXPOSURERS_EXTERNAL_ENABLING          0x01000000 
#define DAQ_NEXPOSURERS_EXTERNAL_ENABLING_POLARITY 0x02000000
#define DAQ_NEXPOSURERS_EXTERNAL_TRIGGER_POLARITY  0x04000000

#define DAQ_NEXPOSURERS_INTERNAL_ACQUISITION          0x00000000 //internally controlled
#define DAQ_NEXPOSURERS_EXTERNAL_ACQUISITION_START    0x08000000 //external acquisition start
#define DAQ_NEXPOSURERS_EXTERNAL_IMAGE_START          0x10000000 //external image start
#define DAQ_NEXPOSURERS_EXTERNAL_IMAGE_START_AND_STOP 0x18000000 //externally controlly, external image start and stop

#define DAQ_NEXPOSURERS_ACTIVATE_AUTO_SUBIMAGING      0x20000000 
#define DAQ_NEXPOSURERS_ACTIVATE_RATE_CORRECTION      0x40000000 

//#define DAQ_MASTER_HALF_MODULE                        0x80000000 currently not used


//chips static bits 
#define DAQ_STATIC_BIT_PROGRAM      0x00000001
#define DAQ_STATIC_BIT_M4           0x00000002 //these are the status bits, not bit mode
#define DAQ_STATIC_BIT_M8           0x00000004 //these are the status bits, not bit mode
#define DAQ_STATIC_BIT_M12          0x00000000 //these are the status bits, not bit mode, ie. "00" is 12 bit mode
#define DAQ_STATIC_BIT_CHIP_TEST    0x00000008
#define DAQ_STATIC_BIT_ROTEST       0x00000010
#define DAQ_CS_BAR_LEFT             0x00000020
#define DAQ_CS_BAR_RIGHT            0x00000040


//status flags
#define DAQ_STATUS_DAQ_RUNNING        0x01
#define DAQ_DATA_COLLISION_ERROR      0x02


#define DAQ_STATUS_CURRENT_M4         0x04 
#define DAQ_STATUS_CURRENT_M8         0x08
#define DAQ_STATUS_CURRENT_M12        0x00 //in 12 bit mode both are cleared
#define DAQ_STATUS_CURRENT_TESTMODE   0x10
#define DAQ_STATUS_TOKEN_OUT          0x20 
#define DAQ_STATUS_SERIAL_OUT         0x40 
#define DAQ_STATUS_PIXELS_ARE_ENABLED 0x80 
#define DAQ_STATUS_DAQ_RUN_TOGGLE     0x200

//data delay registers
#define CHIP_DATA_OUT_DELAY_REG_CTRL       1
#define CHIP_DATA_OUT_DELAY_REG2           2
#define CHIP_DATA_OUT_DELAY_REG3           3
#define CHIP_DATA_OUT_DELAY_REG4           4
#define CHIP_DATA_OUT_DELAY_SET            0x20000000

//module configuration
#define TOP_BIT_MASK				0x00f
#define MASTER_BIT_MASK				0x200
#define NORMAL_MODULE_BIT_MASK		0x400

// Master Slave Top Bottom Definition
#define MODULE_CONFIGURATION_MASK 0x84
//Software Configuration
#define MASTERCONFIG_OFFSET			0x160		//0x20 * 11 (P11)
#define MASTER_BIT					0x1
#define	OVERWRITE_HARDWARE_BIT		0x2
#define DEACTIVATE_BIT				0x4

#define FPGA_TEMP_OFFSET			0x200

#define TXM_DELAY_LEFT_OFFSET		0x180
#define TXM_DELAY_RIGHT_OFFSET		0x1A0
#define TXM_DELAY_FRAME_OFFSET		0x1C0
#define TXM_FLOW_CONTROL_10G		0x140

//command memory
#define LEFT_OFFSET					0x0
#define RIGHT_OFFSET				0x100

#define FIRST_CMD_PART1_OFFSET		0x8
#define FIRST_CMD_PART2_OFFSET		0xc
#define SECOND_CMD_PART1_OFFSET		0x10
#define SECOND_CMD_PART2_OFFSET		0x14
#define COMMAND_COUNTER_OFFSET		0x18
#define STOP_ACQ_OFFSET 			0x1c
#define STOP_ACQ_BIT				0x40000000
#define TWO_REQUESTS_OFFSET			0x1c
#define TWO_REQUESTS_BIT			0x80000000

//version
#define FIRMWARE_VERSION_OFFSET		0x4
#define FIRMWARESOFTWARE_API_OFFSET 0x0

#define FRAME_NUM_RESET_OFFSET		0xA0

//temp so far
#define FEB_REG_STATUS              0xa

//1g counters
#define ONE_GIGA_LEFT_INDEX_LSB_COUNTER		0x04
#define ONE_GIGA_LEFT_INDEX_MSB_COUNTER		0x24

#define ONE_GIGA_LEFT_TXN_DELAY_COUNTER		0x104
#define ONE_GIGA_LEFT_FRAME_DELAY_COUNTER	0x124

#define ONE_GIGA_RIGHT_INDEX_LSB_COUNTER	0x44
#define ONE_GIGA_RIGHT_INDEX_MSB_COUNTER	0x64

#define ONE_GIGA_RIGHT_TXN_DELAY_COUNTER	0x144
#define ONE_GIGA_RIGHT_FRAME_DELAY_COUNTER	0x164

//10g counters
#define TEN_GIGA_LEFT_INDEX_LSB_COUNTER		0x84
#define TEN_GIGA_LEFT_INDEX_MSB_COUNTER		0xa4

#define TEN_GIGA_LEFT_TXN_DELAY_COUNTER		0x184
#define TEN_GIGA_LEFT_FRAME_DELAY_COUNTER	0x1a4

#define TEN_GIGA_RIGHT_INDEX_LSB_COUNTER	0xc4
#define TEN_GIGA_RIGHT_INDEX_MSB_COUNTER	0xe4

#define TEN_GIGA_RIGHT_TXN_DELAY_COUNTER	0x1c4
#define TEN_GIGA_RIGHT_FRAME_DELAY_COUNTER	0x1e4

// udp header (position, id)
#define UDP_HEADER_A_OFST					0x00C0
#define UDP_HEADER_B_OFST					0x00E0
#define UDP_HEADER_X_OFST					(0)
#define UDP_HEADER_X_MSK					(0xFFFF << UDP_HEADER_X_OFST)
#define UDP_HEADER_ID_OFST					(16)
#define UDP_HEADER_ID_MSK					(0xFFFF << UDP_HEADER_ID_OFST)
#define UDP_HEADER_Z_OFST					(0)
#define UDP_HEADER_Z_MSK					(0xFFFF << UDP_HEADER_Z_OFST)
#define UDP_HEADER_Y_OFST					(16)
#define UDP_HEADER_Y_MSK					(0xFFFF << UDP_HEADER_Y_OFST)


















