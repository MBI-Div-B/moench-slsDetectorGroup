#ifndef SLSDETECTORSERVER_DEFS_H
#define SLSDETECTORSERVER_DEFS_H

#include "sls_detector_defs.h" 	//default dynamicgain in settings
#include "RegisterDefs.h"
#include <stdint.h>



#define GOODBYE 					(-200)
#define CTRL_SRVR_INIT_TIME_US		(300 * 1000)
//#define REQUIRED_FIRMWARE_VERSION 16


/* Struct Definitions */
typedef struct ip_header_struct {
	u_int16_t     ip_len;
	u_int8_t      ip_tos;
	u_int8_t      ip_ihl:4 ,ip_ver:4;
	u_int16_t     ip_offset:13,ip_flag:3;
	u_int16_t     ip_ident;
	u_int16_t     ip_chksum;
	u_int8_t      ip_protocol;
	u_int8_t      ip_ttl;
	u_int32_t     ip_sourceip;
	u_int32_t     ip_destip;
} ip_header;

/* Enums */
enum CLK_SPEED_INDEX		{FULL_SPEED, HALF_SPEED, QUARTER_SPEED};
enum ADC_INDEX				{TEMP_FPGA, TEMP_ADC};
enum DAC_INDEX				{VB_COMP, VDD_PROT, VIN_COM, VREF_PRECH, VB_PIXBUF, VB_DS, VREF_DS, VREF_COMP };
#define DEFAULT_DAC_VALS   	{ 	1220,	/* VB_COMP */		\
								3000,	/* VDD_PROT */		\
								1053,	/* VIN_COM */		\
								1450,	/* VREF_PRECH */	\
								750,	/* VB_PIXBUF */		\
								1000,	/* VB_DS */			\
								480,	/* VREF_DS */		\
								420		/* VREF_COMP */		\
							};

#define NUM_SETTINGS		6
#define DEFAULT_SETT_INDX	{DYNAMICGAIN, DYNAMICHG0, FIXGAIN1, FIXGAIN2, FORCESWITCHG1, FORCESWITCHG2};
#define DEFAULT_SETT_VALS	{	0x0f00,		/* DYNAMICGAIN	 	*/	\
 	  	  	  	  	  	  	  	0x0f01,		/* DYNAMICHG0		*/	\
  	  	  	  	  	  	  	  	0x0f02,		/* FIXGAIN1			*/	\
  	  	  	  	  	  	  	  	0x0f06,		/* FIXGAIN2			*/	\
  	  	  	  	  	  	  	  	0x1f00,		/* FORCESWITCHG1	*/	\
 	  	  	  	  	  	  	  	0x3f00		/* FORCESWITCHG2	*/	\
 	  	  	  	  	  	  	 };


/* Hardware Definitions */
#define NMAXMOD 					(1)
#define NMOD 						(1)
#define NCHAN 						(256 * 256)
#define NCHIP 						(8)
#define NADC						(0)
#define NDAC 						(8)
#define DYNAMIC_RANGE				(16)
#define NUM_BITS_PER_PIXEL			(DYNAMIC_RANGE / 8)
#define DATA_BYTES					(NCHIP * NCHAN * NUM_BITS_PER_PIXEL)
#define IP_PACKETSIZE				(0x2052)
#define CLK_RUN						(40)	/* MHz */
#define CLK_SYNC					(20)	/* MHz */


/** Default Parameters */
#define DEFAULT_NUM_FRAMES			(100*1000*1000)
#define DEFAULT_NUM_CYCLES			(1)
#define DEFAULT_EXPTIME				(10*1000)		//ns
#define DEFAULT_PERIOD				(2*1000*1000)	//ns
#define DEFAULT_DELAY				(0)
#define DEFAULT_HIGH_VOLTAGE		(0)
#define DEFAULT_TIMING_MODE			(AUTO_TIMING)
#define DEFAULT_SETTINGS			(DYNAMICGAIN)
#define DEFAULT_TX_UDP_PORT			(0x7e9a)

/* Defines in the Firmware */
#define FIX_PATT_VAL    			(0xACDC2014)
#define ADC_PORT_INVERT_VAL   		(0x453b2a9c)


#define SAMPLE_ADC_HALF_SPEED	 	(SAMPLE_DECMT_FACTOR_2_VAL + SAMPLE_DGTL_SAMPLE_0_VAL + SAMPLE_ADC_DECMT_FACTOR_0_VAL + SAMPLE_ADC_SAMPLE_0_VAL)	/* 0x1000 */
#define SAMPLE_ADC_QUARTER_SPEED 	(SAMPLE_DECMT_FACTOR_4_VAL + SAMPLE_DGTL_SAMPLE_8_VAL + SAMPLE_ADC_DECMT_FACTOR_1_VAL + SAMPLE_ADC_SAMPLE_0_VAL)	/* 0x2810 */
#define CONFIG_HALF_SPEED			(CONFIG_TDMA_TIMESLOT_0_VAL + CONFIG_TDMA_DISABLE_VAL + CONFIG_HALF_SPEED_20MHZ_VAL + CONFIG_MODE_1_X_10GBE_VAL)
#define CONFIG_QUARTER_SPEED		(CONFIG_TDMA_TIMESLOT_0_VAL + CONFIG_TDMA_DISABLE_VAL + CONFIG_QUARTER_SPEED_10MHZ_VAL + CONFIG_MODE_1_X_10GBE_VAL)
#define ADC_OFST_HALF_SPEED_VAL		(0x20) //adc pipeline
#define ADC_OFST_QUARTER_SPEED_VAL	(0x0f)
#define ADC_PHASE_HALF_SPEED 		(0x41)
#define ADC_PHASE_QUARTER_SPEED 	(0x19)

/* Maybe not required for jungfrau */
#define NTRIMBITS 					(6)
#define NCOUNTBITS 					(24)
#define NCHIPS_PER_ADC				(2)
#define TRIM_DR 					(((int)pow(2,NTRIMBITS))-1)
#define COUNT_DR 					(((int)pow(2,NCOUNTBITS))-1)
#define ALLMOD 						(0xffff)
#define ALLFIFO 					(0xffff)

/* LTC2620 DAC DEFINES */
#define LTC2620_DAC_CMD_OFST		(20)
#define LTC2620_DAC_CMD_MSK			(0x0000000F << LTC2620_DAC_CMD_OFST)
#define LTC2620_DAC_ADDR_OFST		(16)
#define LTC2620_DAC_ADDR_MSK		(0x0000000F << LTC2620_DAC_ADDR_OFST)
#define LTC2620_DAC_DATA_OFST		(4)
#define LTC2620_DAC_DATA_MSK		(0x00000FFF << LTC2620_DAC_DATA_OFST)

#define LTC2620_DAC_CMD_WRITE		(0x00000000 << LTC2620_DAC_CMD_OFST)
#define LTC2620_DAC_CMD_SET			(0x00000003 << LTC2620_DAC_CMD_OFST)
#define LTC2620_DAC_CMD_POWER_DOWN	(0x00000004 << LTC2620_DAC_CMD_OFST)
#define LTC2620_DAC_NUMBITS			(24)


/* MAX1932 HV DEFINES */
#define MAX1932_HV_NUMBITS			(8)
#define MAX1932_HV_DATA_OFST		(0)
#define MAX1932_HV_DATA_MSK			(0x000000FF << MAX1932_HV_DATA_OFST)

/* AD9257 ADC DEFINES */
#define AD9257_ADC_NUMBITS			(24)

#define AD9257_DEV_IND_2_REG		(0x04)
#define AD9257_CHAN_H_OFST			(0)
#define AD9257_CHAN_H_MSK			(0x00000001 << AD9257_CHAN_H_OFST)
#define AD9257_CHAN_G_OFST			(1)
#define AD9257_CHAN_G_MSK			(0x00000001 << AD9257_CHAN_G_OFST)
#define AD9257_CHAN_F_OFST			(2)
#define AD9257_CHAN_F_MSK			(0x00000001 << AD9257_CHAN_F_OFST)
#define AD9257_CHAN_E_OFST			(3)
#define AD9257_CHAN_E_MSK			(0x00000001 << AD9257_CHAN_E_OFST)

#define AD9257_DEV_IND_1_REG		(0x05)
#define AD9257_CHAN_D_OFST			(0)
#define AD9257_CHAN_D_MSK			(0x00000001 << AD9257_CHAN_D_OFST)
#define AD9257_CHAN_C_OFST			(1)
#define AD9257_CHAN_C_MSK			(0x00000001 << AD9257_CHAN_C_OFST)
#define AD9257_CHAN_B_OFST			(2)
#define AD9257_CHAN_B_MSK			(0x00000001 << AD9257_CHAN_B_OFST)
#define AD9257_CHAN_A_OFST			(3)
#define AD9257_CHAN_A_MSK			(0x00000001 << AD9257_CHAN_A_OFST)
#define AD9257_CLK_CH_DCO_OFST		(4)
#define AD9257_CLK_CH_DCO_MSK		(0x00000001 << AD9257_CLK_CH_DCO_OFST)
#define AD9257_CLK_CH_IFCO_OFST		(5)
#define AD9257_CLK_CH_IFCO_MSK		(0x00000001 << AD9257_CLK_CH_IFCO_OFST)

#define AD9257_POWER_MODE_REG		(0x08)
#define AD9257_POWER_INTERNAL_OFST	(0)
#define AD9257_POWER_INTERNAL_MSK	(0x00000003 << AD9257_POWER_INTERNAL_OFST)
#define AD9257_INT_RESET_VAL		(0x3)
#define AD9257_INT_CHIP_RUN_VAL		(0x0)
#define AD9257_POWER_EXTERNAL_OFST	(5)
#define AD9257_POWER_EXTERNAL_MSK	(0x00000001 << AD9257_POWER_EXTERNAL_OFST)
#define AD9257_EXT_FULL_POWER_VAL	(0x0)
#define AD9257_EXT_STANDBY_VAL		(0x1)

#define AD9257_OUT_MODE_REG			(0x14)
#define AD9257_OUT_FORMAT_OFST		(0)
#define AD9257_OUT_FORMAT_MSK		(0x00000001 << AD9257_OUT_FORMAT_OFST)
#define AD9257_OUT_BINARY_OFST_VAL	(0)
#define AD9257_OUT_TWOS_COMPL_VAL	(1)
#define AD9257_OUT_LVDS_OPT_OFST	(6)
#define AD9257_OUT_LVDS_OPT_MSK		(0x00000001 << AD9257_OUT_LVDS_OPT_OFST)
#define AD9257_OUT_LVDS_ANSI_VAL	(0)
#define AD9257_OUT_LVDS_IEEE_VAL	(1)

#define AD9257_OUT_PHASE_REG		(0x16)
#define AD9257_OUT_CLK_OFST			(0)
#define AD9257_OUT_CLK_MSK			(0x0000000F << AD9257_OUT_CLK_OFST)
#define AD9257_OUT_CLK_60_VAL		(0x1)
#define AD9257_IN_CLK_OFST			(4)
#define AD9257_IN_CLK_MSK			(0x00000007 << AD9257_IN_CLK_OFST)
#define AD9257_IN_CLK_0_VAL			(0x0)

#define AD9257_VREF_REG				(0x18)
#define AD9257_VREF_OFST			(0)
#define AD9257_VREF_MSK				(0x00000003 << AD9257_VREF_OFST)
#define AD9257_VREF_1_33_VAL		(0x2)

#define AD9257_TEST_MODE_REG		(0x0D)
#define AD9257_OUT_TEST_OFST		(0)
#define AD9257_OUT_TEST_MSK			(0x0000000F << AD9257_OUT_TEST_OFST)
#define AD9257_NONE_VAL				(0x0)
#define AD9257_MIXED_BIT_FREQ_VAL	(0xC)
#define AD9257_TEST_RESET_SHORT_GEN	(4)
#define AD9257_TEST_RESET_LONG_GEN	(5)
#define AD9257_USER_IN_MODE_OFST	(6)
#define AD9257_USER_IN_MODE_MSK		(0x00000003 << AD9257_USER_IN_MODE_OFST)

/** PLL Reconfiguration Registers */
//https://www.altera.com/documentation/mcn1424769382940.html
#define PLL_MODE_REG 				(0x00)
#define PLL_STATUS_REG 				(0x01)
#define PLL_START_REG 				(0x02)
#define PLL_N_COUNTER_REG 			(0x03)
#define PLL_M_COUNTER_REG 			(0x04)
#define PLL_C_COUNTER_REG 			(0x05)
#define PLL_PHASE_SHIFT_REG			(0x06)

#define PLL_SHIFT_NUM_SHIFTS_OFST	(0)
#define PLL_SHIFT_NUM_SHIFTS_MSK	(0x0000FFFF << PLL_SHIFT_NUM_SHIFTS_OFST)

#define PLL_SHIFT_CNT_SELECT_OFST	(16)
#define PLL_SHIFT_CNT_SELECT_MSK	(0x0000001F << PLL_SHIFT_CNT_SELECT_OFST)
#define PLL_SHIFT_CNT_SLCT_C0_VAL	(0x0 << PLL_SHIFT_CNT_SELECT_OFST) & PLL_SHIFT_CNT_SELECT_MSK)
#define PLL_SHIFT_CNT_SLCT_C1_VAL	(0x1 << PLL_SHIFT_CNT_SELECT_OFST) & PLL_SHIFT_CNT_SELECT_MSK)
#define PLL_SHIFT_CNT_SLCT_C2_VAL	(0x2 << PLL_SHIFT_CNT_SELECT_OFST) & PLL_SHIFT_CNT_SELECT_MSK)
#define PLL_SHIFT_CNT_SLCT_C3_VAL	(0x3 << PLL_SHIFT_CNT_SELECT_OFST) & PLL_SHIFT_CNT_SELECT_MSK)
#define PLL_SHIFT_CNT_SLCT_C4_VAL	(0x4 << PLL_SHIFT_CNT_SELECT_OFST) & PLL_SHIFT_CNT_SELECT_MSK)
#define PLL_SHIFT_CNT_SLCT_C5_VAL	(0x5 << PLL_SHIFT_CNT_SELECT_OFST) & PLL_SHIFT_CNT_SELECT_MSK)
#define PLL_SHIFT_CNT_SLCT_C6_VAL	(0x6 << PLL_SHIFT_CNT_SELECT_OFST) & PLL_SHIFT_CNT_SELECT_MSK)
#define PLL_SHIFT_CNT_SLCT_C7_VAL	(0x7 << PLL_SHIFT_CNT_SELECT_OFST) & PLL_SHIFT_CNT_SELECT_MSK)
#define PLL_SHIFT_CNT_SLCT_C8_VAL	(0x8 << PLL_SHIFT_CNT_SELECT_OFST) & PLL_SHIFT_CNT_SELECT_MSK)
#define PLL_SHIFT_CNT_SLCT_C9_VAL	(0x9 << PLL_SHIFT_CNT_SELECT_OFST) & PLL_SHIFT_CNT_SELECT_MSK)
#define PLL_SHIFT_CNT_SLCT_C10_VAL	(0x10 << PLL_SHIFT_CNT_SELECT_OFST) & PLL_SHIFT_CNT_SELECT_MSK)
#define PLL_SHIFT_CNT_SLCT_C11_VAL	(0x11 << PLL_SHIFT_CNT_SELECT_OFST) & PLL_SHIFT_CNT_SELECT_MSK)
#define PLL_SHIFT_CNT_SLCT_C12_VAL	(0x12 << PLL_SHIFT_CNT_SELECT_OFST) & PLL_SHIFT_CNT_SELECT_MSK)
#define PLL_SHIFT_CNT_SLCT_C13_VAL	(0x13 << PLL_SHIFT_CNT_SELECT_OFST) & PLL_SHIFT_CNT_SELECT_MSK)
#define PLL_SHIFT_CNT_SLCT_C14_VAL	(0x14 << PLL_SHIFT_CNT_SELECT_OFST) & PLL_SHIFT_CNT_SELECT_MSK)
#define PLL_SHIFT_CNT_SLCT_C15_VAL	(0x15 << PLL_SHIFT_CNT_SELECT_OFST) & PLL_SHIFT_CNT_SELECT_MSK)
#define PLL_SHIFT_CNT_SLCT_C16_VAL	(0x16 << PLL_SHIFT_CNT_SELECT_OFST) & PLL_SHIFT_CNT_SELECT_MSK)
#define PLL_SHIFT_CNT_SLCT_C17_VAL	(0x17 << PLL_SHIFT_CNT_SELECT_OFST) & PLL_SHIFT_CNT_SELECT_MSK)

#define PLL_SHIFT_UP_DOWN_OFST		(21)
#define PLL_SHIFT_UP_DOWN_MSK		(0x00000001 << PLL_SHIFT_UP_DOWN_OFST)
#define PLL_SHIFT_UP_DOWN_NEG_VAL	((0x0 << PLL_SHIFT_UP_DOWN_OFST) & PLL_SHIFT_UP_DOWN_MSK)
#define PLL_SHIFT_UP_DOWN_POS_VAL	((0x1 << PLL_SHIFT_UP_DOWN_OFST) & PLL_SHIFT_UP_DOWN_MSK)

#define PLL_K_COUNTER_REG 			(0x07)
#define PLL_BANDWIDTH_REG 			(0x08)
#define PLL_CHARGEPUMP_REG 			(0x09)
#define PLL_VCO_DIV_REG 			(0x1c)
#define PLL_MIF_REG 				(0x1f)


#endif /* SLSDETECTORSERVER_DEFS_H */
