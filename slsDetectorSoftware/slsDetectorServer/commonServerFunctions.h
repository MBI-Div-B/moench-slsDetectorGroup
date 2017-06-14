#ifndef COMMON_SERVER_FUNCTIONS_H
#define COMMON_SERVER_FUNCTIONS_H

#include "blackfin.h"

/* global variables */
void serializeToSPI(u_int32_t addr, u_int32_t val, u_int16_t csmask, int numbitstosend, u_int16_t clkmask, u_int16_t digoutmask, int digofset) {
#ifdef VERBOSE
	if (numbitstosend == 16)
		printf("Writing to SPI Register: 0x%04x\n",val);
	else
		printf("Writing to SPI Register: 0x%08x\n", val);
#endif

	u_int16_t valw;

	// start point
	valw = 0xffff; 		/**todo testwith old board 0xff for adc_spi */			// old board compatibility (not using specific bits)
	bus_w16 (addr, valw);

	// chip sel bar down
	valw &= ~csmask; /* todo with test: done a bit different, not with previous value */
	bus_w16 (addr, valw);

	{
		int i = 0;
		for (i = 0; i < numbitstosend; ++i) {

			// clk down
			valw &= ~clkmask;
			bus_w16 (addr, valw);

			// write data (i)
			valw = ((valw & ~digoutmask) + 										// unset bit
					(((val >> (numbitstosend - 1 - i)) & 0x1) << digofset)); 	// each bit from val starting from msb
			bus_w16 (addr, valw);

			// clk up
			valw |= clkmask ;
			bus_w16 (addr, valw);
		}
	}

	// chip sel bar up
	valw |= csmask; /* todo with test: not done for spi */
	bus_w16 (addr, valw);

	//clk down
	valw &= ~clkmask;
	bus_w16 (addr, valw);

	// stop point = start point of course
	valw = 0xffff; 		/**todo testwith old board 0xff for adc_spi */			// old board compatibility (not using specific bits)
	bus_w16 (addr, valw);
}


#endif	//COMMON_SERVER_FUNCTIONS_H
