#ifndef EIGERMODULEDATA_H
#define  EIGERMODULEDATA_H
#include "slsReceiverData.h"



class eigerHalfModuleData : public slsReceiverData<uint32_t> {
public:




	/**
     Implements the slsReceiverData structure for the eiger prototype read out by a half module i.e. using the slsReceiver
     (256*256 pixels, 512 packets for 16 bit mode, 256 for 8, 128 for 4, 1024 for 32, 1040  etc.)
     \param d dynamic range
     \param c crosstalk parameter for the output buffer

	 */


	eigerHalfModuleData(bool t, bool l, int dr, int tg, int psize, int dsize, int npf, int x, int y, double c=0):
		slsReceiverData<uint32_t>(x, y, npf, psize),
		top(t), left(l),
		dynamicRange(dr), tenGiga(tg),
		packetSize(psize), dataSize(dsize), numberofPacketsPerFrame(npf),
		xpixels(x),ypixels(y),
		xtalk(c),
		header_t(0), footer_t(0){


		int **dMap;
		uint32_t **dMask;

		dMap=new int*[ypixels];
		dMask=new uint32_t*[ypixels];


		for (int i = 0; i < ypixels; i++) {
			dMap[i] = new int[xpixels];
			dMask[i] = new uint32_t[xpixels];
		}

		//Map
		int totalNumberOfBytes = numberofPacketsPerFrame * packetSize;
		int iPacket = 8;
		int iData = 0;
		int increment = (dynamicRange/8);
		int ic_increment = 1;
		if (dynamicRange == 4) {
			increment = 1;
			ic_increment = 2;
		}

		if(top){
			for (int ir=0; ir<ypixels; ir++) {
				for (int ic=0; ic<xpixels; ic = ic + ic_increment) {
					dMap[ir][ic]  = iPacket;
					iPacket += increment;
					iData += increment;
					//increment header
					if(iData >= dataSize){
						iPacket += 16;
						iData = 0;
					}

				}
			}
		}

		//bottom
		else{
			iData = 0;
			int numbytesperline;
			switch(dynamicRange){
				case 4: numbytesperline = 256; break;
				case 8:	numbytesperline = 512; break;
				case 16:numbytesperline = 1024; break;
				case 32:numbytesperline = 2048; break;
			}
			iPacket = totalNumberOfBytes - numbytesperline - 8;
			if((dynamicRange == 32) && (!tenGiga))
					iPacket -= 16;

			for (int ir=0; ir<ypixels; ir++) {
				for (int ic=0; ic<xpixels; ic = ic + ic_increment) {
					dMap[ir][ic]  = iPacket;
					iPacket += increment;
					iData += increment;
					//--------------------32 bit 1giga -------------------
					if((dynamicRange == 32) && (!tenGiga)){
						if(iData == numbytesperline){
							iPacket -= (numbytesperline*2 + 16*3);
							iData = 0;
						}
						if(iData == dataSize){
							iPacket += 16;
						}
					}//------------end of 32 bit -------------------------
					else if((iData % numbytesperline) == 0){
						iPacket -= (numbytesperline*2);
						if(iData == dataSize){
							iPacket -= 16;
							iData = 0;
						}
					}
					//---------------------------------------------------
				}
			}
		}




		//Mask
		for(int ir=0; ir<ypixels; ++ir)
			for(int ic=0; ic<xpixels; ++ic)
				dMask[ir][ic] = 0x0;

		setDataMap(dMap);
		setDataMask(dMask);

	};





	/** Returns the frame number for the given dataset.
	     \param buff pointer to the dataset
	     \returns frame number
	 */
	int getFrameNumber(char *buff){
		footer_t = (eiger_packet_footer_t*)(buff + dataSize + sizeof(eiger_packet_header_t));
		return ((uint32_t)(*( (uint64_t*) footer_t)));
	};





	/** gets the packets number
     	 \param buff pointer to the memory
     	 \returns packet number
	 */
	int getPacketNumber(char *buff){
		footer_t = (eiger_packet_footer_t*)(buff + dataSize + sizeof(eiger_packet_header_t));
		return(*( (uint16_t*) footer_t->packetnum));
	};



	/**
    returns the pixel value as double correcting for the output buffer crosstalk
     \param data pointer to the memory
     \param ix coordinate in the x direction
     \param iy coordinate in the y direction
     \returns channel value as double

	 */
	double getValue(char *data, int ix, int iy=0) {
		//  cout << "##" << (void*)data << " " << ix << " " <<iy << endl;
		if (xtalk==0)
			return getChannelwithMissingPackets(data, ix, iy);
		else
			return getChannelwithMissingPackets(data, ix, iy)-xtalk * getChannelwithMissingPackets(data, ix-1, iy);
	};


	/**

	     Returns the value of the selected channel for the given dataset. Virtual function, can be overloaded.
	     \param data pointer to the dataset (including headers etc)
	     \param ix pixel number in the x direction
	     \param iy pixel number in the y direction
	     \returns data for the selected channel, with inversion if required

	 */

	virtual int getChannelwithMissingPackets(char *data, int ix, int iy) {
		uint32_t m=0, n = 0;
		int linesperpacket,newix, newiy,origX;


		//cout <<"ix:"<<ix<<" nx:"<<nx<<" iy:"<<iy<<" ny:"<<ny<<" datamap[iy][ix]:"<< dataMap[iy][ix] <<" datasize:"<< dataSize <<endl;
		if (ix>=0 && ix<nx && iy>=0 && iy<ny && dataMap[iy][ix]>=0 && dataMap[iy][ix]<dataSize) {
			m=dataMask[iy][ix];


			//pixelpos1d = (nx * iy + ix);

			switch(dynamicRange){
			case 4: 	if(tenGiga) linesperpacket=16;	else linesperpacket=4;break;
			case 8: 	if(tenGiga) linesperpacket=8; 	else linesperpacket=2;break;
			case 16: 	if(tenGiga) linesperpacket=4; 	else linesperpacket=1;break;
			case 32:	if(tenGiga) linesperpacket=2; 	else linesperpacket=1;break;
			}



			//each byte is shared by 2 pixels for 4 bit mode
			origX = ix;
			if((dynamicRange == 4) && (ix%2))
				ix--;



			// ------check if missing packet, get to pixel at start of packet-----------------

			//to get the starting of a packet (except 1g 32 bit)
			newix = 0;
			// 0.5 Lines per packet for 1g 32 bit
			if(dynamicRange == 32 && !tenGiga)
				newix = ix - (ix%256);

			//iy divided by linesperpacket depending on bitmode
			if(!(iy%linesperpacket))
				newiy = iy;
			else
				newiy = (iy - (iy%linesperpacket));

			header_t = (eiger_packet_header_t*)((char*)(data +(dataMap[newiy][newix]-8)));
			uint16_t identifier = (uint16_t)*( (uint16_t*) header_t->missingpacket);

			if(identifier==missingPacketValue){
			  //	cprintf(RED,"missing packet\n");
			  return -1;
			}
			else if(identifier==deactivatedPacketValue){
				  //	cprintf(RED,"deactivated packet\n");
				  return -2;
				}
			// -----END OF CHECK -------------------------------------------------------------


		}else{
			cprintf(RED,"outside limits\n");
			return -99;
		}

		//get proper data
		n = ((uint32_t)(*((uint32_t*)(((char*)data)+(dataMap[iy][ix])))));

		//each byte is shared by 2 pixels for 4 bit mode
		if(dynamicRange == 4){
			if(ix != origX)
				return ((n & 0xf0)>>4)^m;
			return (n & 0xf)^m;
		}
		else if(dynamicRange == 8)	return (n & 0xff)^m;
		else if(dynamicRange == 16)	return (n & 0xffff)^m;
		else				return (n & 0xffffffff)^m;


	};


	/** sets the output buffer crosstalk correction parameter
      \param c output buffer crosstalk correction parameter to be set
      \returns current value for the output buffer crosstalk correction parameter

	 */
	double setXTalk(double c) {xtalk=c; return xtalk;}


	/** gets the output buffer crosstalk parameter
      \returns current value for the output buffer crosstalk correction parameter
	 */
	double getXTalk() {return xtalk;}


	//have to make the frame here itself as 1 frame will have less packets than trying to read

	char *readNextFrame(ifstream &filebin, int& fnum) {

		char *data=new char[packetSize*numberofPacketsPerFrame];
		char *retval=0;
		int np=0, nd;
		fnum = -1;
		int pnum=-1;


		if (filebin.is_open()) {

			while (filebin.read(data+np*packetSize,packetSize)) {

				fnum = getFrameNumber(data); cout << "fnum:"<<fnum<<endl;
				pnum = getPacketNumber(data); cout << "pnum:"<<pnum<<endl;

				if (np==(numberofPacketsPerFrame-1)) {cout<<"last packet"<<endl;

					fnum=getFrameNumber(data); //cout << "fnum:"<<fnum<<endl;
					retval=findNextFrame(data,nd,packetSize*numberofPacketsPerFrame);
					np=nd/packetSize;
					cout << "got this mnany packets:"<<np << endl;


					if (retval==data && np==numberofPacketsPerFrame) {
						//      cout << "-" << endl;
						return data;

					} else if (np>numberofPacketsPerFrame) {
						cout << "too many packets!!!!!!!!!!" << endl;
						delete [] data;
						return NULL;
					} else if (retval!=NULL) {
						//  cout << "+" << endl;;
						for (int ip=0; ip<np; ip++)
							memcpy(data+ip*packetSize,retval+ip*packetSize,packetSize);
					}

				} else if (np>numberofPacketsPerFrame) {
					cout << "*******too many packets!!!!!!!!!!" << endl;
					delete [] data;
					return NULL;
				} else {
					//  cout << "." << endl;;
					np++;
					cout<<"np:"<<np<<endl;
				}
			}
		}
		delete [] data;
		return NULL;
	}




private:
	/** Missing Packet identifier value */
	const static uint16_t missingPacketValue = 0xFFFF;
	const static uint16_t deactivatedPacketValue = 0xFEFE;


	const bool top;
	const bool left;
	const int dynamicRange;
	const bool tenGiga;
	const int packetSize;
	const int dataSize;
	const int numberofPacketsPerFrame;
	const int xpixels;
	const int ypixels;
	double xtalk; /**<output buffer crosstalk correction parameter */


	/** structure of an eiger packet*/
	typedef struct{
		unsigned char subframenum[4];
		unsigned char missingpacket[2];
		unsigned char portnum[1];
		unsigned char dynamicrange[1];
	} eiger_packet_header_t;

	typedef struct{
		unsigned char framenum[6];
		unsigned char packetnum[2];
	} eiger_packet_footer_t;

	eiger_packet_header_t* header_t;
	eiger_packet_footer_t* footer_t;
};




#endif
