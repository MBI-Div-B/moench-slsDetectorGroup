#include "energyConversion.h"

#include <fstream>
#include <sstream>
#include <iostream>
#include <cstring>

#include "fileIOStatic.h"

using namespace std;




int energyConversion::readCalibrationFile(string fname, double &gain, double &offset){

  string str;
  ifstream infile;
#ifdef VERBOSE
  std::cout<< "Opening file "<< fname << std::endl;
#endif
  infile.open(fname.c_str(), ios_base::in);
  if (infile.is_open()) {
    getline(infile,str);
#ifdef VERBOSE
    std::cout<< str << std::endl;
#endif
    istringstream ssstr(str);
    ssstr >> offset >> gain;
    infile.close();
	cout << "Calibration file loaded: " << fname << endl;
  } else {
    std::cout<< "Could not open calibration file "<< fname << std::endl;
    gain=0.;
    offset=0.;
#ifndef MYROOT
    return FAIL;
#endif
    return -1;
  }
#ifndef MYROOT
  return OK;
#endif
  return 0;
};

int energyConversion::writeCalibrationFile(string fname, double gain, double offset){
  //std::cout<< "Function not yet implemented " << std::endl;
  ofstream outfile;

  outfile.open (fname.c_str());

  // >> i/o operations here <<
  if (outfile.is_open()) {
    outfile << offset << " " << gain << std::endl;
  } else {
    std::cout<< "Could not open calibration file "<< fname << " for writing" << std::endl;
#ifndef MYROOT
    return FAIL;
#endif
    return -1;
  }

  outfile.close();

#ifndef MYROOT
  return OK;
#endif
  return 0;
};


#ifndef MYROOT

/* I/O */


slsDetectorDefs::sls_detector_module* energyConversion::readSettingsFile(string fname,  detectorType myDetectorType, sls_detector_module *myMod){
	int nflag=0;


	if (myMod==NULL) {
		myMod=createModule(myDetectorType);
		nflag=1;
	}

	int id=0,i;
	string names[100];
	string myfname;
	string str;
	ifstream infile;
	ostringstream oss;
	int iline=0;
	string sargname;
	int ival;
	int ichan=0, ichip=0, idac=0;
	int nch=((myMod->nchan)/(myMod->nchip));

	//ascii settings/trim file
	switch (myDetectorType) {
	case MYTHEN:
		break;
	case MOENCH:
		names[id++]="Vdac0";
		names[id++]="Vdac1";
		names[id++]="Vdac2";
		names[id++]="Vdac3";
		names[id++]="Vdac4";
		names[id++]="Vdac5";
		names[id++]="Vdac6";
		names[id++]="Vdac7";
		break;
	case GOTTHARD:
		names[id++]="Vref";
		names[id++]="VcascN";
		names[id++]="VcascP";
		names[id++]="Vout";
		names[id++]="Vcasc";
		names[id++]="Vin";
		names[id++]="Vref_comp";
		names[id++]="Vib_test";
		break;
	case EIGER:
		break;
	default:
		cout << "Unknown detector type - unknown format for settings file" << endl;
		return NULL;
	}

#ifdef VERBOSE
	std::cout<<   "reading settings file for module number "<< myMod->module << std::endl;
#endif
	myfname=fname;
#ifdef VERBOSE
	std::cout<< "file name is "<< myfname <<   std::endl;
#endif

	switch (myDetectorType) {

	case MYTHEN:
		infile.open(myfname.c_str(), ios_base::in);
		if (infile.is_open()) {
			for (int iarg=0; iarg<myMod->ndac; iarg++) {
				getline(infile,str);
				iline++;
				istringstream ssstr(str);
				ssstr >> sargname >> ival;
#ifdef VERBOSE
				std::cout<< sargname << " dac nr. " << idac << " is " << ival << std::endl;
#endif
				myMod->dacs[idac]=ival;
				idac++;
			}
			for (ichip=0; ichip<myMod->nchip; ichip++) {
				getline(infile,str);
				iline++;
#ifdef VERYVERBOSE
				std::cout<< str << std::endl;
#endif
				istringstream ssstr(str);
				ssstr >> sargname >> ival;
#ifdef VERYVERBOSE
				std::cout<< "chip " << ichip << " " << sargname << " is " << ival << std::endl;
#endif

				myMod->chipregs[ichip]=ival;
				for (ichan=0; ichan<nch; ichan++) {
					getline(infile,str);
#ifdef VERYVERBOSE
					std::cout<< str << std::endl;
#endif
					istringstream ssstr(str);

#ifdef VERYVERBOSE
					std::cout<< "channel " << ichan+ichip*thisDetector->nChans <<" iline " << iline<< std::endl;
#endif
					iline++;
					myMod->chanregs[ichip*nch+ichan]=0;
					for (int iarg=0; iarg<6 ; iarg++) {
						ssstr >>  ival;
						//if (ssstr.good()) {
						switch (iarg) {
						case 0:
#ifdef VERYVERBOSE
							std::cout<< "trimbits " << ival ;
#endif
							myMod->chanregs[ichip*nch+ichan]|=ival&TRIMBITMASK;
							break;
						case 1:
#ifdef VERYVERBOSE
							std::cout<< " compen " << ival ;
#endif
							myMod->chanregs[ichip*nch+ichan]|=ival<<9;
							break;
						case 2:
#ifdef VERYVERBOSE
							std::cout<< " anen " << ival ;
#endif
							myMod->chanregs[ichip*nch+ichan]|=ival<<8;
							break;
						case 3:
#ifdef VERYVERBOSE
							std::cout<< " calen " << ival  ;
#endif
							myMod->chanregs[ichip*nch+ichan]|=ival<<7;
							break;
						case 4:
#ifdef VERBOSE
							std::cout<< " outcomp " << ival  ;
#endif
							myMod->chanregs[ichip*nch+ichan]|=ival<<10;
							break;
						case 5:
#ifdef VERBOSE
							std::cout<< " counts " << ival  << std::endl;
#endif
							myMod->chanregs[ichip*nch+ichan]|=ival<<11;
							break;
						default:
							std::cout<< " too many columns" << std::endl;
							break;
						}
					}
				}
				//	}
			}
#ifdef VERBOSE
			std::cout<< "read " << ichan*ichip << " channels" <<std::endl;
#endif

			infile.close();
			strcpy(settingsFile,fname.c_str());
			cout << "Settings file loaded: " << settingsFile << endl;
			return myMod;

		}


		break;

	case EIGER:
		infile.open(myfname.c_str(),ifstream::binary);
		if (infile.is_open()) {
			infile.read((char*) myMod->dacs,sizeof(int)*(myMod->ndac));
			infile.read((char*) myMod->chanregs,sizeof(int)*(myMod->nchan));
#ifdef VERBOSE
			for(int i=0;i<myMod->ndac;i++)
				std::cout << "dac " << i << ":" << myMod->dacs[i] << std::endl;
#endif
			if(infile.eof()){
				cout<<endl<<"Error, could not load trimbits end of file, "<<myfname<<", reached."<<endl<<endl;
				if (nflag)
					deleteModule(myMod);

				return NULL;
			}
			infile.close();
			strcpy(settingsFile,fname.c_str());
			cout << "Settings file loaded: " << settingsFile << endl;
			return myMod;

		}

		break;

	case MOENCH:
	case GOTTHARD:
		//---------------dacs---------------
		infile.open(myfname.c_str(), ios_base::in);
		if (infile.is_open()) {
			while(infile.good()) {
				getline(infile,str);
				iline++;
#ifdef VERBOSE
				std::cout<< str << std::endl;
#endif
				istringstream ssstr(str);
				ssstr >> sargname >> ival;
				for (i=0;i<id;i++){
					if (!strcasecmp(sargname.c_str(),names[i].c_str())){
						myMod->dacs[i]=ival;
						idac++;
#ifdef VERBOSE
						std::cout<< sargname << " dac nr. " << idac << " is " << ival << std::endl;
#endif
						break;
					}
				}
			}
			if (i < id) {
#ifdef VERBOSE
				std::cout<< sargname << " dac nr. " << idac << " is " << ival << std::endl;
#endif
			}else
				std::cout<< "Unknown dac " << sargname << std::endl;

			infile.close();
			strcpy(settingsFile,fname.c_str());
			cout << "Settings file loaded: " << settingsFile << endl;
			return myMod;

		}

		//----------------------------------
		break;

	default:
		std::cout<< "Unknown detector type - don't know how to read file" <<  myfname << std::endl;
		infile.close();
		deleteModule(myMod);
		return NULL;

	}

	std::cout<< "Error: Could not open settings file " <<  myfname << std::endl;
	if (nflag)
		deleteModule(myMod);

	return NULL;



};


int energyConversion::writeSettingsFile(string fname, detectorType myDetectorType, sls_detector_module mod){

	ofstream outfile;

	int nch=((mod.nchan)/(mod.nchip));

	string names[100];
	int id=0;
	switch (myDetectorType) {
	case MYTHEN:
		names[id++]="Vtrim";
		names[id++]="Vthresh";
		names[id++]="Rgsh1";
		names[id++]="Rgsh2";
		names[id++]="Rgpr";
		names[id++]="Vcal";
		names[id++]="outBuffEnable";
		break;
	case MOENCH:
		names[id++]="Vdac0";
		names[id++]="Vdac1";
		names[id++]="Vdac2";
		names[id++]="Vdac3";
		names[id++]="Vdac4";
		names[id++]="Vdac5";
		names[id++]="Vdac6";
		names[id++]="Vdac7";
		break;
	case GOTTHARD:
		names[id++]="Vref";
		names[id++]="VcascN";
		names[id++]="VcascP";
		names[id++]="Vout";
		names[id++]="Vcasc";
		names[id++]="Vin";
		names[id++]="Vref_comp";
		names[id++]="Vib_test";
		break;
	case EIGER:
		break;
	default:
		cout << "Unknown detector type - unknown format for settings file" << endl;
		return FAIL;
	}

	int iv, ichan, ichip;
	int iv1, idac;
	int nb;

	switch (myDetectorType) {
	case EIGER:
		outfile.open(fname.c_str(), ofstream::binary);
		if (outfile.is_open()) {
			iv = 1150;
			outfile.write((char*)mod.dacs, sizeof(int)*(mod.ndac));
			outfile.write((char*)mod.chanregs, sizeof(int)*(mod.nchan));
			outfile.close();
			return slsDetectorDefs::OK;
		}

		std::cout<< "could not open SETTINGS file " << fname << std::endl;
		return slsDetectorDefs::FAIL;
	default:


		outfile.open(fname.c_str(), ios_base::out);

		if (outfile.is_open()) {
			for (idac=0; idac<mod.ndac; idac++) {
				iv=(int)mod.dacs[idac];
				outfile << names[idac] << " " << iv << std::endl;
			}

			if(myDetectorType==MYTHEN){
				for (ichip=0; ichip<mod.nchip; ichip++) {
					iv1=mod.chipregs[ichip]&1;
					outfile << names[idac] << " " << iv1 << std::endl;
					for (ichan=0; ichan<nch; ichan++) {
						iv=mod.chanregs[ichip*nch+ichan];
						iv1= (iv&TRIMBITMASK);
						outfile <<iv1 << " ";
						nb=9;
						iv1=((iv&(1<<nb))>>nb);
						outfile << iv1 << " ";
						nb=8;
						iv1=((iv&(1<<nb))>>nb);
						outfile << iv1 << " ";
						nb=7;
						iv1=((iv&(1<<nb))>>nb);
						outfile <<iv1  << " ";
						nb=10;
						iv1=((iv&(1<<nb))>>nb);
						outfile << iv1 << " ";
						nb=11;
						iv1= ((iv&0xfffff800)>>nb);
						outfile << iv1  << std::endl;
					}
				}
			}
			outfile.close();
			return OK;
		}
		std::cout<< "could not open SETTINGS file " << fname << std::endl;
		return FAIL;

	}


};



#endif
