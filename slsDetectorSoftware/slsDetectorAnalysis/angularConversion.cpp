#include "angularConversion.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h>
#include "usersFunctions.h"

using namespace std;

angularConversion::angularConversion():  currentPosition(0), 
					  currentPositionIndex(0) 
{
 

}

angularConversion::~angularConversion(){

}

// int angularConversion::setAngularConversionPointer(angleConversionConstant *p, int *nm, int nch, int idet) {
//   if (p) {
//     angOff[idet]=p;
//     nMods[idet]=nm;
//     nCh[idet]=nch;
//   } else {
//     angOff[idet]=NULL;
//     nMods[idet]=NULL;
//   } 
    
//   return OK;
// }




float* angularConversion::convertAngles(float pos) {
  int imod;
  float    *ang=new float[getTotalNumberOfChannels()]; 
  float enc=pos;
  angleConversionConstant *p=NULL;
  for (int ip=0; ip<getTotalNumberOfChannels(); ip++) {
    imod=ip/(getChansPerMod(0)); // always for module 0?????
    p=getAngularConversionPointer(imod);

    if (getMoveFlag(imod)==0)
      enc=0;
    else
      enc=pos;

    if (p)
      ang[ip]=angle(ip%(getChansPerMod()),		\
		    enc,				\
		    (*fineOffset)+(*globalOffset),	\
		    p->r_conversion,			\
		    p->center,				\
		    p->offset,				\
		    p->tilt,				\
		    *angDirection		  );
  }
  return ang;
}



//static!
int angularConversion::readAngularConversion(string fname, int nmod, angleConversionConstant *angOff) {

  ifstream infile;
  string ss;

#ifdef VERBOSE
  std::cout<< "Opening file "<< fname << std::endl;
#endif
  infile.open(fname.c_str(), ios_base::in);
  if (infile.is_open()) {
    readAngularConversion(infile, nmod, angOff);
    infile.close();
  } else {
    std::cout<< "Could not open calibration file "<< fname << std::endl;
    return -1;
  }
  return 0;
}


//static
int angularConversion::readAngularConversion( ifstream& infile, int nmod, angleConversionConstant *angOff) {
  string str;
  int mod;
  float center, ecenter;
  float r_conv, er_conv;
  float off, eoff;
  string ss;
  int interrupt=0;
  int nm=0;
  //" module %i center %E +- %E conversion %E +- %E offset %f +- %f \n"
  while (infile.good() and interrupt==0) {
    getline(infile,str);
#ifdef VERBOSE
    cout << "** mod " << nm << " " ;
    std::cout<< str << std::endl;
#endif
    istringstream ssstr(str);
    ssstr >> ss >> mod;
    ssstr >> ss >> center;
    ssstr >> ss >> ecenter;
    ssstr >> ss >> r_conv;
    ssstr >> ss >> er_conv;
    ssstr >> ss >> off;
    ssstr >> ss >> eoff;
    if (nm<nmod && nm>=0 ) {
      angOff[nm].center=center;
      angOff[nm].r_conversion=r_conv;
      angOff[nm].offset=off;
      angOff[nm].ecenter=ecenter;
      angOff[nm].er_conversion=er_conv;
      angOff[nm].eoffset=eoff;
    } else
      break;
    //cout << nm<<"  " << angOff[nm].offset << endl;
    nm++;
    if (nm>=nmod)
      break;

    


  }
  return nm;
 }

//static
int angularConversion:: writeAngularConversion(string fname, int nmod, angleConversionConstant *angOff) {

  ofstream outfile;
  outfile.open (fname.c_str(),ios_base::out);
  if (outfile.is_open())
  {  
    writeAngularConversion(outfile, nmod, angOff);
    outfile.close();
  } else {
    std::cout<< "Could not open file " << fname << "for writing"<< std::endl;
    return -1;
  }
  //" module %i center %E +- %E conversion %E +- %E offset %f +- %f \n"
  return 0;
}



//static 
int angularConversion:: writeAngularConversion(ofstream& outfile, int nmod, angleConversionConstant *angOff) {

  for (int imod=0; imod<nmod; imod++) {
      outfile << " module " << imod << " center "<< angOff[imod].center<<"  +- "<< angOff[imod].ecenter<<" conversion "<< angOff[imod].r_conversion << " +- "<< angOff[imod].er_conversion <<  " offset "<< angOff[imod].offset << " +- "<< angOff[imod].eoffset << std::endl;
  }
  return 0;
}


//static
int angularConversion::resetMerging(float *mp, float *mv, float *me, int *mm, int nb) {
  
		
#ifdef VERBOSE
  cout << "creating merging arrays "<<  nb << endl;
#endif
  

  for (int ibin=0; ibin<nb; ibin++) {
    mp[ibin]=0;
    mv[ibin]=0;
    me[ibin]=0;
    mm[ibin]=0;
  }
  return OK;
}


//static
int angularConversion::finalizeMerging(float *mp, float *mv, float *me, int *mm,int nb) {
   int np=0;
   for (int ibin=0; ibin<nb; ibin++) {
     if (mm[ibin]>0) {
       mp[np]=mp[ibin]/mm[ibin];
       mv[np]=mv[ibin]/mm[ibin];
      me[np]=me[ibin]/mm[ibin];
      me[np]=sqrt(me[ibin]);
      mm[np]=mm[ibin];
      np++;
    }
  }
  return np;
}

//static
int  angularConversion::addToMerging(float *p1, float *v1, float *e1, float *mp, float *mv,float *me, int *mm, int nchans, float binsize,int nbins, int *badChanMask ) {

  float binmi=-180.;
  int ibin=0;

  if (p1==NULL)
    return 0;
  if (v1==NULL)
    return FAIL;

  if (mp==NULL) //can be changed if we want to use a fixed bin algorithm!
    return FAIL;

  if (mv==NULL)
    return FAIL;
  if (me==NULL)
    return FAIL;
  if (mm==NULL)
    return FAIL;
  if (nchans==0)
    return FAIL;
  
  if (binsize==0)
    return FAIL;
  if (nbins==0)
    return FAIL;
  

  for (int ip=0; ip<nchans; ip++) {
    if (badChanMask) {
      if (badChanMask[ip]) {
#ifdef VERBOSE
	cout << "channel " << ip << " is bad " << endl;
#endif
	continue;
      }
    }
    ibin=(int)((p1[ip]-binmi)/binsize);
   
    if (ibin<nbins && ibin>=0) {
      mp[ibin]+=p1[ip];
      mv[ibin]+=v1[ip];
      if (e1)
	me[ibin]+=(e1[ip]*e1[ip]);
      else
	me[ibin]+=v1[ip];
      mm[ibin]++;
    } else
      return FAIL;
  }
  

  return OK;
  
}

int angularConversion::deleteMerging() {

  if (mergingBins)
    delete [] mergingBins;

  if (mergingCounts)
    delete [] mergingCounts;

  if (mergingErrors)
    delete [] mergingErrors;

}


int angularConversion::resetMerging() {

  mergingBins=new float[nBins];
  

  mergingCounts=new float[nBins];


  mergingErrors=new float[nBins];

  
  mergingMultiplicity=new int[nBins];

  return resetMerging(mergingBins, mergingCounts, mergingErrors, mergingMultiplicity);

}

int angularConversion::resetMerging(float *mp, float *mv, float *me, int *mm) {
  if (nBins)
    return resetMerging(mp, mv, me, mm,nBins);
  else 
    return FAIL;
}




int angularConversion::finalizeMerging() {
  int np=finalizeMerging(mergingBins, mergingCounts, mergingErrors, mergingMultiplicity);

  if (mergingMultiplicity)
    delete [] mergingMultiplicity;

  return np;

}




int angularConversion::finalizeMerging(float *mp, float *mv, float *me, int *mm) {
 
  if (nBins)
    return finalizeMerging(mp, mv, me, mm, nBins);
  else
    return FAIL;
}

int  angularConversion::addToMerging(float *p1, float *v1, float *e1, int *badChanMask ) {
  
  return addToMerging(p1,v1,e1,mergingBins,mergingCounts, mergingErrors, mergingMultiplicity, badChanMask);


}


int  angularConversion::addToMerging(float *p1, float *v1, float *e1, float *mp, float *mv,float *me, int *mm, int *badChanMask ) {


  int del=0;

  if (*binSize==0)
    return FAIL;
  
  if (nBins==0)
    return FAIL;
  
  if (p1==NULL) {
    del=1;
    p1=convertAngles();
  }
  
  int ret=addToMerging(p1, v1, e1, mp, mv,me, mm,getTotalNumberOfChannels(), *binSize,nBins, badChanMask );
  
  
  if (del) {
    delete [] p1;
    p1=NULL;
  }
  return ret;
}



   /**
     sets the value of s angular conversion parameter
     \param c can be ANGULAR_DIRECTION, GLOBAL_OFFSET, FINE_OFFSET, BIN_SIZE
     \param v the value to be set
     \returns the actual value
  */

float angularConversion::setAngularConversionParameter(angleConversionParameter c, float v){
  

  switch (c) {
  case ANGULAR_DIRECTION:
    if (v<0)
      *angDirection=-1;
    else
      *angDirection=1;
    return *angDirection;
  case GLOBAL_OFFSET:
    *globalOffset=v;
    return *globalOffset;
  case FINE_OFFSET:
    *fineOffset=v;
    return *fineOffset;
  case BIN_SIZE:
    if (v>0) {
      *binSize=v;
      nBins=360./(*binSize);
    }
    return *binSize;
  case MOVE_FLAG:
    if (moveFlag) {
      if (v>0)
	*moveFlag=1;
      else if (v==0)
	*moveFlag=0;
      return *moveFlag;
    }
    return -1;
  default:
    return 0;
  }
}

  /**
     returns the value of an angular conversion parameter
     \param c can be ANGULAR_DIRECTION, GLOBAL_OFFSET, FINE_OFFSET, BIN_SIZE
     \returns the actual value

  */

float angularConversion::getAngularConversionParameter(angleConversionParameter c) {

  switch (c) {
  case ANGULAR_DIRECTION:
    return *angDirection;
  case GLOBAL_OFFSET:
    return *globalOffset;
  case FINE_OFFSET:
    return *fineOffset;
  case BIN_SIZE:
    return *binSize;
  case MOVE_FLAG:
    return *moveFlag;
  default:
    return 0;
  }
}




int angularConversion::setAngularConversionFile(string fname) {
  if (fname=="") {
    setAngularCorrectionMask(0);
#ifdef VERBOSE
    std::cout << "Unsetting angular conversion" <<  std::endl;
#endif
  } else {
    if (fname=="default") {
      fname=string(angConvFile);
    }
    
#ifdef VERBOSE
    std::cout << "Setting angular conversion to " << fname << std:: endl;
#endif
    if (readAngularConversionFile(fname)>=0) {
      setAngularCorrectionMask(1);
      strcpy(angConvFile,fname.c_str());
    }
  }
  return setAngularCorrectionMask();
}





  /*
      set  positions for the acquisition
      \param nPos number of positions
      \param pos array with the encoder positions
      \returns number of positions
  */
int angularConversion::setPositions(int nPos, float *pos){
  if (nPos>=0)
    *numberOfPositions=nPos; 
  for (int ip=0; ip<nPos; ip++) 
    detPositions[ip]=pos[ip]; 

  return *numberOfPositions;
}
/* 
   get  positions for the acquisition
   \param pos array which will contain the encoder positions
   \returns number of positions
*/
int angularConversion::getPositions(float *pos){ 
  if (pos) {
    for (int ip=0; ip<(*numberOfPositions); ip++) 
      pos[ip]=detPositions[ip];
  } 

  return  *numberOfPositions;
}
  




