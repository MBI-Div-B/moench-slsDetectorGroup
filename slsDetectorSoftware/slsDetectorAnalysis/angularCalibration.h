
#ifndef ANGULARCALIBRATION_H
#define ANGULARCALIBRATION_H


//#include "usersFunctions.h"


#ifdef ROOT
#include <TROOT.h>
#include <TF1.h>
class TH1;
#endif

  //float angle(int ichan, float encoder, float totalOffset, float conv_r, float center, float offset, float tilt, int direction)



/**
   angular conversion constant for a module
 */
typedef struct  {
  float center;  /**< center of the module (channel at which the radius is perpendicular to the module surface) */
  float ecenter; /**< error in the center determination */
  float r_conversion;  /**<  detector pixel size (or strip pitch) divided by the diffractometer radius */
  float er_conversion;  /**< error in the r_conversion determination */
  float offset; /**< the module offset i.e. the position of channel 0 with respect to the diffractometer 0 */
  float eoffset; /**< error in the offset determination */
  float tilt; /**< ossible tilt in the orthogonal direction (unused)*/
  float etilt; /**< error in the tilt determination */
} angleConversionConstant;



class angularCalibration {

 public:
  angularCalibration(int nm=48);
  ~angularCalibration();
  
  /**
     sets the angular direction of the detector
     \par d 1 or -1 set the angular direction, other valuse simply get
     \returns the angular direction of the detector
  */
  int setDirection(int d=0){if (d==-1 || d==1) direction=d; return direction;};

  /** 
      sets the encoder position
      \param f encoder position to be set
      \returns current encoder position
  */
  float setEncoder(float f) {encoder=f; return encoder;};  

  /** 
      gets the encoder position
      \returns encoder position
  */
  float getEncoder() {return encoder;};

  /** 
      sets the totalOffset of the detector
      \param f total offset to be set
      \returns current total offset
  */
  float setTotalOffset(float f) {totalOffset=f; return totalOffset;};  

  /** 
      gets the encoder position
      \returns encoder position
  */
  float getTotalOffset() {return totalOffset;};




  /**
     sets the angular range for peak fitting
     \param mi minimum of the angular range
     \param ma maximum of the angular range
  */
  void setAngularRange(float mi, float ma){ang_min=mi; ang_max=ma;};


  /**
     gets the angular range for peak fitting
     \param mi reference to the minimum of the angular range
     \param ma reference to the maximum of the angular range
  */
  void getAngularRange(float &mi, float &ma){mi=ang_min; ma=ang_max;};


  /** sets and returns the number of modules
      \param nm number of modules to be set (<0 gets)
      \return current number of modules
  */
  int setNumberOfModules(int nm=-1) {if (nm>=0) nmod=nm; return nmod;};

  /** sets and returns the number of channels per module
      \param n number of channels per module to be set (<0 gets)
      \return current number of channels per module
  */
  int setChannelsPerModule(int n=-1) {if (n>0) nchmod=n; return nchmod;};

  angleConversionConstant *getAngularConversionConstant(int imod=0);
  angleConversionConstant *setAngularConversionConstant(angleConversionConstant *a, int imod=0);


#ifdef ROOT

  /**
     Gaussian with pedestal describing a peak
     par[0] is the heigh of the pean
     par[1] is the peak position
     par[2] is the peak width
     par[3] is the background offset
     par[4] is the background slope
  */
  Double_t peakFunction(Double_t *x, Double_t *par); 


  /**
     Angular conversion function 
     par[0] is the module center
     par[1] is the conversion radius (pitch/radius)
     par[2] is the module offset
  */
  Double_t angleFunction(Double_t *x, Double_t *par);

  /**
     Fits a peak for the angular calibration 
     \param h histogram channels versus intensity
     \returns fitted function or NULL if fit failed
  */
  TF1 *fitPeak(TH1 *h);

#endif
  

 private:
  
  int direction; /**< angular direction of the detector -can be +1 or -1 */

#ifdef ROOT
  TF1 *fpeak; /**< Root function based on function peakFunction */
  
  TF1 *fangle; /**< Root function based on function angleFunction */
  
#endif
  float encoder; /**< position of the detector encoder */
  float totalOffset; /**< total offset of the detector */ 
  float ang_min; /**< minimum of the angular range for peak fitting*/
  float ang_max; /**< maximum of the angular range for peak fitting */
   
  int nmod;
  int nchmod;

  angleConversionConstant *angConv;

  



/* void fitangle(char fname[80],char extension[10], int start, int stop, float startangle, float stopangle);  //fits all datasets and extracts the constants */
/* int fitpeak(char fname[80],char extension[10], int nr, float minang, float maxang); // fits a peak from a pattern using nominal calibration constant */


};

#endif
