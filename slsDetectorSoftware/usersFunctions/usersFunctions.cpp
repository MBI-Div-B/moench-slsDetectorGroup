#include "usersFunctions.h"
#include <math.h>
#include <stdio.h>


double pos;
double i0=0;
#ifdef EPICS

#include <cadef.h>
#include <epicsEvent.h>

static double timeout = 3.0; 

chid ch_pos,ch_i0, ch_getpos;




/* connect to a PV */
int connect_channel(const char *name,  chid *ch_id) {
    int status = ECA_NORMAL; 
    status = ca_create_channel(name, NULL, NULL, CA_PRIORITY_DEFAULT, ch_id);
    if (status != ECA_NORMAL)
      return status; 

    status = ca_pend_io(timeout); 
    return status; 
}

/* disconnect to a PV */
int disconnect_channel(chid ch_id)
{
    ca_clear_channel(ch_id); 
    ca_flush_io(); 
}

int caget(chid ch_id,  double *value) {

    int status = ECA_NORMAL; 

    status = ca_get(DBR_DOUBLE, ch_id, value); 
    if (status  !=  ECA_NORMAL) {
        return status; 
    }

    status = ca_pend_io(timeout); 
    if (status  !=  ECA_NORMAL) {
        return status; 
    }

    return status; 
}

int caputq(chid ch_id,  double value) {
  // does not wait!
    int status = ECA_NORMAL; 

    status = ca_put(DBR_DOUBLE, ch_id, &value); 
    if (status  !=  ECA_NORMAL)
        return status; 

    status = ca_pend_io(timeout); 
    if (status  !=  ECA_NORMAL) {
        return status; 
    }

    return status; 
}

void put_callback(struct event_handler_args args)
{
    epicsEventId eid = (epicsEventId)args.usr; 
    epicsEventSignal(eid); 
}

int caput(chid ch_id,  double value) {

  // waits!
    int status = ECA_NORMAL; 
    epicsEventId eid = epicsEventCreate(epicsEventEmpty); 

    status = ca_put_callback(DBR_DOUBLE, 
            ch_id, 
            &value, 
            put_callback,
            eid); 
    status = ca_pend_io(timeout); 
    if (status  != ECA_NORMAL) 
        return status; 

    if (epicsEventWait(eid) != epicsEventWaitOK) 
        status = ECA_TIMEOUT; 

    return status; 
}

//int main(int argc,  char *argv[]) {



#endif






/* reads the encoder and returns the position */

double defaultGetPosition(void *d) {
#ifdef VERBOSE
  printf("Getting motor position \n");
#endif
  // caget X04SA-ES2-TH2:RO.RBV 

#ifdef EPICS
    int status;
    
    double value = 256;
    if (ch_getpos<0) return -1;

/*     /\* caget *\/ */
    if (caget(ch_getpos, &value) == ECA_NORMAL) {
#ifdef VERBOSE
        printf("caget: %f\n",  value);
#endif
	pos=value;
    } else
        printf(ca_message(status));
#else
    printf("\nmotor position is %f\n",pos);
#endif



  return pos;
}


/* moves the encoder to position p */


int defaultGoToPosition(double p,void *d) {
#ifdef VERBOSE
  printf("Setting  motor position \n");
#endif

#ifdef EPICS
    int status;
    if (ch_pos<0) return -1;
 /*    /\* caput and wait until done *\/ */
    if ((status = caput(ch_pos, p)) == ECA_NORMAL) {
      ;
#ifdef VERBOSE
        printf("caput: success\n");
#endif

    } else
        printf(ca_message(status));
#else
    pos=p;
    printf("\nmoving motor to position %f\n",p);
#endif
  //"caputq X04SA-ES2-TH2:RO p"
  //cawait -nounit  -timeout 3600 X04SA-ES2-TH2:RO.DMOV '==1'



    return (int)p;
}

/* moves the encoder to position p without waiting */

int defaultGoToPositionNoWait(double p,void *d) {
#ifdef VERBOSE
  printf("Setting  motor position no wait \n");
#endif


#ifdef EPICS
    int status;
    if (ch_pos<0) return -1;
 /*    /\* caput and wait until done *\/ */
    if ((status = caputq(ch_pos, p)) == ECA_NORMAL) {
      ;
#ifdef VERBOSE
        printf("caputq: success\n");
#endif
    } else
        printf(ca_message(status));
#else
    pos=p;
#endif
  //"caputq X04SA-ES2-TH2:RO p"

    return (int)p;

  pos=p;
  return (int)pos;
}


/* reads I0 and returns the intensity */

double defaultGetI0(int t,void *d) {
#ifdef VERBOSE
  printf("Getting I0 readout \n");
#endif

#ifdef EPICS
    int status;
    
    double value = 256;
/*     /\* caget *\/ */
    if (ch_i0<0) return -1;
    if (caget(ch_i0, &value) == ECA_NORMAL) {
#ifdef VERBOSE
        printf("caget: %f\n",  value);
#endif

	
	if (t==0)
	  i0=value;
	else
	  i0=value-i0;
	


    } else
        printf(ca_message(status));
#else
    i0++;
#endif

  //"ca_get X04SA-ES2-SC:CH6"
  return i0;

}
  

int defaultConnectChannels(void*d) {
#ifdef EPICS
  //double value = 256;
    /* channel name */
    //const char *name = "ARIDI-PCT:CURRENT";
    /* channel id */
    /* status code */
    int status;

    printf("starting...\n");

    /* init channel access context before any caget/put */
    ca_context_create(ca_enable_preemptive_callback);

    printf("context created\n");

    //"caputq X04SA-ES2-TH2:RO p"
   
    //"ca_get X04SA-ES2-SC:CH6"

    /* open the channel by name and return ch_id */
    status = connect_channel("X04SA-ES2-SC:CH6",  &ch_i0);
    if (status  ==  ECA_NORMAL)
        printf("I0 channel connected \n");
    else {
        printf(ca_message(status));
        //ch_i0=-1;;
    }
    status = connect_channel("X04SA-ES2-TH2:RO",  &ch_pos);
    if (status  ==  ECA_NORMAL)
        printf("Detector position channel connected \n");
    else {
        printf(ca_message(status));
        //ch_i0=-1;;
    }
        status = connect_channel("X04SA-ES2-TH2:RO.RBV",  &ch_getpos);
    if (status  ==  ECA_NORMAL)
        printf("Detector get position channel connected \n");
    else {
        printf(ca_message(status));
        //ch_getpos=-1;;
    }
   
  // caget X04SA-ES2-TH2:RO.RBV 

 //cawait -nounit  -timeout 3600 X04SA-ES2-TH2:RO.DMOV '==1'
#endif
    return 0;
}

int defaultDisconnectChannels(void *d) {  
#ifdef EPICS
    /* close channel connect */
    disconnect_channel(ch_i0);
    disconnect_channel(ch_pos);
    disconnect_channel(ch_getpos);

    /* delete channel access context before program exits */
    ca_context_destroy();
#endif
    return 0;
}



int defaultDataReadyFunc(detectorData* d,  void* p) {
#ifdef VERBOSE
  printf("UUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUU Data received \n");
  printf("Filename %s\n",d->fileName);
  printf("Number of points %d (%d)\n",d->npoints,d->npy);
  for(int i=0;i<d->npoints;i++){
	  printf("%d-[%f]:%f\n",i,*(d->angles+(i*sizeof(double))),*(d->values+(i*sizeof(double))));
	  fflush(stdout);
  }
  printf("UUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUU Finished \n");
#endif
  return 0;
}




