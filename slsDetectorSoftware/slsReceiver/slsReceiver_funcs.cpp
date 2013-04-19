/********************************************//**
 * @file slsReceiver_funcs.h
 * @short interface between receiver and client
 ***********************************************/

#include "slsReceiver_funcs.h"
#include "slsReceiverFunctionList.h"
#include "svnInfoReceiver.h"

#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <stdlib.h>
using namespace std;


int slsReceiverFuncs::file_des(-1);
int slsReceiverFuncs::socketDescriptor(-1);


slsReceiverFuncs::~slsReceiverFuncs() {

	slsReceiverFuncs::closeFile(0);
	cout << "Goodbye!" << endl;
	if(socket) delete socket;

}


slsReceiverFuncs::slsReceiverFuncs(int argc, char *argv[], int &success):
		myDetectorType(GOTTHARD),
		socket(NULL),
		ret(OK),
		lockStatus(0),
		shortFrame(-1),
		packetsPerFrame(GOTTHARD_PACKETS_PER_FRAME){

	int port_no = DEFAULT_PORTNO+2;
	ifstream infile;
	string sLine,sargname;
	int iline = 0;

	success=OK;
	string fname = "";

	//parse command line for config
	for(int iarg=1;iarg<argc;iarg++){
		if((!strcasecmp(argv[iarg],"-config"))||(!strcasecmp(argv[iarg],"-f"))){
			if(iarg+1==argc){
				cout << "no config file name given. Exiting." << endl;
				success=FAIL;
			}else
				fname.assign(argv[iarg+1]);
		}
		if(!strcasecmp(argv[iarg],"-type")){
			if(iarg+1==argc){
				cout << "no detector type given after -type in command line. Exiting." << endl;
				success=FAIL;
			}else{
				if(!strcasecmp(argv[iarg+1],"gotthard"));
				else if(!strcasecmp(argv[iarg+1],"moench"))
					slsReceiverFuncs::myDetectorType = MOENCH;
				else{
					cout << "could not decode detector type in command line. Exiting." << endl;
					success=FAIL;
				}
			}
		}
	}

	if((!fname.empty()) && (success == OK)){
#ifdef VERBOSE
		std::cout<< "config file name "<< fname << std::endl;
#endif
		infile.open(fname.c_str(), ios_base::in);
		if (infile.is_open()) {
			while(infile.good()){
				getline(infile,sLine);
				iline++;
#ifdef VERBOSE
				cout <<  str << endl;
#endif
				if(sLine.find('#')!=string::npos){
#ifdef VERBOSE
					cout << "Line is a comment " << endl;
#endif
					continue;
				}else if(sLine.length()<2){
#ifdef VERBOSE
					cout << "Empty line " << endl;
#endif
					continue;
				}else{
					istringstream sstr(sLine);

					//parameter name
					if(sstr.good())
						sstr >> sargname;



					//value
					if(sargname=="rx_tcpport"){
						if(sstr.good()) {
							sstr >> sargname;
							if(sscanf(sargname.c_str(),"%d",&port_no))
								cout<<"dataport:"<<port_no<<endl;
							else{
								cout << "could not decode port in config file. Exiting." << endl;
								success=FAIL;
							}
						}
					}

					else if(sargname=="type"){
						if(sstr.good()) {
							sstr >> sargname;
							if(!strcasecmp(sargname.c_str(),"gotthard"))
								slsReceiverFuncs::myDetectorType = GOTTHARD;
							else if(!strcasecmp(sargname.c_str(),"moench"))
								slsReceiverFuncs::myDetectorType = MOENCH;
							else{
								cout << "could not decode detector type in config file. Exiting." << endl;
								success=FAIL;
							}
						}
					}



				}
			}
			infile.close();
		}else {
			cout << "Error opening configuration file " << fname << endl;
			success = FAIL;
		}
#ifdef VERBOSE
		cout << "Read configuration file of " << iline << " lines" << endl;
#endif

	}

	if(success == OK){
		//create socket
		socket = new MySocketTCP(port_no);
		if (socket->getErrorStatus()) {
			success = FAIL;
			delete socket;
			socket=NULL;
		}	else {
			//initialize variables
			strcpy(socket->lastClientIP,"none");
			strcpy(socket->thisClientIP,"none1");
			strcpy(mess,"dummy message");

			function_table();
			slsReceiverList =  new slsReceiverFunctionList(myDetectorType);

#ifdef VERBOSE
	cout << "Function table assigned." << endl;
#endif

			file_des=socket->getFileDes();
			socketDescriptor=socket->getsocketDescriptor();

			//success = OK;
		}
	}

	switch(myDetectorType){
	case GOTTHARD: 	cout << "This is a GOTTHARD Receiver" << endl; 	break;
	case MOENCH:	cout << "This is a MOENCH Receiver" << endl; 	break;
	default:		cout << "Unknown Receiver" << endl;success=FAIL;break;
	}
}

void slsReceiverFuncs::start(){

  int v=slsDetectorDefs::OK;

	while(v!=GOODBYE) {
#ifdef VERBOSE
		cout<< endl;
#endif
#ifdef VERY_VERBOSE
		cout << "Waiting for client call" << endl;
#endif
		if(socket->Connect()>=0){
#ifdef VERY_VERBOSE
		cout << "Conenction accepted" << endl;
#endif
		v = decode_function();
#ifdef VERY_VERBOSE
		cout << "function executed" << endl;
#endif
		socket->Disconnect();
#ifdef VERY_VERBOSE
		cout << "connection closed" << endl;
#endif
		}
	}
	


}


int slsReceiverFuncs::function_table(){

	for (int i=0;i<numberOfFunctions;i++)
		flist[i]=&slsReceiverFuncs::M_nofunc;

	flist[F_SET_FILE_NAME]			=	&slsReceiverFuncs::set_file_name;
	flist[F_SET_FILE_PATH]			=	&slsReceiverFuncs::set_file_dir;
	flist[F_SET_FILE_INDEX]		 	=	&slsReceiverFuncs::set_file_index;
	flist[F_SET_FRAME_INDEX]	 	=	&slsReceiverFuncs::set_frame_index;
	flist[F_SETUP_UDP]				=	&slsReceiverFuncs::setup_udp;
	flist[F_START_RECEIVER]			=	&slsReceiverFuncs::start_receiver;
	flist[F_STOP_RECEIVER]			=	&slsReceiverFuncs::stop_receiver;
	flist[F_GET_RECEIVER_STATUS]	=	&slsReceiverFuncs::get_status;
	flist[F_GET_FRAMES_CAUGHT]		=	&slsReceiverFuncs::get_frames_caught;
	flist[F_GET_FRAME_INDEX]		=	&slsReceiverFuncs::get_frame_index;
	flist[F_RESET_FRAMES_CAUGHT]	=	&slsReceiverFuncs::reset_frames_caught;
	flist[F_READ_FRAME]				=	&slsReceiverFuncs::read_frame;
	flist[F_READ_ALL]				=  	&slsReceiverFuncs::read_all;
	flist[F_READ_RECEIVER_FREQUENCY]= 	&slsReceiverFuncs::set_read_frequency;
	flist[F_ENABLE_FILE_WRITE]		=	&slsReceiverFuncs::enable_file_write;
	flist[F_GET_ID]					=	&slsReceiverFuncs::get_version;
	flist[F_CONFIGURE_MAC]			=	&slsReceiverFuncs::set_short_frame;

	//General Functions
	flist[F_LOCK_SERVER]			=	&slsReceiverFuncs::lock_receiver;
	flist[F_SET_PORT]				=	&slsReceiverFuncs::set_port;
	flist[F_GET_LAST_CLIENT_IP]		=	&slsReceiverFuncs::get_last_client_ip;
	flist[F_UPDATE_CLIENT]			=	&slsReceiverFuncs::update_client;
	flist[F_EXIT_SERVER]			=	&slsReceiverFuncs::exit_server;		//not implemented in client
	flist[F_EXEC_COMMAND]			=	&slsReceiverFuncs::exec_command;	//not implemented in client


#ifdef VERBOSE
	for (i=0;i<numberOfFunctions;i++)
		cout << "function " << i << "located at " << flist[i] << endl;
#endif
	return OK;

}





int slsReceiverFuncs::decode_function(){
	ret = FAIL;
	int n,fnum;
#ifdef VERBOSE
	cout <<  "receive data" << endl;
#endif
	n = socket->ReceiveDataOnly(&fnum,sizeof(fnum));
	if (n <= 0) {
#ifdef VERBOSE
		cout << "ERROR reading from socket " << n << ", " << fnum << endl;
#endif
		return FAIL;
	}
#ifdef VERBOSE
	else
		cout << "size of data received " << n <<endl;
#endif

#ifdef VERBOSE
	cout <<  "calling function fnum = "<< fnum << hex << ":"<< flist[fnum] << endl;
#endif

	if (fnum<0 || fnum>numberOfFunctions-1)
		fnum = numberOfFunctions-1;
	//calling function
	(this->*flist[fnum])();
	if (ret==FAIL)
		cout <<  "Error executing the function = " << fnum << endl;

	return ret;
}






int slsReceiverFuncs::M_nofunc(){

	ret=FAIL;
	sprintf(mess,"Unrecognized Function\n");
	cout << mess << endl;

	socket->SendDataOnly(&ret,sizeof(ret));
	socket->SendDataOnly(mess,sizeof(mess));

	return GOODBYE;
}




void slsReceiverFuncs::closeFile(int p){


	if(slsReceiverFunctionList::listening_thread_running)
		fclose(slsReceiverFunctionList::sfilefd);


	close(file_des);
	//shutdown(socketDescriptor,SHUT_RDWR);
	close(socketDescriptor);

}




int slsReceiverFuncs::set_file_name() {
	ret=OK;
	char retval[MAX_STR_LENGTH]="";
	char fName[MAX_STR_LENGTH];
	strcpy(mess,"Could not set file name");

	// receive arguments
	if(socket->ReceiveDataOnly(fName,MAX_STR_LENGTH) < 0 ){
		strcpy(mess,"Error reading from socket\n");
		ret = FAIL;
	}

	// execute action if the arguments correctly arrived
#ifdef SLS_RECEIVER_FUNCTION_LIST
	if (ret==OK) {

		if (lockStatus==1 && socket->differentClients==1){
			sprintf(mess,"Receiver locked by %s\n", socket->lastClientIP);
			ret=FAIL;
		}
 	 	 else
			strcpy(retval,slsReceiverList->setFileName(fName));
	}
#ifdef VERBOSE
	if(ret!=FAIL)
		cout << "file name:" << retval << endl;
	else
		cout << mess << endl;
#endif
#endif

	if(ret==OK && socket->differentClients){
		cout << "Force update" << endl;
		ret=FORCE_UPDATE;
	}

	// send answer
	socket->SendDataOnly(&ret,sizeof(ret));
	if(ret==FAIL)
		socket->SendDataOnly(mess,sizeof(mess));
	socket->SendDataOnly(retval,MAX_STR_LENGTH);

	//return ok/fail
	return ret;
}






int slsReceiverFuncs::set_file_dir() {
	ret=OK;
	char retval[MAX_STR_LENGTH]="";
	char fPath[MAX_STR_LENGTH];
	strcpy(mess,"Could not set file path\n");

	// receive arguments
	if(socket->ReceiveDataOnly(fPath,MAX_STR_LENGTH) < 0 ){
		strcpy(mess,"Error reading from socket\n");
		ret = FAIL;
	}

	// execute action if the arguments correctly arrived
#ifdef SLS_RECEIVER_FUNCTION_LIST
	if (ret==OK) {
		if (lockStatus==1 && socket->differentClients==1){
			sprintf(mess,"Receiver locked by %s\n", socket->lastClientIP);
			ret=FAIL;
		}
		else if((strlen(fPath))&&(slsReceiverList->getStatus()==RUNNING)){
			strcpy(mess,"Can not set file path while receiver running\n");
			ret = FAIL;
		}
		else{
			strcpy(retval,slsReceiverList->setFilePath(fPath));
			// if file path doesnt exist
			if(strlen(fPath))
				if (strcmp(retval,fPath)){
					strcpy(mess,"receiver file path does not exist\n");
					ret=FAIL;
				}
		}

	}
#ifdef VERBOSE
	if(ret!=FAIL)
		cout << "file path:" << retval << endl;
	else
		cout << mess << endl;
#endif
#endif

	if(ret==OK && socket->differentClients){
		cout << "Force update" << endl;
		ret=FORCE_UPDATE;
	}

	// send answer
	socket->SendDataOnly(&ret,sizeof(ret));
	if(ret==FAIL)
		socket->SendDataOnly(mess,sizeof(mess));
	socket->SendDataOnly(retval,MAX_STR_LENGTH);

	//return ok/fail
	return ret;
}






int slsReceiverFuncs::set_file_index() {
	ret=OK;
	int retval=-1;
	int index;
	strcpy(mess,"Could not set file index\n");


	// receive arguments
	if(socket->ReceiveDataOnly(&index,sizeof(index)) < 0 ){
		strcpy(mess,"Error reading from socket\n");
		ret = FAIL;
	}

	// execute action if the arguments correctly arrived
#ifdef SLS_RECEIVER_FUNCTION_LIST
	if (ret==OK) {
		if (lockStatus==1 && socket->differentClients==1){//necessary???
			sprintf(mess,"Receiver locked by %s\n", socket->lastClientIP);
			ret=FAIL;
		}
		else
			retval=slsReceiverList->setFileIndex(index);
	}
#ifdef VERBOSE
	if(ret!=FAIL)
		cout << "file index:" << retval << endl;
	else
		cout << mess << endl;
#endif
#endif

	if(ret==OK && socket->differentClients){
		cout << "Force update" << endl;
		ret=FORCE_UPDATE;
	}

	// send answer
	socket->SendDataOnly(&ret,sizeof(ret));
	if(ret==FAIL)
		socket->SendDataOnly(mess,sizeof(mess));
	socket->SendDataOnly(&retval,sizeof(retval));

	//return ok/fail
	return ret;
}








int slsReceiverFuncs::set_frame_index() {
	ret=OK;
	int retval=-1;
	int index;
	strcpy(mess,"Could not set frame index\n");


	// receive arguments
	if(socket->ReceiveDataOnly(&index,sizeof(index)) < 0 ){
		strcpy(mess,"Error reading from socket\n");
		ret = FAIL;
	}

	// execute action if the arguments correctly arrived
#ifdef SLS_RECEIVER_FUNCTION_LIST
	if (ret==OK) {
		if (lockStatus==1 && socket->differentClients==1){//necessary???
			sprintf(mess,"Receiver locked by %s\n", socket->lastClientIP);
			ret=FAIL;
		}
		else
			retval=slsReceiverList->setFrameIndexNeeded(index);
	}
#ifdef VERBOSE
	if(ret!=FAIL)
		cout << "frame index:" << retval << endl;
	else
		cout << mess << endl;
#endif
#endif

	if(ret==OK && socket->differentClients){
		cout << "Force update" << endl;
		ret=FORCE_UPDATE;
	}

	// send answer
	socket->SendDataOnly(&ret,sizeof(ret));
	if(ret==FAIL)
		socket->SendDataOnly(mess,sizeof(mess));
	socket->SendDataOnly(&retval,sizeof(retval));

	//return ok/fail
	return ret;
}







int slsReceiverFuncs::setup_udp(){
	ret=OK;
	strcpy(mess,"could not set up udp connection");
	char retval[MAX_STR_LENGTH]="";
	char args[2][MAX_STR_LENGTH];

	string temp;
	int udpport;
	char eth[MAX_STR_LENGTH];


	// receive arguments

	if(socket->ReceiveDataOnly(args,sizeof(args)) < 0 ){
		strcpy(mess,"Error reading from socket\n");
		ret = FAIL;
	}

	// execute action if the arguments correctly arrived
#ifdef SLS_RECEIVER_FUNCTION_LIST
	if (ret==OK) {
		if (lockStatus==1 && socket->differentClients==1){//necessary???
			sprintf(mess,"Receiver locked by %s\n", socket->lastClientIP);
			ret=FAIL;
		}
		else if(slsReceiverList->getStatus()==RUNNING){
			ret = FAIL;
			strcpy(mess,"cannot set up udp when receiver is running\n");
		}
		else{
			//set up udp port
			 sscanf(args[1],"%d",&udpport);
			 slsReceiverList->setUDPPortNo(udpport);

			 //setup udpip
			 //get ethernet interface or IP to listen to
			 temp = genericSocket::ipToName(args[0]);
			 if(temp=="none"){
				 ret = FAIL;
				 strcpy(mess, "failed to get ethernet interface or IP to listen to\n");
			 }
			 else{
				 strcpy(eth,temp.c_str());
				 if (strchr(eth,'.')!=NULL) {
					 strcpy(eth,"");
					 ret = FAIL;
				 }
				 cout<<"eth:"<<eth<<endl;
				 slsReceiverList->setEthernetInterface(eth);

				 //get mac address from ethernet interface
				 if (ret != FAIL)
					temp = genericSocket::nameToMac(eth);


				 if ((temp=="00:00:00:00:00:00") || (ret == FAIL)){
					 ret = FAIL;
					 strcpy(mess,"failed to get mac adddress to listen to\n");
					 cout << "mess:" << mess << endl;
				 }
				 else{
					 strcpy(retval,temp.c_str());
					 cout<<"mac:"<<retval<<endl;
				 }
			 }
		}
	}
#endif

	if(ret==OK && socket->differentClients){
		cout << "Force update" << endl;
		ret=FORCE_UPDATE;
	}

	// send answer
	socket->SendDataOnly(&ret,sizeof(ret));
	if(ret==FAIL)
		socket->SendDataOnly(mess,sizeof(mess));
	socket->SendDataOnly(retval,MAX_STR_LENGTH);

	//return ok/fail
	return ret;
}






int slsReceiverFuncs::start_receiver(){
	ret=OK;

	strcpy(mess,"Could not start receiver\n");

	// execute action if the arguments correctly arrived
#ifdef SLS_RECEIVER_FUNCTION_LIST
	if (lockStatus==1 && socket->differentClients==1){
		sprintf(mess,"Receiver locked by %s\n", socket->lastClientIP);
		ret=FAIL;
	}
	else if(!strlen(slsReceiverList->getFilePath())){
		strcpy(mess,"receiver not set up. set receiver ip again.\n");
		ret = FAIL;
	}
	else if(slsReceiverList->getStatus()!=RUNNING)
		ret=slsReceiverList->startReceiver();
#endif

	if(ret==OK && socket->differentClients){
		cout << "Force update" << endl;
		ret=FORCE_UPDATE;
	}

	// send answer
	socket->SendDataOnly(&ret,sizeof(ret));
	if(ret==FAIL)
		socket->SendDataOnly(mess,sizeof(mess));
	//return ok/fail
	return ret;


}


int slsReceiverFuncs::stop_receiver(){
	ret=OK;

	strcpy(mess,"Could not stop receiver\n");

	// execute action if the arguments correctly arrived
#ifdef SLS_RECEIVER_FUNCTION_LIST
	if (lockStatus==1 && socket->differentClients==1){
		sprintf(mess,"Receiver locked by %s\n", socket->lastClientIP);
		ret=FAIL;
	}
	else if(slsReceiverList->getStatus()!=IDLE)
		ret=slsReceiverList->stopReceiver();
#endif

	if(ret==OK && socket->differentClients){
		cout << "Force update" << endl;
		ret=FORCE_UPDATE;
	}

	// send answer
	socket->SendDataOnly(&ret,sizeof(ret));
	if(ret==FAIL)
		socket->SendDataOnly(mess,sizeof(mess));
	//return ok/fail
	return ret;


}


int	slsReceiverFuncs::get_status(){
	ret=OK;
	enum runStatus retval;

	// execute action if the arguments correctly arrived
#ifdef SLS_RECEIVER_FUNCTION_LIST
	retval=slsReceiverList->getStatus();
#endif

	if(socket->differentClients){
		cout << "Force update" << endl;
		ret=FORCE_UPDATE;
	}

	// send answer
	socket->SendDataOnly(&ret,sizeof(ret));
	socket->SendDataOnly(&retval,sizeof(retval));
	//return ok/fail
	return ret;


}


int	slsReceiverFuncs::get_frames_caught(){
	ret=OK;
	int retval=-1;

	// execute action if the arguments correctly arrived
#ifdef SLS_RECEIVER_FUNCTION_LIST
	retval=slsReceiverList->getTotalFramesCaught();
#endif
	if(socket->differentClients){
		cout << "Force update" << endl;
		ret=FORCE_UPDATE;
	}

	// send answer
	socket->SendDataOnly(&ret,sizeof(ret));
	socket->SendDataOnly(&retval,sizeof(retval));
	//return ok/fail
	return ret;


}


int	slsReceiverFuncs::get_frame_index(){
	ret=OK;
	int retval=-1;

	// execute action if the arguments correctly arrived
#ifdef SLS_RECEIVER_FUNCTION_LIST
	retval=slsReceiverList->getAcquisitionIndex();
#endif

	if(socket->differentClients){
		cout << "Force update" << endl;
		ret=FORCE_UPDATE;
	}

	// send answer
	socket->SendDataOnly(&ret,sizeof(ret));
	socket->SendDataOnly(&retval,sizeof(retval));
	//return ok/fail
	return ret;


}


int	slsReceiverFuncs::reset_frames_caught(){
	ret=OK;

	strcpy(mess,"Could not reset frames caught\n");


	// execute action if the arguments correctly arrived
#ifdef SLS_RECEIVER_FUNCTION_LIST
	if (ret==OK) {
		if (lockStatus==1 && socket->differentClients==1){
			sprintf(mess,"Receiver locked by %s\n", socket->lastClientIP);
			ret=FAIL;
		}
		else
			slsReceiverList->resetTotalFramesCaught();
	}
#endif

	if(ret==OK && socket->differentClients){
		cout << "Force update" << endl;
		ret=FORCE_UPDATE;
	}

	// send answer
	socket->SendDataOnly(&ret,sizeof(ret));
	if(ret==FAIL)
		socket->SendDataOnly(mess,sizeof(mess));

	//return ok/fail
	return ret;


}






int slsReceiverFuncs::set_short_frame() {
	ret=OK;
	int index=0;
	int retval=-100;
	strcpy(mess,"Could not set/reset short frame for receiver\n");

	//does not exist for moench
	if(myDetectorType==MOENCH){
		strcpy(mess,"can not set short frame for moench\n");
		ret = FAIL;
	}

	// receive arguments
	if(socket->ReceiveDataOnly(&index,sizeof(index)) < 0 ){
		strcpy(mess,"Error reading from socket\n");
		ret = FAIL;
	}


	// execute action if the arguments correctly arrived
#ifdef SLS_RECEIVER_FUNCTION_LIST
	if (ret==OK) {
		if (lockStatus==1 && socket->differentClients==1){//necessary???
			sprintf(mess,"Receiver locked by %s\n", socket->lastClientIP);
			ret=FAIL;
		}
		else if(slsReceiverList->getStatus()==RUNNING){
			strcpy(mess,"Cannot set short frame while status is running\n");
			ret=FAIL;
		}
		else{
			retval=slsReceiverList->setShortFrame(index);
			shortFrame = retval;
			if(shortFrame==-1)
				packetsPerFrame=GOTTHARD_PACKETS_PER_FRAME;
			else
				packetsPerFrame=GOTTHARD_SHORT_PACKETS_PER_FRAME;
		}
	}
#endif

	if(ret==OK && socket->differentClients){
		cout << "Force update" << endl;
		ret=FORCE_UPDATE;
	}

	// send answer
	socket->SendDataOnly(&ret,sizeof(ret));
	if(ret==FAIL)
		socket->SendDataOnly(mess,sizeof(mess));
	socket->SendDataOnly(&retval,sizeof(retval));

	//return ok/fail
	return ret;
}




int	slsReceiverFuncs::read_frame(){
	switch(myDetectorType){
	case MOENCH:
		return moench_read_frame();
	default:
		return gotthard_read_frame();
	}
}



int	slsReceiverFuncs::moench_read_frame(){
	ret=OK;
	char fName[MAX_STR_LENGTH]="";
	int arg = -1,i,j,x,y;


	int bufferSize = MOENCH_BUFFER_SIZE;
	char* raw 		= new char[bufferSize];

	int rnel 		= bufferSize/(sizeof(int));
	int* retval 	= new int[rnel];
	int* origVal 	= new int[rnel];
	//all initialized to 0
	for(i=0;i<rnel;i++)	retval[i]=0;
	for(i=0;i<rnel;i++)	origVal[i]=0;

	int onebuffersize = bufferSize/MOENCH_PACKETS_PER_FRAME;
	int onedatasize = MOENCH_DATA_BYTES/MOENCH_PACKETS_PER_FRAME;

	int index=-1;//,index2=-1;
	int startIndex=-1;
	int count=0;

	int offset=0;
	int partsPerFrame;

	strcpy(mess,"Could not read frame\n");


	// execute action if the arguments correctly arrived
#ifdef SLS_RECEIVER_FUNCTION_LIST
	ret=FAIL;
	strcpy(mess,"did not start index\n");
	for(i=0;i<10;i++){
		startIndex=slsReceiverList->getStartFrameIndex();
		if(startIndex!=-1){
			ret=OK;
			break;
		}else
			usleep(500000);
	}


	if(ret==FAIL)
		cout<<"failed to start"<<endl;
	else{
		//get frame
		slsReceiverList->readFrame(fName,&raw);
		if (raw == NULL){
			index = startIndex;
			cout << "didnt get data. Gui will try again" << endl;
		}else{
			//upto 40 indices, look at just index1 and index2 for now
			index=(int)(*(int*)raw);
			memcpy(origVal,raw,bufferSize);
			raw=NULL;
		}

/*
			offset = 4;j=0;
			partsPerFrame = onedatasize/MOENCH_BYTES_PER_ADC;
			//filling up in y direction and then in x direcction
			for(x=0;x<MOENCH_BYTES_IN_ONE_DIMENSION/MOENCH_BYTES_PER_ADC;x++){
				for(y=0;j<MOENCH_PIXELS_IN_ONE_DIMENSION;j++){

					memcpy((((char*)retval) + y*MOENCH_BYTES_IN_ONE_DIMENSION + x*MOENCH_BYTES_PER_ADC),
							(((char*) origVal) + offset + j*MOENCH_BYTES_PER_ADC) ,
							MOENCH_BYTES_PER_ADC);

					j++;
					//after 1280 bytes(16 parts of 80 bytes), add 6 bytes to offset
					if (!(j%partsPerFrame))
						offset+=6;
				}
			}
*/

		for(i=0;i<MOENCH_PACKETS_PER_FRAME;i=i+2){
			memcpy((((char*)retval)+ onedatasize*i),		(((char*) origVal)+4+    			onebuffersize*i) , 	 onedatasize);
			memcpy((((char*)retval)+ onedatasize*(i+1)), 	(((char*) origVal)+10+onedatasize+	onebuffersize*i),onedatasize);
		}


		arg=((index - startIndex)/MOENCH_PACKETS_PER_FRAME)-1;
#ifdef VERBOSE
		cout << "\nstartIndex:" << startIndex << endl;
		cout << "fName:" << fName << endl;
		cout << "index:" << arg << endl;
#endif

	}

#endif

	if(ret==OK && socket->differentClients){
		cout << "Force update" << endl;
		ret=FORCE_UPDATE;
	}

	// send answer
	socket->SendDataOnly(&ret,sizeof(ret));
	if(ret==FAIL){
		cout << "mess:" << mess << endl;
		socket->SendDataOnly(mess,sizeof(mess));
	}
	else{
		socket->SendDataOnly(fName,MAX_STR_LENGTH);
		socket->SendDataOnly(&arg,sizeof(arg));
		socket->SendDataOnly(retval,MOENCH_DATA_BYTES);
	}
	//return ok/fail

    ///ADDED BY ANNA?!?!?!?   //?????????FOR GOTTHARD AS WELL?????

	delete [] origVal;
 	delete [] retval;


	return ret;

}




int	slsReceiverFuncs::gotthard_read_frame(){
	ret=OK;
	char fName[MAX_STR_LENGTH]="";
	int arg = -1,i;


	//retval is a full frame
	int rnel 		= GOTTHARD_BUFFER_SIZE/(sizeof(int));
	int* retval 	= new int[rnel];
	//all initialized to 0
	for(i=0;i<rnel;i++)	retval[i]=0;

	//only for full frames
	int onebuffersize = GOTTHARD_BUFFER_SIZE/GOTTHARD_PACKETS_PER_FRAME;
	int onedatasize = GOTTHARD_DATA_BYTES/GOTTHARD_PACKETS_PER_FRAME;


	//depending on shortframe or not
	int bufferSize = GOTTHARD_BUFFER_SIZE;
	if(shortFrame!=-1)
		bufferSize=GOTTHARD_SHORT_BUFFER_SIZE;
	char* raw 		= new char[bufferSize];
	int* origVal 	= new int[rnel];

	//int* emptys 	= new int[rnel]();

	int index=-1,index2=-1;
	int startIndex=-1;
	int count=0;

	strcpy(mess,"Could not read frame\n");


	// execute action if the arguments correctly arrived
#ifdef SLS_RECEIVER_FUNCTION_LIST

	//wait till you get first frame 1. to get index(from startindex 2. filename corresponds to buffer value
	ret=FAIL;
	strcpy(mess,"did not start index\n");
	for(int i=0;i<10;i++){
		startIndex=slsReceiverList->getStartFrameIndex();
		if(startIndex!=-1){
			ret=OK;
			break;
		}else
			usleep(500000);
	}


	//got atleast first frame, read buffer
	if(ret==OK){
		if(shortFrame!=-1){
			//get frame
			slsReceiverList->readFrame(fName,&raw);
			//sending garbage values with index = -1 to try again
			if (raw == NULL){
				index = startIndex;
				memcpy((((char*)retval)+(GOTTHARD_SHORT_DATABYTES*shortFrame)),((char*) origVal)+4, GOTTHARD_SHORT_DATABYTES);
			}else{
				index=(int)(*(int*)raw);
				memcpy(origVal,raw,bufferSize);
				raw=NULL;

				//copy only the adc part
				memcpy((((char*)retval)+(GOTTHARD_SHORT_DATABYTES*shortFrame)),((char*) origVal)+4, GOTTHARD_SHORT_DATABYTES);
			}

		}else{
			//for full frames
			while(count<10){
				//get frame
				slsReceiverList->readFrame(fName,&raw);
				if (raw == NULL){
					count = -1;
					break;
				}

				index=(int)(*(int*)raw);
				index2= (int)(*((int*)((char*)(raw+onebuffersize))));
				memcpy(origVal,raw,bufferSize);
				raw=NULL;
				//cout<<"funcs\tindex:"<<index<<"\tindex2:"<<index2<<endl;


				//1 odd, 1 even
				if((index%2)!=index2%2){
					//ideal situation (should be odd, even(index+1))
					if(index%2){
						memcpy(retval,((char*) origVal)+4, onedatasize);
						memcpy((((char*)retval)+onedatasize), ((char*) origVal)+10+onedatasize, onedatasize);
						break;
					}

					//swap to even,odd
					if(index2%2){
						memcpy((((char*)retval)+onedatasize),((char*) origVal)+4, onedatasize);
						memcpy(retval, ((char*) origVal)+10+onedatasize, onedatasize);
						index=index2;
						break;
					}

				}
				strcpy(mess,"could not read frame due to more than 20 mismatched indices\n");

				usleep(100000);
				count++;
			}

			if ((count==10) || (count == -1)){
				if (count == -10)
					cout << "same type: index:" << index << "\tindex2:" << index2 << endl;
				else
					cout << "no data to read for gui" << endl;
				/**send garbage with -1 index to try again*/
				index = startIndex;
				memcpy(retval,((char*) origVal)+4, onedatasize);
				memcpy((((char*)retval)+onedatasize), ((char*) origVal)+10+onedatasize, onedatasize);
			}
		}

		arg=((index - startIndex)/packetsPerFrame)-1;

#ifdef VERBOSE
		cout << "\nstartIndex:" << startIndex << endl;
		cout << "fName:" << fName << endl;
		cout << "index:" << arg << endl;
#endif

	}else cout<<"failed to start"<<endl;


#endif

	if(ret==OK && socket->differentClients){
		cout << "Force update" << endl;
		ret=FORCE_UPDATE;
	}

	// send answer
	socket->SendDataOnly(&ret,sizeof(ret));
	if(ret==FAIL){
		cout << "mess:" << mess << endl;
		socket->SendDataOnly(mess,sizeof(mess));
	}
	else{
		socket->SendDataOnly(fName,MAX_STR_LENGTH);
		socket->SendDataOnly(&arg,sizeof(arg));
		socket->SendDataOnly(retval,GOTTHARD_DATA_BYTES);
	}
	//return ok/fail

	///ADDED BY ANNA?!?!?!?
	delete [] retval;
	delete [] origVal;


	return ret;
}




int slsReceiverFuncs::set_read_frequency(){
	ret=OK;
	int retval=-1;
	int index;
	strcpy(mess,"Could not set receiver read frequency\n");


	// receive arguments
	if(socket->ReceiveDataOnly(&index,sizeof(index)) < 0 ){
		strcpy(mess,"Error reading from socket\n");
		ret = FAIL;
	}

	// execute action if the arguments correctly arrived
#ifdef SLS_RECEIVER_FUNCTION_LIST
	if (ret==OK) {
		if (lockStatus==1 && socket->differentClients==1){//necessary???
			sprintf(mess,"Receiver locked by %s\n", socket->lastClientIP);
			ret=FAIL;
		}
		else
			retval=slsReceiverList->setNFrameToGui(index);
	}

#endif

	if(ret==OK && socket->differentClients){
		cout << "Force update" << endl;
		ret=FORCE_UPDATE;
	}

	// send answer
	socket->SendDataOnly(&ret,sizeof(ret));
	if(ret==FAIL)
		socket->SendDataOnly(mess,sizeof(mess));
	socket->SendDataOnly(&retval,sizeof(retval));

	//return ok/fail
	return ret;
}




/**needs to be implemented */
int slsReceiverFuncs::read_all(){
	return ret;
}





int slsReceiverFuncs::enable_file_write(){
	ret=OK;
	int retval=-1;
	int enable;
	strcpy(mess,"Could not set/get enable file write\n");


	// receive arguments
	if(socket->ReceiveDataOnly(&enable,sizeof(enable)) < 0 ){
		strcpy(mess,"Error reading from socket\n");
		ret = FAIL;
	}

	// execute action if the arguments correctly arrived
#ifdef SLS_RECEIVER_FUNCTION_LIST
	if (ret==OK) {
		if (lockStatus==1 && socket->differentClients==1){//necessary???
			sprintf(mess,"Receiver locked by %s\n", socket->lastClientIP);
			ret=FAIL;
		}
		else{
			retval=slsReceiverList->setEnableFileWrite(enable);
			if((enable!=-1)&&(enable!=retval))
				ret=FAIL;
		}
	}
#endif

	if(ret==OK && socket->differentClients){
		cout << "Force update" << endl;
		ret=FORCE_UPDATE;
	}

	// send answer
	socket->SendDataOnly(&ret,sizeof(ret));
	if(ret==FAIL)
		socket->SendDataOnly(mess,sizeof(mess));
	socket->SendDataOnly(&retval,sizeof(retval));

	//return ok/fail
	return ret;
}



int slsReceiverFuncs::get_version(){
	ret=OK;
	int64_t retval=-1;

	// execute action if the arguments correctly arrived
#ifdef SLS_RECEIVER_FUNCTION_LIST
	retval= SVNREV;
	retval= (retval <<32) | SVNDATE;
#endif

	if(socket->differentClients){
		cout << "Force update" << endl;
		ret=FORCE_UPDATE;
	}

	// send answer
	socket->SendDataOnly(&ret,sizeof(ret));
	socket->SendDataOnly(&retval,sizeof(retval));

	//return ok/fail
	return ret;
}



































int slsReceiverFuncs::lock_receiver() {
	ret=OK;
	int lock;

	// receive arguments
	if(socket->ReceiveDataOnly(&lock,sizeof(lock)) < 0 ){
		sprintf(mess,"Error reading from socket\n");
		cout << "Error reading from socket (lock)" << endl;
		ret=FAIL;
	}
	// execute action if the arguments correctly arrived
	if(ret==OK){
		if (lock>=0) {
			if (lockStatus==0 || strcmp(socket->lastClientIP,socket->thisClientIP)==0 || strcmp(socket->lastClientIP,"none")==0) {
				lockStatus=lock;
				strcpy(socket->lastClientIP,socket->thisClientIP);
			}   else {
				ret=FAIL;
				sprintf(mess,"Receiver already locked by %s\n", socket->lastClientIP);
			}
		}
	}

	if (socket->differentClients && ret==OK)
		ret=FORCE_UPDATE;

	// send answer
	socket->SendDataOnly(&ret,sizeof(ret));
	if (ret==FAIL)
		socket->SendDataOnly(mess,sizeof(mess));
	else
		socket->SendDataOnly(&lockStatus,sizeof(lockStatus));

	//return ok/fail
	return ret;
}







int slsReceiverFuncs::set_port() {
	ret=OK;
	MySocketTCP* mySocket=NULL;
	int sd=-1;
	enum portType p_type; /** data? control? stop? Unused! */
	int p_number; /** new port number */

	// receive arguments
	if(socket->ReceiveDataOnly(&p_type,sizeof(p_type)) < 0 ){
		strcpy(mess,"Error reading from socket\n");
		cout << mess << endl;
		ret=FAIL;
	}

	if(socket->ReceiveDataOnly(&p_number,sizeof(p_number)) < 0 ){
		strcpy(mess,"Error reading from socket\n");
		cout << mess << endl;
		ret=FAIL;
	}

	// execute action if the arguments correctly arrived
	if (ret==OK) {
		if (socket->differentClients==1 && lockStatus==1 ) {
			ret=FAIL;
			sprintf(mess,"Detector locked by %s\n",socket->lastClientIP);
		}
		else {
			if (p_number<1024) {
				sprintf(mess,"Too low port number %d\n", p_number);
				cout << mess << endl;
				ret=FAIL;
			}
			cout << "set port " << p_type << " to " << p_number <<endl;
			mySocket = new MySocketTCP(p_number);
		}
		if(mySocket){
			sd = mySocket->getErrorStatus();
			if (!sd){
				ret=OK;
				if (mySocket->differentClients)
					ret=FORCE_UPDATE;
			} else {
				ret=FAIL;
				sprintf(mess,"Could not bind port %d\n", p_number);
				cout << mess << endl;
				if (sd==-10) {
					sprintf(mess,"Port %d already set\n", p_number);
					cout << mess << endl;
				}
			}
		}
	}

	// send answer
	socket->SendDataOnly(&ret,sizeof(ret));
	if (ret==FAIL) {
		socket->SendDataOnly(mess,sizeof(mess));
	} else {
		socket->SendDataOnly(&p_number,sizeof(p_number));
		if(sd>=0){
			socket->Disconnect();
			delete socket;
			socket = mySocket;
			file_des=socket->getFileDes();
		}
	}

	//return ok/fail
	return ret;
}






int slsReceiverFuncs::get_last_client_ip() {
	ret=OK;

	if (socket->differentClients )
		ret=FORCE_UPDATE;

	socket->SendDataOnly(&ret,sizeof(ret));
	socket->SendDataOnly(socket->lastClientIP,sizeof(socket->lastClientIP));

	return ret;
}







int slsReceiverFuncs::send_update() {
	ret=OK;
	int ind;
	char path[MAX_STR_LENGTH];

	socket->SendDataOnly(socket->lastClientIP,sizeof(socket->lastClientIP));

	//index
#ifdef SLS_RECEIVER_FUNCTION_LIST
	ind=slsReceiverList->getFileIndex();
#endif
	socket->SendDataOnly(&ind,sizeof(ind));


	//filepath
#ifdef SLS_RECEIVER_FUNCTION_LIST
	strcpy(path,slsReceiverList->getFilePath());
#endif
	socket->SendDataOnly(path,MAX_STR_LENGTH);


	//filename
#ifdef SLS_RECEIVER_FUNCTION_LIST
	strcpy(path,slsReceiverList->getFileName());
#endif
	socket->SendDataOnly(path,MAX_STR_LENGTH);


	if (lockStatus==0) {
		strcpy(socket->lastClientIP,socket->thisClientIP);
	}

	return ret;


}






int slsReceiverFuncs::update_client() {
	ret=OK;
	socket->SendDataOnly(&ret,sizeof(ret));

	return send_update();
}







int slsReceiverFuncs::exit_server() {
	ret=GOODBYE;
	socket->SendDataOnly(&ret,sizeof(ret));
	strcpy(mess,"closing server");
	socket->SendDataOnly(mess,sizeof(mess));
	cout << mess << endl;
	return ret;
}





int slsReceiverFuncs::exec_command() {
	ret = OK;
	char cmd[MAX_STR_LENGTH];
	char answer[MAX_STR_LENGTH];
	int sysret=0;

	// receive arguments
	if(socket->ReceiveDataOnly(cmd,MAX_STR_LENGTH) < 0 ){
		sprintf(mess,"Error reading from socket\n");
		ret=FAIL;
	}

	// execute action if the arguments correctly arrived
	if (ret==OK) {
#ifdef VERBOSE
		cout << "executing command " << cmd << endl;
#endif
		if (lockStatus==0 || socket->differentClients==0)
			sysret=system(cmd);

		//should be replaced by popen
		if (sysret==0) {
			strcpy(answer,"Succeeded\n");
			if (lockStatus==1 && socket->differentClients==1)
				sprintf(answer,"Detector locked by %s\n", socket->lastClientIP);
		} else {
			strcpy(answer,"Failed\n");
			ret=FAIL;
		}
	} else
		strcpy(answer,"Could not receive the command\n");


	// send answer
	socket->SendDataOnly(&ret,sizeof(ret));
	if(socket->SendDataOnly(answer,MAX_STR_LENGTH) < 0){
		strcpy(mess,"Error writing to socket");
		ret=FAIL;
	}

	//return ok/fail
	return ret;
}





