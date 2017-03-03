#pragma once
/************************************************
 * @file zmqSocket.h
 * @short functions to open/close zmq sockets
 ***********************************************/
/**
 *@short functions to open/close zmq sockets
 */

#include "ansi.h"

#include <zmq.h>
#include <errno.h>
#include <netdb.h>				//gethostbyname()
#include <arpa/inet.h>			//inet_ntoa
#include <rapidjson/document.h> //json header in zmq stream
using namespace rapidjson;

#define DEFAULT_ZMQ_PORTNO 	70001
#define DUMMY_MSG_SIZE 		3
#define DUMMY_MSG			"end"

class ZmqSocket {

public:

	//Socket Options for optimization
	//ZMQ_LINGER default is already -1 means no messages discarded. use this options if optimizing required
	//ZMQ_SNDHWM default is 0 means no limit. use this to optimize if optimizing required
	// eg. int value = -1;
	// if (zmq_setsockopt(socketDescriptor, ZMQ_LINGER, &value,sizeof(value))) {
	//	Close();
	// }

	/**
	 * Constructor for a client
	 * Creates socket, context and connects to server
	 * @param hostname hostname or ip of server
	 * @param portnumber port number
	 */
	ZmqSocket (const char* const hostname_or_ip, const uint32_t  portnumber):
		portno (portnumber),
		server (false),
		contextDescriptor (NULL),
		socketDescriptor (NULL)
	{
		char ip[MAX_STR_LENGTH] = "";
		strcpy(ip, hostname_or_ip);

		// construct address
		if (strchr (hostname_or_ip, '.') != NULL) {
			// convert hostname to ip
			char* ptr = ConvertHostnameToIp (hostname_or_ip);
			if (ptr == NULL)
				return;
			strcpy(ip, ptr);
			delete ptr;
		}
		sprintf (serverAddress, "tcp://%s:%d", ip, portno);

		// create context
		contextDescriptor = zmq_ctx_new();
		if (contextDescriptor == NULL)
			return;

		// create publisher
		socketDescriptor = zmq_socket (contextDescriptor, ZMQ_PULL);
		if (socketDescriptor == NULL) {
			PrintError ();
			Close ();
		}

		//Socket Options provided above

		//connect socket
		if (zmq_connect(socketDescriptor, serverAddress) < 0) {
			PrintError ();
			Close ();
		}
	};

	/**
	 * Constructor for a server
	 * Creates socket, context and connects to server
	 * @param hostname hostname or ip of server
	 * @param portnumber port number
	 */
	ZmqSocket (const uint32_t portnumber):
		portno (portnumber),
		server (true),
		contextDescriptor (NULL),
		socketDescriptor (NULL)
	{
		// create context
		contextDescriptor = zmq_ctx_new();
		if (contextDescriptor == NULL)
			return;
		// create publisher
		socketDescriptor = zmq_socket (contextDescriptor, ZMQ_PUSH);
		if (socketDescriptor == NULL) {
			PrintError ();
			Close ();
		}

		//Socket Options provided above

		// construct address
		sprintf (serverAddress,"tcp://*:%d", portno);
		// bind address
		if (zmq_bind (socketDescriptor, serverAddress) < 0) {
			PrintError ();
			Close ();
		}
	};

	/**
	 * Destructor
	 */
	~ZmqSocket () {
		Disconnect();
		Close();
	};

	/**
	 * Returns error status
	 * @returns true if error else false
	 */
	bool IsError() { if (socketDescriptor == NULL) return true; return false; };

	/**
	 * Returns Server Address
	 * @returns Server Address
	 */
	char* GetZmqServerAddress () { return serverAddress; };

	/**
	 * Returns Port Number
	 * @returns Port Number
	 */
	uint32_t GetPortNumber () { return portno; };

	/**
	 * Returns Socket Descriptor
	 * @reutns Socket descriptor
	 */

	void* GetsocketDescriptor () { return socketDescriptor; };

	/**
	 * Unbinds the Socket
	 */
	void Disconnect () {
		if (server)
			zmq_unbind (socketDescriptor, serverAddress);
		else
			zmq_disconnect (socketDescriptor, serverAddress);
	};

	/**
	 * Close Socket and destroy Context
	 */
	void Close () {
		if (socketDescriptor != NULL) {
			zmq_close (socketDescriptor);
			socketDescriptor = NULL;
		}
		if (contextDescriptor != NULL) {
			zmq_ctx_destroy (contextDescriptor);
			contextDescriptor = NULL;
		}
	};

	/**
	 * Convert Hostname to ip
	 * @param hostname hostname
	 * @returns string with ip or NULL if error
	 */
	char* ConvertHostnameToIp (const char* const hostname) {
		struct hostent *he = gethostbyname (hostname);
		if (he == NULL){
			cprintf (RED,"Error: Could not convert hostname to ip (zmq)\n");
			return NULL;
		}
		return inet_ntoa (*(struct in_addr*)he->h_addr);
	};

	/**
	 * Send Message Header
	 * @returns 0 if error, else 1
	 */
	int SendHeaderData (char* buf, int length) {
		if(zmq_send (socketDescriptor, buf, length, ZMQ_SNDMORE) < 0) {
			PrintError ();
			return 0;
		}
		return 1;
	};

	/**
	 * Send Message Body
	 * @returns 0 if error, else 1
	 */
	int SendData (char* buf, int length) {
		if(zmq_send (socketDescriptor, buf, length, 0) < 0) {
			PrintError ();
			return 0;
		}
		return 1;
	};


	/**
	 * Receive Message
	 * @param index self index for debugging
	 * @returns length of message, -1 if error
	 */
	int ReceiveMessage(const int index) {
		int length = zmq_msg_recv (&message, socketDescriptor, 0);
		if (length == -1) {
			PrintError ();
			cprintf (BG_RED,"Error: Could not read header for socket %d\n",index);
		}
		return length;
	};


	/**
	 * Receive Header
	 * @param index self index for debugging
	 * @param acqIndex address of acquisition index
	 * @param frameIndex address of frame index
	 * @param subframeIndex address of subframe index
	 * @param filename address of file name
	 * @returns 0 if error, else 1
	 */
	int ReceiveHeader(const int index, uint64_t &acqIndex,
			uint64_t &frameIndex, uint32_t &subframeIndex, string &filename)
	{
		zmq_msg_init (&message);
		if (ReceiveMessage(index) > 0) {
			if (ParseHeader(index, acqIndex, frameIndex, subframeIndex, filename)) {
				zmq_msg_close(&message);
#ifdef VERBOSE
				cprintf(BLUE,"%d header rxd\n",index);
#endif
				return 1;
			}
		}
		zmq_msg_close(&message);
		return 0;
	};

	/**
	 * Receive Data
	 * @param index self index for debugging
	 * @param buf buffer to copy image data to
	 * @param size size of image
	 * @returns 0 if error, else 1
	 */
	int ReceiveData(const int index, int* buf, const int size)
	{
		zmq_msg_init (&message);
		int length = ReceiveMessage(index);

		//dummy
		if (length == DUMMY_MSG_SIZE) {
#ifdef VERBOSE
			cprintf(RED,"%d Received end of acquisition\n", index);
#endif
			zmq_msg_close(&message);
			return 0;
		}

		//actual data
		if (length == size) {
#ifdef VERBOSE
			cprintf(BLUE,"%d actual data\n", index);
#endif
			memcpy((char*)buf, (char*)zmq_msg_data(&message), size);
		}

		//incorrect size
		else {
			cprintf(RED,"Error: Received weird packet size %d for socket %d\n", length, index);
			memset((char*)buf,0xFF,size);
		}

		zmq_msg_close(&message);
		return 1;
	};


	/**
	 * Parse Header
	 * @param index self index for debugging
	 * @param acqIndex address of acquisition index
	 * @param frameIndex address of frame index
	 * @param subframeIndex address of subframe index
	 * @param filename address of file name
	 */
	int ParseHeader(const int index, uint64_t &acqIndex,
			uint64_t &frameIndex, uint32_t &subframeIndex, string &filename)
	{
		Document d;
		if (d.Parse( (char*)zmq_msg_data(&message), zmq_msg_size(&message)).HasParseError()) {
			cprintf (RED,"Error: Could not parse header for socket %d\n",index);
			return 0;
		}
#ifdef VERYVERBOSE
		printf("version:%.1f\n", d["version"].GetDouble());

		// shape is an array of ints
		rapidjson::Value::Array shape = d["shape"].GetArray();
		printf("%d: shape: ", index);
		for (int i = 0; i < shape.Size(); i++)
			printf("%d: %d ", index, shape[i].GetInt());
		printf("\n");

		printf("%d: type: %s\n", index, d["type"].GetString());
#endif

		if(d["acqIndex"].GetUint64()!=-1){
			acqIndex 		= d["acqIndex"].GetUint64();
			frameIndex 		= d["fIndex"].GetUint64();
			subframeIndex 	= d["subfnum"].GetUint();
			filename 		= d["fname"].GetString();
#ifdef VERYVERBOSE
			cout << "Acquisition index: " << acqIndex << endl;
			cout << "Frame index: " << frameIndex << endl;
			cout << "Subframe index: " << subframeIndex << endl;
			cout << "File name: " << filename << endl;
#endif
		}
		return 1;
	};


	/**
	 * Print error
	 */
	void PrintError () {
		switch (errno) {
		case EINVAL:
			cprintf(RED, "Error: The socket type/option or value/endpoint supplied is invalid (zmq)\n");
			break;
		case EAGAIN:
			cprintf(RED, "Error: Non-blocking mode was requested and the message cannot be sent/available at the moment (zmq)\n");
			break;
		case ENOTSUP:
			cprintf(RED, "Error: The zmq_send()/zmq_msg_recv() operation is not supported by this socket type (zmq)\n");
			break;
		case EFSM:
			cprintf(RED, "Error: The zmq_send()/zmq_msg_recv() unavailable now as socket in inappropriate state (eg. ZMQ_REP). Look up messaging patterns (zmq)\n");
			break;
		case EFAULT:
			cprintf(RED, "Error: The provided context/message is invalid (zmq)\n");
			break;
		case EMFILE:
			cprintf(RED, "Error: The limit on the total number of open ØMQ sockets has been reached (zmq)\n");
			break;
		case EPROTONOSUPPORT:
			cprintf(RED, "Error: The requested transport protocol is not supported (zmq)\n");
			break;
		case ENOCOMPATPROTO:
			cprintf(RED, "Error: The requested transport protocol is not compatible with the socket type (zmq)\n");
			break;
		case EADDRINUSE:
			cprintf(RED, "Error: The requested address is already in use (zmq)\n");
			break;
		case EADDRNOTAVAIL:
			cprintf(RED, "Error: The requested address was not local (zmq)\n");
			break;
		case ENODEV:
			cprintf(RED, "Error: The requested address specifies a nonexistent interface (zmq)\n");
			break;
		case ETERM:
			cprintf(RED, "Error: The ØMQ context associated with the specified socket was terminated (zmq)\n");
			break;
		case ENOTSOCK:
			cprintf(RED, "Error: The provided socket was invalid (zmq)\n");
			break;
		case EINTR:
			cprintf(RED, "Error: The operation was interrupted by delivery of a signal (zmq)\n");
			break;
		case EMTHREAD:
			cprintf(RED, "Error: No I/O thread is available to accomplish the task (zmq)\n");
			break;
		default:
			cprintf(RED, "Error: Unknown socket error (zmq)\n");
			break;
		}
	};


private:
	/** Port Number */
	uint32_t portno;

	/** true if server, else false */
	bool server;

	/** Context Descriptor */
	void* contextDescriptor;

	/** Socket Descriptor */
	void* socketDescriptor;

	/** Server Address */
	char serverAddress[1000];

	/** Zmq Message */
	zmq_msg_t message;
};
