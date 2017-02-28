#pragma once
/************************************************
 * @file DataStreamer.h
 * @short streams data from receiver via ZMQ
 ***********************************************/
/**
 *@short creates & manages a data streamer thread each
 */

#include "ThreadObject.h"

class GeneralData;
class Fifo;
class DataStreamer;
class ZmqSocket;

class DataStreamer : private virtual slsReceiverDefs, public ThreadObject {
	
 public:
	/**
	 * Constructor
	 * Calls Base Class CreateThread(), sets ErrorMask if error and increments NumberofDataStreamers
	 * @param f address of Fifo pointer
	 * @param dr pointer to dynamic range
	 * @param freq poiner to streaming frequency
	 * @param timer poiner to timer if streaming frequency is random
	 */
	DataStreamer(Fifo*& f, uint32_t* dr, uint32_t* freq, uint32_t* timer);

	/**
	 * Destructor
	 * Calls Base Class DestroyThread() and decrements NumberofDataStreamers
	 */
	~DataStreamer();


	//*** static functions ***
	/**
	 * Get RunningMask
	 * @return RunningMask
	 */
	static uint64_t GetErrorMask();

	/**
	 * Get RunningMask
	 * @return RunningMask
	 */
	static uint64_t GetRunningMask();

	/**
	 * Reset RunningMask
	 */
	static void ResetRunningMask();

	//*** non static functions ***
	//*** getters ***



	//*** setters ***
	/**
	 * Set bit in RunningMask to allow thread to run
	 */
	void StartRunning();

	/**
	 * Reset bit in RunningMask to prevent thread from running
	 */
	void StopRunning();

	/**
	 * Set Fifo pointer to the one given
	 * @param f address of Fifo pointer
	 */
	void SetFifo(Fifo*& f);

	/**
	 * Reset parameters for new acquisition (including all scans)
	 */
	void ResetParametersforNewAcquisition();

	/**
	 * Reset parameters for new measurement (eg. for each scan)
	 */
	void ResetParametersforNewMeasurement();

	/**
	 * Set GeneralData pointer to the one given
	 * @param g address of GeneralData (Detector Data) pointer
	 */
	void SetGeneralData(GeneralData* g);

	/**
	 * Set thread priority
	 * @priority priority
	 * @returns OK or FAIL
	 */
	int SetThreadPriority(int priority);

	/**
	 * Creates Zmq Sockets
	 * @param dindex pointer to detector index
	 * @param nunits pointer to number of theads/ units per detector
	 * @return OK or FAIL
	 */
	int CreateZmqSockets(int* dindex, int* nunits);

	/**
	 * Shuts down and deletes Zmq Sockets
	 */
	void CloseZmqSocket();

 private:

	/**
	 * Get Type
	 * @return type
	 */
	std::string GetType();

	/**
	 * Returns if the thread is currently running
	 * @returns true if thread is running, else false
	 */
	bool IsRunning();

	/**
	 * Create Part1 of Json Header which includes common attributes in an acquisition
	 */
	void CreateHeaderPart1();

	/**
	 * Record First Indices (firstAcquisitionIndex, firstMeasurementIndex)
	 * @param fnum frame index to record
	 */
	void RecordFirstIndices(uint64_t fnum);

	/**
	 * Thread Exeution for DataStreamer Class
	 * Stream an image via zmq
	 */
	void ThreadExecution();

	/**
	 * Frees dummy buffer,
	 * reset running mask by calling StopRunning()
	 * @param buf address of pointer
	 */
	void StopProcessing(char* buf);

	/**
	 * Process an image popped from fifo,
	 * write to file if fw enabled & update parameters
	 * @param buffer
	 */
	void ProcessAnImage(char* buf);

	/**
	 * This function should be called only in random frequency mode
	 * Checks if timer is done and ready to send data
	 * @returns true if ready to send data, else false
	 */
	bool CheckTimer();

	/**
	 * This function should be called only in non random frequency mode
	 * Checks if count is done and ready to send data
	 * @returns true if ready to send data, else false
	 */
	bool CheckCount();

	/**
	 * Create and send Json Header
	 * @param fnum frame number
	 * @param snum sub frame number
	 * @param dummy true if its a dummy header
	 * @returns 0 if error, else 1
	 */
	int SendHeader(uint64_t fnum, uint32_t snum, bool dummy = false);

	/** type of thread */
	static const std::string TypeName;

	/** Total Number of DataStreamer Objects */
	static int NumberofDataStreamers;

	/** Mask of errors on any object eg.thread creation */
	static uint64_t ErrorMask;

	/** Mask of all listener objects running */
	static uint64_t RunningMask;

	/** mutex to update static items among objects (threads)*/
	static pthread_mutex_t Mutex;

	/** Json Header Format for each measurement part */
	static const char *jsonHeaderFormat_part1;

	/** Json Header Format */
	static const char *jsonHeaderFormat;

	/** GeneralData (Detector Data) object */
	const GeneralData* generalData;

	/** Fifo structure */
	Fifo* fifo;

	/** ZMQ Socket - Receiver to Client */
	ZmqSocket* zmqSocket;

	/** Pointer to dynamic range */
	uint32_t* dynamicRange;

	/** Pointer to Streaming frequency, if 0, sending random images with a timer */
	uint32_t* streamingFrequency;

	/** Pointer to the timer if Streaming frequency is random */
	uint32_t* streamingTimerInMs;

	/** Current frequency count */
	uint32_t currentFreqCount;

	/** timer beginning stamp for random streaming */
	struct timespec timerBegin;

	/** Current Json Header  prefix*/
	char* currentHeader;

	/** Aquisition Started flag */
	bool acquisitionStartedFlag;

	/** Measurement Started flag */
	bool measurementStartedFlag;

	/** Frame Number of First Frame of an entire Acquisition (including all scans) */
	uint64_t firstAcquisitionIndex;

	/** Frame Number of First Frame for each real time acquisition (eg. for each scan) */
	uint64_t firstMeasurementIndex;
};

