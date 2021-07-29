#pragma once

#include "SharedMemory.h"
#include "sls/Result.h"
#include "sls/logger.h"
#include "sls/sls_detector_defs.h"

class ZmqSocket;
class detectorData;

#include <memory>
#include <mutex>
#include <semaphore.h>
#include <string>
#include <thread>
#include <vector>

#define MULTI_SHMAPIVERSION 0x190809
#define MULTI_SHMVERSION    0x201007
#define SHORT_STRING_LENGTH 50

#include <future>
#include <numeric>

namespace sls {

class Module;

/**
 * @short structure allocated in shared memory to store detector settings
 * for IPC and cache
 */
struct sharedMultiSlsDetector {

    /* FIXED PATTERN FOR STATIC FUNCTIONS. DO NOT CHANGE, ONLY APPEND
     * ------*/

    /** shared memory version */
    int shmversion;

    /** last process id accessing the shared memory */
    pid_t lastPID;

    /** last user name accessing the shared memory */
    char lastUser[SHORT_STRING_LENGTH];

    /** last time stamp when accessing the shared memory */
    char lastDate[SHORT_STRING_LENGTH];

    int numberOfDetectors;
    slsDetectorDefs::detectorType multiDetectorType;

    /** END OF FIXED PATTERN
     * -----------------------------------------------*/

    /** Number of detectors operated at once */
    slsDetectorDefs::xy numberOfDetector;

    /**  max number of channels for complete detector*/
    slsDetectorDefs::xy numberOfChannels;

    bool acquiringFlag;
    bool initialChecks;
    bool gapPixels;
    /** high water mark of listening tcp port (only data) */
    int zmqHwm;
};

class DetectorImpl : public virtual slsDetectorDefs {
  public:
    /**
     * Constructor
     * @param multi_id multi detector id
     * @param verify true to verify if shared memory version matches existing
     * one
     * @param update true to update last user pid, date etc
     */
    explicit DetectorImpl(int multi_id = 0, bool verify = true,
                          bool update = true);

    /**
     * Destructor
     */
    virtual ~DetectorImpl();

    template <class CT> struct NonDeduced { using type = CT; };
    template <typename RT, typename... CT>
    sls::Result<RT> Parallel(RT (sls::Module::*somefunc)(CT...),
                             std::vector<int> positions,
                             typename NonDeduced<CT>::type... Args) {

        if (detectors.empty())
            throw sls::RuntimeError("No detectors added");
        if (positions.empty() ||
            (positions.size() == 1 && positions[0] == -1)) {
            positions.resize(detectors.size());
            std::iota(begin(positions), end(positions), 0);
        }
        std::vector<std::future<RT>> futures;
        futures.reserve(positions.size());
        for (size_t i : positions) {
            if (i >= detectors.size())
                throw sls::RuntimeError("Detector out of range");
            futures.push_back(std::async(std::launch::async, somefunc,
                                         detectors[i].get(), Args...));
        }
        sls::Result<RT> result;
        result.reserve(positions.size());
        for (auto &i : futures) {
            result.push_back(i.get());
        }
        return result;
    }

    template <typename RT, typename... CT>
    sls::Result<RT> Parallel(RT (sls::Module::*somefunc)(CT...) const,
                             std::vector<int> positions,
                             typename NonDeduced<CT>::type... Args) const {

        if (detectors.empty())
            throw sls::RuntimeError("No detectors added");
        if (positions.empty() ||
            (positions.size() == 1 && positions[0] == -1)) {
            positions.resize(detectors.size());
            std::iota(begin(positions), end(positions), 0);
        }
        std::vector<std::future<RT>> futures;
        futures.reserve(positions.size());
        for (size_t i : positions) {
            if (i >= detectors.size())
                throw sls::RuntimeError("Detector out of range");
            futures.push_back(std::async(std::launch::async, somefunc,
                                         detectors[i].get(), Args...));
        }
        sls::Result<RT> result;
        result.reserve(positions.size());
        for (auto &i : futures) {
            result.push_back(i.get());
        }
        return result;
    }

    template <typename... CT>
    void Parallel(void (sls::Module::*somefunc)(CT...),
                  std::vector<int> positions,
                  typename NonDeduced<CT>::type... Args) {

        if (detectors.empty())
            throw sls::RuntimeError("No detectors added");
        if (positions.empty() ||
            (positions.size() == 1 && positions[0] == -1)) {
            positions.resize(detectors.size());
            std::iota(begin(positions), end(positions), 0);
        }
        std::vector<std::future<void>> futures;
        futures.reserve(positions.size());
        for (size_t i : positions) {
            if (i >= detectors.size())
                throw sls::RuntimeError("Detector out of range");
            futures.push_back(std::async(std::launch::async, somefunc,
                                         detectors[i].get(), Args...));
        }
        for (auto &i : futures) {
            i.get();
        }
    }

    template <typename... CT>
    void Parallel(void (sls::Module::*somefunc)(CT...) const,
                  std::vector<int> positions,
                  typename NonDeduced<CT>::type... Args) const {

        if (detectors.empty())
            throw sls::RuntimeError("No detectors added");
        if (positions.empty() ||
            (positions.size() == 1 && positions[0] == -1)) {
            positions.resize(detectors.size());
            std::iota(begin(positions), end(positions), 0);
        }
        std::vector<std::future<void>> futures;
        futures.reserve(positions.size());
        for (size_t i : positions) {
            if (i >= detectors.size())
                throw sls::RuntimeError("Detector out of range");
            futures.push_back(std::async(std::launch::async, somefunc,
                                         detectors[i].get(), Args...));
        }
        for (auto &i : futures) {
            i.get();
        }
    }

    /** set acquiring flag in shared memory */
    void setAcquiringFlag(bool flag);

    /** return multi detector shared memory ID */
    int getMultiId() const;

    /** Free specific shared memory from the command line without creating
     * object */
    static void freeSharedMemory(int multiId, int detPos = -1);

    /** Free all modules from current multi Id shared memory and delete members
     */
    void freeSharedMemory();

    /** Get user details of shared memory */
    std::string getUserDetails();

    bool getInitialChecks() const;

    /** initial compaibility and other server start up checks
     * default enabled */
    void setInitialChecks(const bool value);

    /**
     * Connect to Virtual Detector Servers at local host
     * @param numdet number of detectors
     * @param port starting port number
     */
    void setVirtualDetectorServers(const int numdet, const int port);

    /** Sets the hostname of all sls detectors in shared memory and updates
     * local cache */
    void setHostname(const std::vector<std::string> &name);

    /** Gets the total number of detectors */
    int size() const;

    slsDetectorDefs::xy getNumberOfDetectors() const;

    slsDetectorDefs::xy getNumberOfChannels() const;

    /** Must be set before setting hostname
     * Sets maximum number of channels of all sls detectors */
    void setNumberOfChannels(const slsDetectorDefs::xy c);

    /** [Eiger][Jungfrau] */
    bool getGapPixelsinCallback() const;
    /** [Eiger][Jungfrau] */
    void setGapPixelsinCallback(const bool enable);

    bool getDataStreamingToClient();
    void setDataStreamingToClient(bool enable);
    int getClientStreamingHwm() const;
    void setClientStreamingHwm(const int limit);

    /**
     * register callback for accessing acquisition final data
     * @param func function to be called at the end of the acquisition.
     * gets detector status and progress index as arguments
     * @param pArg argument
     */
    void registerAcquisitionFinishedCallback(void (*func)(double, int, void *),
                                             void *pArg);

    /**
     * register calbback for accessing detector final data,
     * also enables data streaming in client and receiver
     * @param userCallback function for plotting/analyzing the data.
     * Its arguments are
     * the data structure d and the frame number f,
     * s is for subframe number for eiger for 32 bit mode
     * @param pArg argument
     */
    void registerDataCallback(void (*userCallback)(detectorData *, uint64_t,
                                                   uint32_t, void *),
                              void *pArg);

    /**
     * Performs a complete acquisition
     * resets frames caught in receiver, starts receiver, starts detector,
     * blocks till detector finished acquisition, stop receiver, increments file
     * index, loops for measurements, calls required call backs.
     * @returns OK or FAIL depending on if it already started
     */
    int acquire();

    /**
     * Combines data from all readouts and gives it to the gui
     * or just gives progress of acquisition by polling receivers
     */
    void processData(bool receiver);

    /**
     * Convert raw file
     * [Jungfrau][Ctb] from pof file
     * [Mythen3][Gotthard2] from rbf file
     * @param fname name of pof/rbf file
     * @returns binary of the program
     */
    std::vector<char> readProgrammingFile(const std::string &fname);

    sls::Result<int> getNumberofUDPInterfaces(Positions pos) const;
    void setNumberofUDPInterfaces(int n, Positions pos);
    sls::Result<int> getDefaultDac(defs::dacIndex index,
                                   defs::detectorSettings sett,
                                   Positions pos = {});
    void setDefaultDac(defs::dacIndex index, int defaultValue,
                       defs::detectorSettings sett, Positions pos);

  private:
    /**
     * Creates/open shared memory, initializes detector structure and members
     * Called by constructor/ set hostname / read config file
     * @param verify true to verify if shared memory version matches existing
     * one
     * @param update true to update last user pid, date etc
     */
    void setupMultiDetector(bool verify = true, bool update = true);

    /**
     * Creates shm and initializes shm structure OR
     * Open shm and maps to structure
     * @param verify true to verify if shm size matches existing one
     */
    void initSharedMemory(bool verify = true);

    /** Initialize detector structure for the shared memory just created */
    void initializeDetectorStructure();

    /** Initialize members (eg. slsDetectors from shm, zmqsockets)
     * @param verify true to verify if shm size matches existing one
     */
    void initializeMembers(bool verify = true);

    /** Update in shm */
    void updateUserdetails();

    bool isAcquireReady();

    /** Execute command in terminal and return result */
    std::string exec(const char *cmd);

    void addSlsDetector(const std::string &hostname);

    void updateDetectorSize();

    int destroyReceivingDataSockets();
    int createReceivingDataSockets();

    /**
     * Reads frames from receiver through a constant socket
     * Called during acquire() when call back registered or when using gui
     */
    void readFrameFromReceiver();

    /** [Eiger][Jungfrau]
     * add gap pixels to the imag
     * @param image pointer to image without gap pixels
     * @param gpImage poiner to image with gap pixels, if NULL, allocated
     * @param quadEnable quad enabled
     * @param dr dynamic range
     * @param nPixelsx number of pixels in X axis (updated)
     * @param nPixelsy number of pixels in Y axis (updated)
     * @returns total data bytes for updated image
     */
    int InsertGapPixels(char *image, char *&gpImage, bool quadEnable, int dr,
                        int &nPixelsx, int &nPixelsy);

    void printProgress(double progress);

    void startProcessingThread(bool receiver);

    /**
     * Check if processing thread is ready to join main thread
     * @returns true if ready, else false
     */
    bool getJoinThreadFlag() const;

    /**
     * Main thread sets if the processing thread should join it
     * @param value true if it should join, else false
     */
    void setJoinThreadFlag(bool value);

    /**
     * Listen to key event to stop acquiring
     * when using acquire command
     */
    int kbhit();

    /** Multi detector Id */
    const int multiId{0};

    /** Shared Memory object */
    sls::SharedMemory<sharedMultiSlsDetector> multi_shm{0, -1};

    /** pointers to the Module structures */
    std::vector<std::unique_ptr<sls::Module>> detectors;

    /** data streaming (down stream) enabled in client (zmq sckets created) */
    bool client_downstream{false};

    /** ZMQ Socket - Receiver to Client */
    std::vector<std::unique_ptr<ZmqSocket>> zmqSocket;

    /** number of zmq sockets running currently */
    volatile int numZmqRunning{0};

    /** mutex to synchronize main and data processing threads */
    mutable std::mutex mp;

    /** sets when the acquisition is finished */
    bool jointhread{false};

    /** the data processing thread */
    std::thread dataProcessingThread;

    /** detector data packed for the gui */
    detectorData *thisData{nullptr};

    void (*acquisition_finished)(double, int, void *){nullptr};
    void *acqFinished_p{nullptr};

    void (*dataReady)(detectorData *, uint64_t, uint32_t, void *){nullptr};
    void *pCallbackArg{nullptr};
};

} // namespace sls