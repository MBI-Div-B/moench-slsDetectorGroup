// SPDX-License-Identifier: LGPL-3.0-or-other
// Copyright (C) 2021 Contributors to the SLS Detector Package
#include "ClientInterface.h"
#include "sls/ServerSocket.h"
#include "sls/StaticVector.h"
#include "sls/ToString.h"

#include "sls/sls_detector_exceptions.h"
#include "sls/string_utils.h"
#include "sls/versionAPI.h"

#include <array>
#include <chrono>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <limits.h>
#include <map>
#include <memory>
#include <sstream>
#include <string>
#include <sys/syscall.h>
#include <unistd.h>
#include <vector>

using ns = std::chrono::nanoseconds;
using sls::RuntimeError;
using sls::SocketError;
using Interface = sls::ServerInterface;

ClientInterface::~ClientInterface() {
    killTcpThread = true;
    LOG(logINFO) << "Shutting down TCP Socket on port " << portNumber;
    server.shutdown();
    LOG(logDEBUG) << "TCP Socket closed on port " << portNumber;
    tcpThread->join();
}

ClientInterface::ClientInterface(int portNumber)
    : detType(GOTTHARD),
      portNumber(portNumber > 0 ? portNumber : DEFAULT_PORTNO + 2),
      server(portNumber) {
    functionTable();
    parentThreadId = syscall(SYS_gettid);
    tcpThread =
        sls::make_unique<std::thread>(&ClientInterface::startTCPServer, this);
}

int64_t ClientInterface::getReceiverVersion() { return APIRECEIVER; }

/***callback functions***/
void ClientInterface::registerCallBackStartAcquisition(
    int (*func)(std::string, std::string, uint64_t, uint32_t, void *),
    void *arg) {
    startAcquisitionCallBack = func;
    pStartAcquisition = arg;
}

void ClientInterface::registerCallBackAcquisitionFinished(void (*func)(uint64_t,
                                                                       void *),
                                                          void *arg) {
    acquisitionFinishedCallBack = func;
    pAcquisitionFinished = arg;
}

void ClientInterface::registerCallBackRawDataReady(
    void (*func)(char *, char *, uint32_t, void *), void *arg) {
    rawDataReadyCallBack = func;
    pRawDataReady = arg;
}

void ClientInterface::registerCallBackRawDataModifyReady(
    void (*func)(char *, char *, uint32_t &, void *), void *arg) {
    rawDataModifyReadyCallBack = func;
    pRawDataReady = arg;
}

void ClientInterface::startTCPServer() {
    tcpThreadId = syscall(SYS_gettid);
    LOG(logINFOBLUE) << "Created [ TCP server Tid: " << tcpThreadId << "]";
    LOG(logINFO) << "SLS Receiver starting TCP Server on port " << portNumber
                 << '\n';

    while (!killTcpThread) {
        LOG(logDEBUG1) << "Start accept loop";
        try {
            auto socket = server.accept();
            try {
                verifyLock(); // lock should be checked only for set (not get),
                              // Move it back?
                ret = decodeFunction(socket);
            } catch (const RuntimeError &e) {
                // We had an error needs to be sent to client
                char mess[MAX_STR_LENGTH]{};
                sls::strcpy_safe(mess, e.what());
                socket.Send(FAIL);
                socket.Send(mess);
            }
            // if tcp command was to exit server
            if (ret == GOODBYE) {
                break;
            }
        } catch (const RuntimeError &e) {
            LOG(logERROR) << "Accept failed";
        }
    }

    if (receiver) {
        receiver->shutDownUDPSockets();
    }
    LOG(logINFOBLUE) << "Exiting [ TCP server Tid: " << tcpThreadId << "]";
}

// clang-format off
int ClientInterface::functionTable(){
	flist[F_LOCK_RECEIVER]					=	&ClientInterface::lock_receiver;
	flist[F_GET_LAST_RECEIVER_CLIENT_IP]	=	&ClientInterface::get_last_client_ip;
	flist[F_GET_RECEIVER_VERSION]			=	&ClientInterface::get_version;
	flist[F_SETUP_RECEIVER]				    =	&ClientInterface::setup_receiver;
	flist[F_RECEIVER_SET_ROI]				=	&ClientInterface::set_roi;
	flist[F_RECEIVER_SET_NUM_FRAMES]        =   &ClientInterface::set_num_frames;  
	flist[F_SET_RECEIVER_NUM_TRIGGERS]      =   &ClientInterface::set_num_triggers;           
	flist[F_SET_RECEIVER_NUM_BURSTS]        =   &ClientInterface::set_num_bursts;         
	flist[F_SET_RECEIVER_NUM_ADD_STORAGE_CELLS] = &ClientInterface::set_num_add_storage_cells;                 
	flist[F_SET_RECEIVER_TIMING_MODE]       =   &ClientInterface::set_timing_mode;          
	flist[F_SET_RECEIVER_BURST_MODE]        =   &ClientInterface::set_burst_mode;  
	flist[F_RECEIVER_SET_NUM_ANALOG_SAMPLES]=   &ClientInterface::set_num_analog_samples;            
	flist[F_RECEIVER_SET_NUM_DIGITAL_SAMPLES]=  &ClientInterface::set_num_digital_samples;           
	flist[F_RECEIVER_SET_EXPTIME]           =   &ClientInterface::set_exptime;
	flist[F_RECEIVER_SET_PERIOD]            =   &ClientInterface::set_period;
	flist[F_RECEIVER_SET_SUB_EXPTIME]       =   &ClientInterface::set_subexptime;    
	flist[F_RECEIVER_SET_SUB_DEADTIME]      =   &ClientInterface::set_subdeadtime;    
	flist[F_SET_RECEIVER_DYNAMIC_RANGE]		= 	&ClientInterface::set_dynamic_range;
	flist[F_SET_RECEIVER_STREAMING_FREQUENCY] = 	&ClientInterface::set_streaming_frequency;
	flist[F_GET_RECEIVER_STREAMING_FREQUENCY] = 	&ClientInterface::get_streaming_frequency;
	flist[F_GET_RECEIVER_STATUS]			=	&ClientInterface::get_status;
	flist[F_START_RECEIVER]					=	&ClientInterface::start_receiver;
	flist[F_STOP_RECEIVER]					=	&ClientInterface::stop_receiver;
	flist[F_SET_RECEIVER_FILE_PATH]			=	&ClientInterface::set_file_dir;
	flist[F_GET_RECEIVER_FILE_PATH]			=	&ClientInterface::get_file_dir;
	flist[F_SET_RECEIVER_FILE_NAME]			=	&ClientInterface::set_file_name;
	flist[F_GET_RECEIVER_FILE_NAME]			=	&ClientInterface::get_file_name;
	flist[F_SET_RECEIVER_FILE_INDEX]		=	&ClientInterface::set_file_index;
	flist[F_GET_RECEIVER_FILE_INDEX]		=	&ClientInterface::get_file_index;
	flist[F_GET_RECEIVER_FRAME_INDEX]		=	&ClientInterface::get_frame_index;    
	flist[F_GET_RECEIVER_FRAMES_CAUGHT]		=	&ClientInterface::get_frames_caught;
    flist[F_GET_NUM_MISSING_PACKETS]		=	&ClientInterface::get_missing_packets;
	flist[F_SET_RECEIVER_FILE_WRITE]		=	&ClientInterface::set_file_write;
	flist[F_GET_RECEIVER_FILE_WRITE]		=	&ClientInterface::get_file_write;
	flist[F_SET_RECEIVER_MASTER_FILE_WRITE]	=	&ClientInterface::set_master_file_write;
	flist[F_GET_RECEIVER_MASTER_FILE_WRITE]	=	&ClientInterface::get_master_file_write;
	flist[F_SET_RECEIVER_OVERWRITE]		    = 	&ClientInterface::set_overwrite;
	flist[F_GET_RECEIVER_OVERWRITE]		    = 	&ClientInterface::get_overwrite;
	flist[F_ENABLE_RECEIVER_TEN_GIGA]		= 	&ClientInterface::enable_tengiga;
	flist[F_SET_RECEIVER_FIFO_DEPTH]		= 	&ClientInterface::set_fifo_depth;
	flist[F_RECEIVER_ACTIVATE]				= 	&ClientInterface::set_activate;
	flist[F_SET_RECEIVER_STREAMING]		    = 	&ClientInterface::set_streaming;
	flist[F_GET_RECEIVER_STREAMING]		    = 	&ClientInterface::get_streaming;
	flist[F_RECEIVER_STREAMING_TIMER]		= 	&ClientInterface::set_streaming_timer;
	flist[F_GET_FLIP_ROWS_RECEIVER]		    = 	&ClientInterface::get_flip_rows;
	flist[F_SET_FLIP_ROWS_RECEIVER]		    = 	&ClientInterface::set_flip_rows;
	flist[F_SET_RECEIVER_FILE_FORMAT]		= 	&ClientInterface::set_file_format;
	flist[F_GET_RECEIVER_FILE_FORMAT]		= 	&ClientInterface::get_file_format;
	flist[F_SET_RECEIVER_STREAMING_PORT]	= 	&ClientInterface::set_streaming_port;
	flist[F_GET_RECEIVER_STREAMING_PORT]	= 	&ClientInterface::get_streaming_port;
	flist[F_SET_RECEIVER_STREAMING_SRC_IP]	= 	&ClientInterface::set_streaming_source_ip;
	flist[F_GET_RECEIVER_STREAMING_SRC_IP]	= 	&ClientInterface::get_streaming_source_ip;
	flist[F_SET_RECEIVER_SILENT_MODE]		= 	&ClientInterface::set_silent_mode;
	flist[F_GET_RECEIVER_SILENT_MODE]		= 	&ClientInterface::get_silent_mode;
	flist[F_RESTREAM_STOP_FROM_RECEIVER]	= 	&ClientInterface::restream_stop;
	flist[F_SET_ADDITIONAL_JSON_HEADER]     =   &ClientInterface::set_additional_json_header;
	flist[F_GET_ADDITIONAL_JSON_HEADER]     =   &ClientInterface::get_additional_json_header;
    flist[F_RECEIVER_UDP_SOCK_BUF_SIZE]     =   &ClientInterface::set_udp_socket_buffer_size;
    flist[F_RECEIVER_REAL_UDP_SOCK_BUF_SIZE]=   &ClientInterface::get_real_udp_socket_buffer_size;
    flist[F_SET_RECEIVER_FRAMES_PER_FILE]	=   &ClientInterface::set_frames_per_file;
    flist[F_GET_RECEIVER_FRAMES_PER_FILE]	=   &ClientInterface::get_frames_per_file;
    flist[F_RECEIVER_CHECK_VERSION]			=   &ClientInterface::check_version_compatibility;
    flist[F_SET_RECEIVER_DISCARD_POLICY]	=   &ClientInterface::set_discard_policy;
    flist[F_GET_RECEIVER_DISCARD_POLICY]	=   &ClientInterface::get_discard_policy;
	flist[F_SET_RECEIVER_PADDING]		    =   &ClientInterface::set_padding_enable;
	flist[F_GET_RECEIVER_PADDING]		    =   &ClientInterface::get_padding_enable;
	flist[F_RECEIVER_SET_READOUT_MODE] 	    = 	&ClientInterface::set_readout_mode;
	flist[F_RECEIVER_SET_ADC_MASK]			=	&ClientInterface::set_adc_mask;
	flist[F_SET_RECEIVER_DBIT_LIST]			=	&ClientInterface::set_dbit_list;
	flist[F_GET_RECEIVER_DBIT_LIST]			=	&ClientInterface::get_dbit_list;
	flist[F_SET_RECEIVER_DBIT_OFFSET]		= 	&ClientInterface::set_dbit_offset;
	flist[F_GET_RECEIVER_DBIT_OFFSET]		= 	&ClientInterface::get_dbit_offset;
    flist[F_SET_RECEIVER_QUAD]			    = 	&ClientInterface::set_quad_type;
    flist[F_SET_RECEIVER_READ_N_ROWS]       =   &ClientInterface::set_read_n_rows;
    flist[F_SET_RECEIVER_UDP_IP]            =   &ClientInterface::set_udp_ip;
	flist[F_SET_RECEIVER_UDP_IP2]           =   &ClientInterface::set_udp_ip2;
	flist[F_SET_RECEIVER_UDP_PORT]          =   &ClientInterface::set_udp_port;
	flist[F_SET_RECEIVER_UDP_PORT2]         =   &ClientInterface::set_udp_port2;
	flist[F_SET_RECEIVER_NUM_INTERFACES]    =   &ClientInterface::set_num_interfaces;
	flist[F_RECEIVER_SET_ADC_MASK_10G]		=	&ClientInterface::set_adc_mask_10g;
    flist[F_RECEIVER_SET_COUNTER_MASK]      =   &ClientInterface::set_counter_mask;
    flist[F_INCREMENT_FILE_INDEX]           =   &ClientInterface::increment_file_index;
    flist[F_SET_ADDITIONAL_JSON_PARAMETER]  =   &ClientInterface::set_additional_json_parameter;
	flist[F_GET_ADDITIONAL_JSON_PARAMETER]  =   &ClientInterface::get_additional_json_parameter;
	flist[F_GET_RECEIVER_PROGRESS]          =   &ClientInterface::get_progress;
    flist[F_SET_RECEIVER_NUM_GATES]         =   &ClientInterface::set_num_gates;    
    flist[F_SET_RECEIVER_GATE_DELAY]        =   &ClientInterface::set_gate_delay;        
    flist[F_GET_RECEIVER_THREAD_IDS]        =   &ClientInterface::get_thread_ids;
    flist[F_GET_RECEIVER_STREAMING_START_FNUM] = &ClientInterface::get_streaming_start_fnum;
    flist[F_SET_RECEIVER_STREAMING_START_FNUM] = &ClientInterface::set_streaming_start_fnum;
    flist[F_SET_RECEIVER_RATE_CORRECT]      =   &ClientInterface::set_rate_correct;
    flist[F_SET_RECEIVER_SCAN]              =   &ClientInterface::set_scan;
    flist[F_RECEIVER_SET_THRESHOLD]         =   &ClientInterface::set_threshold;
    flist[F_GET_RECEIVER_STREAMING_HWM]     =   &ClientInterface::get_streaming_hwm;
    flist[F_SET_RECEIVER_STREAMING_HWM]     =   &ClientInterface::set_streaming_hwm;
    flist[F_RECEIVER_SET_ALL_THRESHOLD]     =   &ClientInterface::set_all_threshold;
    flist[F_RECEIVER_SET_DATASTREAM]        =   &ClientInterface::set_detector_datastream;
    

	for (int i = NUM_DET_FUNCTIONS + 1; i < NUM_REC_FUNCTIONS ; i++) {
		LOG(logDEBUG1) << "function fnum: " << i << " (" <<
				getFunctionNameFromEnum((enum detFuncs)i) << ") located at " << flist[i];
	}

	return OK;
}
// clang-format on
int ClientInterface::decodeFunction(Interface &socket) {
    ret = FAIL;
    socket.Receive(fnum);
    socket.setFnum(fnum);
    if (fnum <= NUM_DET_FUNCTIONS || fnum >= NUM_REC_FUNCTIONS) {
        throw RuntimeError("Unrecognized Function enum " +
                           std::to_string(fnum) + "\n");
    } else {
        LOG(logDEBUG1) << "calling function fnum: " << fnum << " ("
                       << getFunctionNameFromEnum((enum detFuncs)fnum) << ")";
        ret = (this->*flist[fnum])(socket);
        LOG(logDEBUG1) << "Function "
                       << getFunctionNameFromEnum((enum detFuncs)fnum)
                       << " finished";
    }
    return ret;
}

void ClientInterface::functionNotImplemented() {
    std::ostringstream os;
    os << "Function: " << getFunctionNameFromEnum((enum detFuncs)fnum)
       << ", is is not implemented for this detector";
    throw RuntimeError(os.str());
}

void ClientInterface::modeNotImplemented(const std::string &modename,
                                         int mode) {
    std::ostringstream os;
    os << modename << " (" << mode << ") is not implemented for this detector";
    throw RuntimeError(os.str());
}

template <typename T>
void ClientInterface::validate(T arg, T retval, const std::string &modename,
                               numberMode hex) {
    if (ret == OK && arg != -1 && retval != arg) {
        auto format = (hex == HEX) ? std::hex : std::dec;
        auto prefix = (hex == HEX) ? "0x" : "";
        std::ostringstream os;
        os << "Could not " << modename << ". Set " << prefix << format << arg
           << ", but read " << prefix << retval << '\n';
        throw RuntimeError(os.str());
    }
}

void ClientInterface::verifyLock() {
    if (lockedByClient && server.getThisClient() != server.getLockedBy()) {
        throw sls::SocketError("Receiver locked\n");
    }
}

void ClientInterface::verifyIdle(Interface &socket) {
    if (impl()->getStatus() != IDLE) {
        std::ostringstream oss;
        oss << "Can not execute "
            << getFunctionNameFromEnum((enum detFuncs)fnum)
            << " when receiver is not idle";
        throw sls::SocketError(oss.str());
    }
}

int ClientInterface::lock_receiver(Interface &socket) {
    auto lock = socket.Receive<int>();
    LOG(logDEBUG1) << "Locking Server to " << lock;
    if (lock >= 0) {
        if (!lockedByClient ||
            (server.getLockedBy() == server.getThisClient())) {
            lockedByClient = lock;
            lock ? server.setLockedBy(server.getThisClient())
                 : server.setLockedBy(sls::IpAddr{});
        } else {
            throw RuntimeError("Receiver locked\n");
        }
    }
    return socket.sendResult(lockedByClient);
}

int ClientInterface::get_last_client_ip(Interface &socket) {
    return socket.sendResult(server.getLastClient());
}

int ClientInterface::get_version(Interface &socket) {
    return socket.sendResult(getReceiverVersion());
}

int ClientInterface::setup_receiver(Interface &socket) {
    auto arg = socket.Receive<rxParameters>();
    LOG(logDEBUG) << sls::ToString(arg);

    // if object exists, verify unlocked and idle, else only verify lock
    // (connecting first time)
    if (receiver != nullptr) {
        verifyIdle(socket);
    }

    // basic setup
    setDetectorType(arg.detType);
    impl()->setDetectorSize(arg.numberOfModule);
    impl()->setModulePositionId(arg.moduleIndex);
    impl()->setDetectorHostname(arg.hostname);

    // udp setup
    // update retvals only if detmac is not the same as in detector
    sls::MacAddr retvals[2];
    if (arg.udp_dstip != 0) {
        sls::MacAddr r = setUdpIp(sls::IpAddr(arg.udp_dstip));
        sls::MacAddr detMac{arg.udp_dstmac};
        if (detMac != r) {
            retvals[0] = r;
        }
    }
    if (arg.udp_dstip2 != 0) {
        sls::MacAddr r = setUdpIp2(sls::IpAddr(arg.udp_dstip2));
        sls::MacAddr detMac{arg.udp_dstmac2};
        if (detMac != r) {
            retvals[1] = r;
        }
    }
    impl()->setUDPPortNumber(arg.udp_dstport);
    impl()->setUDPPortNumber2(arg.udp_dstport2);
    if (detType == JUNGFRAU || detType == GOTTHARD2) {
        try {
            impl()->setNumberofUDPInterfaces(arg.udpInterfaces);
        } catch (const RuntimeError &e) {
            throw RuntimeError("Failed to set number of interfaces to " +
                               std::to_string(arg.udpInterfaces));
        }
    }
    impl()->setUDPSocketBufferSize(0);

    // acquisition parameters
    impl()->setNumberOfFrames(arg.frames);
    impl()->setNumberOfTriggers(arg.triggers);
    if (detType == GOTTHARD2) {
        impl()->setNumberOfBursts(arg.bursts);
    }
    if (detType == JUNGFRAU) {
        impl()->setNumberOfAdditionalStorageCells(arg.additionalStorageCells);
    }
    if (detType == MOENCH || detType == CHIPTESTBOARD) {
        try {
            impl()->setNumberofAnalogSamples(arg.analogSamples);
        } catch (const RuntimeError &e) {
            throw RuntimeError("Could not set num analog samples to " +
                               std::to_string(arg.analogSamples) +
                               " due to fifo structure memory allocation.");
        }
    }
    if (detType == CHIPTESTBOARD) {
        try {
            impl()->setNumberofDigitalSamples(arg.digitalSamples);
        } catch (const RuntimeError &e) {
            throw RuntimeError("Could not set num digital samples to " +
                               std::to_string(arg.analogSamples) +
                               " due to fifo structure memory allocation.");
        }
    }
    if (detType != MYTHEN3) {
        impl()->setAcquisitionTime(std::chrono::nanoseconds(arg.expTimeNs));
    }
    impl()->setAcquisitionPeriod(std::chrono::nanoseconds(arg.periodNs));
    if (detType == EIGER) {
        impl()->setSubExpTime(std::chrono::nanoseconds(arg.subExpTimeNs));
        impl()->setSubPeriod(std::chrono::nanoseconds(arg.subExpTimeNs) +
                             std::chrono::nanoseconds(arg.subDeadTimeNs));
        impl()->setActivate(static_cast<bool>(arg.activate));
        impl()->setDetectorDataStream(LEFT, arg.dataStreamLeft);
        impl()->setDetectorDataStream(RIGHT, arg.dataStreamRight);
        try {
            impl()->setQuad(arg.quad == 0 ? false : true);
        } catch (const RuntimeError &e) {
            throw RuntimeError("Could not set quad to " +
                               std::to_string(arg.quad) +
                               " due to fifo strucutre memory allocation");
        }
        impl()->setThresholdEnergy(arg.thresholdEnergyeV[0]);
    }
    if (detType == EIGER || detType == JUNGFRAU) {
        impl()->setReadNRows(arg.readNRows);
    }
    if (detType == MYTHEN3) {
        std::array<int, 3> val;
        for (int i = 0; i < 3; ++i) {
            val[i] = arg.thresholdEnergyeV[i];
        }
        impl()->setThresholdEnergy(val);
    }
    if (detType == EIGER || detType == MYTHEN3) {
        try {
            impl()->setDynamicRange(arg.dynamicRange);
        } catch (const RuntimeError &e) {
            throw RuntimeError(
                "Could not set dynamic range. Could not allocate "
                "memory for fifo or could not start listening/writing threads");
        }
    }
    impl()->setTimingMode(arg.timMode);
    if (detType == EIGER || detType == MOENCH || detType == CHIPTESTBOARD ||
        detType == MYTHEN3) {
        try {
            impl()->setTenGigaEnable(arg.tenGiga);
        } catch (const RuntimeError &e) {
            throw RuntimeError("Could not set 10GbE.");
        }
    }
    if (detType == CHIPTESTBOARD) {
        try {
            impl()->setReadoutMode(arg.roMode);
        } catch (const RuntimeError &e) {
            throw RuntimeError("Could not set read out mode "
                               "due to fifo memory allocation.");
        }
    }
    if (detType == CHIPTESTBOARD || detType == MOENCH) {
        try {
            impl()->setADCEnableMask(arg.adcMask);
        } catch (const RuntimeError &e) {
            throw RuntimeError("Could not set adc enable mask "
                               "due to fifo memory allcoation");
        }
        try {
            impl()->setTenGigaADCEnableMask(arg.adc10gMask);
        } catch (const RuntimeError &e) {
            throw RuntimeError("Could not set 10Gb adc enable mask "
                               "due to fifo memory allcoation");
        }
    }
    if (detType == GOTTHARD) {
        try {
            impl()->setROI(arg.roi);
        } catch (const RuntimeError &e) {
            throw RuntimeError("Could not set ROI");
        }
    }
    if (detType == MYTHEN3) {
        impl()->setCounterMask(arg.countermask);
        impl()->setAcquisitionTime1(std::chrono::nanoseconds(arg.expTime1Ns));
        impl()->setAcquisitionTime2(std::chrono::nanoseconds(arg.expTime2Ns));
        impl()->setAcquisitionTime3(std::chrono::nanoseconds(arg.expTime3Ns));
        impl()->setGateDelay1(std::chrono::nanoseconds(arg.gateDelay1Ns));
        impl()->setGateDelay2(std::chrono::nanoseconds(arg.gateDelay2Ns));
        impl()->setGateDelay3(std::chrono::nanoseconds(arg.gateDelay3Ns));
        impl()->setNumberOfGates(arg.gates);
    }
    if (detType == GOTTHARD2) {
        impl()->setBurstMode(arg.burstType);
    }
    impl()->setScan(arg.scanParams);

    return socket.sendResult(retvals);
}

void ClientInterface::setDetectorType(detectorType arg) {
    switch (arg) {
    case GOTTHARD:
    case EIGER:
    case CHIPTESTBOARD:
    case MOENCH:
    case JUNGFRAU:
    case MYTHEN3:
    case GOTTHARD2:
        break;
    default:
        throw RuntimeError("Unknown detector type: " + std::to_string(arg));
        break;
    }

    try {
        detType = GENERIC;
        receiver = sls::make_unique<Implementation>(arg);
        detType = arg;
    } catch (...) {
        throw RuntimeError("Could not set detector type");
    }

    // callbacks after (in setdetectortype, the object is reinitialized)
    if (startAcquisitionCallBack != nullptr)
        impl()->registerCallBackStartAcquisition(startAcquisitionCallBack,
                                                 pStartAcquisition);
    if (acquisitionFinishedCallBack != nullptr)
        impl()->registerCallBackAcquisitionFinished(acquisitionFinishedCallBack,
                                                    pAcquisitionFinished);
    if (rawDataReadyCallBack != nullptr)
        impl()->registerCallBackRawDataReady(rawDataReadyCallBack,
                                             pRawDataReady);
    if (rawDataModifyReadyCallBack != nullptr)
        impl()->registerCallBackRawDataModifyReady(rawDataModifyReadyCallBack,
                                                   pRawDataReady);

    impl()->setThreadIds(parentThreadId, tcpThreadId);
}

int ClientInterface::set_roi(Interface &socket) {
    auto arg = socket.Receive<ROI>();
    LOG(logDEBUG1) << "Set ROI: [" << arg.xmin << ", " << arg.xmax << "]";

    if (detType != GOTTHARD)
        functionNotImplemented();

    verifyIdle(socket);
    try {
        impl()->setROI(arg);
    } catch (const RuntimeError &e) {
        throw RuntimeError("Could not set ROI");
    }
    return socket.Send(OK);
}

int ClientInterface::set_num_frames(Interface &socket) {
    auto value = socket.Receive<int64_t>();
    if (value <= 0) {
        throw RuntimeError("Invalid number of frames " + std::to_string(value));
    }
    verifyIdle(socket);
    LOG(logDEBUG1) << "Setting num frames to " << value;
    impl()->setNumberOfFrames(value);
    return socket.Send(OK);
}

int ClientInterface::set_num_triggers(Interface &socket) {
    auto value = socket.Receive<int64_t>();
    if (value <= 0) {
        throw RuntimeError("Invalid number of triggers " +
                           std::to_string(value));
    }
    verifyIdle(socket);
    impl()->setNumberOfTriggers(value);
    return socket.Send(OK);
}

int ClientInterface::set_num_bursts(Interface &socket) {
    auto value = socket.Receive<int64_t>();
    if (value <= 0) {
        throw RuntimeError("Invalid number of bursts " + std::to_string(value));
    }
    verifyIdle(socket);
    impl()->setNumberOfBursts(value);
    return socket.Send(OK);
}

int ClientInterface::set_num_add_storage_cells(Interface &socket) {
    auto value = socket.Receive<int>();
    if (value < 0) {
        throw RuntimeError("Invalid number of additional storage cells " +
                           std::to_string(value));
    }
    verifyIdle(socket);
    LOG(logDEBUG1) << "Setting num additional storage cells to " << value;
    impl()->setNumberOfAdditionalStorageCells(value);
    return socket.Send(OK);
}

int ClientInterface::set_timing_mode(Interface &socket) {
    auto value = socket.Receive<int>();
    if (value < 0 || value >= NUM_TIMING_MODES) {
        throw RuntimeError("Invalid timing mode " + std::to_string(value));
    }
    verifyIdle(socket);
    LOG(logDEBUG1) << "Setting timing mode to " << value;
    impl()->setTimingMode(static_cast<timingMode>(value));
    return socket.Send(OK);
}

int ClientInterface::set_burst_mode(Interface &socket) {
    auto value = socket.Receive<int>();
    if (value < 0 || value >= NUM_BURST_MODES) {
        throw RuntimeError("Invalid burst mode " + std::to_string(value));
    }
    verifyIdle(socket);
    LOG(logDEBUG1) << "Setting burst mode to " << value;
    impl()->setBurstMode(static_cast<burstMode>(value));
    return socket.Send(OK);
}

int ClientInterface::set_num_analog_samples(Interface &socket) {
    auto value = socket.Receive<int>();
    LOG(logDEBUG1) << "Setting num analog samples to " << value;
    if (detType != CHIPTESTBOARD && detType != MOENCH) {
        functionNotImplemented();
    }
    try {
        impl()->setNumberofAnalogSamples(value);
    } catch (const RuntimeError &e) {
        throw RuntimeError("Could not set num analog samples to " +
                           std::to_string(value) +
                           " due to fifo structure memory allocation.");
    }
    return socket.Send(OK);
}

int ClientInterface::set_num_digital_samples(Interface &socket) {
    auto value = socket.Receive<int>();
    LOG(logDEBUG1) << "Setting num digital samples to " << value;
    if (detType != CHIPTESTBOARD) {
        functionNotImplemented();
    }
    try {
        impl()->setNumberofDigitalSamples(value);
    } catch (const RuntimeError &e) {
        throw RuntimeError("Could not set num digital samples to " +
                           std::to_string(value) +
                           " due to fifo structure memory allocation.");
    }
    return socket.Send(OK);
}

int ClientInterface::set_exptime(Interface &socket) {
    int64_t args[2]{-1, -1};
    socket.Receive(args);
    int gateIndex = static_cast<int>(args[0]);
    ns value = std::chrono::nanoseconds(args[1]);
    LOG(logDEBUG1) << "Setting exptime to " << sls::ToString(value)
                   << " (gateIndex: " << gateIndex << ")";
    switch (gateIndex) {
    case -1:
        if (detType == MYTHEN3) {
            impl()->setAcquisitionTime1(value);
            impl()->setAcquisitionTime2(value);
            impl()->setAcquisitionTime3(value);
        } else {
            impl()->setAcquisitionTime(value);
        }
        break;
    case 0:
        if (detType != MYTHEN3) {
            functionNotImplemented();
        }
        impl()->setAcquisitionTime1(value);
        break;
    case 1:
        if (detType != MYTHEN3) {
            functionNotImplemented();
        }
        impl()->setAcquisitionTime2(value);
        break;
    case 2:
        if (detType != MYTHEN3) {
            functionNotImplemented();
        }
        impl()->setAcquisitionTime3(value);
        break;
    default:
        throw RuntimeError("Unknown gate index for exptime " +
                           std::to_string(gateIndex));
    }
    return socket.Send(OK);
}

int ClientInterface::set_period(Interface &socket) {
    auto value = std::chrono::nanoseconds(socket.Receive<int64_t>());
    LOG(logDEBUG1) << "Setting period to " << sls::ToString(value);
    impl()->setAcquisitionPeriod(value);
    return socket.Send(OK);
}

int ClientInterface::set_subexptime(Interface &socket) {
    auto value = std::chrono::nanoseconds(socket.Receive<int64_t>());
    LOG(logDEBUG1) << "Setting period to " << sls::ToString(value);
    ns subdeadtime = impl()->getSubPeriod() - impl()->getSubExpTime();
    impl()->setSubExpTime(value);
    impl()->setSubPeriod(impl()->getSubExpTime() + subdeadtime);
    return socket.Send(OK);
}

int ClientInterface::set_subdeadtime(Interface &socket) {
    auto value = std::chrono::nanoseconds(socket.Receive<int64_t>());
    LOG(logDEBUG1) << "Setting sub deadtime to " << sls::ToString(value);
    impl()->setSubPeriod(value + impl()->getSubExpTime());
    LOG(logDEBUG1) << "Setting sub period to "
                   << sls::ToString(impl()->getSubPeriod());
    return socket.Send(OK);
}

int ClientInterface::set_dynamic_range(Interface &socket) {
    auto dr = socket.Receive<int>();
    if (dr >= 0) {
        verifyIdle(socket);
        LOG(logDEBUG1) << "Setting dynamic range: " << dr;
        bool exists = false;
        switch (dr) {
        case 16:
            exists = true;
            break;
        /*case 1: //TODO: Not yet implemented in firmware
            if (detType == MYTHEN3) {
                exists = true;
            }
            break;
        */
        case 4:
            if (detType == EIGER) {
                exists = true;
            }
            break;
        case 8:
        case 32:
            if (detType == EIGER || detType == MYTHEN3) {
                exists = true;
            }
            break;
        default:
            break;
        }
        if (!exists) {
            modeNotImplemented("Dynamic range", dr);
        } else {
            try {
                impl()->setDynamicRange(dr);
            } catch (const RuntimeError &e) {
                throw RuntimeError("Could not allocate memory for fifo or "
                                   "could not start listening/writing threads");
            }
        }
    }
    int retval = impl()->getDynamicRange();
    validate(dr, retval, "set dynamic range", DEC);
    LOG(logDEBUG1) << "dynamic range: " << retval;
    return socket.sendResult(retval);
}

int ClientInterface::set_streaming_frequency(Interface &socket) {
    auto index = socket.Receive<int>();
    if (index < 0) {
        throw RuntimeError("Invalid streaming frequency: " +
                           std::to_string(index));
    }
    verifyIdle(socket);
    impl()->setStreamingFrequency(index);
    return socket.Send(OK);
}

int ClientInterface::get_streaming_frequency(Interface &socket) {
    int retval = impl()->getStreamingFrequency();
    LOG(logDEBUG1) << "streaming freq:" << retval;
    return socket.sendResult(retval);
}

int ClientInterface::get_status(Interface &socket) {
    auto retval = impl()->getStatus();
    LOG(logDEBUG1) << "Status:" << sls::ToString(retval);
    return socket.sendResult(retval);
}

int ClientInterface::start_receiver(Interface &socket) {
    if (impl()->getStatus() == IDLE) {
        LOG(logDEBUG1) << "Starting Receiver";
        impl()->startReceiver();
    }
    return socket.Send(OK);
}

int ClientInterface::stop_receiver(Interface &socket) {
    auto arg = socket.Receive<int>();
    if (impl()->getStatus() == RUNNING) {
        LOG(logDEBUG1) << "Stopping Receiver";
        impl()->setStoppedFlag(static_cast<bool>(arg));
        impl()->stopReceiver();
    }
    auto s = impl()->getStatus();
    if (s != IDLE)
        throw RuntimeError("Could not stop receiver. It as it is: " +
                           sls::ToString(s));

    return socket.Send(OK);
}

int ClientInterface::set_file_dir(Interface &socket) {
    std::string fpath = socket.Receive(MAX_STR_LENGTH);

    if (fpath.empty()) {
        throw RuntimeError("Cannot set empty file path");
    }
    if (fpath[0] != '/')
        throw RuntimeError("Receiver path needs to be absolute path");

    LOG(logDEBUG1) << "Setting file path: " << fpath;
    impl()->setFilePath(fpath);
    return socket.Send(OK);
}

int ClientInterface::get_file_dir(Interface &socket) {
    auto fpath = impl()->getFilePath();
    LOG(logDEBUG1) << "file path:" << fpath;
    fpath.resize(MAX_STR_LENGTH);
    return socket.sendResult(fpath);
}

int ClientInterface::set_file_name(Interface &socket) {
    std::string fname = socket.Receive(MAX_STR_LENGTH);
    if (fname.empty()) {
        throw RuntimeError("Cannot set empty file name");
    }
    LOG(logDEBUG1) << "Setting file name: " << fname;
    impl()->setFileName(fname);
    return socket.Send(OK);
}

int ClientInterface::get_file_name(Interface &socket) {
    auto fname = impl()->getFileName();
    LOG(logDEBUG1) << "file name:" << fname;
    fname.resize(MAX_STR_LENGTH);
    return socket.sendResult(fname);
}

int ClientInterface::set_file_index(Interface &socket) {
    auto index = socket.Receive<int64_t>();
    if (index < 0) {
        throw RuntimeError("Invalid file index: " + std::to_string(index));
    }
    verifyIdle(socket);
    impl()->setFileIndex(index);
    return socket.Send(OK);
}

int ClientInterface::get_file_index(Interface &socket) {
    int64_t retval = impl()->getFileIndex();
    LOG(logDEBUG1) << "file index:" << retval;
    return socket.sendResult(retval);
}

int ClientInterface::get_frame_index(Interface &socket) {
    uint64_t retval = impl()->getCurrentFrameIndex();
    LOG(logDEBUG1) << "frame index:" << retval;
    return socket.sendResult(retval);
}

int ClientInterface::get_missing_packets(Interface &socket) {
    auto missing_packets = impl()->getNumMissingPackets();
    LOG(logDEBUG1) << "missing packets:" << sls::ToString(missing_packets);
    auto size = static_cast<int>(missing_packets.size());
    socket.Send(OK);
    socket.Send(size);
    socket.Send(missing_packets);
    return OK;
}

int ClientInterface::get_frames_caught(Interface &socket) {
    int64_t retval = impl()->getFramesCaught();
    LOG(logDEBUG1) << "frames caught:" << retval;
    return socket.sendResult(retval);
}

int ClientInterface::set_file_write(Interface &socket) {
    auto enable = socket.Receive<int>();
    if (enable < 0) {
        throw RuntimeError("Invalid file write enable: " +
                           std::to_string(enable));
    }
    verifyIdle(socket);
    LOG(logDEBUG1) << "Setting File write enable:" << enable;
    impl()->setFileWriteEnable(enable);
    return socket.Send(OK);
}

int ClientInterface::get_file_write(Interface &socket) {
    int retval = impl()->getFileWriteEnable();
    LOG(logDEBUG1) << "file write enable:" << retval;
    return socket.sendResult(retval);
}

int ClientInterface::set_master_file_write(Interface &socket) {
    auto enable = socket.Receive<int>();
    if (enable < 0) {
        throw RuntimeError("Invalid master file write enable: " +
                           std::to_string(enable));
    }
    verifyIdle(socket);
    impl()->setMasterFileWriteEnable(enable);
    return socket.Send(OK);
}

int ClientInterface::get_master_file_write(Interface &socket) {
    int retval = impl()->getMasterFileWriteEnable();
    LOG(logDEBUG1) << "master file write enable:" << retval;
    return socket.sendResult(retval);
}

int ClientInterface::set_overwrite(Interface &socket) {
    auto index = socket.Receive<int>();
    if (index < 0) {
        throw RuntimeError("Invalid over write enable: " +
                           std::to_string(index));
    }
    verifyIdle(socket);
    impl()->setOverwriteEnable(index);
    return socket.Send(OK);
}

int ClientInterface::get_overwrite(Interface &socket) {
    int retval = impl()->getOverwriteEnable();
    LOG(logDEBUG1) << "file overwrite enable:" << retval;
    return socket.sendResult(retval);
}

int ClientInterface::enable_tengiga(Interface &socket) {
    auto val = socket.Receive<int>();
    if (detType != EIGER && detType != CHIPTESTBOARD && detType != MOENCH &&
        detType != MYTHEN3)
        functionNotImplemented();

    if (val >= 0) {
        verifyIdle(socket);
        LOG(logDEBUG1) << "Setting 10GbE:" << val;
        try {
            impl()->setTenGigaEnable(val);
        } catch (const RuntimeError &e) {
            throw RuntimeError("Could not set 10GbE.");
        }
    }
    int retval = impl()->getTenGigaEnable();
    validate(val, retval, "set 10GbE", DEC);
    LOG(logDEBUG1) << "10Gbe:" << retval;
    return socket.sendResult(retval);
}

int ClientInterface::set_fifo_depth(Interface &socket) {
    auto value = socket.Receive<int>();
    if (value >= 0) {
        verifyIdle(socket);
        LOG(logDEBUG1) << "Setting fifo depth:" << value;
        try {
            impl()->setFifoDepth(value);
        } catch (const RuntimeError &e) {
            throw RuntimeError("Could not set fifo depth due to fifo structure "
                               "memory allocation.");
        }
    }
    int retval = impl()->getFifoDepth();
    validate(value, retval, std::string("set fifo depth"), DEC);
    LOG(logDEBUG1) << "fifo depth:" << retval;
    return socket.sendResult(retval);
}

int ClientInterface::set_activate(Interface &socket) {
    auto enable = socket.Receive<int>();
    if (detType != EIGER)
        functionNotImplemented();

    if (enable >= 0) {
        verifyIdle(socket);
        LOG(logDEBUG1) << "Setting activate:" << enable;
        impl()->setActivate(static_cast<bool>(enable));
    }
    auto retval = static_cast<int>(impl()->getActivate());
    validate(enable, retval, "set activate", DEC);
    LOG(logDEBUG1) << "Activate: " << retval;
    return socket.sendResult(retval);
}

int ClientInterface::set_streaming(Interface &socket) {
    auto index = socket.Receive<int>();
    if (index < 0) {
        throw RuntimeError("Invalid streaming enable: " +
                           std::to_string(index));
    }
    verifyIdle(socket);
    LOG(logDEBUG1) << "Setting data stream enable:" << index;
    try {
        impl()->setDataStreamEnable(index);
    } catch (const RuntimeError &e) {
        throw RuntimeError("Could not set data stream enable to " +
                           std::to_string(index));
    }
    return socket.Send(OK);
}

int ClientInterface::get_streaming(Interface &socket) {
    auto retval = static_cast<int>(impl()->getDataStreamEnable());
    LOG(logDEBUG1) << "data streaming enable:" << retval;
    return socket.sendResult(retval);
}

int ClientInterface::set_streaming_timer(Interface &socket) {
    auto index = socket.Receive<int>();
    if (index >= 0) {
        verifyIdle(socket);
        LOG(logDEBUG1) << "Setting streaming timer:" << index;
        impl()->setStreamingTimer(index);
    }
    int retval = impl()->getStreamingTimer();
    validate(index, retval, "set data stream timer", DEC);
    LOG(logDEBUG1) << "Streaming timer:" << retval;
    return socket.sendResult(retval);
}

int ClientInterface::get_flip_rows(Interface &socket) {
    if (detType != EIGER)
        functionNotImplemented();

    int retval = impl()->getFlipRows();
    LOG(logDEBUG1) << "Flip rows:" << retval;
    return socket.sendResult(retval);
}

int ClientInterface::set_flip_rows(Interface &socket) {
    auto arg = socket.Receive<int>();

    if (detType != EIGER)
        functionNotImplemented();

    if (arg != 0 && arg != 1) {
        throw RuntimeError("Could not set flip rows. Invalid argument: " +
                           std::to_string(arg));
    }
    verifyIdle(socket);
    LOG(logDEBUG1) << "Setting flip rows:" << arg;
    impl()->setFlipRows(static_cast<bool>(arg));

    int retval = impl()->getFlipRows();
    validate(arg, retval, std::string("set flip rows"), DEC);
    LOG(logDEBUG1) << "Flip rows:" << retval;
    return socket.sendResult(retval);
}

int ClientInterface::set_file_format(Interface &socket) {
    auto f = socket.Receive<fileFormat>();
    if (f < 0 || f > NUM_FILE_FORMATS) {
        throw RuntimeError("Invalid file format: " + std::to_string(f));
    }
    verifyIdle(socket);
    LOG(logDEBUG1) << "Setting file format:" << f;
    impl()->setFileFormat(f);

    auto retval = impl()->getFileFormat();
    validate(f, retval, "set file format", DEC);
    LOG(logDEBUG1) << "File Format: " << retval;
    return socket.Send(OK);
}

int ClientInterface::get_file_format(Interface &socket) {
    auto retval = impl()->getFileFormat();
    LOG(logDEBUG1) << "File Format: " << retval;
    return socket.sendResult(retval);
}

int ClientInterface::set_streaming_port(Interface &socket) {
    auto port = socket.Receive<int>();
    if (port < 0) {
        throw RuntimeError("Invalid zmq port " + std::to_string(port));
    }
    verifyIdle(socket);
    impl()->setStreamingPort(port);
    return socket.Send(OK);
}

int ClientInterface::get_streaming_port(Interface &socket) {
    int retval = impl()->getStreamingPort();
    LOG(logDEBUG1) << "streaming port:" << retval;
    return socket.sendResult(retval);
}

int ClientInterface::set_streaming_source_ip(Interface &socket) {
    auto ip = socket.Receive<sls::IpAddr>();
    if (ip == 0)
        throw RuntimeError("Invalid zmq ip " + ip.str());
    verifyIdle(socket);
    impl()->setStreamingSourceIP(ip);
    return socket.Send(OK);
}

int ClientInterface::get_streaming_source_ip(Interface &socket) {
    sls::IpAddr retval = impl()->getStreamingSourceIP();
    LOG(logDEBUG1) << "streaming IP:" << retval;
    return socket.sendResult(retval);
}

int ClientInterface::set_silent_mode(Interface &socket) {
    auto value = socket.Receive<int>();
    if (value < 0) {
        throw RuntimeError("Invalid silent mode: " + std::to_string(value));
    }
    verifyIdle(socket);
    LOG(logDEBUG1) << "Setting silent mode:" << value;
    impl()->setSilentMode(value);
    return socket.Send(OK);
}

int ClientInterface::get_silent_mode(Interface &socket) {
    auto retval = static_cast<int>(impl()->getSilentMode());
    LOG(logDEBUG1) << "silent mode:" << retval;
    return socket.sendResult(retval);
}

int ClientInterface::restream_stop(Interface &socket) {
    verifyIdle(socket);
    if (!impl()->getDataStreamEnable()) {
        throw RuntimeError(
            "Could not restream stop packet as data Streaming is disabled");
    } else {
        LOG(logDEBUG1) << "Restreaming stop";
        impl()->restreamStop();
    }
    return socket.Send(OK);
}

int ClientInterface::set_additional_json_header(Interface &socket) {
    std::map<std::string, std::string> json;
    auto size = socket.Receive<int>();
    if (size > 0) {
        std::string buff(size, '\0');
        socket.Receive(&buff[0], buff.size());
        std::istringstream iss(buff);
        std::string key, value;
        while (iss >> key) {
            iss >> value;
            json[key] = value;
        }
    }
    // verifyIdle(socket); allowing it to be set on the fly
    LOG(logDEBUG1) << "Setting additional json header: " << sls::ToString(json);
    impl()->setAdditionalJsonHeader(json);
    return socket.Send(OK);
}

int ClientInterface::get_additional_json_header(Interface &socket) {
    std::map<std::string, std::string> json = impl()->getAdditionalJsonHeader();
    LOG(logDEBUG1) << "additional json header:" << sls::ToString(json);
    std::ostringstream oss;
    for (auto &it : json) {
        oss << it.first << ' ' << it.second << ' ';
    }
    auto buff = oss.str();
    auto size = static_cast<int>(buff.size());
    socket.sendResult(size);
    if (size > 0)
        socket.Send(buff);
    return OK;
}

int ClientInterface::set_udp_socket_buffer_size(Interface &socket) {
    auto size = socket.Receive<int>();
    if (size == 0) {
        throw RuntimeError("Receiver socket buffer size must be > 0.");
    }
    if (size > 0) {
        verifyIdle(socket);
        if (size > INT_MAX / 2) {
            throw RuntimeError(
                "Receiver socket buffer size exceeded max (INT_MAX/2)");
        }
        LOG(logDEBUG1) << "Setting UDP Socket Buffer size: " << size;
        impl()->setUDPSocketBufferSize(size);
    }
    int retval = impl()->getUDPSocketBufferSize();
    if (size != 0)
        validate(size, retval,
                 "set udp socket buffer size (No CAP_NET_ADMIN privileges?)",
                 DEC);
    LOG(logDEBUG1) << "UDP Socket Buffer Size:" << retval;
    return socket.sendResult(retval);
}

int ClientInterface::get_real_udp_socket_buffer_size(Interface &socket) {
    auto size = impl()->getActualUDPSocketBufferSize();
    LOG(logDEBUG1) << "Actual UDP socket size :" << size;
    return socket.sendResult(size);
}

int ClientInterface::set_frames_per_file(Interface &socket) {
    auto index = socket.Receive<int>();
    if (index < 0) {
        throw RuntimeError("Invalid frames per file: " + std::to_string(index));
    }
    verifyIdle(socket);
    LOG(logDEBUG1) << "Setting frames per file: " << index;
    impl()->setFramesPerFile(index);
    return socket.Send(OK);
}

int ClientInterface::get_frames_per_file(Interface &socket) {
    auto retval = static_cast<int>(impl()->getFramesPerFile());
    LOG(logDEBUG1) << "frames per file:" << retval;
    return socket.sendResult(retval);
}

int ClientInterface::check_version_compatibility(Interface &socket) {
    auto arg = socket.Receive<int64_t>();
    LOG(logDEBUG1) << "Checking versioning compatibility with value " << arg;
    int64_t client_requiredVersion = arg;
    int64_t rx_apiVersion = APIRECEIVER;
    int64_t rx_version = getReceiverVersion();

    if (rx_apiVersion > client_requiredVersion) {
        std::ostringstream os;
        os << "Incompatible versions.\n Client's receiver API Version: (0x"
           << std::hex << client_requiredVersion
           << "). Receiver API Version: (0x" << std::hex
           << ").\n Please update the client!\n";
        throw RuntimeError(os.str());
    } else if (client_requiredVersion > rx_version) {
        std::ostringstream os;
        os << "This receiver is incompatible.\n Receiver Version: (0x"
           << std::hex << rx_version << "). Client's receiver API Version: (0x"
           << std::hex << client_requiredVersion
           << ").\n Please update the receiver";
        throw RuntimeError(os.str());
    } else {
        LOG(logINFO) << "Compatibility with Client: Successful";
    }
    return socket.Send(OK);
}

int ClientInterface::set_discard_policy(Interface &socket) {
    auto index = socket.Receive<int>();
    if (index < 0 || index > NUM_DISCARD_POLICIES) {
        throw RuntimeError("Invalid discard policy " + std::to_string(index));
    }
    verifyIdle(socket);
    LOG(logDEBUG1) << "Setting frames discard policy: " << index;
    impl()->setFrameDiscardPolicy(static_cast<frameDiscardPolicy>(index));
    return socket.Send(OK);
}

int ClientInterface::get_discard_policy(Interface &socket) {
    int retval = impl()->getFrameDiscardPolicy();
    LOG(logDEBUG1) << "frame discard policy:" << retval;
    return socket.sendResult(retval);
}

int ClientInterface::set_padding_enable(Interface &socket) {
    auto index = socket.Receive<int>();
    if (index < 0) {
        throw RuntimeError("Invalid padding enable: " + std::to_string(index));
    }
    verifyIdle(socket);
    LOG(logDEBUG1) << "Setting frames padding enable: " << index;
    impl()->setFramePaddingEnable(static_cast<bool>(index));
    return socket.Send(OK);
}

int ClientInterface::get_padding_enable(Interface &socket) {
    auto retval = static_cast<int>(impl()->getFramePaddingEnable());
    LOG(logDEBUG1) << "Frame Padding Enable:" << retval;
    return socket.sendResult(retval);
}

int ClientInterface::set_readout_mode(Interface &socket) {
    auto arg = socket.Receive<readoutMode>();

    if (detType != CHIPTESTBOARD)
        functionNotImplemented();

    if (arg >= 0) {
        verifyIdle(socket);
        LOG(logDEBUG1) << "Setting readout mode: " << arg;
        try {
            impl()->setReadoutMode(arg);
        } catch (const RuntimeError &e) {
            throw RuntimeError(
                "Could not set read out mode due to fifo memory allocation.");
        }
    }
    auto retval = impl()->getReadoutMode();
    validate(static_cast<int>(arg), static_cast<int>(retval),
             "set readout mode", DEC);
    LOG(logDEBUG1) << "Readout mode: " << retval;
    return socket.sendResult(retval);
}

int ClientInterface::set_adc_mask(Interface &socket) {
    auto arg = socket.Receive<uint32_t>();
    verifyIdle(socket);
    LOG(logDEBUG1) << "Setting 1Gb ADC enable mask: " << arg;
    try {
        impl()->setADCEnableMask(arg);
    } catch (const RuntimeError &e) {
        throw RuntimeError(
            "Could not set adc enable mask due to fifo memory allcoation");
    }
    auto retval = impl()->getADCEnableMask();
    if (retval != arg) {
        std::ostringstream os;
        os << "Could not set 1Gb ADC enable mask. Set 0x" << std::hex << arg
           << " but read 0x" << std::hex << retval;
        throw RuntimeError(os.str());
    }
    LOG(logDEBUG1) << "1Gb ADC enable mask retval: " << retval;
    return socket.sendResult(retval);
}

int ClientInterface::set_dbit_list(Interface &socket) {
    sls::StaticVector<int, MAX_RX_DBIT> args;
    socket.Receive(args);
    if (detType != CHIPTESTBOARD)
        functionNotImplemented();
    LOG(logDEBUG1) << "Setting DBIT list";
    for (auto &it : args) {
        LOG(logDEBUG1) << it << " ";
    }
    LOG(logDEBUG1) << '\n';
    verifyIdle(socket);
    impl()->setDbitList(args);
    return socket.Send(OK);
}

int ClientInterface::get_dbit_list(Interface &socket) {
    if (detType != CHIPTESTBOARD)
        functionNotImplemented();
    sls::StaticVector<int, MAX_RX_DBIT> retval;
    retval = impl()->getDbitList();
    LOG(logDEBUG1) << "Dbit list size retval:" << retval.size();
    return socket.sendResult(retval);
}

int ClientInterface::set_dbit_offset(Interface &socket) {
    auto arg = socket.Receive<int>();
    if (detType != CHIPTESTBOARD)
        functionNotImplemented();
    if (arg < 0) {
        throw RuntimeError("Invalid dbit offset: " + std::to_string(arg));
    }
    verifyIdle(socket);
    LOG(logDEBUG1) << "Setting Dbit offset: " << arg;
    impl()->setDbitOffset(arg);
    return socket.Send(OK);
}

int ClientInterface::get_dbit_offset(Interface &socket) {
    if (detType != CHIPTESTBOARD)
        functionNotImplemented();
    int retval = impl()->getDbitOffset();
    LOG(logDEBUG1) << "Dbit offset retval: " << retval;
    return socket.sendResult(retval);
}

int ClientInterface::set_quad_type(Interface &socket) {
    auto quadEnable = socket.Receive<int>();
    if (quadEnable >= 0) {
        verifyIdle(socket);
        LOG(logDEBUG1) << "Setting quad:" << quadEnable;
        try {
            impl()->setQuad(quadEnable == 0 ? false : true);
        } catch (const RuntimeError &e) {
            throw RuntimeError("Could not set quad to " +
                               std::to_string(quadEnable) +
                               " due to fifo strucutre memory allocation");
        }
    }
    int retval = impl()->getQuad() ? 1 : 0;
    validate(quadEnable, retval, "set quad", DEC);
    LOG(logDEBUG1) << "quad retval:" << retval;
    return socket.Send(OK);
}

int ClientInterface::set_read_n_rows(Interface &socket) {
    auto arg = socket.Receive<int>();
    if (arg >= 0) {
        verifyIdle(socket);
        if (detType != EIGER && detType != JUNGFRAU) {
            throw RuntimeError("Could not set number of rows. Not implemented "
                               "for this detector");
        }
        LOG(logDEBUG1) << "Setting number of rows:" << arg;
        impl()->setReadNRows(arg);
    }
    int retval = impl()->getReadNRows();
    validate(arg, retval, "set number of rows", DEC);
    LOG(logDEBUG1) << "read number of rows:" << retval;
    return socket.Send(OK);
}

sls::MacAddr ClientInterface::setUdpIp(sls::IpAddr arg) {
    LOG(logINFO) << "Received UDP IP: " << arg;
    // getting eth
    std::string eth = sls::IpToInterfaceName(arg.str());
    if (eth == "none") {
        throw RuntimeError("Failed to get udp ethernet interface from IP " +
                           arg.str());
    }
    if (eth.find('.') != std::string::npos) {
        eth = "";
        LOG(logERROR) << "Failed to get udp ethernet interface from IP " << arg
                      << ". Got " << eth;
    }
    impl()->setEthernetInterface(eth);
    if (detType == EIGER) {
        impl()->setEthernetInterface2(eth);
    }
    // get mac address
    auto retval = sls::InterfaceNameToMac(eth);
    if (retval == 0) {
        throw RuntimeError("Failed to get udp mac adddress to listen to (eth:" +
                           eth + ", ip:" + arg.str() + ")\n");
    }
    LOG(logINFO) << "Receiver MAC Address: " << retval;
    return retval;
}

int ClientInterface::set_udp_ip(Interface &socket) {
    auto arg = socket.Receive<sls::IpAddr>();
    verifyIdle(socket);
    auto retval = setUdpIp(arg);
    return socket.sendResult(retval);
}

sls::MacAddr ClientInterface::setUdpIp2(sls::IpAddr arg) {
    LOG(logINFO) << "Received UDP IP2: " << arg;
    // getting eth
    std::string eth = sls::IpToInterfaceName(arg.str());
    if (eth == "none") {
        throw RuntimeError("Failed to get udp ethernet interface2 from IP " +
                           arg.str());
    }
    if (eth.find('.') != std::string::npos) {
        eth = "";
        LOG(logERROR) << "Failed to get udp ethernet interface2 from IP " << arg
                      << ". Got " << eth;
    }
    impl()->setEthernetInterface2(eth);

    // get mac address
    auto retval = sls::InterfaceNameToMac(eth);
    if (retval == 0) {
        throw RuntimeError(
            "Failed to get udp mac adddress2 to listen to (eth:" + eth +
            ", ip:" + arg.str() + ")\n");
    }
    LOG(logINFO) << "Receiver MAC Address2: " << retval;
    return retval;
}

int ClientInterface::set_udp_ip2(Interface &socket) {
    auto arg = socket.Receive<sls::IpAddr>();
    verifyIdle(socket);
    if (detType != JUNGFRAU && detType != GOTTHARD2) {
        throw RuntimeError(
            "UDP Destination IP2 not implemented for this detector");
    }
    auto retval = setUdpIp2(arg);
    return socket.sendResult(retval);
}

int ClientInterface::set_udp_port(Interface &socket) {
    auto arg = socket.Receive<int>();
    verifyIdle(socket);
    LOG(logDEBUG1) << "Setting UDP Port:" << arg;
    impl()->setUDPPortNumber(arg);
    return socket.Send(OK);
}

int ClientInterface::set_udp_port2(Interface &socket) {
    auto arg = socket.Receive<int>();
    verifyIdle(socket);
    if (detType != JUNGFRAU && detType != EIGER && detType != GOTTHARD2) {
        throw RuntimeError(
            "UDP Destination Port2 not implemented for this detector");
    }
    LOG(logDEBUG1) << "Setting UDP Port:" << arg;
    impl()->setUDPPortNumber2(arg);
    return socket.Send(OK);
}

int ClientInterface::set_num_interfaces(Interface &socket) {
    auto arg = socket.Receive<int>();
    arg = (arg > 1 ? 2 : 1);
    verifyIdle(socket);
    if (detType != JUNGFRAU && detType != GOTTHARD2) {
        throw RuntimeError(
            "Number of interfaces not implemented for this detector");
    }
    LOG(logDEBUG1) << "Setting Number of UDP Interfaces:" << arg;
    try {
        impl()->setNumberofUDPInterfaces(arg);
    } catch (const RuntimeError &e) {
        throw RuntimeError("Failed to set number of interfaces to " +
                           std::to_string(arg));
    }
    return socket.Send(OK);
}

int ClientInterface::set_adc_mask_10g(Interface &socket) {
    auto arg = socket.Receive<uint32_t>();
    verifyIdle(socket);
    LOG(logDEBUG1) << "Setting 10Gb ADC enable mask: " << arg;
    try {
        impl()->setTenGigaADCEnableMask(arg);
    } catch (const RuntimeError &e) {
        throw RuntimeError(
            "Could not set 10Gb adc enable mask due to fifo memory allcoation");
    }
    auto retval = impl()->getTenGigaADCEnableMask();
    if (retval != arg) {
        std::ostringstream os;
        os << "Could not 10gb ADC enable mask. Set 0x" << std::hex << arg
           << " but read 0x" << std::hex << retval;
        throw RuntimeError(os.str());
    }
    LOG(logDEBUG1) << "10Gb ADC enable mask retval: " << retval;
    return socket.sendResult(retval);
}

int ClientInterface::set_counter_mask(Interface &socket) {
    auto arg = socket.Receive<uint32_t>();
    verifyIdle(socket);
    LOG(logDEBUG1) << "Setting counters: " << arg;
    impl()->setCounterMask(arg);
    return socket.Send(OK);
}

int ClientInterface::increment_file_index(Interface &socket) {
    verifyIdle(socket);
    if (impl()->getFileWriteEnable()) {
        LOG(logDEBUG1) << "Incrementing file index";
        impl()->setFileIndex(impl()->getFileIndex() + 1);
    }
    return socket.Send(OK);
}

int ClientInterface::set_additional_json_parameter(Interface &socket) {
    char args[2][SHORT_STR_LENGTH]{};
    socket.Receive(args);
    // verifyIdle(socket); allowing it to be set on the fly
    LOG(logDEBUG1) << "Setting additional json parameter (" << args[0]
                   << "): " << args[1];
    impl()->setAdditionalJsonParameter(args[0], args[1]);
    return socket.Send(OK);
}

int ClientInterface::get_additional_json_parameter(Interface &socket) {
    std::string key = socket.Receive(SHORT_STR_LENGTH);
    std::string value = impl()->getAdditionalJsonParameter(key);
    value.resize(SHORT_STR_LENGTH);
    return socket.sendResult(value);
}

int ClientInterface::get_progress(Interface &socket) {
    double retval = impl()->getProgress();
    LOG(logDEBUG1) << "progress retval: " << retval;
    return socket.sendResult(retval);
}

int ClientInterface::set_num_gates(Interface &socket) {
    auto value = socket.Receive<int>();
    LOG(logDEBUG1) << "Setting num gates to " << value;
    if (detType != MYTHEN3) {
        functionNotImplemented();
    }
    impl()->setNumberOfGates(value);
    return socket.Send(OK);
}

int ClientInterface::set_gate_delay(Interface &socket) {
    int64_t args[2]{-1, -1};
    socket.Receive(args);
    int gateIndex = static_cast<int>(args[0]);
    auto value = std::chrono::nanoseconds(args[1]);
    LOG(logDEBUG1) << "Setting gate delay to " << sls::ToString(value)
                   << " (gateIndex: " << gateIndex << ")";
    if (detType != MYTHEN3) {
        functionNotImplemented();
    }
    switch (gateIndex) {
    case -1:
        impl()->setGateDelay1(value);
        impl()->setGateDelay2(value);
        impl()->setGateDelay3(value);
        break;
    case 0:
        impl()->setGateDelay1(value);
        break;
    case 1:
        impl()->setGateDelay2(value);
        break;
    case 2:
        impl()->setGateDelay3(value);
        break;
    default:
        throw RuntimeError("Unknown gate index for gate delay " +
                           std::to_string(gateIndex));
    }
    return socket.Send(OK);
}

int ClientInterface::get_thread_ids(Interface &socket) {
    auto retval = impl()->getThreadIds();
    LOG(logDEBUG1) << "thread ids retval: " << sls::ToString(retval);
    return socket.sendResult(retval);
}

int ClientInterface::get_streaming_start_fnum(Interface &socket) {
    int retval = impl()->getStreamingStartingFrameNumber();
    LOG(logDEBUG1) << "streaming start fnum:" << retval;
    return socket.sendResult(retval);
}

int ClientInterface::set_streaming_start_fnum(Interface &socket) {
    auto index = socket.Receive<int>();
    if (index < 0) {
        throw RuntimeError("Invalid streaming start frame number: " +
                           std::to_string(index));
    }
    verifyIdle(socket);
    LOG(logDEBUG1) << "Setting streaming start fnum: " << index;
    impl()->setStreamingStartingFrameNumber(index);
    return socket.Send(OK);
}

int ClientInterface::set_rate_correct(Interface &socket) {
    auto index = socket.Receive<int>();
    if (index <= 0) {
        throw RuntimeError("Invalid number of rate correction values: " +
                           std::to_string(index));
    }
    LOG(logDEBUG) << "Number of detectors for rate correction: " << index;
    std::vector<int64_t> t(index);
    socket.Receive(t);
    verifyIdle(socket);
    LOG(logINFO) << "Setting rate corrections[" << index << ']';
    impl()->setRateCorrections(t);
    return socket.Send(OK);
}

int ClientInterface::set_scan(Interface &socket) {
    auto arg = socket.Receive<scanParameters>();
    LOG(logDEBUG) << "Scan Mode: " << sls::ToString(arg);
    verifyIdle(socket);
    impl()->setScan(arg);
    return socket.Send(OK);
}

int ClientInterface::set_threshold(Interface &socket) {
    auto arg = socket.Receive<int>();
    LOG(logDEBUG) << "Threshold: " << arg << " eV";
    if (detType != EIGER)
        functionNotImplemented();
    verifyIdle(socket);
    impl()->setThresholdEnergy(arg);
    return socket.Send(OK);
}

int ClientInterface::get_streaming_hwm(Interface &socket) {
    int retval = impl()->getStreamingHwm();
    LOG(logDEBUG1) << "zmq send hwm limit:" << retval;
    return socket.sendResult(retval);
}

int ClientInterface::set_streaming_hwm(Interface &socket) {
    auto limit = socket.Receive<int>();
    if (limit < -1) {
        throw RuntimeError("Invalid zmq send hwm limit " +
                           std::to_string(limit));
    }
    verifyIdle(socket);
    impl()->setStreamingHwm(limit);
    return socket.Send(OK);
}

int ClientInterface::set_all_threshold(Interface &socket) {
    auto eVs = socket.Receive<std::array<int, 3>>();
    LOG(logDEBUG) << "Threshold:" << sls::ToString(eVs);
    if (detType != MYTHEN3)
        functionNotImplemented();
    verifyIdle(socket);
    impl()->setThresholdEnergy(eVs);
    return socket.Send(OK);
}

int ClientInterface::set_detector_datastream(Interface &socket) {
    int args[2]{-1, -1};
    socket.Receive(args);
    portPosition port = static_cast<portPosition>(args[0]);
    switch (port) {
    case LEFT:
    case RIGHT:
        break;
    default:
        throw RuntimeError("Invalid port type");
    }
    bool enable = static_cast<int>(args[1]);
    LOG(logDEBUG1) << "Setting datastream (" << sls::ToString(port) << ") to "
                   << sls::ToString(enable);
    if (detType != EIGER)
        functionNotImplemented();
    verifyIdle(socket);
    impl()->setDetectorDataStream(port, enable);
    return socket.Send(OK);
}
