/* WARINING This file is auto generated any edits might be overwritten without
 * warning */

#include <pybind11/chrono.h>
#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "sls_detector_defs.h"
namespace py = pybind11;
void init_enums(py::module &m) {
    py::class_<slsDetectorDefs> Defs(m, "slsDetectorDefs");
    py::class_<slsDetectorDefs::xy> xy(m, "xy");
    xy.def_readwrite("x", &slsDetectorDefs::xy::x);
    xy.def_readwrite("y", &slsDetectorDefs::xy::y);

    py::enum_<slsDetectorDefs::detectorType>(Defs, "detectorType")
        .value("GET_DETECTOR_TYPE",
               slsDetectorDefs::detectorType::GET_DETECTOR_TYPE)
        .value("GENERIC", slsDetectorDefs::detectorType::GENERIC)
        .value("EIGER", slsDetectorDefs::detectorType::EIGER)
        .value("GOTTHARD", slsDetectorDefs::detectorType::GOTTHARD)
        .value("JUNGFRAU", slsDetectorDefs::detectorType::JUNGFRAU)
        .value("CHIPTESTBOARD", slsDetectorDefs::detectorType::CHIPTESTBOARD)
        .value("MOENCH", slsDetectorDefs::detectorType::MOENCH)
        .value("MYTHEN3", slsDetectorDefs::detectorType::MYTHEN3)
        .value("GOTTHARD2", slsDetectorDefs::detectorType::GOTTHARD2)
        .export_values();

    py::enum_<slsDetectorDefs::timerIndex>(Defs, "timerIndex")
        .value("FRAME_NUMBER", slsDetectorDefs::timerIndex::FRAME_NUMBER)
        .value("ACQUISITION_TIME",
               slsDetectorDefs::timerIndex::ACQUISITION_TIME)
        .value("FRAME_PERIOD", slsDetectorDefs::timerIndex::FRAME_PERIOD)
        .value("DELAY_AFTER_TRIGGER",
               slsDetectorDefs::timerIndex::DELAY_AFTER_TRIGGER)
        .value("CYCLES_NUMBER", slsDetectorDefs::timerIndex::CYCLES_NUMBER)
        .value("ACTUAL_TIME", slsDetectorDefs::timerIndex::ACTUAL_TIME)
        .value("MEASUREMENT_TIME",
               slsDetectorDefs::timerIndex::MEASUREMENT_TIME)
        .value("PROGRESS", slsDetectorDefs::timerIndex::PROGRESS)
        .value("FRAMES_FROM_START",
               slsDetectorDefs::timerIndex::FRAMES_FROM_START)
        .value("FRAMES_FROM_START_PG",
               slsDetectorDefs::timerIndex::FRAMES_FROM_START_PG)
        .value("ANALOG_SAMPLES", slsDetectorDefs::timerIndex::ANALOG_SAMPLES)
        .value("DIGITAL_SAMPLES", slsDetectorDefs::timerIndex::DIGITAL_SAMPLES)
        .value("SUBFRAME_ACQUISITION_TIME",
               slsDetectorDefs::timerIndex::SUBFRAME_ACQUISITION_TIME)
        .value("STORAGE_CELL_NUMBER",
               slsDetectorDefs::timerIndex::STORAGE_CELL_NUMBER)
        .value("SUBFRAME_DEADTIME",
               slsDetectorDefs::timerIndex::SUBFRAME_DEADTIME)
        .value("MEASURED_PERIOD", slsDetectorDefs::timerIndex::MEASURED_PERIOD)
        .value("MEASURED_SUBPERIOD",
               slsDetectorDefs::timerIndex::MEASURED_SUBPERIOD)
        .value("STORAGE_CELL_DELAY",
               slsDetectorDefs::timerIndex::STORAGE_CELL_DELAY)
        .value("MAX_TIMERS", slsDetectorDefs::timerIndex::MAX_TIMERS)
        .export_values();

    py::enum_<slsDetectorDefs::runStatus>(Defs, "runStatus")
        .value("IDLE", slsDetectorDefs::runStatus::IDLE)
        .value("ERROR", slsDetectorDefs::runStatus::ERROR)
        .value("WAITING", slsDetectorDefs::runStatus::WAITING)
        .value("RUN_FINISHED", slsDetectorDefs::runStatus::RUN_FINISHED)
        .value("TRANSMITTING", slsDetectorDefs::runStatus::TRANSMITTING)
        .value("RUNNING", slsDetectorDefs::runStatus::RUNNING)
        .value("STOPPED", slsDetectorDefs::runStatus::STOPPED)
        .export_values();

    py::enum_<slsDetectorDefs::frameDiscardPolicy>(Defs, "frameDiscardPolicy")
        .value("GET_FRAME_DISCARD_POLICY",
               slsDetectorDefs::frameDiscardPolicy::GET_FRAME_DISCARD_POLICY)
        .value("NO_DISCARD", slsDetectorDefs::frameDiscardPolicy::NO_DISCARD)
        .value("DISCARD_EMPTY_FRAMES",
               slsDetectorDefs::frameDiscardPolicy::DISCARD_EMPTY_FRAMES)
        .value("DISCARD_PARTIAL_FRAMES",
               slsDetectorDefs::frameDiscardPolicy::DISCARD_PARTIAL_FRAMES)
        .value("NUM_DISCARD_POLICIES",
               slsDetectorDefs::frameDiscardPolicy::NUM_DISCARD_POLICIES)
        .export_values();

    py::enum_<slsDetectorDefs::fileFormat>(Defs, "fileFormat")
        .value("GET_FILE_FORMAT", slsDetectorDefs::fileFormat::GET_FILE_FORMAT)
        .value("BINARY", slsDetectorDefs::fileFormat::BINARY)
        .value("HDF5", slsDetectorDefs::fileFormat::HDF5)
        .value("NUM_FILE_FORMATS",
               slsDetectorDefs::fileFormat::NUM_FILE_FORMATS)
        .export_values();

    py::enum_<slsDetectorDefs::networkParameter>(Defs, "networkParameter")
        .value("DETECTOR_MAC", slsDetectorDefs::networkParameter::DETECTOR_MAC)
        .value("DETECTOR_IP", slsDetectorDefs::networkParameter::DETECTOR_IP)
        .value("RECEIVER_HOSTNAME",
               slsDetectorDefs::networkParameter::RECEIVER_HOSTNAME)
        .value("RECEIVER_UDP_IP",
               slsDetectorDefs::networkParameter::RECEIVER_UDP_IP)
        .value("RECEIVER_UDP_PORT",
               slsDetectorDefs::networkParameter::RECEIVER_UDP_PORT)
        .value("RECEIVER_UDP_MAC",
               slsDetectorDefs::networkParameter::RECEIVER_UDP_MAC)
        .value("RECEIVER_UDP_PORT2",
               slsDetectorDefs::networkParameter::RECEIVER_UDP_PORT2)
        .value("DETECTOR_TXN_DELAY_LEFT",
               slsDetectorDefs::networkParameter::DETECTOR_TXN_DELAY_LEFT)
        .value("DETECTOR_TXN_DELAY_RIGHT",
               slsDetectorDefs::networkParameter::DETECTOR_TXN_DELAY_RIGHT)
        .value("DETECTOR_TXN_DELAY_FRAME",
               slsDetectorDefs::networkParameter::DETECTOR_TXN_DELAY_FRAME)
        .value("FLOW_CONTROL_10G",
               slsDetectorDefs::networkParameter::FLOW_CONTROL_10G)
        .value("FLOW_CONTROL_WR_PTR",
               slsDetectorDefs::networkParameter::FLOW_CONTROL_WR_PTR)
        .value("FLOW_CONTROL_RD_PTR",
               slsDetectorDefs::networkParameter::FLOW_CONTROL_RD_PTR)
        .value("RECEIVER_STREAMING_PORT",
               slsDetectorDefs::networkParameter::RECEIVER_STREAMING_PORT)
        .value("CLIENT_STREAMING_PORT",
               slsDetectorDefs::networkParameter::CLIENT_STREAMING_PORT)
        .value("RECEIVER_STREAMING_SRC_IP",
               slsDetectorDefs::networkParameter::RECEIVER_STREAMING_SRC_IP)
        .value("CLIENT_STREAMING_SRC_IP",
               slsDetectorDefs::networkParameter::CLIENT_STREAMING_SRC_IP)
        .value("ADDITIONAL_JSON_HEADER",
               slsDetectorDefs::networkParameter::ADDITIONAL_JSON_HEADER)
        .value("RECEIVER_UDP_SCKT_BUF_SIZE",
               slsDetectorDefs::networkParameter::RECEIVER_UDP_SCKT_BUF_SIZE)
        .value(
            "RECEIVER_REAL_UDP_SCKT_BUF_SIZE",
            slsDetectorDefs::networkParameter::RECEIVER_REAL_UDP_SCKT_BUF_SIZE)
        .export_values();

    py::enum_<slsDetectorDefs::numberOf>(Defs, "numberOf")
        .value("MAXMODX", slsDetectorDefs::numberOf::MAXMODX)
        .value("MAXMODY", slsDetectorDefs::numberOf::MAXMODY)
        .value("NMODX", slsDetectorDefs::numberOf::NMODX)
        .value("NMODY", slsDetectorDefs::numberOf::NMODY)
        .value("NCHANSX", slsDetectorDefs::numberOf::NCHANSX)
        .value("NCHANSY", slsDetectorDefs::numberOf::NCHANSY)
        .value("NCHIPSX", slsDetectorDefs::numberOf::NCHIPSX)
        .value("NCHIPSY", slsDetectorDefs::numberOf::NCHIPSY)
        .export_values();

    py::enum_<slsDetectorDefs::dimension>(Defs, "dimension")
        .value("X", slsDetectorDefs::dimension::X)
        .value("Y", slsDetectorDefs::dimension::Y)
        .export_values();

    py::enum_<slsDetectorDefs::externalSignalFlag>(Defs, "externalSignalFlag")
        .value("GET_EXTERNAL_SIGNAL_FLAG",
               slsDetectorDefs::externalSignalFlag::GET_EXTERNAL_SIGNAL_FLAG)
        .value("SIGNAL_OFF", slsDetectorDefs::externalSignalFlag::SIGNAL_OFF)
        .value("GATE_IN_ACTIVE_HIGH",
               slsDetectorDefs::externalSignalFlag::GATE_IN_ACTIVE_HIGH)
        .value("GATE_IN_ACTIVE_LOW",
               slsDetectorDefs::externalSignalFlag::GATE_IN_ACTIVE_LOW)
        .value("TRIGGER_IN_RISING_EDGE",
               slsDetectorDefs::externalSignalFlag::TRIGGER_IN_RISING_EDGE)
        .value("TRIGGER_IN_FALLING_EDGE",
               slsDetectorDefs::externalSignalFlag::TRIGGER_IN_FALLING_EDGE)
        .value("RO_TRIGGER_IN_RISING_EDGE",
               slsDetectorDefs::externalSignalFlag::RO_TRIGGER_IN_RISING_EDGE)
        .value("RO_TRIGGER_IN_FALLING_EDGE",
               slsDetectorDefs::externalSignalFlag::RO_TRIGGER_IN_FALLING_EDGE)
        .value("GATE_OUT_ACTIVE_HIGH",
               slsDetectorDefs::externalSignalFlag::GATE_OUT_ACTIVE_HIGH)
        .value("GATE_OUT_ACTIVE_LOW",
               slsDetectorDefs::externalSignalFlag::GATE_OUT_ACTIVE_LOW)
        .value("TRIGGER_OUT_RISING_EDGE",
               slsDetectorDefs::externalSignalFlag::TRIGGER_OUT_RISING_EDGE)
        .value("TRIGGER_OUT_FALLING_EDGE",
               slsDetectorDefs::externalSignalFlag::TRIGGER_OUT_FALLING_EDGE)
        .value("RO_TRIGGER_OUT_RISING_EDGE",
               slsDetectorDefs::externalSignalFlag::RO_TRIGGER_OUT_RISING_EDGE)
        .value("RO_TRIGGER_OUT_FALLING_EDGE",
               slsDetectorDefs::externalSignalFlag::RO_TRIGGER_OUT_FALLING_EDGE)
        .value("OUTPUT_LOW", slsDetectorDefs::externalSignalFlag::OUTPUT_LOW)
        .value("OUTPUT_HIGH", slsDetectorDefs::externalSignalFlag::OUTPUT_HIGH)
        .value(
            "MASTER_SLAVE_SYNCHRONIZATION",
            slsDetectorDefs::externalSignalFlag::MASTER_SLAVE_SYNCHRONIZATION)
        .export_values();

    py::enum_<slsDetectorDefs::timingMode>(Defs, "timingMode")
        .value("GET_TIMING_MODE", slsDetectorDefs::timingMode::GET_TIMING_MODE)
        .value("AUTO_TIMING", slsDetectorDefs::timingMode::AUTO_TIMING)
        .value("TRIGGER_EXPOSURE",
               slsDetectorDefs::timingMode::TRIGGER_EXPOSURE)
        .value("GATED", slsDetectorDefs::timingMode::GATED)
        .value("BURST_TRIGGER", slsDetectorDefs::timingMode::BURST_TRIGGER)
        .export_values();

    py::enum_<slsDetectorDefs::idMode>(Defs, "idMode")
        .value("DETECTOR_SERIAL_NUMBER",
               slsDetectorDefs::idMode::DETECTOR_SERIAL_NUMBER)
        .value("DETECTOR_FIRMWARE_VERSION",
               slsDetectorDefs::idMode::DETECTOR_FIRMWARE_VERSION)
        .value("DETECTOR_SOFTWARE_VERSION",
               slsDetectorDefs::idMode::DETECTOR_SOFTWARE_VERSION)
        .value("THIS_SOFTWARE_VERSION",
               slsDetectorDefs::idMode::THIS_SOFTWARE_VERSION)
        .value("RECEIVER_VERSION", slsDetectorDefs::idMode::RECEIVER_VERSION)
        .value("SOFTWARE_FIRMWARE_API_VERSION",
               slsDetectorDefs::idMode::SOFTWARE_FIRMWARE_API_VERSION)
        .value("CLIENT_SOFTWARE_API_VERSION",
               slsDetectorDefs::idMode::CLIENT_SOFTWARE_API_VERSION)
        .value("CLIENT_RECEIVER_API_VERSION",
               slsDetectorDefs::idMode::CLIENT_RECEIVER_API_VERSION)
        .export_values();

    py::enum_<slsDetectorDefs::digitalTestMode>(Defs, "digitalTestMode")
        .value("DETECTOR_FIRMWARE_TEST",
               slsDetectorDefs::digitalTestMode::DETECTOR_FIRMWARE_TEST)
        .value("DETECTOR_BUS_TEST",
               slsDetectorDefs::digitalTestMode::DETECTOR_BUS_TEST)
        .value("IMAGE_TEST", slsDetectorDefs::digitalTestMode::IMAGE_TEST)
        .export_values();

    py::enum_<slsDetectorDefs::dacIndex>(Defs, "dacIndex")
        .value("THRESHOLD", slsDetectorDefs::dacIndex::THRESHOLD)
        .value("CALIBRATION_PULSE",
               slsDetectorDefs::dacIndex::CALIBRATION_PULSE)
        .value("TRIMBIT_SIZE", slsDetectorDefs::dacIndex::TRIMBIT_SIZE)
        .value("PREAMP", slsDetectorDefs::dacIndex::PREAMP)
        .value("SHAPER1", slsDetectorDefs::dacIndex::SHAPER1)
        .value("SHAPER2", slsDetectorDefs::dacIndex::SHAPER2)
        .value("TEMPERATURE_ADC", slsDetectorDefs::dacIndex::TEMPERATURE_ADC)
        .value("TEMPERATURE_FPGA", slsDetectorDefs::dacIndex::TEMPERATURE_FPGA)
        .value("HUMIDITY", slsDetectorDefs::dacIndex::HUMIDITY)
        .value("DETECTOR_BIAS", slsDetectorDefs::dacIndex::DETECTOR_BIAS)
        .value("VA_POT", slsDetectorDefs::dacIndex::VA_POT)
        .value("VDD_POT", slsDetectorDefs::dacIndex::VDD_POT)
        .value("VSH_POT", slsDetectorDefs::dacIndex::VSH_POT)
        .value("VIO_POT", slsDetectorDefs::dacIndex::VIO_POT)
        .value("G_VREF_DS", slsDetectorDefs::dacIndex::G_VREF_DS)
        .value("G_VCASCN_PB", slsDetectorDefs::dacIndex::G_VCASCN_PB)
        .value("G_VCASCP_PB", slsDetectorDefs::dacIndex::G_VCASCP_PB)
        .value("G_VOUT_CM", slsDetectorDefs::dacIndex::G_VOUT_CM)
        .value("G_VCASC_OUT", slsDetectorDefs::dacIndex::G_VCASC_OUT)
        .value("G_VIN_CM", slsDetectorDefs::dacIndex::G_VIN_CM)
        .value("G_VREF_COMP", slsDetectorDefs::dacIndex::G_VREF_COMP)
        .value("G_IB_TESTC", slsDetectorDefs::dacIndex::G_IB_TESTC)
        .value("E_SvP", slsDetectorDefs::dacIndex::E_SvP)
        .value("E_SvN", slsDetectorDefs::dacIndex::E_SvN)
        .value("E_Vtr", slsDetectorDefs::dacIndex::E_Vtr)
        .value("E_Vrf", slsDetectorDefs::dacIndex::E_Vrf)
        .value("E_Vrs", slsDetectorDefs::dacIndex::E_Vrs)
        .value("E_Vtgstv", slsDetectorDefs::dacIndex::E_Vtgstv)
        .value("E_Vcmp_ll", slsDetectorDefs::dacIndex::E_Vcmp_ll)
        .value("E_Vcmp_lr", slsDetectorDefs::dacIndex::E_Vcmp_lr)
        .value("E_cal", slsDetectorDefs::dacIndex::E_cal)
        .value("E_Vcmp_rl", slsDetectorDefs::dacIndex::E_Vcmp_rl)
        .value("E_Vcmp_rr", slsDetectorDefs::dacIndex::E_Vcmp_rr)
        .value("E_rxb_rb", slsDetectorDefs::dacIndex::E_rxb_rb)
        .value("E_rxb_lb", slsDetectorDefs::dacIndex::E_rxb_lb)
        .value("E_Vcp", slsDetectorDefs::dacIndex::E_Vcp)
        .value("E_Vcn", slsDetectorDefs::dacIndex::E_Vcn)
        .value("E_Vis", slsDetectorDefs::dacIndex::E_Vis)
        .value("IO_DELAY", slsDetectorDefs::dacIndex::IO_DELAY)
        .value("ADC_VPP", slsDetectorDefs::dacIndex::ADC_VPP)
        .value("HIGH_VOLTAGE", slsDetectorDefs::dacIndex::HIGH_VOLTAGE)
        .value("TEMPERATURE_FPGAEXT",
               slsDetectorDefs::dacIndex::TEMPERATURE_FPGAEXT)
        .value("TEMPERATURE_10GE", slsDetectorDefs::dacIndex::TEMPERATURE_10GE)
        .value("TEMPERATURE_DCDC", slsDetectorDefs::dacIndex::TEMPERATURE_DCDC)
        .value("TEMPERATURE_SODL", slsDetectorDefs::dacIndex::TEMPERATURE_SODL)
        .value("TEMPERATURE_SODR", slsDetectorDefs::dacIndex::TEMPERATURE_SODR)
        .value("TEMPERATURE_FPGA2",
               slsDetectorDefs::dacIndex::TEMPERATURE_FPGA2)
        .value("TEMPERATURE_FPGA3",
               slsDetectorDefs::dacIndex::TEMPERATURE_FPGA3)
        .value("M_vIpre", slsDetectorDefs::dacIndex::M_vIpre)
        .value("M_vIbias", slsDetectorDefs::dacIndex::M_vIbias)
        .value("M_vIinSh", slsDetectorDefs::dacIndex::M_vIinSh)
        .value("M_VdcSh", slsDetectorDefs::dacIndex::M_VdcSh)
        .value("M_Vth2", slsDetectorDefs::dacIndex::M_Vth2)
        .value("M_VPL", slsDetectorDefs::dacIndex::M_VPL)
        .value("M_Vth3", slsDetectorDefs::dacIndex::M_Vth3)
        .value("M_casSh", slsDetectorDefs::dacIndex::M_casSh)
        .value("M_cas", slsDetectorDefs::dacIndex::M_cas)
        .value("M_vIbiasSh", slsDetectorDefs::dacIndex::M_vIbiasSh)
        .value("M_vIcin", slsDetectorDefs::dacIndex::M_vIcin)
        .value("M_vIpreOut", slsDetectorDefs::dacIndex::M_vIpreOut)
        .value("V_POWER_A", slsDetectorDefs::dacIndex::V_POWER_A)
        .value("V_POWER_B", slsDetectorDefs::dacIndex::V_POWER_B)
        .value("V_POWER_C", slsDetectorDefs::dacIndex::V_POWER_C)
        .value("V_POWER_D", slsDetectorDefs::dacIndex::V_POWER_D)
        .value("V_POWER_IO", slsDetectorDefs::dacIndex::V_POWER_IO)
        .value("V_POWER_CHIP", slsDetectorDefs::dacIndex::V_POWER_CHIP)
        .value("I_POWER_A", slsDetectorDefs::dacIndex::I_POWER_A)
        .value("I_POWER_B", slsDetectorDefs::dacIndex::I_POWER_B)
        .value("I_POWER_C", slsDetectorDefs::dacIndex::I_POWER_C)
        .value("I_POWER_D", slsDetectorDefs::dacIndex::I_POWER_D)
        .value("I_POWER_IO", slsDetectorDefs::dacIndex::I_POWER_IO)
        .value("V_LIMIT", slsDetectorDefs::dacIndex::V_LIMIT)
        .value("SLOW_ADC0", slsDetectorDefs::dacIndex::SLOW_ADC0)
        .value("SLOW_ADC1", slsDetectorDefs::dacIndex::SLOW_ADC1)
        .value("SLOW_ADC2", slsDetectorDefs::dacIndex::SLOW_ADC2)
        .value("SLOW_ADC3", slsDetectorDefs::dacIndex::SLOW_ADC3)
        .value("SLOW_ADC4", slsDetectorDefs::dacIndex::SLOW_ADC4)
        .value("SLOW_ADC5", slsDetectorDefs::dacIndex::SLOW_ADC5)
        .value("SLOW_ADC6", slsDetectorDefs::dacIndex::SLOW_ADC6)
        .value("SLOW_ADC7", slsDetectorDefs::dacIndex::SLOW_ADC7)
        .value("SLOW_ADC_TEMP", slsDetectorDefs::dacIndex::SLOW_ADC_TEMP)
        .export_values();

    py::enum_<slsDetectorDefs::detectorSettings>(Defs, "detectorSettings")
        .value("GET_SETTINGS", slsDetectorDefs::detectorSettings::GET_SETTINGS)
        .value("STANDARD", slsDetectorDefs::detectorSettings::STANDARD)
        .value("FAST", slsDetectorDefs::detectorSettings::FAST)
        .value("HIGHGAIN", slsDetectorDefs::detectorSettings::HIGHGAIN)
        .value("DYNAMICGAIN", slsDetectorDefs::detectorSettings::DYNAMICGAIN)
        .value("LOWGAIN", slsDetectorDefs::detectorSettings::LOWGAIN)
        .value("MEDIUMGAIN", slsDetectorDefs::detectorSettings::MEDIUMGAIN)
        .value("VERYHIGHGAIN", slsDetectorDefs::detectorSettings::VERYHIGHGAIN)
        .value("DYNAMICHG0", slsDetectorDefs::detectorSettings::DYNAMICHG0)
        .value("FIXGAIN1", slsDetectorDefs::detectorSettings::FIXGAIN1)
        .value("FIXGAIN2", slsDetectorDefs::detectorSettings::FIXGAIN2)
        .value("FORCESWITCHG1",
               slsDetectorDefs::detectorSettings::FORCESWITCHG1)
        .value("FORCESWITCHG2",
               slsDetectorDefs::detectorSettings::FORCESWITCHG2)
        .value("VERYLOWGAIN", slsDetectorDefs::detectorSettings::VERYLOWGAIN)
        .value("UNDEFINED", slsDetectorDefs::detectorSettings::UNDEFINED)
        .value("UNINITIALIZED",
               slsDetectorDefs::detectorSettings::UNINITIALIZED)
        .export_values();

    py::enum_<slsDetectorDefs::speedVariable>(Defs, "speedVariable")
        .value("CLOCK_DIVIDER", slsDetectorDefs::speedVariable::CLOCK_DIVIDER)
        .value("ADC_CLOCK", slsDetectorDefs::speedVariable::ADC_CLOCK)
        .value("ADC_PHASE", slsDetectorDefs::speedVariable::ADC_PHASE)
        .value("ADC_PIPELINE", slsDetectorDefs::speedVariable::ADC_PIPELINE)
        .value("DBIT_CLOCK", slsDetectorDefs::speedVariable::DBIT_CLOCK)
        .value("DBIT_PHASE", slsDetectorDefs::speedVariable::DBIT_PHASE)
        .value("DBIT_PIPELINE", slsDetectorDefs::speedVariable::DBIT_PIPELINE)
        .value("MAX_ADC_PHASE_SHIFT",
               slsDetectorDefs::speedVariable::MAX_ADC_PHASE_SHIFT)
        .value("MAX_DBIT_PHASE_SHIFT",
               slsDetectorDefs::speedVariable::MAX_DBIT_PHASE_SHIFT)
        .value("SYNC_CLOCK", slsDetectorDefs::speedVariable::SYNC_CLOCK)
        .export_values();

    py::enum_<slsDetectorDefs::readoutMode>(Defs, "readoutMode")
        .value("ANALOG_ONLY", slsDetectorDefs::readoutMode::ANALOG_ONLY)
        .value("DIGITAL_ONLY", slsDetectorDefs::readoutMode::DIGITAL_ONLY)
        .value("ANALOG_AND_DIGITAL",
               slsDetectorDefs::readoutMode::ANALOG_AND_DIGITAL)
        .export_values();

    py::enum_<slsDetectorDefs::speedLevel>(Defs, "speedLevel")
        .value("FULL_SPEED", slsDetectorDefs::speedLevel::FULL_SPEED)
        .value("HALF_SPEED", slsDetectorDefs::speedLevel::HALF_SPEED)
        .value("QUARTER_SPEED", slsDetectorDefs::speedLevel::QUARTER_SPEED)
        .export_values();

    py::enum_<slsDetectorDefs::portType>(Defs, "portType")
        .value("CONTROL_PORT", slsDetectorDefs::portType::CONTROL_PORT)
        .value("STOP_PORT", slsDetectorDefs::portType::STOP_PORT)
        .value("DATA_PORT", slsDetectorDefs::portType::DATA_PORT)
        .export_values();

    py::enum_<slsDetectorDefs::masterFlags>(Defs, "masterFlags")
        .value("GET_MASTER", slsDetectorDefs::masterFlags::GET_MASTER)
        .value("NO_MASTER", slsDetectorDefs::masterFlags::NO_MASTER)
        .value("IS_MASTER", slsDetectorDefs::masterFlags::IS_MASTER)
        .value("IS_SLAVE", slsDetectorDefs::masterFlags::IS_SLAVE)
        .export_values();

    py::enum_<slsDetectorDefs::frameModeType>(Defs, "frameModeType")
        .value("GET_FRAME_MODE", slsDetectorDefs::frameModeType::GET_FRAME_MODE)
        .value("PEDESTAL", slsDetectorDefs::frameModeType::PEDESTAL)
        .value("NEW_PEDESTAL", slsDetectorDefs::frameModeType::NEW_PEDESTAL)
        .value("FLATFIELD", slsDetectorDefs::frameModeType::FLATFIELD)
        .value("NEW_FLATFIELD", slsDetectorDefs::frameModeType::NEW_FLATFIELD)
        .export_values();

    py::enum_<slsDetectorDefs::detectorModeType>(Defs, "detectorModeType")
        .value("GET_DETECTOR_MODE",
               slsDetectorDefs::detectorModeType::GET_DETECTOR_MODE)
        .value("COUNTING", slsDetectorDefs::detectorModeType::COUNTING)
        .value("INTERPOLATING",
               slsDetectorDefs::detectorModeType::INTERPOLATING)
        .value("ANALOG", slsDetectorDefs::detectorModeType::ANALOG)
        .export_values();
}
