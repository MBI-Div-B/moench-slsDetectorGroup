// SPDX-License-Identifier: LGPL-3.0-or-other
// Copyright (C) 2021 Contributors to the SLS Detector Package
#pragma once

#include "File.h"

#include <mutex>

class HDF5DataFile : private virtual slsDetectorDefs, public File {

  public:
    HDF5DataFile(const int index, std::mutex *hdf5Lib);
    ~HDF5DataFile();

    std::array<std::string, 2> GetFileAndDatasetName() const override;
    uint32_t GetFilesInAcquisition() const override;
    DataType GetPDataType() const override;
    std::vector<std::string> GetParameterNames() const override;
    std::vector<DataType> GetParameterDataTypes() const override;

    void CloseFile() override;

    void CreateFirstHDF5DataFile(
        const std::string filePath, const std::string fileNamePrefix,
        const uint64_t fileIndex, const bool overWriteEnable,
        const bool silentMode, const int modulePos,
        const int numUnitsPerReadout, const uint32_t udpPortNumber,
        const uint32_t maxFramesPerFile, const uint64_t numImages,
        const uint32_t nPixelsX, const uint32_t nPixelsY,
        const uint32_t dynamicRange) override;

    void WriteToFile(char *buffer, const int buffersize,
                     const uint64_t currentFrameNumber,
                     const uint32_t numPacketsCaught) override;

  private:
    void CreateFile();
    void WriteDataFile(const uint64_t currentFrameNumber, char *buffer);
    void WriteParameterDatasets(const uint64_t currentFrameNumber,
                                sls_receiver_header *rheader);
    void ExtendDataset();

    int index_;
    std::mutex *hdf5Lib_;
    H5File *fd_{nullptr};
    std::string fileName_;
    std::string dataSetName_;
    DataSpace *dataSpace_{nullptr};
    DataSet *dataSet_{nullptr};
    DataType dataType_{PredType::STD_U16LE};

    DataSpace *dataSpacePara_{nullptr};
    std::vector<DataSet *> dataSetPara_{nullptr};
    std::vector<std::string> parameterNames_;
    std::vector<DataType> parameterDataTypes_;

    uint32_t subFileIndex_{0};
    uint32_t numFramesInFile_{0};
    uint32_t numFilesInAcquisition_{0};
    uint32_t maxFramesPerFile_{0};
    uint64_t numImages_{0};
    uint64_t extNumImages_{0};
    uint32_t nPixelsX_{0};
    uint32_t nPixelsY_{0};
    uint32_t dynamicRange_{0};

    std::string filePath_;
    std::string fileNamePrefix_;
    uint64_t fileIndex_{0};
    bool overWriteEnable_{false};
    bool silentMode_{false};
    int detIndex_{0};
    int numUnitsPerReadout_{0};
    uint32_t udpPortNumber_{0};
};