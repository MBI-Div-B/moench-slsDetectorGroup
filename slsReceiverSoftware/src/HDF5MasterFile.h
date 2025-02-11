// SPDX-License-Identifier: LGPL-3.0-or-other
// Copyright (C) 2021 Contributors to the SLS Detector Package
#pragma once

#include "File.h"

#include <mutex>

class HDF5MasterFile : private virtual slsDetectorDefs, public File {

  public:
    HDF5MasterFile(std::mutex *hdf5Lib);
    ~HDF5MasterFile();

    void CloseFile() override;
    void LinkDataFile(std::string dataFilename, std::string dataSetname,
                      const std::vector<std::string> parameterNames,
                      const bool silentMode) override;
    void CreateMasterFile(const std::string filePath,
                          const std::string fileNamePrefix,
                          const uint64_t fileIndex, const bool overWriteEnable,
                          const bool silentMode,
                          MasterAttributes *attr) override;
    void UpdateMasterFile(MasterAttributes *attr, bool silentMode) override;

  private:
    std::mutex *hdf5Lib_;
    H5File *fd_{nullptr};
    std::string fileName_;
    std::string attrGroupName_;
};