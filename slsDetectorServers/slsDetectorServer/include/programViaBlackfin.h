// SPDX-License-Identifier: LGPL-3.0-or-other
// Copyright (C) 2021 Contributors to the SLS Detector Package
#pragma once

#include "common.h"

#include <stdint.h>
#include <stdio.h>
#include <sys/types.h>

#define TEMP_PROG_FOLDER_NAME           "/var/tmp/"
#define TEMP_PROG_FOLDER_NAME_ALL_FILES "/var/tmp/*"
#define TEMP_PROG_FILE_NAME             TEMP_PROG_FOLDER_NAME "tmp.rawbin"

int defineGPIOpins(char *mess);
int FPGAdontTouchFlash(char *mess);
int FPGATouchFlash(char *mess);
int resetFPGA(char *mess);

int emptyTempFolder(char *mess);
/**
 * deletes old file
 * verify memory available to copy
 * open file to copy
 */
int preparetoCopyProgram(char *mess, char *functionType, FILE **fd,
                         uint64_t fsize);
int eraseAndWriteToFlash(char *mess, enum PROGRAM_INDEX index,
                         char *functionType, char *clientChecksum,
                         ssize_t fsize);
int getDrive(char *mess, enum PROGRAM_INDEX index);
/** Notify fpga not to touch flash, open src and flash drive to write */
int openFileForFlash(char *mess, FILE **flashfd, FILE **srcfd);
int eraseFlash(char *mess);
/* write from tmp file to flash */
int writeToFlash(char *mess, ssize_t fsize, FILE *flashfd, FILE *srcfd);
/** Notify fpga to pick up firmware from flash and wait for status confirmation
 */
int waitForFPGAtoTouchFlash(char *mess);
int moveBinaryFile(char *mess, char *serverName);
