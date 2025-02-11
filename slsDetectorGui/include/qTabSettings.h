// SPDX-License-Identifier: LGPL-3.0-or-other
// Copyright (C) 2021 Contributors to the SLS Detector Package
#pragma once
#include "sls/Detector.h"
#include "ui_form_tab_settings.h"
#include <QCheckBox>

class qTabSettings : public QWidget, private Ui::TabSettingsObject {
    Q_OBJECT

  public:
    qTabSettings(QWidget *parent, sls::Detector *detector);
    ~qTabSettings();
    void Refresh();
  public slots:
    void SetExportMode(bool exportMode);

  private slots:
    void SetSettings(int index);
    void SetGainMode(int index);
    void SetDynamicRange(int index);
    void SetThresholdEnergy(int index);
    void SetThresholdEnergies();
    void SetCounterMask();

  private:
    void SetupWidgetWindow();
    void SetupDetectorSettings();
    void SetupGainMode();
    void ShowFixG0(bool expertMode);
    void Initialization();

    void GetSettings();
    void GetGainMode();
    void GetDynamicRange();
    void GetThresholdEnergy();
    void GetThresholdEnergies();
    void GetCounterMask();

    sls::Detector *det;
    std::vector<QCheckBox *> counters;

    enum {
        STANDARD,
        FAST,
        HIGHGAIN,
        DYNAMICGAIN,
        LOWGAIN,
        MEDIUMGAIN,
        VERYHIGHGAIN,
        HIGHGAIN0,
        FIXGAIN1,
        FIXGAIN2,
        VERLOWGAIN,
        G1_HIGHGAIN,
        G1_LOWGAIN,
        G2_HIGHCAP_HIGHGAIN,
        G2_HIGHCAP_LOWGAIN,
        G2_LOWCAP_HIGHGAIN,
        G2_LOWCAP_LOWGAIN,
        G4_HIGHGAIN,
        G4_LOWGAIN,
        GAIN0,
        UNDEFINED,
        UNINITIALIZED,
        NUMSETTINGS
    };

    enum { DYNAMIC, FORCE_SWITCH_G1, FORCE_SWITCH_G2, FIX_G1, FIX_G2, FIX_G0 };
    bool isVisibleFixG0{false};

    enum { DYNAMICRANGE_32, DYNAMICRANGE_16, DYNAMICRANGE_8, DYNAMICRANGE_4 };
};
