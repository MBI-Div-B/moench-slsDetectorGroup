#pragma once
#include "sls/Detector.h"
#include "ui_form_tab_dataoutput.h"

class qTabDataOutput : public QWidget, private Ui::TabDataOutputObject {
    Q_OBJECT

  public:
    qTabDataOutput(QWidget *parent, sls::Detector *detector);
    ~qTabDataOutput();
    void Refresh();

  private slots:
    void GetOutputDir();
    void BrowseOutputDir();
    void SetOutputDir(bool force = false);
    void ForceSetOutputDir();
    void SetFileFormat(int format);
    void SetOverwriteEnable(bool enable);
    void SetTenGigaEnable(bool enable);
    void EnableRateCorrection();
    void SetRateCorrection();
    void SetSpeed(int speed);
    void SetParallel(bool enable);
    void SetCounterMask();

  private:
    void SetupWidgetWindow();
    void Initialization();
    void PopulateDetectors();
    void EnableBrowse();
    void GetFileWrite();
    void GetFileName();
    void GetFileFormat();
    void GetFileOverwrite();
    void GetTenGigaEnable();
    void GetRateCorrection();
    void GetSpeed();
    void GetParallel();
    void GetCounterMask();

    sls::Detector *det;
    // Button group for radiobuttons for rate
    QButtonGroup *btnGroupRate;
};
