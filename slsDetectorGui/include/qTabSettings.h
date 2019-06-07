#pragma once

#include "qDefs.h"

#include "ui_form_tab_settings.h"

class multiSlsDetector;

/**
 *@short sets up the Settings parameters
 */
class qTabSettings: public QWidget, private Ui::TabSettingsObject{
	Q_OBJECT

public:
	/**
	 * The constructor
	 * @param parent is the parent tab widget
	 * @param detector is the detector returned from the detector tab
	 */
	qTabSettings(QWidget *parent, multiSlsDetector* detector);

	/**
	 * Destructor
	 */
	~qTabSettings();

	/**
	 * Refresh and update widgets
	 */
	void Refresh();


private slots:
	/**
	 * Set settings according to selection
	 * @param index index of selection
	 */
	void SetSettings(int index);

	/**
	 * Set dynamic range if possible
	 * @param index selection
	 */
	void SetDynamicRange(int index);

	/**
	 * Set threshold energy
	 * @param index selection
	 */
	void SetThresholdEnergy(int index);


private:

	/**
	 * Sets up the widget
	 */
	void SetupWidgetWindow();

	/**
	 * Sets up the detector settings
	 */
	void SetupDetectorSettings();

	/**
	 * Sets up all the slots and signals
	 */
	void Initialization();

	/**
	 * Get Settings
	 */
	void GetSettings();

	/**
	 * Gets the dynamic range and sets it on the gui
	 */
	void GetDynamicRange();

	/**
	 * Gets the threshold energy and update widget
	 */
	void GetThresholdEnergy();

	/** The sls detector object */
	multiSlsDetector *myDet;

	enum {
		STANDARD, 
		FAST, 
		HIGHGAIN, 
		DYNAMICGAIN, 
		LOWGAIN, 
		MEDIUMGAIN, 
		VERYHIGHGAIN, 
		DYNAMICHG0, 
		FIXGAIN1, 
		FIXGAIN2, 
		FORCESWITCHG1, 
		FORCESWITCHG2, 
		VERLOWGAIN,
		UNDEFINED, 
		UNINITIALIZED, 
		NUMSETTINGS
	};
	
	enum {
		DYNAMICRANGE_32,
		DYNAMICRANGE_16,
		DYNAMICRANGE_8,
		DYNAMICRANGE_4
	};
};
