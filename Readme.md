# A set of tools to process data collected for the Treadmill Oscillation Walking (TOW) project, including a script to compensate for the inertial force induced by accelerating an instrumented treadmill with visual quality check, a GUI to select and inspect EMG signal from MVIC trials for walking trial normalization, a pipeline to transform CoP collected from a moving instrumented treadmill to the Lab coordinate system where kinematics were represented in, and extract participant bodyweight from the static calibration trial.

**Note: Using the tool requires installation of [ezc3d package](https://github.com/pyomeca/ezc3d)**


## File:
1. **c3d EMG_filtering**: Get MVIC files and extract MVIC with manual selection, filter EMG from gait trials and normalize to MVIC

	- quality check completed 4/5/2023, merged from folder EMG_filtering_test

2. **Inertial_Compensation**: Get Unload_OS trials, match R_tread1X marker and subtract force from ..OS0.. trials, lowpass filter all trials

	- quality check completed 4/6/2023

3. **FP_LAB_Match**: use 4 treadmill markers to compute *(1) X offset of 4 markers between static trial in the current folder and the standard trial (s02 origin)*, *(2) Tread coordinate system defined using L1,3 and R1 in current trial, and FP origin location in tread CS in standard trial.* 
	- Combined to get FP origin in current trial and see how it deviates from (0,0,0)
	- Output from this script determines the FP segment origin landmark that should be input into Visual3D static model.

	- quality check completed 4/5/2023

4. **vGRFCadenceTest_validation**: validate the cadence and vGRF asymmetry calculated in D-flow and save a pdf for vGRF
	- Note: the force is filtered in this script, while not in D-flow

	- quality check completed 4/5/2023

	- 4/12/2023 found figure off, bug in the code. Will fix later. Skip for now and manually check during data collection.

5. **Data_processing_tool_runner**: Run files 1-4 and save command window output.

 	- quality check completed 4/6/2023

6. **Static_BW_validation**: using Fz from static trail to validate if BW is off. Visualize and display values and their differences.


## How to use:
1. Copy all 5 files to each folder containing original .c3d files from Vicon. 
2. Open "Data_processing_tool_runner" 
3. Check the filename, filter parameter, and EMG sensor number
4. Run "Data_processing_tool_runner"
5. In the Force processing session, each "OS" file will spit out compensation result. Close figure(3) will continue the for loop.
6. In the EMG prosession session, a MVIC manual checking GUI will pop up. 
	Select a new MVIC if needed
	click "store value" when finished selecting to store the value and the graph.
	click "Done selecting" to close the figure and resume the code.
7. The vGRFCadenceTest_validation will show cadence and vGRF. Compare with Google sheet (from D-flow).
8. The FP_LAB_Match will show FP origin offset. Put the value into Visual 3D FP segment origin lanmark if largely deviated from zero.
9. Use EMGfiltered...c3d to create the final .cmz file


## Note:
4/6/2023: Decide to still compensate all six channels of force
8/19/2023: some post-training session has no weight scale measurement, **added Static_BW_validation module** using Fz from static trail to validate if BW is off.
