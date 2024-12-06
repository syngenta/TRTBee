# TRTBee
Scripts and application for the BYOM based TRT calculator for bee risk assessment

TRTBee installation guide
To install the Standalone software
1.	Navigate to https://github.com/syngenta/TRTBee/releases/ and select the relevant version of the app. Unless replicating a previous use, this will be the latest release.
2.	Download TRTBee_installer.7z
3.	Extract the .exe installation file. To extract the installer requires a program that can extract 7z files. We recommend the free, open-source program 7-zip. Available from https://www.7-zip.org/
4.	Open the .exe file and follow the instructions in the installer.
5.	You will then have the program installed and be able to run TRT assessments.

Using the software
1.	Open the installed TRTBee software, or the Matlab app version
2.	Click the (Re)start button to prepare the software
3.	Click “1 - Load in mortality data” and navigate to the mortality data that you wish to use. Ensure that the data is correctly formatted
4.	Once you can see the data populating the main table, click “2 – Run TRT assessment”. This step may take several minutes, especially for large datasets or examples with time-variable exposure concentrations
5.	Once complete, a summary of the results will be presented in the dialog box under the main table. The “Save as” button can now be clicked, selecting a file location will produce a .csv of the results in that folder.
There is an option to Save output figures. These figures show the comparison of the GUTS model fits to the data and likelihood plots. If the option is turned on, these figures will appear in the same location as the data that was loaded into the tool.

To install and run the Matlab based version:
1.	The tool cannot be guaranteed to work with versions of Matlab older than R2016b. Older versions may need to be updated.
2.	Ensure that you have the BYOM package downloaded and in the Matlab path. BYOM can be downloaded from: https://debtox.info/byom.html
3.	Download the source code folder from https://github.com/syngenta/TRTBee/releases/, extract and place within the BYOM folder.
4.	Open and run the “run_fit_BYOM.m” script and follow the instructions.
