# Hayes-Kaye_SELECT2.0_manuscript
Data and code for the SELECT 2.0 manuscript by Hayes and Kaye et al., 

Data and code are deposited here and at Zenodo DOI: 
The preprint DOI is here: 

## Outline 
This repository contains the data and code for the manuscript 'SELECT 2.0: Refined and open access SELection Endpoints in Communities of bacTeria (SELECT) method to determine concentrations of antibiotics that may select for antimicrobial resistance in the environment' by Hayes & Kaye et al., XX. The code here can recreate all analyses, figures and tables. 

## Information about the files

### Antibiotic OD Files 
The OD files needed for the code to run are  deposited here. Each csv file is titled by antibiotic. These are processed with _**01_select2.0**_
- conc - concentration of antibiotic (mg/L)
- replicate - in plate replicate per concentration
- columns labelled 0-72 (minutes) - OD readings per ten minute interval

### For summary figures and comparison between SELECT 1.0 and 2.0
- _**results.csv**_: this file contains four columns that are then used in script  _**02_comparison_between_methods**_ to generate PNECRs and to compare between the SELECT 1.0 and 2.0 methods
  - antibiotic - all antibiotics in the study
  - class - class of antibiotic
  - select1.0 - SELECT 1.0 LOEC
  - select2.0 - SELECT 2.0 EC1   

### Environmental Risk Assessments
- _**uba_mecs_transposed.csv**_: contains all measured environmental concentrations used in the global environmental risk assessments
  - antibiotic - antibiotic
  - water_type - wastewater influent or wastewater effluent
  - all other columns contain μg/L concentrations (unlabelled)
   
- _**cip3_mecs**_: contains all measured environmental concentrations used in the England and Wales environmental risk assessments
  - SampleLocationName - wastewater type
  - SampleValue - μg/L concentrations
  - NameDeterminandName - antibioitc
  - BelowMinReading - whether below the minimum reading for quantification (Y/N)     

- _**naoh_control.txt**_: contains OD readings for NaOH solvent control
  - Kinetic_read - time
  - A1-H12 - OD readings for each well 

### Code
Each file is named for the analyses within
- _**01_select2.0**_: this script will generate SELECT 1.0 LOECs and PNECRs, and 2.0 EC1s and PNECRs for all antibiotics tested in this study
- _**02_comparison_between_methods**_: this script determines the difference between the SELECT 1.0 and 2.0 PNECRs for the antibiotics in this study. It also allows you to perform a Bland-Altman analysis.
- _**03_environment_risk_assessment**_: this script runs the global and England and Wales environmental risk assessments
- _**04_solvent_controls**_: this script provides the analysis for the solvent controls needed for this study (acetic acid and 0.5M NaOH)
