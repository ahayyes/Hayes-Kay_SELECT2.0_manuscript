# Hayes-Kaye_SELECT2.0_manuscript
Data and code for the SELECT 2.0 manuscript by Hayes and Kaye et al., 

Data and code are deposited here and at Zenodo DOI: 
The preprint DOI is here: 

## Outline 
This repository contains the data and code for the manuscript 'SELECT 2.0: Refined and open access SELection Endpoints in Communities of bacTeria (SELECT) method to determine concentrations of antibiotics that may select for antimicrobial resistance in the environment' by Hayes & Kaye et al., XX. The code here can recreate all analyses, figures and tables. 

## Information about the files

### OD Files
The OD files needed for the code to run are also deposited here. Each csv file is titled by antibiotic. Here, the columns read the concentration, replicate number, and subsequent columns are time points. 

### Environmental Risk Assessments
For the global and England and Wales environmental risk assessments, the CSV files .......

**Code**
Each file is named for the analyses within
- _**01_select2.0**_: this script will generate SELECT 1.0 LOECs and PNECRs, and 2.0 EC1s and PNECRs for all antibiotics tested in this study
- _**02_comparison_between_methods**_: this script determines the difference between the SELECT 1.0 and 2.0 PNECRs for the antibiotics in this study. It also allows you to perform a Bland-Altman analysis.
- _**03_environment_risk_assessment**_: this script runs the global and England and Wales environmental risk assessments
- _**04_solvent_controls**_: this script provides the analysis for the solvent controls needed for this study (acetic acid and 0.5M NaOH)
