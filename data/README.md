### This folder contains all of the data files to conduct the different analyses.
__________________________________________________________________________________________________________________________________________

## Data
Our project uses the following three datasets:

1) transect_survey_effort.csv       
Description: This file contains a list of all standardized transect survey events conducted during the two time periods: pre-Bd (Dec 1997 - Dec 2004) and post-Bd (Sept 2006 - July 2012). Each survey event (row) was assigned a unique ID (surveyID), ordered by year (Year), month (Mo), day (Day), transect name (Transect), and time of day (AM or PM; ampm). Also included is the length surveyed (Meters), number of people who participated in the survey event (People), and the number or person minutes surveyed (Pminutes). This file was used in both the community occurrence model as well as the composition analysis.

2) snake_body_condition_data.csv         
Description: This file contains a list of every snake encountered during both standardized surveys (i.e., those in the effort file) as well as those observed opportunistically during the two time periods of our analysis: pre-Bd invasion (Dec 1997 - Dec 2004) and post-Bd (Sept 2006 - July 2012). For each snake observation (row), we include the genus and species of the snake and the year, month, day, transect, and time of day (ampm and/or time) that the individual was observed. When it was possible to capture snakes, their snout-to-vent lenght (SVL; cm) and mass (g) was recorded. Note that not every individual encountered has these data available and that efforts to capture snakes to collect this information increased post-epizootic. This data file was used for the body condition analysis.

3) snake_occurrence_data_from_transects.csv       
Description: This file contains a list of each snake encountered during the standardized surveys (i.e., the transect survey events listed in the effort file) during the two time periods of our analysis (1997-2004 and 2006-2012). Each row represents a detected individual and includes information on the species and when it was observed (surveyID, date, and transect). Snout-to-vent lenght (SVL; cm) and mass (g) are included when it was possible to obtain those measurements. This file was used in both the community occurrence model as well as the composition analysis. 
