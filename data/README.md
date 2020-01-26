### This folder contains all of the data files to conduct the different analyses.
__________________________________________________________________________________________________________________________________________

## Data
Our project uses the following three datasets:

1) transect_survey_effort.csv       
Description: This file contains information on each transect that was surveyed (listed by surveyID, a unique survey event identification number) during the two time periods of our analysis (1997-2004 and 2006-2012). The file contains the date and time that each transect was surveyed as well as if the survey occurred during the am or pm (ampm), the length surveyed (Meters), number of people who participated (People), and the number or person minutes surveyed (Pminutes). This file was used in both the community occurrence model as well as the composition analysis.

2) snake_body_condition_data.csv         
Description: This file contains a list of every snake encountered during both standardized surveys (i.e., those in the effort file) as well as those observed opportunistically during the two time periods of our analysis (1997-2004 and 2006-2012). When it was possible to capture snakes, their snout-to-vent lenght (SVL; cm) and mass (g) was recorded. Note that not every individual encountered has these data available and that efforts to capture snakes to collect this information increased post-epizootic. This data file was used for the body condition analysis.

3) snake_occurrence_data_from_transects.csv       
Description: This file contains a list of each snake encountered during the standardized surveys (i.e., the surveys listed in the effort file by surveyID) during the two time periods of our analysis (1997-2004 and 2006-2012). Each row represents a detected individual and includes information on the species and surveyID (including the date and transect) on which it was observed. Snout-to-vent lenght (SVL; cm) and mass (g) are provided when they were recorded. This file was used in both the community occurrence model as well as the composition analysis. 
