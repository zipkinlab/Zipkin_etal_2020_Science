# [Tropical snake diversity collapses after widespread amphibian loss](https://xxx)

### Elise F. Zipkin, Graziella V. DiRenzo, Julie M. Ray, Sam Rossman & Karen R. Lips

### In press

### Code/Data DOI:    TBD

### Please contact the first author for questions about the code or data: Elise Zipkin (ezipkin@msu.edu)
__________________________________________________________________________________________________________________________________________

## Abstract:
Biodiversity is declining at unprecedented rates worldwide. Yet, cascading effects of biodiversity loss on other taxa are largely unknown because baseline data are often unavailable. We document the collapse of a Neotropical snake community after the invasive fungal pathogen Batrachochytrium dendrobatidis caused a chytridiomycosis epizootic leading to the catastrophic loss of amphibians, a food source for snakes. Following mass mortality of amphibians, the snake community contained fewer species and was more homogeneous across the study site with several species in poorer body condition, despite no other systematic changes in the environment. The demise of the snake community following amphibian loss demonstrates the repercussive and often unnoticed consequences of the biodiversity crisis and calls attention to the invisible declines of rare and data deficient species.

## Code 
1. [body_condition_analysis](./body_condition_analysis/): This folder contains the code to run the body condition analyses for the six species with a sufficient number of samples pre- and post-epizootic.

2. [community_occurrence_model](./community_occurrence_model/): This folder contains the code to run the community occurrence model. The wrapper includes code to calcuate species richness, individual species occurrence rates, and detection rates for both observed and unobserved species.

3. [composition_analysis](./composition_analysis/): This folder contains the code to run the composition analysis using data from observed species.


## Data
This project uses three datasets (all found in the [data folder](./data)).


1) transect_survey_effort.csv       
Description: This file contains a list of all standardized transect survey events conducted during the two time periods: pre-Bd (Dec 1997 - Dec 2004) and post-Bd (Sept 2006 - July 2012). Each survey event (row) was assigned a unique ID (surveyID), ordered by year (Year), month (Mo), day (Day), transect name (Transect), and time of day (AM or PM; ampm). Also included is the length surveyed (Meters), number of people who participated in the survey event (People), and the number or person minutes surveyed (Pminutes). This file was used in both the community occurrence model as well as the composition analysis.

2) snake_body_condition_data.csv         
Description: This file contains a list of every snake encountered during both standardized surveys (i.e., those in the effort file) as well as those observed opportunistically during the two time periods of our analysis: pre-Bd invasion (Dec 1997 - Dec 2004) and post-Bd (Sept 2006 - July 2012). For each snake observation (row), we include the genus and species of the snake and the year, month, day, transect, and time of day (ampm and/or time) that the individual was observed. When it was possible to capture snakes, their snout-to-vent lenght (SVL; cm) and mass (g) was recorded. Note that not every individual encountered has these data available and that efforts to capture snakes to collect this information increased post-epizootic. This data file was used for the body condition analysis.

3) snake_occurrence_data_from_transects.csv       
Description: This file contains a list of each snake encountered during the standardized surveys (i.e., the transect survey events listed in the effort file) during the two time periods of our analysis (1997-2004 and 2006-2012). Each row represents a detected individual and includes information on the species and when it was observed (surveyID, date, and transect). Snout-to-vent lenght (SVL; cm) and mass (g) are included when it was possible to obtain those measurements. This file was used in both the community occurrence model as well as the composition analysis. 
