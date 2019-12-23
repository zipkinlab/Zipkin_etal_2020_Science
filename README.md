# [Tropical snake diversity collapses following widespread amphibian loss](https://xxx)

### Elise F. Zipkin, Graziella V. DiRenzo, Julie M. Ray, Sam Rossman & Karen R. Lips

### In review

### Code/Data DOI: TBD

### Please contact the first author for questions about the code or data: Elise F. Zipkin (ezipkin@msu.edu)
__________________________________________________________________________________________________________________________________________

## Abstract:
Biodiversity is declining at unprecedented rates worldwide. Yet, cascading effects of biodiversity loss on other taxa are largely unknown because baseline data are often unavailable. We document the collapse of a Neotropical snake community after the invasive fungal pathogen Batrachochytrium dendrobatidis caused a chytridiomycosis epizootic leading to the catastrophic loss of amphibians, a food source for snakes. Following mass mortality of amphibians, the snake community contained fewer species and was more homogeneous across the study site, despite no other systematic changes in the environment. The demise of the snake community following amphibian loss demonstrates the repercussive and often unnoticed consequences of the biodiversity crisis and calls attention to the invisible declines of rare and data deficient species.

## Code 
1. [community occurrence model](https://github.com/ezipkin/snake_community_model/community_occurrence_model/): Data management and manipulation code to run the community occurrence model. Also contains the code to process results. 


## Data
This project uses three datasets (which can all be found in the [data folder](https://github.com/ezipkin/snake_community_model/tree/master/data)).

1) effort 1997 to 2012 transect surveys for model_June2018.csv       
Description: This file contains information on each transect that was surveyed (listed by surveyID) durig the two time periods of our analysis (1997-2004 and 2006-2012). The file contains the date and time that each transect was surveyed as well as if the survey occured during the am or pm (ampm), the length surveyed (Meters), number of people who participated (People), and the number or person minutes (Pminutes). This file was used in both the community occurrence model as well as the composition analysis.

2) snake 1997 2012 dataset for body condition analysis_June2018.csv         
Description: This file contains a list of every snake encountered during both standardized surveys (i.e., those in the effort file) as well as those observed opportunistically during the two time periods of our analysis (1997-2004 and 2006-2012). When it was possible to capture snakes, their snout-to-vent lenght (SVL; cm) and mass (g) was recorded. Note that not every sample has these data available and that efforts to capture snakes to collect this information increased post-epizootic. This data file was used for the body condition analysis.

3) snake 1997 2012 dataset for body condition analysis_June2018.csv       
Description: This file contains a list of each snake encountered during the standardized surveys (i.e., the surveys listed in the effort file by surveyID) during the two time periods of our analysis (1997-2004 and 2006-2012). Each row represents a detected individual, its species and the date, time, and transect on which it was observed. This file was used in both the community occurrence model as well as the composition analysis.
