#This R script contains all of the R code to do the (hopefully this time
# for sure) final analyses for Andre's bacterial transmission manuscript.
# The manuscript is being reformulated so that it is centered around a
# paper by Alizon that suggests that four relationships need to be
# described in order to understand the direction in which the virulence
# of a pathogen should evolve. These four relationships are (x versus y):
#  (1) pathogen replication rate vs virulence
#  (2) pathogen replication rate vs transmission rate
#  (3) virulence vs transmission rate
#  (4) recovery rate vs virulence
#
# For various reasons, Andre also wants two more regressions:
#  (5) virulence vs recovery rate
#  (6) virulence vs fitness (where fitness = duration of infection X transmission
#       rate)
# In all cases, we need to test whether there are linear or quadratic
# relationships between predictor and response.
#
#Andre is proposing to use the following measured responses to 
# experimental infection in order to fit the four models:
#  (1) replication rate: summed bacterial load to d18 post-inoculation 
#       (or the closest measurement day to d18, which could be d17 or d19
#       as well) [I need to create estimates from the raw data]
#  (2) virulence: summed eye scores to day 18 post-inoculation [I need to
#       create estimates from the raw data]
#  (3) transmission rate: isolate-specific values estimated by Steve Ellner
#  (4) recovery rate: duration of infection based on eye score [I need to
#       create estimates from the raw data]
#
#This script was written by Wesley Hochachka, starting on 9 October 2023.



## Set Up, Read Data into R, and Initial Formatting

#Set R up by reading in libraries and specify the working directory.
#
library(tidyverse)
library(readxl)
#
#This is the office computer working directory.
WorkingDir <- "D:/users/Wes/Projects/CLOnonCitSci/HOFI_related/Dhondt_transmission_2017/FourRelationshipsAnalyses_Oct2023"
#This is the laptop working directory.
WorkingDir <- "C:/Users/Wes/November2023Work/Dhondt_analyses/FourRelationshipsAnalyses_Oct2023"
#
#setwd(WorkingDir)


#Read in raw data files. The second of these files will require substantial 
# manipulation before it can be used.
TransmissionRate <- read_xlsx("beta.xlsx",
                              sheet = "Sheet1")
#
#Commented out are the changes needed for the previous version of the second
# data table, left for the sake of documentation in case there is need to
# ever return to this file.
#OtherRawDataInput <- read_csv(file = "disease_trajectory_with_day_0_ added D_updated.csv") %>%
#  rename(LOGqPCR = 'log(qPCR)') %>%
#  rename(EyeScore = 'eyescore(LR)') %>%
#  rename(dayPI = 'day PI') %>%
#  mutate(isolate = replace(isolate, isolate == "CA 06", "CA06"),
#         region = replace(region, region == "E" & isolate == "CA06", "W"))
#
#Read in the most recent version of the second table of raw information, and
# then do a number of manipulations in order to turn the information into a
# useful R data object:
#  - rename all variables that have spaces in their names, to remove the spaces
#  - rename a few other variables to names that for me are more explanatory
#  - rename the "E-W" column to "E_W" because R seems to interpret the hyphen as
#     as a subtraction sign
#  - remove the RPA column because it appears irrelevant and also so that I do
#     not want to fix the use of "nd" to replace it with NA values and turn the
#     column into a numeric variable (upon reading in, it is character)
#  - add the sex for bird "W.G index" from set D with the CA2008 innoculum
#     on its day 4 record; hard code this missing sex to "F"
#  - remove some blank rows at the bottom that Excel things are information but
#     are totally blank lines with which R shouldn't be made to deal
#  - deal with the issue of encoding the first sampling period on which the a
#     bird was dead by writing "DEAD" or "dead" in the right-eye eye-score
#     column; some sort of first-day-dead column needs to be created
#There's still one issue that Andre needs to clarify: 
#  - how is it possible to have a qPCR value of "??", and yet have a log-qPCR 
#     numeric value?
#  - what is the single missing "dayPI" value for bird "YB._ (index F)" from
#     CA2010 set B experiment; actually what values should be in all of the 
#     columns?  It looks like the bird was known dead on the D71 sampling
#     period, so presumably this row can be deleted?  Actually, there are at
#     least 2 birds with almost entirely blank rows for sampling periods after
#     they died; the second one that I found was BirdID "R.WW (index)"
#  - there are apparently 2 separate birds banded "W.G" used in separate sets
#     (experiment sets "C" and "D"). One is named "W.G (index)" and the other
#     "W.G index", used with different innocula. Double check that these are 
#     are separate birds (they almost certainly are) and then name sure to
#     never summarize data only based on the band combination, but include
#     other grouping variables like the isolate of the inoculum
#First we do all of the corrections changes that each require a single line of
# R code.
OtherRawDataInput <- read_xlsx("index birds to DPI around 60.xlsx",
                              sheet = "all") %>%
  rename(isolate = "inoculum") %>%
  rename(BirdID = "Bird ID") %>%
  rename(Eye_R = "Eye R") %>%
  rename(Eye_L = "Eye L") %>%
  rename(LOGqPCR = 'log MGC2') %>%
  rename(SumEyeScore = 'eye_sum') %>%
  rename(dayPI = 'DPI') %>%
  rename(region = "E-W") %>%
  select(-RPA) %>%
  filter(!is.na(region)) %>%
  mutate(Sex = replace(Sex, is.na(Sex), "F"))
#Now add a binary variable, IsAlive, to indicate whether a bird is alive or dead
# during each of the sampling periods.  "1" means that the bird is alive during
# the sampling period, and "0" means that the bird is dead. The creation of
# this variable depends on the text string either "DEAD" or "dead" to have been
# placed into the Eye_R column on the first sampling periods on which the bird
# was known to be dead.
#Also do a bit more filtering to remove birds that died early in the experiments.
# I have talked with Andre, and for the bird that died before the d18 sampling
# periods we will just remove all of that bird's data. For the birds that died
# after the d18 sampling period, I will use the bird's data to estimate
# replication rate and virulence, but *not* duration of infection or recovery
# rate. I will need to tweak the data file to nuke the pre-d18 dying bird, and
# then I think create a boolean flag to indicate whether a bird's data should
# be used to calculating duration of infection and recovery rate.
#Add the IsAlive boolean.
FirstDayDead <- OtherRawDataInput %>% 
  filter(Eye_R == "DEAD" | Eye_R == "dead") %>%
  select(Set, BirdID, isolate, Eye_R, dayPI) %>%
  rename(FirstDayDead = "dayPI") %>%
  mutate(IsAlive = 0)
OtherRawDataWorking <- left_join(x = OtherRawDataInput,
                                 y = FirstDayDead,
                                 by = c("Set", "BirdID", "isolate")) %>%
  select(-Eye_R.y) %>%
  mutate(IsAlive = case_when(!is.na(FirstDayDead) & dayPI < FirstDayDead ~ 1,
                             !is.na(FirstDayDead) & dayPI >= FirstDayDead ~ 0,
                             is.na(FirstDayDead) ~ 1)) %>%
  rename(Eye_R = Eye_R.x)
#Add the DiedEarly boolean.
DiedEarly = OtherRawDataWorking %>%
  filter((Eye_R == "DEAD" | Eye_R == "dead") & dayPI < 24 & IsAlive == 0) %>%
  select(Set, BirdID, isolate) %>%
  mutate(DiedEarly = 1)
OtherRawDataWorking <- left_join(x = OtherRawDataWorking,
                                 y = DiedEarly,
                                 by = c("Set", "BirdID", "isolate")) %>%
  mutate(DiedEarly = replace(DiedEarly, is.na(DiedEarly), 0))
#
#Tidy up.
rm(DiedEarly, FirstDayDead)
                                 
 
  



## Summarize Raw Data

#Now that we have the raw data in a usable form, it is time to summarize down
# these raw data, in which each observation day is a separate record, into 
# summary statistics for each bird. We need to have the following 4 summaries
# from each bird's data:
#  (1) summed bacterial load to d18 post-inoculation (replication rate index)
#  (2) summed eye scores to day 18 post-inoculation (virulence index)
#  (3) duration of infection based on eye score (recovery rate index)
#
#
#Create the index of (1), replication rate. This is a summation of the 
# bacterial load values up to "day-18" (the actual last day may be either
# day-17 or day-19 too). There is one bird (I think just one) that died
# before this final sampling day, and all data from this bird have already
# been removed before this summarization process...but just in case this 
# chunk of code will remove any birds that I might have missed. Note that I am
#  also removing the value for day-0, because this represents inoculation and
# not replication.
#
LastSampleDate <- OtherRawDataWorking %>%
  group_by(isolate, region, BirdID) %>%
  summarise(LastSamplePeriod = max(dayPI[!is.na(LOGqPCR)])) %>%
  ungroup
ReplRateIndex <- left_join(x = OtherRawDataWorking,
                           y = LastSampleDate,
                           by = c("isolate", "region", "BirdID")) %>%
  filter(dayPI > 0 & dayPI <= 19,
         LastSamplePeriod >= 17) %>%
  group_by(isolate, region, BirdID) %>%
  summarise(sumqPCR = sum(LOGqPCR),
            nSampPeriods = n()) %>%
  ungroup()
rm(LastSampleDate)


#Create (2), the index of virulence, which is the summed eye scores up to and 
# including the "day-18" sampling periods (i.e. up to day-19 to include 
# 3 sampling periods for all birds.)
LastSampleDate <- OtherRawDataWorking %>%
  group_by(isolate, region, BirdID) %>%
  summarise(LastSamplePeriod = max(dayPI[!is.na(LOGqPCR)])) %>%
  ungroup
VirulenceIndex <- left_join(x = OtherRawDataWorking,
                           y = LastSampleDate,
                           by = c("isolate", "region", "BirdID")) %>%
  filter(dayPI > 0 & dayPI <= 19,
         LastSamplePeriod >= 17) %>%
  group_by(isolate, region, BirdID) %>%
  summarise(SumEyeScoreToD18 = sum(SumEyeScore),
            nSampPeriods = n()) %>%
  ungroup()
#
#Tidy up.
rm(LastSampleDate)



#Create (3) duration in infection (i.e. the "recovery rate" index). Recovery is
# defined based on having an eye score of 0, following some period of non-zero
# eye score (regardless of qPCR detection of bacteria). "Death" is defined as
# having occurred on the day after an eye score of 6 (i.e. both eyes highly
# swollen and the bird being essentially blind). All birds that died before the
# end of the experiment are removed from the set of data used.
# I STILL NEED TO MAKE SURE THAT BIRDS THAT NEVER HAD EYE SCORES > 0 ARE 
# INCLUDED IN THE OUTPUT DATA!
#Determine the first day of disease (non-zero eye score) for every bird
# that ever showed signs of disease.
FirstDiseaseDay <- OtherRawDataWorking %>%
  filter(DiedEarly == 0,
         dayPI > 0,
         SumEyeScore > 0) %>%
  group_by(region, isolate, BirdID) %>%
  summarise(FirstDiseaseDay = min(dayPI)) %>%
  ungroup()
#
#Determine the first day of an eye score of six, when this occurs, and set
# the next day as the day of "death".
FirstDayEye6 <- OtherRawDataWorking %>%
  filter(dayPI > 0,
         SumEyeScore == 6) %>%
  group_by(region, isolate, BirdID) %>%
  summarise(DeathDay = min(dayPI) + 1) %>%  #"death" on day after eye score of 6
  ungroup()
#
#Identify the last day of disease so that we can calculate the interval of
# disease.
LastDiseaseDay <- OtherRawDataWorking %>%
  filter(DiedEarly == 0,
         dayPI > 0,
         SumEyeScore > 0) %>%
  group_by(region, isolate, BirdID) %>%
  summarise(LastDiseaseDay = max(dayPI)) %>%
  ungroup()
#
#Identify the birds that never showed any signs of disease.  After extracting
# these rows, add some columns as fillers to match the columns found for
# the birds that had non-zero eyescores on at least one day.  We will need
# this filler information to concatinate the data from the never-diseased birds
# with the data from the diseased birds.
NoDiseaseBirds <- OtherRawDataWorking %>%
  filter(DiedEarly == 0,
         dayPI > 0) %>%
  group_by(region, isolate, BirdID) %>%
  summarise(NoDiseaseDay = max(SumEyeScore)) %>%
  ungroup() %>%
  filter(NoDiseaseDay == 0) %>%
  mutate(FirstDiseaseDay = NA,
         LastDiseaseDay = NA,
         DeathDay = NA,
         DisDuration = 0) %>%
  select(-NoDiseaseDay)
#
#Now combine all of this information to calculate  the duration of disease.
# The logic to the end of the first "case_when" statement does not work when there
# was only a single  sampling day on which disease was observed, which causes
# the calculations to return a duration of disease of zero days. In these few
# cases we need to manually an arbitrary half of the interval between the
# single day on which disease was observed and the next sampling period on which
# disease was not observed; the addition is of 3 days (it could be 4...or any
# other reasonable but arbitrary number). Finally, add back in the data  from
# the birds that  were never diseased (i.e. that have a disease duration of
# zero).
TimeToRecoverDie <- full_join(x = FirstDiseaseDay,
                              y = LastDiseaseDay,
                              by = c("region", "isolate", "BirdID")) %>%
  left_join(x = .,
            y = FirstDayEye6,
            by = c("region", "isolate", "BirdID")) %>%
  mutate(DisDuration = case_when(LastDiseaseDay < DeathDay ~ LastDiseaseDay - FirstDiseaseDay,
                                 LastDiseaseDay > DeathDay ~ DeathDay - FirstDiseaseDay,
                                 is.na(DeathDay) ~ LastDiseaseDay - FirstDiseaseDay)) %>%
  mutate(DisDuration = case_when(DisDuration == 0 ~ FirstDiseaseDay + 3,
                                 TRUE ~ DisDuration))
#Concatenate the columns of data from the non-diseased birds to the end.
TimeToRecoverDie <- bind_rows(TimeToRecoverDie, NoDiseaseBirds)
#
#Tidy up.
rm(FirstDiseaseDay, LastDiseaseDay, NoDiseaseBirds, FirstDayEye6)




#Now that all of the necessary tables of summary data are available
# (TransmissionRate, ReplRateIndex, TimeToRecoverDie, and VirulenceIndex)
# the relationships can be examined.



  