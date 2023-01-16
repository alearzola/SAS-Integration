# NHANES - PAD
# 
# National Health and Nutrition Examination Survey (NHANES) - Peripheral Artery Disease (PAD)    
# 
# Data collection:
# Household screener, interview, and physical examination
# 
# Objectives:    
# Understand the survey data and create a predictive model to identify the
# main factors that are related to the disease. The model can also be     
# useful to prioritize the physical exams and to support the diagnostics. 
# 
# Activities:                                                              
# - Start the session                                                     
# - Prepare the data for Modelling
# - Data Partition (Training and Validation) 
# - Feature Engineering (add additional features)                                                   
# - Modelling                                                             
# - Scoring

####### SET UP THE CONFIGURATION FOR RStudio IN ORDER TO WORK IN CAS FROM R

## Run the following code on R to install following packages before running this:
# 1) install.packages("pkgbuild")
# 2) install.packages('jsonlite')
# 3) install.packages("tidyverse")

## Download and install RTools from:
# 4) https://cran.r-project.org/bin/windows/Rtools/

## Install the SWAT vX.X.X package as indicated here:
# 5) https://github.com/sassoftware/R-swat/releases

## Start the session and prepare the environment
library(swat)
library(ggplot2)

# Connect to CAS (this depends on the environment configuration)
your_host <- ' '
your_port <- 
your_username <- ' '
your_password <- ' '

conn <- swat::CAS(your_host, your_port, protocol='auto', username=your_username, password=your_password)

# CAS Server connection details
# out <- cas.builtins.serverStatus(conn)
# print(out)

### Import action sets
cas.builtins.loadActionSet(conn,actionSet="dataStep")
cas.builtins.loadActionSet(conn,actionSet="dataPreprocess")
cas.builtins.loadActionSet(conn,actionSet="cardinality")
cas.builtins.loadActionSet(conn,actionSet="sampling")
cas.builtins.loadActionSet(conn,actionSet="decisionTree")
cas.builtins.loadActionSet(conn,actionSet="astore")
cas.builtins.loadActionSet(conn,actionSet="percentile")

## Prepare the data for Modelling and assign variable to the table
path_to_data <- ' /nhanes_nof.sas7bdat'
cas.upload(conn,path_to_data,
           casOut=list(name="NHANES_NOF", caslib="CASUSER(alarzo)", replace=TRUE))

## Check the columns
# cas.table.columnInfo(conn,table="NHANES_NOF")

## Create the target variable
cas.dataStep.runCode(conn, code="
    data CASUSER.NHANES_PAD1 promote; 
          set CASUSER.NHANES_NOF;
            if LEXRABPI = . then LEXRABPI = LEXLABPI;
            if ((LEXLABPI < 0.9) OR (LEXRABPI< 0.9 )) then PAD_Target = 1;
                else PAD_Target = 0;
    run;"  )

## Data Partition (Training and Validation)
cas.sampling.srs(conn,
                 table="NHANES_PAD1",
                 samppct=30,
                 partind=TRUE,
                 output=list(casout = list(name="NHANES_PAD_PART", replace=TRUE),
                             copyvars="ALL")
                 )

## Present the partitions on a frequancy table and a bar chart
tablesize <- 6929
# With "CAS.TABLE.FETCH"-ACTION we extract a column from a in-memory table into
# R-Studio local-memory. Since it is a very long column, SAS pull it as a list
# with many elements. To get it as a unique vector, we use "UNLIST".
# Finally, we use only the numeric elements. For that reason we use "AS.NUMERIC"
pad_target_column <- as.numeric(unlist(
                        cas.table.fetch(conn, 
                                        table="NHANES_PAD_PART",
                                        fetchVars="PAD_TARGET",
                                        index=FALSE,
                                        to=tablesize) 
                        ))

# To produce frequency table of partition elements
freq_partition <- data.frame(table(pad_target_column))
freq_partition

# To produce bar chart of partition elements
ggplot(data=freq_partition, aes(x=pad_target_column, y=Freq, fill=pad_target_column)) + 
  geom_bar(stat="identity")+theme_minimal()

## Feature Engineering (add additional features)
cas.dataStep.runCode(conn,
                       code="data CASUSER.NHANES_PAD1(replace=yes); 
                              set CASUSER.NHANES_PAD_PART; 
                              PulsePreassure = BPXSAR - BPXDAR;
                              TC_HDL = LBXTC / LBDHDL;
                              IF ((DIQ010 In ('Yes','Borderline')) OR (DIQ050 In ('Yes')) OR (LBXGH > 6.5))
                                 then Diabetes = 1;
                                  else Diabetes = 0;
                              IF ( BPXSAR >= 140 OR BPXDAR >= 90 ) 
                                 then Hypertension = 1; 
                                  else Hypertension = 0;
                               run; ")

## Check the columns
# cas.table.columnInfo(conn,table="NHANES_PAD1")

## Modelling

# Specify the data set inputs and target 
interval_inputs <- c('RIDAGEMN_Recode', 'PulsePreassure', 'BMXBMI', 'TC_HDL', 'LBXGH', 'Diabetes', 'Hypertension')
class_inputs <- c('INDHHINC', 'DMDEDUC2', 'RIDRETH1', 'DIQ150', 'DIQ110', 'SMQ040', 'ALQ100', 'RIAGENDR')
class_vars   <- c('INDHHINC', 'DMDEDUC2', 'RIDRETH1', 'DIQ150', 'DIQ110', 'SMQ040', 'ALQ100', 'RIAGENDR', 'PAD_Target')
target       <- 'PAD_Target'

# Specify a generic cut-off
Gen_cutoff <- 0.5

# Train the model
cas.decisionTree.gbtreeTrain(conn,
                              table=list(name="NHANES_PAD1", where="strip(put(_PartInd_, best.))='0'"),
                              target=target,
                              inputs=c(class_inputs,interval_inputs),
                              nominals=class_vars,
                              nTree=150,  m=7,  lasso=0.777,  learningrate=1,  subsamplerate=0.883,  ridge=6.03,  seed=1634211770,
                              leafsize=5,  maxbranch=2,  binorder=TRUE,  encodename=TRUE,  mergebin=TRUE,  nBins=20,  maxLevel=6,
                              varImp=TRUE,  missing="USEINSEARCH",
                              casOut=list(name="gb_model", replace=TRUE)
)

## Score the model
cas.decisionTree.gbtreeScore(conn,
                             table=list(name="NHANES_PAD1"),
                             modelTable=list(name="gb_model"),
                             casOut=list(name='scored_gb', replace=TRUE),
                             copyVars=list("PAD_Target", "_PartInd_"),
                             encodename = TRUE,
                             assessonerow = TRUE
)

## Assess and compare model
assessed <- cas.percentile.assess(conn,
                                  table = 'scored_gb',
                                  inputs = 'P_PAD_Target1',
                                  casout = list(name = 'assessed', replace = TRUE),
                                  response = target,
                                  event = '1'
)

## We get 2 produced tables from the previous step that we can use to plot ROC and LIFT
assessed[["OutputCasTables"]]
cas.table.columnInfo(conn,table="assessed")
cas.table.columnInfo(conn,table="assessed_ROC")

## Plot ROC-Curve
specifity_ROC <- as.numeric(unlist(
  cas.table.fetch(conn, 
                  table="assessed_ROC",
                  fetchVars="_Specificity_",
                  index=FALSE,
                  to=tablesize) 
))
False_Positive_Rate <- 1-specifity_ROC

sensitivity_ROC <- as.numeric(unlist(
  cas.table.fetch(conn, 
                  table="assessed_ROC",
                  fetchVars="_Sensitivity_",
                  index=FALSE,
                  to=tablesize) 
))
True_Positive_Rate <- sensitivity_ROC

df <- data.frame(False_Positive_Rate, True_Positive_Rate)

ROC_Curve <- ggplot(df, aes(x = False_Positive_Rate, y = True_Positive_Rate)) + 
             geom_path(color='red') + geom_point(size = 2, color='red') + 
             geom_abline(slope = 1,linetype="dotted") +
             labs(title="ROC Curve for Gradient Boosting", x="False Positive Rate (1-Specificity)", y="True Positive Rate") 
ROC_Curve

## Plot Cumulative Lift Curve
depth_LIFT <- as.numeric(unlist(
  cas.table.fetch(conn, 
                  table="assessed",
                  fetchVars="_Depth_",
                  index=FALSE,
                  to=tablesize) 
))

cumulative_LIFT <- as.numeric(unlist(
  cas.table.fetch(conn, 
                  table="assessed",
                  fetchVars="_CumLift_",
                  index=FALSE,
                  to=tablesize) 
))

df2 <- data.frame(depth_LIFT, cumulative_LIFT)

LIFT_Curve <- ggplot(df2, aes(x = depth_LIFT, y = cumulative_LIFT)) + 
  geom_path(color='red') + geom_point(size = 2, color='red') +
  labs(title="Cumulative Lift Curve for Gradient Boosting", x="Depth", y="Cumulative Lift") 
LIFT_Curve

cas.terminate(conn)
