# This code is to calibrate SWAT model with GLUE using RSWAT
args <- commandArgs(trailingOnly=TRUE)
#                              Start user input                                #
# -----------------------------------------------------------------------------#

# 1. Main directory (containning all other folders)
#HPCPath <- "/nesi/nobackup/niwa03253/Thanh/rswat/"
HPCPath <- "/nesi/nobackup/pletzera/rswat2/"

# 2. RSWAT source code
RSWATsourceFile <- paste0(HPCPath,"R-SWAT/R")

# 3. TRUE = SWAT, FALSE = SWAT+
SWATproject <- TRUE

# 3. Working folder (all SWAT simulations are saved here
workingFolder <- paste0(HPCPath,"workingFolder")
if (length(args) >= 1) {
  # let the user specify a work directory
  workingFolder <- args[1]
}
print(sprintf("working folder: %s", workingFolder))

# 4. TxtInOut folder (SWAT setup folder)
TxtInOutFolder <- paste0(HPCPath,"TxtInOut_Hauraki_simple")

# 5. SWAT (or SWAT+) executable file
SWATexeFile <- paste0(HPCPath,"SWAT/build/src/swat2012.685.ifort.rel")

# 6. SWAT (or SWAT+) parameter file (only modify swatParam.txt if you want to add more parameters)
SWATParamFile <- paste0(HPCPath,"TxtInOut_Toenepi/swatParam.txt")

# 7. Select parameter for calibration or sensitivity analysis with SWAT
paraSelection <- data.frame(
  Parameter = c("CN2.mgt" , "ALPHA_BF.gw", "CH_N2.rte", "CH_K2.rte", "ESCO.hru", "EPCO.hru", "SURLAG.bsn", "TDRAIN.mgt", "GDRAIN.mgt", "SOL_AWC.sol", "SOL_Z.sol", "SOL_K.sol", "OV_N.hru"),
  Change =    c("relative", "replace"    , "replace"  , "replace"  , "replace" , "replace" , "replace"   , "replace"   ,  "replace"  , "relative"   , "relative" , "relative" , "replace"),
  Min =       c(-0.25     , 0            ,  0.001     ,  0         ,  0        ,  0        ,  0.05       ,  12         ,  12         , -0.25        , -0.25      , -0.25      ,  0.01),
  Max =       c(0.25      , 1            ,  0.2       ,  150       ,  1        ,  1        ,  24         ,  48         ,  48         ,  0.25        , 0.25       , 0.25       ,  30),
  Subbasin =  c("All"     , "All"        ,  "All"     ,  "All"     ,  "All"    ,  "All"    ,  "All"      ,  "All"      ,  "All"      ,  "All"       , "All"      , "All"      ,  "All"),
  Landuse =   c("All"     , "All"        ,  "All"     ,  "All"     ,  "All"    ,  "All"    ,  "All"      ,  "All"      ,  "All"      ,  "All"       , "All"      , "All"      ,  "All"),
  Soil =      c("All"     , "All"        ,  "All"     ,  "All"     ,  "All"    ,  "All"    ,  "All"      ,  "All"      ,  "All"      ,  "All"       , "All"      , "All"      ,  "All"),
  Slope =     c("All"     , "All"        ,  "All"     ,  "All"     ,  "All"    ,  "All"    ,  "All"      ,  "All"      ,  "All"      ,  "All"       , "All"      , "All"      ,  "All")
)

# 8. Select parameter for calibration or sensitivity analysis with SWAT+
# paraSelection <- data.frame(
#   Parameter = c("cn2.hru"   , "canmx.hru"),
#   Change =    c("relative"  , "replace"),
#   Min =       c(-0.25       , 1.0),
#   Max =       c(0.25        , 10.0),
#   Object  =   c("All"       , "All"),
#   Conditions = c("All"      , "All")
# )

# 9. Parameter calibration/sensitivity analysis approach
samplingApproach <- "Cali_(Generalized_Likelihood_Uncertainty_Estimation)"

# 10. Additional information about parameter calibration sensitivity analysis
sensCaliCommand <- 10 

# 11. Output extraction for SWAT
outputExtraction <- data.frame(
  FileType = c("output.rch"),    # if for two files: = c("watout.dat", "output.rch"),
  FileName = c("output.rch"),    #                   = c("watout.dat", "output.rch"),
  Column = c("7"),               #                   = c("4"         , "6"),
  Reach = c("369")                 #                   = c(" "         , "2")
)

# 12. Output extraction for SWAT+
# outputExtraction <- data.frame(
#   FileType = c("channel_sd_day.txt"),
#   FileName = c("channel_sd_day.txt"),
#   Column = c("48"),
#   Reach = c("1")
# )

# 13. Date range for extraction
dateRangeCali <- as.Date(c("2004-01-01", "2012-12-31"), format = "%Y-%m-%d")

# 14. Number of parallel runs
ncores <- 21

# 15. Objective function could be "NSE", "KGE", "RMSE", "R2", "aBIAS" is absolute bias ranging from [0, inf]
objFunction <- "NSE"

# 16. Observed data file(s)
observedDataFile <- c(paste0(HPCPath,"TxtInOut_Hauraki_simple/obs_var_1.txt"))

# 17. Behavioral threshold
behThreshold <- -0.2


# -----------------------------------------------------------------------------#
#                            End user input                                    #
# -----------------------------------------------------------------------------#
# Load R-SWAT functions
# Require packages
requiredPackages <- c('foreach', 'doParallel', 'lhs', 'sensitivity', 'boot', 
                      'optimization', 'nloptr', 'spsComps', 'doMPI') #'hydroPSO',

# Check and install only missing packages
install.packages(setdiff(requiredPackages,rownames(installed.packages())), dependencies = TRUE, repos = "http://cran.us.r-project.org") 

# Load these packages
lapply(requiredPackages, library, character.only = TRUE)

# Load R-SWAT functions
RSWATfiles <- c("executeSWAT.R", "glue.R", "multiRegression.R",
                "objFunction.R", "readSaveOutput.R", "updateTxtInOut.R",
                "userObjFunction.R", "userReadSwatOutput.R", "displayOutput.R")

# Load R-SWAT function
setwd(RSWATsourceFile)
lapply(RSWATfiles, source)

# Generate parameter samples (uniform distribution)
parameterValue <- runifSampling(sensCaliCommand, getParamRange(paraSelection)[,1],getParamRange(paraSelection)[,2])

# Get HRU infor (HRU names, landuse, soil, slope, subbasin)
if (SWATproject){
  HRUinfo <- getHruInfo(TxtInOutFolder)
} else {
  HRUinfo <- read.table( paste(TxtInOutFolder, "/hru-data.hru", sep =""),
                         header = TRUE, skip = 1, sep = "")
}

# Read SWAT parameter
SWATParam <- loadSwatParam(SWATParamFile)

# Get location of parameters in TxtInOut files and load TxtInOut file content
caliParam <- loadParamChangeFileContent(HRUinfo,paraSelection,SWATParam,
                                        TxtInOutFolder)

# Set first run is true so R-SWAT will delete previous simulation
firstRun <- TRUE
copyUnchangeFiles <- TRUE

# Get content of the file.cio file (about simulation time)
fileCioInfo <- getFileCioInfo(TxtInOutFolder)

# Now start to run SWAT
runSWATpar(workingFolder,TxtInOutFolder,outputExtraction,ncores,SWATexeFile,
           parameterValue,paraSelection,caliParam,copyUnchangeFiles,fileCioInfo,
           dateRangeCali,firstRun)

# Number of output variables (from the output extraction data frame)
OutputVar <- getNumberOutputVar(outputExtraction)
nOutputVar <- OutputVar$nOutputVar

# Check if users use their own output extraction function
userReadSwatOutput <- OutputVar$userReadSwatOutput

# Get observed data (first need to sort observed data file)
observedDataFile <- sortObservedDataFile(observedDataFile)

# Now read observed data (as list) from observed files
observedData <- list()
for (i in 1:length(observedDataFile)){
  # Read observed files and save to a dummy variable
  temp <- read.table(observedDataFile[i], skip = 1, sep = "")
  # Get bbserved data from dummy variable
  observedData [[i]] <- data.frame(Date = as.POSIXct(paste(temp[,1], 
                                                           temp[,2], 
                                                           sep = " "), 
                                                     format = "%Y-%m-%d %H:%M", 
                                                     tz = ""),
                                   Value = temp[,3],
                                   Flag = temp[,4])
}

# Calculate objective function (this function goes through the output/TxtInOut_x and reads the simulated data)
obj <- calObjFunction(parameterValue,
                      ncores,
                      nOutputVar,
                      userReadSwatOutput,
                      observedData,
                      workingFolder,
                      objFunction)

# Model performance for calibration (obj$objValueCali) and validation (obj$objValueValid)
objValueCali <- obj$objValueCali
objValueValid <- obj$objValueValid

# Parameter sensitivity using multi-variable regression
# Prepare table with parameter values and objective function
parameterObj <- as.data.frame(cbind(objValueCali, parameterValue[,-c(1)]))
colnames(parameterObj) <- c("objFunction", paraSelection[,1])

# Parameter sensitivity using multivariate regression analysis
ParamSensitivity <- summary(lm(formula = objFunction ~ ., parameterObj))[4]$coefficients[,3:4] 

# The objective function should be maximum "Maximize" or minimum "Minimize"
# e.g., if objective function is NSE, KGE , R2 then it should be "Maximize"
#       if objective function is "RMSE", "aBIAS" then it should be "Minimize"
minOrmax <- "Maxmimize"

# Which output variable number you want to calculate behavioral simulations
varNumber <- 1

# Calculate behavior
shinyCatch(
    behavioral <- behaSimulation(obj$objValueCali,
                             obj$simData,
                             parameterValue,
                             behThreshold,
                             varNumber,
                             objFunction,
                             observedData,
                             minOrmax,
                             samplingApproach)
)

# Save all data + function 
save.image(paste0(workingFolder, "/RSWATproject.RData", sep = ""))

# Load all data function for later work
# load("./workingFolder/RSWATproject.RData")

# Print the best objective function values for model calibration and validation
print(paste0(objFunction," calibration = ",max(objValueCali)))
print(paste0(objFunction," validation = ",max(objValueValid)))
