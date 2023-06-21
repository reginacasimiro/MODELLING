#########
##PACKAGE INSTALLATION
#########
# download the package from CRAN and install it
install.packages(c("epanetReader","epanet2toolkit","forecast","GA"))
library(epanetReader)
library(epanet2toolkit)
library(forecast)
library(GA)

# path to Net1.inp example file included with this package
inp <- file.path( find.package("epanet2toolkit"), "extdata","Net1.inp")
print(inp)

#set wd in same place of Net1.inp
setwd("C:/Users/ReginaCasimiro/Documents/R/win-library/3.4/epanet2toolkit/extdata")

#########
##RUNNING A FULL SIMULATION
#########
ENepanet( inp, "Net1.rpt")# results stored in the report Net1.rpt

# Accessing and changing network properties
ENopen("Net1.inp", "Net1.rpt")
Net1 # to see  the .inp file
ENgetflowunits()
ENgetlinkvalue(2, "EN_LENGTH") # to extract the value from link
ENsetlinkvalue(2, "EN_LENGTH", 6789)
ENgetlinkvalue(2, "EN_LENGTH")
ENclose()


#########
##OBJECTIVE FUNCTION FOR PIPE ROUGHNESS CALIBRATION
#########
# calls the helper function setAllPipeRoughness to change the roughnessvalue of all pipes in the network to the provided value 
calibObj<-function(roughness, obsindex, obs){
  setAllPipeRoughness(roughness)
  ENsolveH() # carries out a full hydraulic simulation
  sse(obsindex, obs)
}

# loops over all links in the network, only updating the roughness if the link if of type pipe as opposed to pump or valve
setAllPipeRoughness <- function( roughness ){
  for(i in 1:ENgetcount("EN_LINKCOUNT")){
    if( ENgetlinktype(i) < 2)
      ENsetlinkvalue(i, "EN_ROUGHNESS", roughness )
  }
}
# computes the sum of squared error between measured and modeled pressure
sse <- function( nodeIndex, measuredPressure){
  N <- length(nodeIndex)
  error <- rep(NA, N)
  for( i in 1:N){
    modeledPressure <- ENgetnodevalue( nodeIndex[i], "EN_PRESSURE")
    error[i] <- modeledPressure - measuredPressure[i]
  }
  sse <- sum( error * error )
}

#########
## PIPE ROUGHNESS CALIBRATION WITH UNIVERATE OPTIMIZATION
## assuming the same roughness value applies to all pipes in the network
#########

# to initializes Epanet and processes network information
ENopen("Net1.inp", "Net1.rpt")
# update the network to the conditions of the example
ENsettimeparam("EN_DURATION", 0) # contitions were typical for the time 00:00
ENsetnodevalue(7, "EN_BASEDEMAND", 2000) # demand of 2000 gallons per minute was induced at node 23 (index 7)
# under that condition, pressure measurements of  112.11, ...psi were collected at nodes with indices 4, 6, and 8
optimize( calibObj,interval = c(50, 150),obsindex = c(4,6,8),obs = c(112.11, 110.87, 110.32) ) # to find a value of pipe ruohnessthat minimizes the sum of squarrel errors between the measured and modeled pressures 
# the roughness value,from the Hazen-Williams formula, is constrained to fall in the interval between 50 and 150
# optimize returns a list with components minimum (or maximum) and objective which give the location of the minimum (or maximum) and the value of the function at that point
# running this example yields a C-value of 131.2724 ($minimum) and a objective value ($objective) of 0.8828811.

ENclose()

#########
## MINIMIZING ENERGY COST
#########
setwd("C:/Users/ReginaCasimiro/Desktop")

updatePumpSchedule <- function(x){
  pump_patterns <- c( "pmp1Apat", "pmp2Apat", "pmp3Apat", "pmp4Bpat",
                      "pmp5Cpat", "pmp6Dpat", "pmp7Fpat")
  labeled_inputs <- split( x, rep(pump_patterns, each = 48 ) )
  # Update the pattern/schedule for each pump
  for( i in 1:length(pump_patterns)){
    pattern_index <- ENgetpatternindex( pump_patterns[i] )
    new_schedule <- labeled_inputs[[i]]
    ENsetpattern(pattern_index, new_schedule)
  }
}
objFun <- function(x){
  ENopen("richmond-modified.inp","temp.rpt","temp.bin")
  updatePumpSchedule(x)
  # penalty for a hydraulic warning or error
  penalty <- tryCatch({
    ENsolveH()
    0 # no penalty if all is well
  }, warning = function(w){
    1000 # penalty for warning
  }, error = function(e){
    5000 # penalty for error
  })
  # solve & save & close
  ENsolveQ(); ENreport(); ENclose()
  # read the result out of the rpt file
  energyUsage <-suppressWarnings( read.rpt("temp.rpt")$energyUsage )
  cost <- sum(energyUsage$dailyCost) + penalty
}

# use existing schedule as starting point
rmod <- read.inp("richmond-modified.inp")
sugg <- c(rmod$Patterns$pmp1Apat, rmod$Patterns$pmp2Apat, rmod$Patterns$pmp3Apat,
          rmod$Patterns$pmp4Bpat, rmod$Patterns$pmp5Cpat, rmod$Patterns$pmp6Dpat,
          rmod$Patterns$pmp7Fpat)
# apply genetic algorithm
myga <- ga( type = "binary", fitness = function(x) -objFun(x),
            nBits = 7*48, maxiter=10, monitor = FALSE, suggestions = sugg)
summary(myga)

# save optimal schedule in new file
ENopen("richmond-modified.inp","temp.rpt","temp.bin")
updatePumpSchedule( myga@solution[1,] )
ENsaveinpfile("richmond-modified-solved.inp")
## NULL
ENclose()



########
##???????????????????????
#######
setwd("C:/Users/ReginaCasimiro/Desktop")
pred_mean <- day_ahead$mean
# standard normal quantile for upper 80% interval
Z_10 <- qnorm( 0.90, mean=0, sd=1)
# calc the standardized error of the prediction
pred_std_err <- (day_ahead$upper[,"80%"] - pred_mean ) / Z_10
# now we can use this to sample
sampled_demands <- function( pred, se ){
  N <- length(pred)
  Z <- rnorm( N, mean=0, sd = 1)
  x <- pred + Z * se
}
numRealizations <- 200
daily_cost <- rep(NA, numRealizations)
for(i in 1:numRealizations){
  dmd <- sampled_demands( pred_mean, pred_std_err)
  update_pattern("richmond-modified-solved.inp", "forecast", dmd)
  ENepanet("richmond-modified-solved.inp", "mc.rpt")
  rm_rpt <- read.rpt("mc.rpt")
  daily_cost[i] <- sum(rm_rpt$energyUsage$dailyCost)
}
hist(daily_cost)
#############################################################################################
##Stochastic Simulation ???
#########
install.packages("astsa")
library(astsa)
data<-read.csv("C:/Users/ReginaCasimiro/Documentos/data.csv",header=TRUE) # the total water demand in the network measured at hourly intervals for a period of one week is available in a comma-separated values file data.csv
# cannot open file 'data.csv': No such file or directory
data <- read.csv("data.csv")
forecast <- sarima.for(data$Measurement, n.ahead = 24, p = 0, d = 1, q = 4,
                         P = 0, D = 1, Q = 1, S = 24)
newpattern <- as.numeric(forecast$pred / mean(forecast$pred))
ENopen("Net1.inp", "Net1.rpt")
ENsettimeparam("EN_PATTERNSTEP", 3600)
ENsetpattern(1, newpattern)
ENsaveinpfile("Net1-forecast.inp")
ENclose()
#########
##Running the package test suite ???????????
#########
# is designed to run automatically to detect changes that alter package behaviour and to document bugs and fixes
install.packages("testthat")
library(testthat)
test_package("epanet2toolkit",reporter=default_reporter())




