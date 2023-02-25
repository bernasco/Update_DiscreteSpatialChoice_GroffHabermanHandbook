##################################################
# Crime Location Choice Models
# Wim Bernasco, 2020
##################################################
library(survival)  # clogit() function
library(mlogit)    # mlogit() function
library(car)       # linearHypothesis() function for coefficient tests
library(here)      # file localization 



####################################################
# The Public Domain (CrimeStat) version of the 
#   The Hague burglary data contain average
#   neighborhood attributes across 1996-2002
####################################################

# Read datafiles
B    <- read.csv(here("TheHagueBurglars.csv"))
N    <- read.csv(here("TheHagueNeighborhoods.csv"))


# Create data structure for discrete location choice
BXN <- merge(B,N)
# Create dependent variable: was neighborhood target or not
BXN$CHOSEN <- as.numeric(BXN$NHOODBUR==BXN$NHOODID)
# Create independent variable: distance home - potential target
BXN$DIST <- sqrt((BXN$XRESID-BXN$X)^2 + (BXN$YRESID-BXN$Y)^2) / 1000
# Correction: Distance within the same neighborhood
BXN$DIST[BXN$DIST==0] <- sqrt(BXN$SURFACE[BXN$DIST==0]) / 2
# Transform distance measure to proximity measure for convenience
BXN$PROXIMITY <- -BXN$DIST

# Bernasco & Nieuwbeerta (2005) Table 2 (p. 308)
Model1_clogit <- clogit(CHOSEN~PROPVAL+SINGFAM+RESMOBIL+ETNHETERO+PROXIMITY+
                 PROXCITY+RESUNITS+strata(CASE),data=BXN, method="exact")
summary(Model1_clogit)

# robust SE estimates
Model1a_clogit<-clogit(CHOSEN~PROPVAL+SINGFAM+RESMOBIL+ETNHETERO+PROXIMITY+
                 PROXCITY+RESUNITS+strata(CASE),data=BXN,
                 method="efron", cluster=PERSONID)
summary(Model1a_clogit)



# sort by CASEID, NHOODID
BXN <- BXN[order(BXN$CASE, BXN$NHOODID),]

# Create data frame for use with mlogit
BXN.mldata <- mlogit.data(data=BXN, 
                          shape="long",
                          choice = "CHOSEN", 
                          alt.var="NHOODID", 
                          chid.var="CASE",
                          id.var="PERSONID")
# estimate conditional logit model
Model1_mlogit <-
  mlogit(formula = CHOSEN~PROPVAL+SINGFAM+RESMOBIL+
           ETNHETERO+PROXIMITY+
           PROXCITY+RESUNITS |0,
         data=BXN.mldata)
summary(Model1_mlogit)


# Bernasco & Nieuwbeerta (2005), Table 3 (p. 309)
BXN$EH_NATIVE  <- BXN$ETNHETERO * BXN$B_NATIVE
BXN$EH_FOREIGN <- BXN$ETNHETERO * BXN$B_FOREIGN
BXN$PR_MINOR   <- BXN$PROXIMITY * BXN$B_MINOR
BXN$PR_ADULT   <- BXN$PROXIMITY * BXN$B_ADULT

Model2_clogit<-clogit(CHOSEN~PROPVAL+SINGFAM+RESMOBIL+PROXCITY+RESUNITS+
                        PR_ADULT+PR_MINOR+EH_NATIVE+EH_FOREIGN+strata(CASE),
                      data=BXN)
summary(Model2_clogit)

# test coefficient equality
linearHypothesis(Model2_clogit, "EH_NATIVE = EH_FOREIGN")
linearHypothesis(Model2_clogit, "PR_MINOR = PR_ADULT")

# mixed logit models
Mixedlogit1_mlogit <-
  mlogit(formula = CHOSEN~PROPVAL+SINGFAM+RESMOBIL+
           ETNHETERO+PROXIMITY+
           PROXCITY+RESUNITS |0,
         data=BXN.mldata, rpar = c(PROXIMITY= 'n'))
summary(Mixedlogit1_mlogit)

# Random sampling from alternatives
# subset of chosen alternatives
BXNChosen    <- BXN[BXN$CHOSEN==1,]
# subset of non-chosen alternatives
BXNNotChosen <- BXN[BXN$CHOSEN==0,]

# Create a list of subsets (by crime_id) of the non-chosen alternatives
#   and apply the sampling function to each subset
df2 <- lapply(split(BXNNotChosen, BXNNotChosen$CASE),
              function(subdf) subdf[sample(1:nrow(subdf), 29, replace=FALSE),]
)
# convert the result (which is a list) back to a dataframe by rbind
df3 <- rbind(do.call('rbind', df2), BXNChosen)
# sort for convenience
BXN_SA <- df3[order(df3$CASE, df3$NHOODID),]
# estimate model on subset of alternatives
Model_SA_clogit <- clogit(CHOSEN~PROPVAL+SINGFAM+RESMOBIL+ETNHETERO+PROXIMITY+
                            PROXCITY+RESUNITS+strata(CASE),data=BXN_SA, method="exact")
summary(Model_SA_clogit)

# Importance sampling from alternatives
# subset of 'important' alternatives + chosen
BXNSelected    <- BXN[BXN$CHOSEN==1 | BXN$PROXCITY > -1.5,]
# subset of non-chosen alternatives
BXNNotSelected <- BXN[BXN$CHOSEN==0 & BXN$PROXCITY <= -1.5,]

# Create a list of subsets (by crime_id) of the non-chosen alternatives
#   and apply the sampling function to each subset
df2 <- lapply(split(BXNNotSelected, BXNNotSelected$CASE),
              function(subdf) subdf[sample(1:nrow(subdf), 30, replace=FALSE),] 
)
# convert list to dataframe by rbind
BXNSampled <- do.call('rbind', df2)
# assign correction factors
BXNSelected$log_INCLUSION <- log(1)
BXNSampled$log_INCLUSION <- 
  log(1/( nrow(BXNSampled) / nrow(BXNNotSelected)))
# merge both parts of dataset 
df3 <- rbind(BXNSampled, BXNSelected)
# sort for convenience
BXN_SAI <- df3[order(df3$CASE, df3$NHOODID),]
# estimate model on subset of alternatives
Model_SAI_clogit <-clogit(CHOSEN~PROPVAL+SINGFAM+RESMOBIL+ETNHETERO+PROXIMITY+
                       PROXCITY+RESUNITS+offset(log_INCLUSION)+strata(CASE),data=BXN_SAI,
                     method="efron", robust=TRUE)
summary(Model_SAI_clogit)
