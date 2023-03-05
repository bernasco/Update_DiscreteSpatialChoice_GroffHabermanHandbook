##################################################
# Crime Location Choice Models (update March 2023)
# Wim Bernasco, 2023
##################################################
library(tidyverse) # new
library(broom)

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
B    <- read_csv(here("TheHagueBurglars.csv"))
N    <- read_csv(here("TheHagueNeighborhoods.csv"))



# Create data structure for discrete location choice
BXN <- cross_join(B,N) |>
  mutate(
    # Create dependent variable: was neighborhood target or not?
    CHOSEN = as.numeric(NHOODBUR==NHOODID),
    # new variable: distance home - potential target
    DIST   = sqrt((XRESID-X)^2 + (YRESID-Y)^2) / 1000,
    # Correction: Distance within the same neighborhood
    DIST   = if_else(DIST==0, sqrt(SURFACE) / 2, DIST),
    # Transform distance measure to proximity measure for convenience
    PROXIMITY = -DIST
    )



# Bernasco & Nieuwbeerta (2005) Table 2 (p. 308)
Model1_clogit <- clogit(CHOSEN~PROPVAL+SINGFAM+RESMOBIL+ETNHETERO+PROXIMITY+
                 PROXCITY+RESUNITS+strata(CASE),data=BXN, method="exact")
# Traditional summary (not a tibble)
summary(Model1_clogit)
# Model-level results (single-row tibble)
glance(Model1_clogit)
# Coeffcient table with estimate, standard error, statistic, and p-value
tidy(Model1_clogit)
# Same with estimate replace with exp(estimate) 
tidy(Model1_clogit, exponentiate=TRUE)

# robust SE estimates
Model1a_clogit<-clogit(CHOSEN~PROPVAL+SINGFAM+RESMOBIL+ETNHETERO+PROXIMITY+
                 PROXCITY+RESUNITS+strata(CASE),data=BXN,
                 method="efron", cluster=PERSONID)
summary(Model1a_clogit)
# You cannot (yet) use tidy to obtain robust estimates in a tibble
tidy(Model1a_clogit, exponentiate=TRUE)


BXN_mlogit <- 
  BXN |>
  dfidx(idx = c("CASE", "NHOODID"))


# # Create data frame for use with mlogit
# BXN.mldata <- mlogit.data(data=BXN, 
#                           shape="long",
#                           choice = "CHOSEN", 
#                           alt.var="NHOODID", 
#                           chid.var="CASE",
#                           id.var="PERSONID")
# estimate conditional logit model
Model1_mlogit <-
  mlogit(formula = CHOSEN~PROPVAL+SINGFAM+RESMOBIL+
           ETNHETERO+PROXIMITY+
           PROXCITY+RESUNITS |0,
         data=BXN_mlogit)
summary(Model1_mlogit)
tidy(Model1_mlogit)


# Bernasco & Nieuwbeerta (2005), Table 3 (p. 309)
BXN <- 
  BXN |>
  mutate(EH_NATIVE  = ETNHETERO * B_NATIVE,
         EH_FOREIGN = ETNHETERO * B_FOREIGN,
         PR_MINOR   = PROXIMITY * B_MINOR,
         PR_ADULT   = PROXIMITY * B_ADULT
         )


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
         data=BXN_mlogit, rpar = c(PROXIMITY= 'n'))
summary(Mixedlogit1_mlogit)

#--------------------------------------------------------------------
# Random sampling from alternatives
BXN_SFA <- 
  bind_rows(
    # first select all non-chosen rows
    BXN |> filter(CHOSEN==0) |>
      group_by(CASE) |>
      # within each crime, sample 29 of the non-chosen rows
      sample_n(size = 29, replace=FALSE),
    # append the result to the chosen rows (which are this sampled with P=1)
    BXN |> filter(CHOSEN==1) 
  ) |>
  # sort for convenience
  arrange(CASE, NHOODID)
# show result 
BXN_SFA 
Model_SFA_clogit <- clogit(CHOSEN~PROPVAL+SINGFAM+RESMOBIL+ETNHETERO+PROXIMITY+
                            PROXCITY+RESUNITS+strata(CASE),data=BXN_SFA, method="exact")
summary(Model_SFA_clogit)


# # Importance sampling from alternatives
# # subset of 'important' alternatives + chosen

BXN_ISFA <- 
  bind_rows(
    # first select all non-chosen rows
    BXN %>% filter(CHOSEN==0 & PROXCITY <= -1.5) |>
      group_by(CASE) |> 
      # count 'leftover' alternatives per crime
      mutate(cases=n(),
             # log inverse sampling probability 
             loginvp = log( cases / 2)) |>
      # remove redundant variable
      select(-cases) |>
      # per crime, sample 29 of the leftover alternatives
      sample_n(size = 29, replace=FALSE),
    # append the result to the chosen rows (which are thus sampled with P=1)
    BXN |> filter(CHOSEN==1 | PROXCITY > -1.5) |>
      # p=1, sp 1/p=1, sp log(1/p)=0 for this subset
      mutate(loginvp = log(1))  
  ) |>
  # sort for convenience
  arrange(CASE, NHOODID)
# show result 
BXN_ISFA 

#  The offset corrects the estimates)
Model_ISFA_clogit <-
  clogit(CHOSEN ~ PROPVAL + SINGFAM + RESMOBIL + ETNHETERO + PROXIMITY +
         PROXCITY + RESUNITS + strata(CASE) + offset(loginvp), 
         data=BXN_ISFA)
summary(Model_ISFA_clogit)



# Extra: Conditional logit and Poisson regression equivalance ------------------
# Estimate a conditional logit model without PROXIMITY variable. This implies 
#   all active variables only vary across neighborhoods, not across offences.
alt_model_clogit <- 
  clogit(CHOSEN ~ PROPVAL + SINGFAM + RESMOBIL + ETNHETERO +
           PROXCITY + RESUNITS + strata(CASE),
         data=BXN, method="exact")

# Aggregate crimes across neighborhoods
aggregated_BXN <-
  BXN |>
  # The next two lines calculate # of crimes per neighborhood,
  #   but they do not aggregate
  group_by(NHOODID) |>
  mutate(BURGLARY_COUNT = sum(CHOSEN)) |>
  # The next line aggregates because all relevant variables are
  #   constants within neighborhood
  filter(row_number() == 1) 

# Estimate a Poisson model on the aggregated data
alt_model_poisson <-
  glm(BURGLARY_COUNT ~  PROPVAL + SINGFAM + RESMOBIL + ETNHETERO +
        PROXCITY + RESUNITS, 
      family = "poisson", data = aggregated_BXN) 

# Compare results
tidy(alt_model_clogit)
tidy(alt_model_poisson)


