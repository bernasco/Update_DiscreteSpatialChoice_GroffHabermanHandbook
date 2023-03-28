# Update_DiscreteSpatialChoice_GroffHabermanHandbook

Materials for supplement to chapter 9.3 on discrete spatial choice in "Understanding Crime and Place: A Methods Handbook" by Groff and Haberman

This repository contains an update that addresses three issues in particular
 1. The new script uses the tidyverse ecosystem (mostly dplyr, readr, tidyr, broom)
 2. The new script uses the 'dfidx::dfidx' function rather than the deprecated 'mlogit::mlogit.data' function 
 3. The new script demonstrates the nested logit model 
 
 In addition, it provides a computational example showing the conditional logit model is equivalent to a Poisson model with choices aggregated across alternatives. This equivalence only holds if all alternatives have the same values for all decision makers. 
