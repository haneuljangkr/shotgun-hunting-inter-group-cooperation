## Load packages and data
library(brms)
library(usdm)
library(dplyr)
library(tidyverse)
library(rethinking)
library(dagitty)

# set working directory
setwd("/Users/haneuljangkr/shotgun-hunting-inter-group-cooperation")

d <- read.csv("hunting_cooperation_data.csv")
ky <- read.csv("Yambe_kin.csv", row.names = 1)
kb <- read.csv("BaYaka_kin.csv", row.names = 1)

##################################################################
######################## Data preparation ########################
##################################################################

# YKINID and BKINID is an alias of YAID and BYID, but name differently for brms to estimate kinship effects
d$YKINID <- d$YAID 
d$BKINID <- d$BYID 

length(unique(d$YAID))
length(unique(d$BYID))
table(d$hunt_bin)

# standardize continuous variables
d$BY_age_z <- as.numeric(scale(d$BY_age))
d$YA_age_z <- as.numeric(scale(d$YA_age))
d$BY_noChild_z <- as.numeric(scale(d$BY_noChild))
d$YA_noChild_z <- as.numeric(scale(d$YA_noChild))
d$NoGuns_z <- as.numeric(scale(d$NoGuns)) 
d$BY_skill_z <- as.numeric(scale(d$SmithsS_BYSkill)) 

# make factors
d$Mother_natal <- as.factor(d$Mother_natal)
d$Wife_natal <- as.factor(d$Wife_natal)
d$Council <- as.factor(d$Council)
d$MetalRoof <- as.factor(d$MetalRoof)

##################################################################
############## Express causal assumptions as a DAG ###############
##################################################################
g <- dagitty(
  "dag{
  BaYaka_Age -> BaYaka_Skill
  BaYaka_Age -> BaYaka_No_Children
  BaYaka_Age -> Cooperate
  
  BaYaka_Skill -> Cooperate
  BaYaka_No_Children -> Cooperate
  
  Yambe_Age -> Wife_Natal
  Yambe_Age -> Mother_Natal
  Yambe_Age -> Cooperate
  Yambe_Age -> Yambe_No_Children
  Yambe_Age -> U
  Yambe_Age -> Yambe_No_Guns
  Yambe_Age -> Council
  Yambe_Age -> Metal_Roof
  
  Wife_Natal -> Cooperate
  Mother_Natal -> Cooperate
  
  U -> Yambe_No_Children
  U -> Yambe_No_Guns
  U -> Council
  U -> Metal_Roof
  
  Yambe_No_Children -> Cooperate
  Yambe_No_Guns -> Cooperate
  Council -> Cooperate
  Metal_Roof -> Cooperate
  }"
)

# Visualize DAG to make sure we set it up as intended
plot(g)

##################################################################
######################### Baseline model #########################
##################################################################

# Set priors
prior_0 <- c(prior(normal(0, 1), class = Intercept),
             prior(exponential(1), class = sd))

# The baseline modelaccounts for multilevel structure of data but includes no covariates
m0 <- brm(hunt_bin ~ (1|YAID) + (1|BYID) + (1|gr(YKINID, cov = ky)) + (1|gr(BKINID, cov = kb)), 
          family = bernoulli(),
          prior = prior_0,
          data = d,
          data2 = list(ky=ky, kb=kb),
          iter = 4000,
          warmup = 2000,
          chains = 4,
          cores = 4,
          control = list(adapt_delta = 0.9),
          init = "0",
          seed = 14)

# Variance decomposition: R^2 for Yambe ID, Yambe Kinship, and BaYaka ID
post <- as_draws_df(m0)

total_var <- post$sd_YAID__Intercept^2 + post$sd_BYID__Intercept^2 + post$sd_YKINID__Intercept^2 + post$sd_BKINID__Intercept^2 + pi^2 / 3

r2_YAID <- (post$sd_YAID__Intercept^2) / total_var
r2_BYID <- (post$sd_BYID__Intercept^2) / total_var
r2_YKINID <- (post$sd_YKINID__Intercept^2) / total_var
r2_BKINID <- (post$sd_BKINID__Intercept^2) / total_var

# Note there's high uncertainty in how much variance in Yambe random intercepts is due to shared kinship vs unique individual variance.
hist(r2_YAID)
hist(r2_YKINID)

mean(r2_YAID)
mean(r2_YKINID)

# However, we can see that their shared variance is substantial
hist(r2_YAID + r2_YKINID)

# Substanial variance across BaYaka IDs as well
hist(r2_BYID)
hist(r2_BKINID)

# Difference in variance: BaYaka ID vs Yambe ID
hist((r2_BYID + r2_BKINID) - (r2_YAID + r2_YKINID))
mean((r2_BYID + r2_BKINID) - (r2_YAID + r2_YKINID))
mean( (r2_BYID + r2_BKINID) - (r2_YAID + r2_YKINID) > 0 )


##################################################################
################## Models for BaYaka attributes ##################
##################################################################

# Use the DAG to get our adjustment sets for covariates. 
# For BaYaka attributes, we only need to adjust for BaYaka age
adjustmentSets(g, exposure = "BaYaka_No_Children", outcome = "Cooperate")
adjustmentSets(g, exposure = "BaYaka_Skill", outcome = "Cooperate")

# Set priors
prior <- c(prior(normal(0, 1), class = Intercept),
           prior(exponential(1), class = sd),
           prior(normal(0, 1), class = b)
)

# Model with the number of children of BaYaka hunters
m_bayaka_children <- brm(hunt_bin ~ (1|YAID) + (1|BYID) + (1|gr(YKINID, cov = ky)) + (1|gr(BKINID, cov = kb)) + 
                         BY_age_z + BY_noChild_z,
                         family = bernoulli(),
                         prior = prior,
                         data = d,
                         data2 = list(ky=ky, kb=kb),
                         iter = 4000,
                         warmup = 2000,
                         chains = 4,
                         cores = 4,
                         control = list(adapt_delta = 0.97),
                         init = "0",
                         seed = 14)

# Model with hunting skill score of BaYaka hunters
m_bayaka_skill <- brm(hunt_bin ~ (1|YAID) + (1|BYID) + (1|gr(YKINID, cov = ky)) + (1|gr(BKINID, cov = kb)) + 
                      BY_age_z + BY_skill_z,
                      family = bernoulli(),
                      prior = prior,
                      data = d,
                      data2 = list(ky=ky, kb=kb),
                      iter = 4000,
                      warmup = 2000,
                      chains = 4,
                      cores = 4,
                      control = list(adapt_delta = 0.97),
                      init = "0",
                      seed = 14)

# Visualize effects
conditional_effects(m_bayaka_children, "BY_noChild_z")
conditional_effects(m_bayaka_skill, "BY_skill_z")

# Model comparison
round(brms::model_weights(m0, m_bayaka_children), 3)
round(brms::model_weights(m0, m_bayaka_skill), 3)


##################################################################
################## Models for Yambe attributes ###################
##################################################################

adjustmentSets(g, exposure = "Yambe_No_Children", outcome = "Cooperate")
adjustmentSets(g, exposure = "Yambe_No_Guns", outcome = "Cooperate")
adjustmentSets(g, exposure = "Council", outcome = "Cooperate")
adjustmentSets(g, exposure = "Metal_Roof", outcome = "Cooperate")
adjustmentSets(g, exposure = "Wife_Natal", outcome = "Cooperate")
adjustmentSets(g, exposure = "Mother_Natal", outcome = "Cooperate")

# Set priors
prior <- c(prior(normal(0, 1), class = Intercept),
           prior(exponential(1), class = sd),
           prior(normal(0, 1), class = b)
)

# Model with the number of children, the number of guns, metal roof presence, council attendance
m_yambe_1 <- brm(hunt_bin ~ (1|YAID) + (1|BYID) + (1|gr(YKINID, cov = ky)) + (1|gr(BKINID, cov = kb)) + 
                 YA_age_z + YA_noChild_z + mo(NoGuns) + MetalRoof + Council,
                 family = bernoulli(),
                 prior = prior,
                 data = d,
                 data2 = list(ky=ky, kb=kb),
                 iter = 4000,
                 warmup = 2000,
                 chains = 4,
                 cores = 4,
                 control = list(adapt_delta = 0.97),
                 init = "0",
                 seed = 14)

# Model with Yambe's wife natal village
m_wife_natal <- brm(hunt_bin ~ (1|YAID) + (1|BYID) + (1|gr(YKINID, cov = ky)) + (1|gr(BKINID, cov = kb)) + 
                    YA_age_z + Wife_natal,
                    family = bernoulli(),
                    prior = prior,
                    data = d,
                    data2 = list(ky=ky, kb=kb),
                    iter = 4000,
                    warmup = 2000,
                    chains = 4,
                    cores = 4,
                    control = list(adapt_delta = 0.97),
                    init = "0",
                    seed = 14)

# Model with Yambe's mother natal village
m_mother_natal <- brm(hunt_bin ~ (1|YAID) + (1|BYID) + (1|gr(YKINID, cov = ky)) + (1|gr(BKINID, cov = kb)) + 
                      YA_age_z + Mother_natal,
                      family = bernoulli(),
                      prior = prior,
                      data = d,
                      data2 = list(ky=ky, kb=kb),
                      iter = 4000,
                      warmup = 2000,
                      chains = 4,
                      cores = 4,
                      control = list(adapt_delta = 0.97),
                      init = "0",
                      seed = 14)

# Visualize effects
conditional_effects(m_yambe_1, effects = "YA_noChild_z")
conditional_effects(m_wife_natal, effects = "Wife_natal")
conditional_effects(m_mother_natal, effects = "Mother_natal")

# Model comparison
round(brms::model_weights(m0, m_yambe_1, weights = "stacking"), 3)
round(brms::model_weights(m0, m_wife_natal), 3)
round(brms::model_weights(m0, m_mother_natal), 3)

