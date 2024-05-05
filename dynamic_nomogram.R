rm(list=ls())
library(this.path);setwd(this.path::this.dir());cat('当前脚本执行路径为',getwd());
library(survival)
library(survminer)
library(lattice)
library(Formula)
library(ggplot2)
library(Hmisc)
library(rms)
library(regplot)

library(pec)
library(dplyr)
library(DynNom)
library(shiny)
library(plotly)
library(stargazer)
library(compare)
library(prediction)

# █████████████████████████████████████████████████████████████████████████████
# load data
# █████████████████████████████████████████████████████████████████████████████
df <- read.csv('dynamic_nomogram_data.csv', header=T,row.names = 1)

# Delete the maximum and minimum values of riskscore
df <- df[!(rownames(df) == "X0928906"), ]  # X0928906 minimum, RiskScore_h0_log = -3.96
df <- df[!(rownames(df) == "X1127693"), ]  # X1127693 minimum, RiskScore_h0_log = -2.1
df <- df[!(rownames(df) == "X0908051"), ] # X0908051 maximum 38.7


# 筛出训练集、测试集数据
TCGA_LGG <- read.csv('groups/TRAIN_tcga_nonGBM_54.txt', header=F)
TCGA_GBM <- read.csv('groups/TRAIN_tcga_gbm_83.txt', header=F)
JSPH_LGG <- read.csv('groups/TRAIN_jsph_nonGBM_39.txt', header=F) 
JSPH_GBM <- read.csv('groups/TRAIN_jsph_gbm_68.txt', header=F)

test_JSPH_LGG <- read.csv('groups/TEST_jsph_nonGBM_40.txt', header=F)
test_JSPH_GBM <- read.csv('groups/TEST_jsph_gbm_68.txt', header=F)

LGG <- rbind(TCGA_LGG,JSPH_LGG,test_JSPH_LGG)
GBM <- rbind(TCGA_GBM,JSPH_GBM,test_JSPH_GBM)

LGG_train <- rbind(TCGA_LGG,JSPH_LGG)
GBM_train <- rbind(TCGA_GBM,JSPH_GBM)
train <- rbind(LGG_train,GBM_train)

LGG_test <- rbind(test_JSPH_LGG)
GBM_test <- rbind(test_JSPH_GBM)
test <- rbind(LGG_test,GBM_test)

all <-  rbind(train,test)

# ←~~~~~~~~~~~~~~~~~~~~~ CHOSE COHORT ←~~~~~~~~~~~~~~~~~~~~~
df = df[train$V1,]


# clean 
rm(
  LGG,LGG_train,LGG_test,
  GBM,GBM_train,GBM_test,
  TCGA_LGG,TCGA_GBM,
  JSPH_LGG,JSPH_GBM,
  test_JSPH_LGG,test_JSPH_GBM,
  train, test
   )

# █████████████████████████████████████████████████████████████████████████████
# data organization
# █████████████████████████████████████████████████████████████████████████████

# Divide the time variable by 12 to change its units from months to adults
df$time <- df$time / 12

# Convert the age column to numerical type
df$age <- as.numeric(df$age)
df$RiskScore_h0 <- as.numeric(df$RiskScore_h0)
df$RiskScore_h0_log <- as.numeric(df$RiskScore_h0_log)

names(df)

features_needed <- c(
  'time',
  'status',

  'RiskScore_h0',
  'RiskScore_h0_log',
  
  'age'
)

df <- df[ ,colnames(df) %in% features_needed]

# Remove all empty rows
df <- na.omit(df)

summary(df)
str(df)

# █████████████████████████████████████████████████████████████████████████████
# Draw a dynamic column chart
# █████████████████████████████████████████████████████████████████████████████

res.cox = coxph(Surv(time, status)~ 
                  # RiskScore_h0+
                  RiskScore_h0_log+
                  # risktype+
                  age,
                
                data =df)

DynNom(res.cox, 
       clevel = 0.55,
       DNtitle = 'Nomogram',
       DNxlab = 'probability',
       DNylab = NULL,
       DNlimits = NULL)




