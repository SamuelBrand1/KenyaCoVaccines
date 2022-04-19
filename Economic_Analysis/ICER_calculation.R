# this script uses the output of the CEA function to calculate the ICERs
# and plot cost effectiveness plane
#load and install required libraries
library(BCEA)
library(plotly)
library(ggplot2)
library(pastecs)
library(dplyr)
library(tidyverse)
#set the working directory
setwd("C:/Users/corlendo/KenyaCoVaccinesPrivate")

#load in all the scenario files
load(file = "model_outputs/societal_persp/CEA_MODEL_OUTPUT_10_scenario_2baseline_scenario_1.rda")
case.21=as.data.frame(CEA_output)
load(file = "model_outputs/societal_persp/CEA_MODEL_OUTPUT_10_scenario_3baseline_scenario_1.rda")
case.31=as.data.frame(CEA_output)
load(file = "model_outputs/societal_persp/CEA_MODEL_OUTPUT_10_scenario_4baseline_scenario_1.rda")
case.41=as.data.frame(CEA_output)
load(file = "model_outputs/societal_persp/CEA_MODEL_OUTPUT_10_scenario_5baseline_scenario_1.rda")
case.51=as.data.frame(CEA_output)
load(file = "model_outputs/societal_persp/CEA_MODEL_OUTPUT_10_scenario_6baseline_scenario_1.rda")
case.61=as.data.frame(CEA_output)
load(file = "model_outputs/societal_persp/CEA_MODEL_OUTPUT_10_scenario_7baseline_scenario_1.rda")
case.71=as.data.frame(CEA_output)

#labels vectors
ints1 <- c("30%","0%")
ints2 <- c("50%","30%")
ints3 <- c("70%","50%")
ints4 <- c("Rapid_30%","Rapid_0%")
ints5 <- c("Rapid_50%","Rapid_30%")
ints6 <- c("Rapid_70%","Rapid_30%")

#EXTRACT OUTCOMES (DALYS)
index.70.effectiveness <-case.71 %>% pull(Total.Daly)
index.60.effectiveness <-case.61 %>% pull(Total.Daly)
index.50.effectiveness <-case.51 %>% pull(Total.Daly)
index.40.effectiveness <-case.41 %>% pull(Total.Daly)
index.30.effectiveness <-case.31 %>% pull(Total.Daly)
index.20.effectiveness <-case.21 %>% pull(Total.Daly)
index.0.effectiveness <-case.21 %>% pull(Total.DALY_baseline)

e1 <- cbind(index.20.effectiveness, index.0.effectiveness)
e2 <- cbind(index.30.effectiveness, index.20.effectiveness)
e3 <- cbind(index.40.effectiveness, index.30.effectiveness)
e4 <- cbind(index.50.effectiveness, index.0.effectiveness)
e5 <- cbind(index.60.effectiveness, index.50.effectiveness)
e6 <- cbind(index.70.effectiveness, index.60.effectiveness)

#extract cases
index.70.cases <-case.71 %>% pull(Total.cases)
index.60.cases <-case.61 %>% pull(Total.cases)
index.50.cases <-case.51 %>% pull(Total.cases)
index.40.cases <-case.41 %>% pull(Total.cases)
index.30.cases <-case.31 %>% pull(Total.cases)
index.20.cases<-case.21 %>% pull(Total.cases)
index.0.cases <-case.21 %>% pull(Total.cases_baseline)

#extract deaths
index.70.deaths <-case.71 %>% pull(Total.deaths)
index.60.deaths <-case.61 %>% pull(Total.deaths)
index.50.deaths <-case.51 %>% pull(Total.deaths)
index.40.deaths <-case.41 %>% pull(Total.deaths)
index.30.deaths <-case.31 %>% pull(Total.deaths)
index.20.deaths<-case.21 %>% pull(Total.deaths)
index.0.deaths <-case.21 %>% pull(Total.deaths.baseline)


#EXTRACT COSTS
index.70.cost <-case.71%>% pull(Total.Cost)
index.60.cost <-case.61%>% pull(Total.Cost)
index.50.cost <-case.51%>% pull(Total.Cost)
index.40.cost <-case.41%>% pull(Total.Cost)
index.30.cost <-case.31%>% pull(Total.Cost)
index.20.cost <-case.21%>% pull(Total.Cost)
index.0.cost <-case.21%>% pull(Total.Cost.baseline)

c1 <- cbind(index.20.cost,index.0.cost)
c2 <- cbind(index.30.cost,index.20.cost)
c3 <- cbind(index.40.cost,index.30.cost)
c4 <- cbind(index.50.cost,index.0.cost)
c5 <- cbind(index.60.cost,index.50.cost)
c6 <- cbind(index.70.cost,index.60.cost)

c1 <- -c1
c2 <- -c2
c3 <- -c3
c4 <- -c4
c5 <- -c5
c6 <- -c6


#ANALYSE THE BCEA FUNCTION
Kmax=25000
m1 <- bcea(e1,c1,ref=2,interventions=ints1,Kmax=Kmax)
m2 <- bcea(e2,c2,ref=2,interventions=ints2,Kmax=Kmax)
m3 <- bcea(e3,c3,ref=2,interventions=ints3,Kmax=Kmax)
m4 <- bcea(e4,c4,ref=2,interventions=ints4,Kmax=Kmax)
m5 <- bcea(e5,c5,ref=2,interventions=ints5,Kmax=Kmax)
m6 <- bcea(e6,c6,ref=2,interventions=ints6,Kmax=Kmax)
#output a summary of the results
summary(m1)


#PLOT ALL PLOTS
#plot health perspective
par(mfrow=c(2,3))

ceplane.plot(m1,wtp=919.105,
             pos="false",
             graph="base",
             ICER_sizes=1,
             xlab ="Effectiveness differential",
             ylab ="Cost differential",
             title ="0% vs 30% ",
             xlim=c(-8e3,1.2e5),
             ylim=c(-3.5e7,1.7e8))
polygon(x=c(-8e3, 1.2e5, 1.2e5, -8e3),
          y=c(-3.5e7, -3.5e7, 1.1e8, 0),
          col=rgb(0.5,0.4,0.4,0.25),border = FALSE)
 polygon(x=c(1e4, 1.15e5, 1.15e5, 1e4),
          y=c(1.45e8, 1.45e8, 1.7e8, 1.7e8),
         col="#FFFFFF",border = FALSE)
text(x = 4e4, y = 1e8,                # Add text element
     "ICER = 1086.5",cex=0.8)
abline(h=1.7e8)

#m2plot
ceplane.plot(m2,wtp=919.105,xlab ="Effectiveness differential",
             ylab ="Cost differential",
             title ="30% vs 50% ",pos="true",graph="base",
             ICER_sizes=2,
             xlim=c(-1.2e3,1.5e4),
             ylim=c(-2.5e7, 1e8))
  polygon(x=c(-1.2e3, 1.5e4, 1.5e4, -1.2e3),
          y=c(-2.5e7, -2.5e7, 1.4e7, -1.5e6),
          col=rgb(0.5,0.4,0.4,0.25),border = FALSE)
  polygon(x=c(2e3, 1.45e4, 1.45e4, 2e3),
          y=c(8e7, 8e7, 1e8, 1e8),
          col="#FFFFFF",border = FALSE)
  text(x = 4.5e3, y = 5e7,                # Add text element
       "ICER = 5951.1",cex=0.8)
  abline(h=1e8)

#m3plot

# ceplane.plot(m3,wtp=919.105,xlab ="Effectiveness differential",
#              ylab ="Cost differential",
#              title ="50% vs 70% ",pos="true",graph="base",
#              ICER_sizes=2,
#              xlim=c(-5e2,8e3),
#              ylim=c(-3e7,1.4e8),
#              cex.axis = 0.5,main="Main title")
#   polygon(x=c(-5e2,8e3, 8e3, -5e2),
#           y=c(-3e7, -3e7, 8e6, -2.3e6),
#           col=rgb(0.5,0.4,0.4,0.25),border = FALSE)
#   polygon(x=c(1e3, 8e3, 8e3, 1e3),
#           y=c(1.15e8, 1.15e8, 1.4e8, 1.4e8),
#           col="#FFFFFF",border = FALSE)
#   text(x = 3e3, y = 7e7,                # Add text element
#        "ICER = 16536.9",cex=0.8)
#   abline(h=1.4e8)
#
#
# #m4plot
# ceplane.plot(m4,wtp=919.105,xlab ="Effectiveness differential",
#              ylab ="Cost differential",
#              title ="0% vs 30% (Rapid) ",pos="true",graph="base",
#              ICER_sizes=2,
#              xlim=c(-1e4,1.5e5),
#              ylim=c(-8e7,3e8))
#   polygon(x=c(-1e4,1.5e5, 1.5e5, -1e4),
#           y=c(-8e7, -8e7, 1.4e8, -1e7),
#           col=rgb(0.5,0.4,0.4,0.25),border = FALSE)
#   polygon(x=c(2e4, 1.4e5, 1.4e5, 2e4),
#           y=c(2e8, 2e8, 3e8, 3e8),
#           col="#FFFFFF",border = FALSE)
#   text(x = 4e4, y = 2e8,                # Add text element
#        "ICER = 633.8",cex=0.8)
#   abline(h=3e8)
#
# #m5plot
#
# ceplane.plot(m5,wtp=919.105,xxlab ="Effectiveness differential",
#              ylab ="Cost differential",
#              title ="30% vs 50% (Rapid) ",label.pos="false",graph="base",
#              ICER_sizes=2,
#              xlim=c(-1e3,1.2e4),
#              ylim=c(-2e7,1e8),
#              cex.axis = 0.5)
#   polygon(x=c(-1e3,1.2e4, 1.2e4, -1e3),
#           y=c(-2e7, -2e7, 1.1e7, -1.2e6),
#           col=rgb(0.5,0.4,0.4,0.25),border = FALSE)
# polygon(x=c(2e3, 1.15e4, 1.15e4, 2e3),
#         y=c(8.5e7, 8.5e7, 1e8, 1e8),
#         col="#FFFFFF",border = FALSE)
# text(x = 4e3, y = 5e7,                # Add text element
#      "ICER = 9131.1",cex=0.8)
# abline(h=1e8)
#
# #m6plot
# ceplane.plot(m6,wtp=919.105,xlab ="Effectiveness differential",
#              ylab ="Cost differential",
#              title ="50% vs 70% (Rapid)",pos="true",graph="base",
#              ICER_sizes=2,
#              xlim=c(-8e2,1e4),
#              ylim=c(-2.5e7,1.4e8))
#   polygon(x=c(-8e2,1e4, 1e4, -8e2),
#           y=c(-2.5e7, -2.5e7, 1e7, -1.7e6),
#           col=rgb(0.5,0.4,0.4,0.25),border = FALSE)
# polygon(x=c(1e3, 1e4, 1e4, 1e3),
#         y=c(1.2e8, 1.2e8, 1.4e8, 1.4e8),
#         col="#FFFFFF",border = FALSE)
# text(x = 4e3, y = 7e7,                # Add text element
#      "ICER = 17640.4",cex=0.8)
# abline(h=1.4e8)


# #plto society perspective graphs
par(mfrow=c(2,3))

ceplane.plot(m1,wtp=919.105,
             pos="false",
             graph="base",
             ICER_sizes=1,
             xlab ="Effectiveness differential",
             ylab ="Cost differential",
             title ="0% vs 30% ",
             xlim=c(-7e3,1.2e5),
             ylim=c(-2e8,1.7e8))
polygon(x=c(-7e3, 1.2e5, 1.2e5, -7e3),
        y=c(-2e8, -2e8, 1.1e8, 0),
        col=rgb(0.5,0.4,0.4,0.25),border = FALSE)
polygon(x=c(1e4, 1.15e5, 1.15e5, 1e4),
        y=c(1.2e8, 1.2e8, 1.7e8, 1.7e8),
        col="#FFFFFF",border = FALSE)
text(x = 4e4, y = -1e8,                # Add text element
     "ICER = -812.4",cex=0.75)
abline(h=1.7e8)

#m2plot
ceplane.plot(m2,wtp=919.105,xlab ="Effectiveness differential",
             ylab ="Cost differential",
             title ="30% vs 50% ",pos="true",graph="base",
             ICER_sizes=2,
             xlim=c(-1e3,1.5e4),
             ylim=c(-2e7,1e8))
polygon(x=c(-1e3, 1.5e4, 1.5e4, -1e3),
        y=c(-2e7, -2e7, 1.4e7, -1.5e6),
        col=rgb(0.5,0.4,0.4,0.25),border = FALSE)
polygon(x=c(2e3, 1.45e4, 1.45e4, 2e3),
        y=c(8e7, 8e7, 1e8, 1e8),
        col="#FFFFFF",border = FALSE)
text(x = 4e3, y = 5e7,                # Add text element
     "ICER = 4033.9",cex=0.75)
abline(h=1e8)

#m3plot

ceplane.plot(m3,wtp=919.105,xlab ="Effectiveness differential",
             ylab ="Cost differential",
             title ="50% vs 70% ",pos="true",graph="base",
             ICER_sizes=2,
             xlim=c(-5e2,8e3),
             ylim=c(-3e7,1.4e8))
polygon(x=c(-1e3,8e3, 8e3, -1e3),
        y=c(-3e7, -3e7, 8e6, -2.3e6),
        col=rgb(0.5,0.4,0.4,0.25),border = FALSE)
polygon(x=c(1e3, 8e3, 8e3, 1e3),
        y=c(1.15e8, 1.15e8, 1.4e8, 1.4e8),
        col="#FFFFFF",border = FALSE)
text(x = 3e3, y = 7e7,                # Add text element
     "ICER = 14615.1",cex=0.75)
abline(h=1.4e8)


#m4plot
ceplane.plot(m4,wtp=919.105,xlab ="Effectiveness differential",
             ylab ="Cost differential",
             title ="0% vs 30% (Rapid) ",pos="true",graph="base",
             ICER_sizes=2,
             xlim=c(-1e4,1.5e5),
             ylim=c(-3e8,3e8))
polygon(x=c(-1e4,1.5e5, 1.5e5, -1e4),
        y=c(-3e8, -3e8, 1.4e8, -1e7),
        col=rgb(0.5,0.4,0.4,0.25),border = FALSE)
polygon(x=c(2e4, 1.45e5, 1.45e5, 2e4),
        y=c(2e8, 2e8, 3e8, 3e8),
        col="#FFFFFF",border = FALSE)
text(x = 4e4, y = -2e8,                # Add text element
     "ICER = -1269",cex=0.75)
abline(h=3e8)

#m5plot

ceplane.plot(m5,wtp=919.105,xxlab ="Effectiveness differential",
             ylab ="Cost differential",
             title ="30% vs 50% (Rapid) ",label.pos="false",graph="base",
             ICER_sizes=2,
             xlim=c(-8e2,1.2e4),
             ylim=c(-2e7,1e8),
             cex.axis = 0.5)
polygon(x=c(-8e2,1.2e4, 1.2e4, -8e2),
        y=c(-2e7, -2e7, 1.1e7, -1.2e6),
        col=rgb(0.5,0.4,0.4,0.25),border = FALSE)
polygon(x=c(2e3, 1.15e4, 1.15e4, 2e3),
        y=c(8.5e7, 8.5e7, 1e8, 1e8),
        col="#FFFFFF",border = FALSE)
text(x = 4e3, y = 5e7,                # Add text element
     "ICER = 7218.6",cex=0.75)
abline(h=1e8)

#m6plot
ceplane.plot(m6,wtp=919.105,xlab ="Effectiveness differential",
             ylab ="Cost differential",
             title ="50% vs 70% (Rapid)",pos="true",graph="base",
             ICER_sizes=2,
             xlim=c(-6e2,1e4),
             ylim=c(-1.9e7,1.4e8))
polygon(x=c(-6e2,1e4, 1e4, -6e2),
        y=c(-1.9e7, -1.9e7, 1e7, -1.7e6),
        col=rgb(0.5,0.4,0.4,0.25),border = FALSE)
polygon(x=c(1e3, 1e4, 1e4, 1e3),
        y=c(1.2e8, 1.2e8, 1.4e8, 1.4e8),
        col="#FFFFFF",border = FALSE)
text(x = 4e3, y = 7e7,                # Add text element
     "ICER = 15731",cex=0.75)
abline(h=1.4e8)

