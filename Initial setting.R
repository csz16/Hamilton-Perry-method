#Load package
library("ggplot2")
library("MASS")
library("tidyverse")
library("tibble")
library("readr")
library("MASS")
library("tidyr")
library("plyr")
library("dplyr")
library("ftsa")
library("reshape2")
library("R2jags")



#Read data
pop_m <- readRDS("pop_m.rds")
pop_f <- readRDS("pop_f.rds")



#initial value
Age <- unique(data$Age)
Year <- unique(data$Year)
n <- length(Age)
m <- length(Year)
m.ratio <- 105.7/(105.7+100) 
f.ratio <- 1-m.ratio
method <- "classical"
