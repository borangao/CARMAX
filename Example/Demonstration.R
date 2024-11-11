rm(list=ls())

devtools::install_github("ZikunY/CARMAX")  

library(CARMAX)
setwd('...')# Set your local directory
dt<-readRDS('./CARMAX/Example/Example.rds') 
dt<-readRDS('/Volumes/My_Passport/1000Genome/Multi_study/data/summ_data/package/Example.rds')
summ <- dt$`Summ data`
head(summ)
# Initialize lists
lambda.power <- 1
z.list <- list()
ld.list <- list()
full.w.list <- list()
lambda.list <- c()
label.list <- list()
input.prior.list <- list()
# Read the Z-scores and compute the LD for each ancestry 
# AFR ancestry
g <- 1
z.list[[g]] <- as.matrix(summ$`Z_AFR`)
lambda.list[[g]] <- lambda.power
ld.list[[g]] <- cor(dt$`AFR genotype`)
# EUR ancestry
g <- 2
z.list[[g]] <- as.matrix(summ$`Z_EUR`)
lambda.list[[g]] <- lambda.power
ld.list[[g]] <- cor(dt$`EUR genotype`)

prior.name='Spike-slab'

set.seed(124)

results <- CARMAX(z.list, ld.list, lambda.list = lambda.list, effect.size.prior = prior.name, outlier.switch = TRUE)


