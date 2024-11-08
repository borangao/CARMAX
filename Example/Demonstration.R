rm(list=ls())

devtools::install_github("ZikunY/CARMAX")  

library(CARMAX)
setwd('...')# Set your local directory
dt<-readRDS('./CARMAX/Example/Example.rds') 

summ<-read.table(paste0(''))

lambda.power<-1
z.list<-list()
ld.list<-list()
full.w.list<-list()
lambda.list<-c()
label.list<-list()
input.prior.list<-list()

summ<-dt[[1]]
p<-nrow(summ)

head(summ)

g=1
z.list[[g]]<-as.matrix(summ$BETA_EUR/summ$SE_EUR)
lambda.list[[g]]<-lambda.power

g=2
z.list[[g]]<-as.matrix(summ$BETA_AMR/summ$SE_AMR)
lambda.list[[g]]<-lambda.power

ld.list[[1]]<-dt[[2]] #######
ld.list[[2]]<-dt[[3]]
prior.name='Spike-slab'
results<-CARMAX(z.list,ld.list,lambda.list=lambda.list,outlier.switch=T,
          effect.size.prior=prior.name,
          epsilon.threshold=1e-3)



