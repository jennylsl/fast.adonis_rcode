# paramter
rm(list = ls())
library(data.table)
library(vegan)
library(readxl)
require(dplyr)
library(plyr)
# load distance matrix for PLCO
D.PLCO <- fread("plo/PLCOdistMatBrayCurtis.txt")[,-1]
phenoPLCO <- read_excel("plo/PLCO.phenotype.R2.Shilan.xlsx")
phenoPLCO <- phenoPLCO[match(colnames(D.PLCO),phenoPLCO$Sample_ID),]
# data preparation
phenoPLCO$agecat_at_collection <- as.factor(phenoPLCO$agecat_at_collection)
phenoPLCO$bmi_curr_cat <- as.factor(phenoPLCO$bmi_curr_cat)
phenoPLCO$race7_new <- as.factor(phenoPLCO$race7_new)
phenoPLCO$sex <- as.factor(phenoPLCO$sex)
phenoPLCO$educat_new <- as.factor(phenoPLCO$educat_new)
phenoPLCO$smk_new <- as.factor(phenoPLCO$smk_new)
phenoPLCO$alcohol_dpd <- as.factor(phenoPLCO$alcohol_dpd )
D<- D.PLCO

boot.times=100
time2<-Sys.time()
formula= ~ alcohol_dpd
X<- phenoPLCO[,c("agecat_at_collection","race7_new",
                 "sex","educat_new","smk_new","alcohol_dpd","bmi_curr_cat")]
X<- phenoPLCO[,c("agecat_at_collection",
                 "smk_new",
                 "alcohol_dpd","bmi_curr_cat","race7_new",
                 "sex","educat_new")]

ind.NA <-!is.na(diag(as.matrix(D))) &!(apply(as.matrix(X), 1, function(x)any(is.na(x)))) 

D.new <- as.matrix(D)[ind.NA,ind.NA]
w.new<-phenoPLCO$weight
boot.sample.size=dim(D)[1]
A <- -0.5*D.new^2
dim(A)
dim(X)
formula <- A ~ alcohol_dpd
permutations=9
boot.times= 10
boot.se = "WCB"
boot.sample.size= NULL
weights <-rep(1,dim(A)[1])
num_orders=2
order_list=NULL
by="terms"

data <- X
parallel = getOption("mc.cores")
num_orders<-2

if(num_orders>0){
  order_list <- list(sample(1:3,size=3),
                     c(1,3,2))
}
order_list

library(vegan)
git1<-adonis(D.new ~alcohol_dpd ,data=X,permutations = 20)
git2<-adonis2(D.new ~alcohol_dpd + race7_new+smk_new,data=X,permutations = 9)

fi1<-fast.adonis(A ~ alcohol_dpd , data = X, boot.times = 10,
                 permutations = 20,order_list = order_list)

fi1$aov.tab
library(vegan)
data(dune)
data(dune.env)
## default test by terms
adonis2(dune ~ Management*A1, data = dune.env)
