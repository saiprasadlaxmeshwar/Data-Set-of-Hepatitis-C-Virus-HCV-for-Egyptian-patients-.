#Importing the Hepatatis c Dataset 
HCV= read.csv("HCV-Egy-Data.csv")
HCV

#Summary
attach(HCV)
summary(HCV)

#Dimensions of the data set
NROW(HCV)
NCOL(HCV)

#Displaying the column names of the dataset
colnames(HCV)

#Another menthod for dimensions
dim(HCV)

#Preprocessing data was done but did'nt find any discrepancies.
na= is.na(HCV)
na
any(is.na(HCV))

#Displaying the first six rows of the datasets
head(HCV)
tail(HCV)

#Differentiating data set based on gender
Gen_male = HCV[HCV$Gender== '1',]
Gen_female = HCV[HCV$Gender=='2',]

#Exploring symptoms
#**********************Male Data Exploration*****************
#Fever
Fev_male = Gen_male[Gen_male$Fever == '2',]
Fev_male
summary(Fev_male)

#Vomiting and Nausea
Nau_male = Gen_male[Gen_male$Nausea.Vomting =='2',]
Nau_male
summary(Nau_male)

#Fatigue
Fat_male = Gen_male[Gen_male$Fatigue...generalized.bone.ache =='2',]
Fat_male
summary(Fat_male)

#Jaundice
Jau_male = Gen_male[Gen_male$Jaundice  =='2',]
Jau_male
summary(Jau_male)

#Stomack pain
sto_male = Gen_female[Gen_female$Epigastric.pain =='2',]
sto_male
summary(sto_male)

#*********************Female Data Exploration*****************
#Fever
Fev_female = Gen_female[Gen_female$Fever == '2',]
Fev_female
summary(Fev_female)

#Vomiting and Nausea
Nau_female = Gen_female[Gen_female$Nausea.Vomting =='2',]
Nau_female
summary(Nau_female)

#Fatigue
Fat_female = Gen_female[Gen_female$Fatigue...generalized.bone.ache =='2',]
Fat_female
summary(Fat_female)

#Jaundice
Jau_female = Gen_female[Gen_female$Jaundice  =='2',]
Jau_female
summary(Jau_female)

#Stomack pain
sto_female = Gen_female[Gen_female$Epigastric.pain =='2',]
sto_female
summary(sto_female)



#CORRELATION, COVARIANCE AND DISTANCE
covariance<-cov(HCV[,c(11:16,23)]) #variamce-covariance matrix created
correlation<-cor(HCV[,c(11:16,23)]) #standardized
#colmeans
cm<-colMeans(HCV[,c(11:16,23)])
distance<-dist(scale(HCV[,c(11:16,23)],center=FALSE))
#Calculating di(generalized distance for all observations of our data)
#before that first extract all numeric variable in a dataframe
x<-HCV[,c(11:16,23)]
d <- apply(x, MARGIN = 1, function(x) + t(x - cm) %*% solve(covariance) %*% (x - cm))



#Exlporation of the data for high chances of HCV Infection
#Here RNA.base value if it is more than 700000 units then virus is detected in high quantity.
#Here ALT.1 if value is greater than 57 then it is not normal.
#we sorted the data on these two components.

library(dplyr)
HCV_male =  HCV %>% filter(Gender == 1 & RNA.Base>= 700000 & ALT.1 >= 57)
HCV_male

HCV_female =  HCV %>% filter(Gender == 2 & RNA.Base>= 700000 & ALT.1 >= 57)
HCV_female


#Box Plot
boxplot(RNA.Base, main="RNA.BASE Box plot",yaxt="n", xlab="RNA", horizontal=TRUE)
boxplot(ALT.1, main="ALT.1 Box plot",yaxt="n", xlab="ALT", horizontal=TRUE)
boxplot(WBC, main="WBC Box plot",yaxt="n", xlab="WBC", horizontal=TRUE)
boxplot(RBC, main="WBC Box plot",yaxt="n", xlab="RBC", horizontal=TRUE)
boxplot(AST.1, main="AST.1 Box plot",yaxt="n", xlab="AST", horizontal=TRUE)


#plotting, Are they in a straight line.  
#Male Plotting of the dataset is done for five different attributes.
qqnorm(HCV_male[,"RNA.Base"], main = "RNA.Base"); qqline(HCV_male[,"RNA.Base"])
qqnorm(HCV_male[,"ALT.1"], main = "ALT.1"); qqline(HCV_male[,"ALT.1"])
qqnorm(HCV_male[,"WBC"], main = "WBC"); qqline(HCV_male[,"WBC"])
qqnorm(HCV_male[,"RBC"], main = "RBC"); qqline(HCV_male[,"RBC"])
qqnorm(HCV_male[,"AST.1"], main = "AST.1"); qqline(HCV_male[,"AST.1"])

#Female, Are they in a straight line.
#FeMale Plotting of the dataset is done for five different attributes.
qqnorm(HCV_female[,"RNA.Base"], main = "RNA.Base"); qqline(HCV_female[,"RNA.Base"])
qqnorm(HCV_female[,"ALT.1"], main = "ALT.1"); qqline(HCV_female[,"ALT.1"])
qqnorm(HCV_female[,"WBC"], main = "WBC"); qqline(HCV_female[,"WBC"])
qqnorm(HCV_female[,"RBC"], main = "RBC"); qqline(HCV_female[,"RBC"])
qqnorm(HCV_female[,"AST.1"], main = "AST.1"); qqline(HCV_female[,"AST.1"])


#Visualisatiom
#Chiplot
library(HSAUR2)
library(tools)
library(MVA)

#Chiplot
#For male data
with(HCV_male, chiplot(RNA.Base, ALT.1))

#For Female Data
with(HCV_female, chiplot(RNA.Base, ALT.1))



library(GGally)

ggpairs(HCV_male, columns=c("AST.1","RNA.EOT","WBC","ALT.1", "RBC"), color="Survivorship")
ggpairs(HCV_female, columns=c("AST.1","RNA.EOT","WBC","ALT.1", "RBC"), color="Survivorship")
summary(lm(data = HCV , RNA.EOT~Age))
summary(lm(data = HCV , RNA.EOT~Gender))
summary(lm(data = HCV , RNA.EOT~WBC))
summary(lm(data = HCV , RNA.EOT~ALT.1))
cor(HCV)


#Pca || T-test || F-test



#Get the Correlations between the measurements
cor(HCV)
# Using prcomp to compute the principal components (eigenvalues and eigenvectors). 
#With scale=TRUE, variable means are set to zero, and variances set to one
x_pca <- prcomp(HCV,scale=TRUE)
x_pca
summary(x_pca)
x_pca$rotation
# Eigenvalues are sdev^2
(eigen_x <- x_pca$sdev^2)
names(eigen_x) <- paste("PC",1:29)
eigen_x
sumlambdas <- sum(eigen_x)
sumlambdas #total sample variance
propvar <- eigen_x/sumlambdas
propvar
cumvar_x <- cumsum(propvar)
cumvar_x
matlambdas <- rbind(eigen_x,propvar,cumvar_x)
matlambdas
rownames(matlambdas) <- c("Eigenvalues","Prop. variance","Cum. prop. variance")
round(matlambdas,4)
# Sample scores stored in x_pca$x
x_pca$x
xtyp_pca <- cbind(data.frame(HCV),x_pca$x)
xtyp_pca
colnames(xtyp_pca)[colnames(xtyp_pca)=="HCV"] <- "p"


#T-test

t.test(xtyp_pca$PC1,xtyp_pca$p,var.equal = TRUE)
t.test(xtyp_pca$PC2,xtyp_pca$p,var.equal = TRUE)
t.test(xtyp_pca$PC3,xtyp_pca$p,var.equal = TRUE)
t.test(xtyp_pca$PC4,xtyp_pca$p,var.equal = TRUE)
t.test(xtyp_pca$PC5,xtyp_pca$p,var.equal = TRUE)
t.test(xtyp_pca$PC6,xtyp_pca$p,var.equal = TRUE)
t.test(xtyp_pca$PC7,xtyp_pca$p,var.equal = TRUE)
t.test(xtyp_pca$PC8,xtyp_pca$p,var.equal = TRUE)
t.test(xtyp_pca$PC9,xtyp_pca$p,var.equal = TRUE)
t.test(xtyp_pca$PC10,xtyp_pca$p,var.equal = TRUE)
t.test(xtyp_pca$PC11,xtyp_pca$p,var.equal = TRUE)
t.test(xtyp_pca$PC12,xtyp_pca$p,var.equal = TRUE)
t.test(xtyp_pca$PC13,xtyp_pca$p,var.equal = TRUE)
t.test(xtyp_pca$PC14,xtyp_pca$p,var.equal = TRUE)
t.test(xtyp_pca$PC15,xtyp_pca$p,var.equal = TRUE)
t.test(xtyp_pca$PC16,xtyp_pca$p,var.equal = TRUE)
t.test(xtyp_pca$PC17,xtyp_pca$p,var.equal = TRUE)
t.test(xtyp_pca$PC18,xtyp_pca$p,var.equal = TRUE)
t.test(xtyp_pca$PC19,xtyp_pca$p,var.equal = TRUE)
t.test(xtyp_pca$PC20,xtyp_pca$p,var.equal = TRUE)
t.test(xtyp_pca$PC21,xtyp_pca$p,var.equal = TRUE)
t.test(xtyp_pca$PC22,xtyp_pca$p,var.equal = TRUE)
t.test(xtyp_pca$PC23,xtyp_pca$p,var.equal = TRUE)
t.test(xtyp_pca$PC24,xtyp_pca$p,var.equal = TRUE)
t.test(xtyp_pca$PC25,xtyp_pca$p,var.equal = TRUE)
t.test(xtyp_pca$PC26,xtyp_pca$p,var.equal = TRUE)
t.test(xtyp_pca$PC27,xtyp_pca$p,var.equal = TRUE)
t.test(xtyp_pca$PC28,xtyp_pca$p,var.equal = TRUE)
t.test(xtyp_pca$PC29,xtyp_pca$p,var.equal = TRUE)




