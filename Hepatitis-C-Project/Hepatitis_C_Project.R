*************************************
#Data Exploration
*************************************  
  
  
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




#********************************
#Principal Component Analysis
#********************************


#Importing the Hepatatis c Dataset 
HCV_NEW= read.csv("HCV-Egy-Data.csv")
HCV<-  HCV_NEW
attach(HCV)
library(dplyr)
HCV$Survivorship <- if_else( RNA.EOT>= 400000 , 'NC','C')
Survivorship
cbind(Survivorship,HCV)
#Summary
summary(HCV)
#Dimensions of the data set
NROW(HCV)
NCOL(HCV)

#Another menthod for dimensions
dim(HCV)

#Preprocessing data was done but did'nt find any discrepancies.
na= is.na(HCV)
na
any(is.na(HCV))

#Displaying the first six rows of the datasets
head(HCV)
tail(HCV)

correlation<-cor(HCV[1:29]) #standardized
correlation

# Using prcomp to compute the principal components (eigenvalues and eigenvectors). With scale=TRUE, variable means are set to zero, and variances set to one
hcv_pca <- prcomp(HCV[1:29],scale=TRUE)
hcv_pca
summary(hcv_pca)
# sample scores stored in sparrows_pca$x
# singular values (square roots of eigenvalues) stored in sparrow_pca$sdev
# loadings (eigenvectors) are stored in sparrows_pca$rotation
# variable means stored in sparrows_pca$center
# variable standard deviations stored in sparrows_pca$scale
# A table containing eigenvalues and %'s accounted, follows
# Eigenvalues are sdev^2
(eigen_hcv <- hcv_pca$sdev^2)
names(eigen_hcv) <- paste("PC",1:29,sep="")
eigen_hcv
sumlambdas <- sum(eigen_hcv)
sumlambdas
propvar <- eigen_hcv/sumlambdas
propvar
cumvar_hcv <- cumsum(propvar)
cumvar_hcv
matlambdas <- rbind(eigen_hcv,propvar,cumvar_hcv)
rownames(matlambdas) <- c("Eigenvalues","Prop. variance","Cum. prop. variance")
round(matlambdas,29)
summary(hcv_pca)
hcv_pca$rotation
print(hcv_pca)
# Sample scores stored in hcv_pca$x
hcv_pca$x
############################################### Identifying the scores by their survival status
hcvtyp_pca <- cbind(data.frame(Survivorship),hcv_pca$x)
hcvtyp_pca
# Means of scores for all the PC's classified by Survival status
tabmeansPC <- aggregate(hcvtyp_pca[,2:30],by=list(Survivorship=HCV$Survivorship),mean)
tabmeansPC
tabmeansPC <- tabmeansPC[rev(order(tabmeansPC$Survivorship)),]
tabmeansPC
tabfmeans <- t(tabmeansPC[,-1])
tabfmeans
colnames(tabfmeans) <- t(as.vector(tabmeansPC[1]))
tabfmeans
# Standard deviations of scores for all the PC's classified by Survival status
tabsdsPC <- aggregate(hcvtyp_pca[,2:30],by=list(Survivorship=HCV$Survivorship),sd)
tabfsds <- t(tabsdsPC[,-1])
colnames(tabfsds) <- t(as.vector(tabsdsPC[1]))
tabfsds

#t-test
t.test(PC1~HCV$Survivorship,data=hcvtyp_pca)
t.test(PC2~HCV$Survivorship,data=hcvtyp_pca)
t.test(PC3~HCV$Survivorship,data=hcvtyp_pca)
t.test(PC4~HCV$Survivorship,data=hcvtyp_pca)
t.test(PC5~HCV$Survivorship,data=hcvtyp_pca)
t.test(PC6~HCV$Survivorship,data=hcvtyp_pca)
t.test(PC7~HCV$Survivorship,data=hcvtyp_pca)
t.test(PC8~HCV$Survivorship,data=hcvtyp_pca)
t.test(PC9~HCV$Survivorship,data=hcvtyp_pca)
t.test(PC10~HCV$Survivorship,data=hcvtyp_pca)
t.test(PC11~HCV$Survivorship,data=hcvtyp_pca)
t.test(PC12~HCV$Survivorship,data=hcvtyp_pca)
t.test(PC13~HCV$Survivorship,data=hcvtyp_pca)
t.test(PC14~HCV$Survivorship,data=hcvtyp_pca)
t.test(PC15~HCV$Survivorship,data=hcvtyp_pca)
t.test(PC16~HCV$Survivorship,data=hcvtyp_pca)
t.test(PC17~HCV$Survivorship,data=hcvtyp_pca)
t.test(PC18~HCV$Survivorship,data=hcvtyp_pca)
t.test(PC19~HCV$Survivorship,data=hcvtyp_pca)
t.test(PC20~HCV$Survivorship,data=hcvtyp_pca)
t.test(PC21~HCV$Survivorship,data=hcvtyp_pca)
t.test(PC22~HCV$Survivorship,data=hcvtyp_pca)
t.test(PC23~HCV$Survivorship,data=hcvtyp_pca)
t.test(PC24~HCV$Survivorship,data=hcvtyp_pca)
t.test(PC25~HCV$Survivorship,data=hcvtyp_pca)
t.test(PC26~HCV$Survivorship,data=hcvtyp_pca)
t.test(PC27~HCV$Survivorship,data=hcvtyp_pca)
t.test(PC28~HCV$Survivorship,data=hcvtyp_pca)
t.test(PC29~HCV$Survivorship,data=hcvtyp_pca)


# F ratio tests
var.test(PC1~HCV$Survivorship,data=hcvtyp_pca)
var.test(PC2~HCV$Survivorship,data=hcvtyp_pca)
var.test(PC3~HCV$Survivorship,data=hcvtyp_pca)
var.test(PC4~HCV$Survivorship,data=hcvtyp_pca)
var.test(PC5~HCV$Survivorship,data=hcvtyp_pca)
var.test(PC6~HCV$Survivorship,data=hcvtyp_pca)
var.test(PC7~HCV$Survivorship,data=hcvtyp_pca)
var.test(PC8~HCV$Survivorship,data=hcvtyp_pca)
var.test(PC9~HCV$Survivorship,data=hcvtyp_pca)
var.test(PC10~HCV$Survivorship,data=hcvtyp_pca)
var.test(PC11~HCV$Survivorship,data=hcvtyp_pca)
var.test(PC12~HCV$Survivorship,data=hcvtyp_pca)
var.test(PC13~HCV$Survivorship,data=hcvtyp_pca)
var.test(PC14~HCV$Survivorship,data=hcvtyp_pca)
var.test(PC15~HCV$Survivorship,data=hcvtyp_pca)
var.test(PC16~HCV$Survivorship,data=hcvtyp_pca)
var.test(PC17~HCV$Survivorship,data=hcvtyp_pca)
var.test(PC18~HCV$Survivorship,data=hcvtyp_pca)
var.test(PC19~HCV$Survivorship,data=hcvtyp_pca)
var.test(PC20~HCV$Survivorship,data=hcvtyp_pca)
var.test(PC21~HCV$Survivorship,data=hcvtyp_pca)
var.test(PC22~HCV$Survivorship,data=hcvtyp_pca)
var.test(PC23~HCV$Survivorship,data=hcvtyp_pca)
var.test(PC24~HCV$Survivorship,data=hcvtyp_pca)
var.test(PC25~HCV$Survivorship,data=hcvtyp_pca)
var.test(PC26~HCV$Survivorship,data=hcvtyp_pca)
var.test(PC27~HCV$Survivorship,data=hcvtyp_pca)
var.test(PC28~HCV$Survivorship,data=hcvtyp_pca)
var.test(PC29~HCV$Survivorship,data=hcvtyp_pca)

# Plotting the scores for the first and second components
plot(hcvtyp_pca$PC1, hcvtyp_pca$PC2,pch=ifelse(hcvtyp_pca$Survivorship == "S",1,16),xlab="PC1", ylab="PC2", main="49 HCV against values for PC1 & PC2")
abline(h=0)
abline(v=0)
legend("bottomleft", legend=c("Cured","Not_Cured"), pch=c(1,16))
plot(eigen_hcv, xlab = "Component number", ylab = "Component variance", type = "l", main = "Scree diagram")
plot(log(eigen_hcv), xlab = "Component number",ylab = "log(Component variance)", type="l",main = "Log(eigenvalue) diagram")
print(summary(hcv_pca))
View(hcv_pca)
diag(cov(hcv_pca$x))
xlim <- range(hcv_pca$x[,1])
hcv_pca$x[,1]
hcv_pca$x
plot(hcv_pca$x,xlim=xlim,ylim=xlim)
hcv_pca$rotation[,1]
hcv_pca$rotation
plot(HCV[,-1])
hcv_pca$x
plot(hcv_pca)
#get the original value of the data based on PCA
center <- hcv_pca$center
scale <- hcv_pca$scale
new_HCV1 <- as.matrix(HCV[,-1])
new_HCV1
drop(scale(new_HCV1,center=center, scale=scale)%*%hcv_pca$rotation[,1])
predict(hcv_pca)[,1]
#The aboved two gives us the same thing. predict is a good function to know.
out <- sapply(1:5, function(i){plot(HCV$Survivorship,hcv_pca$x[,i],xlab=paste("PC",i,sep=""),ylab="Survivorship")})
pairs(hcv_pca$x[,1:5], ylim = c(-6,4),xlim = c(-6,4),panel=function(x,y,...){text(x,y,HCV$Survivorship)})


#******************************************************************
# Factor Analysis
#*******************************************************************

HCV<- read.csv("HCV-Egy-Data.csv")
View(HCV)
attach(HCV)
library(dplyr)
Survivorship = HCV$Survivorship <- if_else( RNA.EOT>= 400000 , 'NC','C')
cbind(data.frame(Survivorship),HCV)

HCV_pca = select(HCV,Age,Gender,Nausea.Vomting,Jaundice,Epigastric.pain,WBC,Plat,AST.1,ALT.1,ALT.after.24.w,RNA.Base,RNA.EOT)

# Computing Correlation Matrix
corr.HCV <- cor(HCV[1:29])
corr.HCV

hcv_pca <- prcomp(HCV_pca, scale=TRUE)
summary(hcv_pca)
plot(hcv_pca)
# A table containing eigenvalues and %'s accounted, follows. Eigenvalues are the sdev^2
(eigen_hcv <- round(hcv_pca$sdev^2,2))
names(eigen_hcv) <- paste("PC",1:12,sep="")
eigen_hcv
sumlambdas <- sum(eigen_hcv)
sumlambdas
propvar <- round(eigen_hcv/sumlambdas,2)
propvar
cumvar_hcv <- cumsum(propvar)
cumvar_hcv
matlambdas <- rbind(eigen_hcv,propvar,cumvar_hcv)
matlambdas
rownames(matlambdas) <- c("Eigenvalues","Prop. variance","Cum. prop. variance")
rownames(matlambdas)
eigvec.hcv <- hcv_pca$rotation
print(hcv_pca)
# Taking the first four PCs to generate linear combinations for all the variables with four factors
pcafactors.hcv <- eigvec.hcv[,1:4]
pcafactors.hcv
# Multiplying each column of the eigenvector's matrix by the square-root of the corresponding eigenvalue in order to get the factor loadings
unrot.fact.hcv <- sweep(pcafactors.hcv,MARGIN=2,hcv_pca$sdev[1:4],`*`)
unrot.fact.hcv
# Computing communalities
communalities.hcv <- rowSums(unrot.fact.hcv^2)
communalities.hcv
# Performing the varimax rotation. The default in the varimax function is norm=TRUE thus, Kaiser normalization is carried out
rot.fact.hcv <- varimax(unrot.fact.hcv)
View(unrot.fact.hcv)
rot.fact.hcv
# The print method of varimax omits loadings less than abs(0.1). In order to display all the loadings, it is necessary to ask explicitly the contents of the object $loadings
fact.load.hcv <- rot.fact.hcv$loadings[1:9,1:4]
fact.load.hcv
# Computing the rotated factor scores for the 30 European Countries. Notice that signs are reversed for factors F2 (PC2), F3 (PC3) and F4 (PC4)
scale.hcv <- scale(HCV[1:29])
scale.hcv


library(psych)
fit.pc <- principal(HCV[1:29], nfactors=4, rotate="varimax")
fit.pc
round(fit.pc$values, 3)
fit.pc$loadings
# Loadings with more digits
for (i in c(1,3,2,4)) { print(fit.pc$loadings[[1,i]])}
# Communalities
fit.pc$communality
# Rotated factor scores, Notice the columns ordering: RC1, RC3, RC2 and RC4
fit.pc$scores
# Play with FA utilities

fa.parallel(HCV[1:29]) # See factor recommendation
fa.plot(fit.pc) # See Correlations within Factors
fa.diagram(fit.pc) # Visualize the relationship
vss(HCV[1:29]) # See Factor recommendations for a simple structure


**********************************************
#Multi-Linear Regression
**********************************************
  
#Importing the Hepatatis c Dataset 
HCV <- read.csv("C:/Users/saiprasad/Desktop/HCV-Egy-Data.csv")
library(dplyr)
attach(HCV)
Survivorship =HCV$Survivorship <- if_else( RNA.EOT>= 400000 , 0,1)
cbind(data.frame(Survivorship),HCV)

HCV_lr = select(HCV,RNA.12,RNA.EF,RNA.EOT,RNA.Base)

#Multiple Regression
View(HCV)
# Performing multiple regression on HCV dataset
fit <- lm(Survivorship~RNA.EF+RNA.EOT+RNA.Base+RNA.12,data=HCV)
#show the results
summary(fit)
#Summary has three sections. Section1: How well does the model fit the data (before Coefficients). Section2: Is the hypothesis supported? (until sifnif codes). Section3: How well does data fit the model (again).
# Useful Helper Functions
coefficients(fit)
#install.packages("GGally", lib="/Library/Frameworks/R.framework/Versions/3.5/Resources/library")
library(GGally)
ggpairs(data=HCV_lr, title="HCV Data")
confint(fit,level=0.95)
# Predicted Values
fitted(fit)
residuals(fit)
#Anova Table
anova(fit)
vcov(fit)
cov2cor(vcov(fit))
temp <- influence.measures(fit)
temp
View(temp)
#diagnostic plots
plot(fit)
# Assessing Outliers
outlierTest(fit)
qqPlot(fit, main="QQ Plot")
leveragePlots(fit) # leverage plots
# Influential Observations
# Cook's D plot
# identify D values > 4/(n-k-1)
cutoff <- 4/((nrow(HCV_lr)-length(fit$coefficients)-2))
plot(fit, which=4, cook.levels=cutoff)
# Influence Plot
influencePlot(fit, id.method="identify", main="Influence Plot", sub="Circle size is proportial to Cook's Distance" )
# Normality of Residuals
# qq plot for studentized resid
qqPlot(fit, main="QQ Plot")
# distribution of studentized residuals
library(MASS)
sresid <- studres(fit)
hist(sresid, freq=FALSE,
     main="Distribution of Studentized Residuals")
xfit<-seq(min(sresid),max(sresid),length=40)
yfit<-dnorm(xfit)
lines(xfit, yfit)
#Non-constant Error Variance
# Evaluate homoscedasticity
# non-constant error variance test
ncvTest(fit)
# plot studentized residuals vs. fitted values
spreadLevelPlot(fit)
#Multi-collinearity
# Evaluate Collinearity
vif(fit) # variance inflation factors
sqrt(vif(fit)) > 2 # problem?
#Nonlinearity
# component + residual plot
crPlots(fit)
# Ceres plots
ceresPlots(fit)
#Non-independence of Errors
# Test for Autocorrelated Errors
durbinWatsonTest(fit)
# Global test of model assumptions
install.packages("gvlma", lib="/Library/Frameworks/R.framework/Versions/3.5/Resources/library")
library(gvlma)
gvmodel <- gvlma(fit)
summary(gvmodel)
fit
summary(fit)
fit1 <- fit
fit2 <- lm(Survivorship~RNA.EF+RNA.EOT+RNA.12, data=HCV)
# compare models
anova(fit1, fit2)
step <- stepAIC(fit, direction="both")
step$anova # display results
install.packages("leaps", lib="/Library/Frameworks/R.framework/Versions/3.5/Resources/library")
library(leaps)
leaps<-regsubsets(Survivorship~RNA.EF+RNA.EOT+RNA.Base+RNA.12, data=HCV_lr,nbest=10)
# view results
summary(leaps)
# plot a table of models showing variables in each model.
# models are ordered by the selection statistic.
plot(leaps)
plot(leaps,scale="r2")
subsets(leaps, statistic="rsq")
# All Subsets Regression
plot(leaps,scale="bic")
summary(leaps)
?regsubsets
View(leaps)
leaps
coef(leaps,1:4)
#prediction of the cured 
predict.lm(fit, data.frame(RNA.12 =585688 ,RNA.EF=582301,RNA.EOT= 744463,RNA.Base=1041941) )



*******************************************
#Logistic Regression
*******************************************
  
  #install.packages("cowplot", lib="/Library/Frameworks/R.framework/Versions/3.5/Resources/library")
  library(cowplot)
library(dplyr)
library(ggplot2)
data <- read.csv("C:/Users/saiprasad/Desktop/Fall 2019/Multi analysis/MVA/Project/Dataset/HCV-EGY-Data.csv")
attach(data)
Survivorship = data$Survivorship <- if_else( RNA.EOT>= 400000 , 'NC','C')
cbind(data.frame(Survivorship),data)

data$Survivorship <- as.factor(data$Survivorship)
#####################################
##
## Reformat the data so that it is
## 1) Easy to use (add nice column names)
## 2) Interpreted correctly by glm()..
##
#####################################
head(data) # you see data, but no column names
str(data)
# this shows that we need to tell R which columns contain factors it also shows us that there are some missing values. There are "?"s
## in the dataset. These are in the "ca" and "thal" columns. First, convert "?"s to NAs...
data[data == "?"] <- NA
## Now add factors for variables that are factors and clean up the factors that had missing data...
data[data$Gender == 1,]$Gender <- "M"
data[data$Gender == 2,]$Gender <- "F"
data$Gender <- as.factor(data$Gender)
data[data$Fever == 1,]$Fever <- "No"
data[data$Fever == 2,]$Fever <- "Yes"
data$Fever <- as.factor(data$Fever)
data[data$Nausea.Vomting == 1,]$Nausea.Vomting <- "No"
data[data$Nausea.Vomting == 2,]$Nausea.Vomting <- "Yes"
data$Nausea.Vomting <- as.factor(data$Nausea.Vomting)
data[data$Headache == 1,]$Headache <- "No"
data[data$Headache == 2,]$Headache <- "Yes"
data$Headache <- as.factor(data$Headache)
data[data$Diarrhea == 1,]$Diarrhea <- "No"
data[data$Diarrhea == 2,]$Diarrhea <- "Yes"
data$Diarrhea <- as.factor(data$Diarrhea)
data[data$Fatigue...generalized.bone.ache == 1,]$Fatigue...generalized.bone.ache <- "No"
data[data$Fatigue...generalized.bone.ache == 2,]$Fatigue...generalized.bone.ache <- "Yes"
data$Fatigue...generalized.bone.ache <- as.factor(data$Fatigue...generalized.bone.ache)
data[data$Jaundice == 1,]$Jaundice <- "No"
data[data$Jaundice == 2,]$Jaundice <- "Yes"
data$Jaundice <- as.factor(data$Jaundice)
data[data$Epigastric.pain == 1,]$Epigastric.pain <- "No"
data[data$Epigastric.pain == 2,]$Epigastric.pain <- "Yes"
data$Epigastric.pain <- as.factor(data$Epigastric.pain)
data[data$Baselinehistological.staging == 1,]$Baselinehistological.staging <- "Portal Fibrosis"
data[data$Baselinehistological.staging == 2,]$Baselinehistological.staging<- "Few Septa"
data[data$Baselinehistological.staging == 3,]$Baselinehistological.staging <- "Many Septa "
data[data$Baselinehistological.staging == 4,]$Baselinehistological.staging <- "Cirrhosis"
data$Baseline.histological.Grading <- as.factor(data$Baseline.histological.Grading)
data$Baselinehistological.staging <- as.factor(data$Baselinehistological.staging)
str(data)

###################################
xtabs(~ Survivorship + Gender, data=data)
xtabs(~ Survivorship + Fever, data=data)
xtabs(~ Survivorship + Nausea.Vomting, data=data)
xtabs(~ Survivorship + Headache, data=data)
xtabs(~ Survivorship + Diarrhea, data=data)
xtabs(~ Survivorship + Fatigue...generalized.bone.ache, data=data)
xtabs(~ Survivorship + Jaundice, data=data)
xtabs(~ Survivorship + Epigastric.pain, data=data)
xtabs(~ Survivorship + Baselinehistological.staging, data=data)
## Now we are ready for some logistic regression. First we'll create a very
## simple model that uses sex to predict heart disease
##
xtabs(~ Survivorship + Gender, data=data)
## Most of the females are healthy and most of the males are unhealthy.
## Being female is likely to decrease the odds in being unhealthy.
##    In other words, if a sample is female, the odds are against it that it
##    will be unhealthy
## Being male is likely to increase the odds in being unhealthy...
##    In other words, if a sample is male, the odds are for it being unhealthy
logistic_simple <- glm(Survivorship ~ Gender, data=data, family="binomial")
summary(logistic_simple)
## The intercept is the log(odds) a female will be unhealthy. This is because
## female is the first factor in "sex" (the factors are ordered,
## alphabetically by default,"female", "male")
## Now let's look at the second coefficient...
##   sexM        1.2737     0.2725   4.674 2.95e-06 ***
##
## sexM is the log(odds ratio) that tells us that if a sample has sex=M, the
## odds of being unhealthy are, on a log scale, 1.27 times greater than if
## a sample has sex=F.
female.log.odds <- log(253 /425)
female.log.odds
# Now you know how these are calculated
male.log.odds.ratio <- log((229 / 478) / (253/425))
male.log.odds.ratio
## Now calculate the overall "Pseudo R-squared" and its p-value
## NOTE: Since we are doing logistic regression...
## Null devaiance = 2*(0 - LogLikelihood(null model))
##               = -2*LogLikihood(null model)
## Residual deviance = 2*(0 - LogLikelihood(proposed model))
##                   = -2*LogLikelihood(proposed model)
ll.null <- logistic_simple$null.deviance/-2
ll.proposed <- logistic_simple$deviance/-2
ll.null
ll.proposed
## McFadden's Pseudo R^2 = [ LL(Null) - LL(Proposed) ] / LL(Null)
(ll.null - ll.proposed) / ll.null
## chi-square value = 2*(LL(Proposed) - LL(Null))
## p-value = 1 - pchisq(chi-square value, df = 2-1)
1 - pchisq(2*(ll.proposed - ll.null), df=1)
1 - pchisq((logistic_simple$null.deviance - logistic_simple$deviance), df=1)
## Lastly, let's  see what this logistic regression predicts, given
## that a patient is either female or male (and no other data about them).
predicted.data <- data.frame(probability.of.Survivorship=logistic_simple$fitted.values,Gender=data$Gender)
predicted.data
## We can plot the data...
ggplot(data=predicted.data, aes(x=Gender, y=probability.of.Survivorship)) +
  geom_point(aes(color=Gender), size=5) +
  xlab("Gender") +
  ylab("Predicted probability of getting HCV Disease")
## Since there are only two probabilities (one for females and one for males),
## we can use a table to summarize the predicted probabilities.
xtabs(~ probability.of.Survivorship + Gender, data=predicted.data)
#####################################
##
## Now we will use all of the data available to predict heart disease. This is not the best way to do this
##
#####################################
logistic <- glm(Survivorship ~ ., data=data, family="binomial")
summary(logistic)
## Now calculate the overall "Pseudo R-squared" and its p-value
ll.null <- logistic$null.deviance/-2
ll.proposed <- logistic$deviance/-2
## McFadden's Pseudo R^2 = [ LL(Null) - LL(Proposed) ] / LL(Null)
(ll.null - ll.proposed) / ll.null
## The p-value for the R^2
1 - pchisq(2*(ll.proposed - ll.null), df=(length(logistic$coefficients)-1))
## now we can plot the data
predicted.data <- data.frame(probability.of.Survivorship=logistic$fitted.values,Survivorship=data$Survivorship)
predicted.data <- predicted.data[order(predicted.data$probability.of.Survivorship, decreasing=FALSE),]
predicted.data$rank <- 1:nrow(predicted.data)
## Lastly, we can plot the predicted probabilities for each sample having
## heart disease and color by whether or not they actually had heart disease
ggplot(data=predicted.data, aes(x=rank, y=probability.of.Survivorship)) +
  geom_point(aes(color=Survivorship), alpha=1, shape=4, stroke=2) +
  xlab("Index") +
  ylab("Predicted probability of getting HCV disease")
# Few packages for confusion matrix. Lets look at them one by one
#install.packages("regclass", lib="/Library/Frameworks/R.framework/Versions/3.5/Resources/library")
library(regclass)
confusion_matrix(logistic)
#install.packages("caret", lib="/Library/Frameworks/R.framework/Versions/3.5/Resources/library")
library(caret)
pdata <- predict(logistic,newdata=data,type="response" )
pdata
data$Survivorship
#pdataF <- as.factor(ifelse(test=as.numeric(pdata>0.5) == 0, yes="Healthy", no="Unhealthy"))
#install.packages("e1071", lib="/Library/Frameworks/R.framework/Versions/3.5/Resources/library")
library(e1071)
#confusionMatrix(pdataF, data$Survivorship)
#install.packages("pROC", lib="/Library/Frameworks/R.framework/Versions/3.5/Resources/library")
library(pROC)
roc(data$Survivorship,logistic$fitted.values,plot=TRUE)
par(pty = "s")
roc(data$Survivorship,logistic$fitted.values,plot=TRUE)
## NOTE: By default, roc() uses specificity on the x-axis and the values range
## from 1 to 0. This makes the graph look like what we would expect, but the
## x-axis itself might induce a headache. To use 1-specificity (i.e. the
## False Positive Rate) on the x-axis, set "legacy.axes" to TRUE.
roc(data$Survivorship,logistic$fitted.values,plot=TRUE, legacy.axes=TRUE)
roc(data$Survivorship,logistic$fitted.values,plot=TRUE, legacy.axes=TRUE, xlab="False Positive Percentage", ylab="True Postive Percentage")
roc(data$Survivorship,logistic$fitted.values,plot=TRUE, legacy.axes=TRUE, xlab="False Positive Percentage", ylab="True Postive Percentage", col="#377eb8", lwd=4)
roc(data$Survivorship,logistic$fitted.values,plot=TRUE, legacy.axes=TRUE, xlab="False Positive Percentage", ylab="True Postive Percentage", col="#377eb8", lwd=4)
## If we want to find out the optimal threshold we can store the
## data used to make the ROC graph in a variable...
roc.info <- roc(data$Survivorship, logistic$fitted.values, legacy.axes=TRUE)
str(roc.info)
roc.df <- data.frame(tpp=roc.info$sensitivities*100, ## tpp = true positive percentage
                     fpp=(1 - roc.info$specificities)*100, ## fpp = false positive precentage
                     thresholds=roc.info$thresholds)
roc.df
head(roc.df) ## head() will show us the values for the upper right-hand corner of the ROC graph, when the threshold is so low
## (negative infinity) that every single sample is called "obese".
## Thus TPP = 100% and FPP = 100%
tail(roc.df) ## tail() will show us the values for the lower left-hand corner
## of the ROC graph, when the threshold is so high (infinity)
## that every single sample is called "not obese".
## Thus, TPP = 0% and FPP = 0%
## now let's look at the thresholds between TPP 60% and 80%
roc.df[roc.df$tpp > 60 & roc.df$tpp < 80,]
roc(data$Survivorship,logistic$fitted.values,plot=TRUE, legacy.axes=TRUE, xlab="False Positive Percentage", ylab="True Postive Percentage", col="#377eb8", lwd=4, percent=TRUE)
roc(data$Survivorship,logistic$fitted.values,plot=TRUE, legacy.axes=TRUE, xlab="False Positive Percentage", ylab="True Postive Percentage", col="#377eb8", lwd=4, percent=TRUE, print.auc=TRUE)
roc(data$Survivorship,logistic$fitted.values,plot=TRUE, legacy.axes=TRUE, xlab="False Positive Percentage", ylab="True Postive Percentage", col="#377eb8", lwd=4, percent=TRUE, print.auc=TRUE, partial.auc=c(100, 90), auc.polygon = TRUE, auc.polygon.col = "#377eb822", print.auc.x=45)
# Lets do two roc plots to understand which model is better
roc(data$Survivorship, logistic_simple$fitted.values, plot=TRUE, legacy.axes=TRUE, percent=TRUE, xlab="False Positive Percentage", ylab="True Postive Percentage", col="#377eb8", lwd=4, print.auc=TRUE)
# Lets add the other graph
plot.roc(data$Survivorship, logistic$fitted.values, percent=TRUE, col="#4daf4a", lwd=4, print.auc=TRUE, add=TRUE, print.auc.y=40)
legend("bottomright", legend=c("Simple", "Non Simple"), col=c("#377eb8", "#4daf4a"), lwd=4) # Make it user friendly


















































