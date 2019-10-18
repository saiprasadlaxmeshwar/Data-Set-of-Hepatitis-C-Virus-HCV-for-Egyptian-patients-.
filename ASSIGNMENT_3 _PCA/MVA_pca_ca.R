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










































