# Factor Analysis
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



