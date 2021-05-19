
library(glmnet)
library(reshape)
library(ggplot2)
library(gbm)
library(pROC)
library(ROCR)
library(grid)
library(gridBase)
library(viridis)


load('~/IvanaJanickova/Cortex(2).rdata')
set.seed(1)
train = sample(1:nrow(cortex), nrow(cortex)*2/3)
test = (-train)
x = model.matrix(cortex$Behavior~., data = cortex[, 1:70])
x.train = x[train,]
x.test = x[test,]
y = cortex$Behavior
y = ifelse(y == "C/S", 1, 0)
y.train = y[train]
y.test = y[test]

correlationMatrix=(cor(cortex[,1:70]))
melted = melt(correlationMatrix)
#How many genes are more that 90% correlated? & How many genes are more than 50% correlated?
cor_0.9 = 0
cor_0.5 = 0
names.cor_0.9 = c()

for(i in 1:nrow(melted)) {
  if (melted[i,1] != melted[i,2]) {
    if (melted[i,3] >= 0.9) {
      cor_0.9 = cor_0.9 + 1
      names.cor_0.9 = c(names.cor_0.9, as.character(melted[i,1]))
    } 
    if (melted[i,3] >= 0.5) {
      cor_0.5 = cor_0.5 + 1
    }
  }
}
names.cor_0.9 = unique(names.cor_0.9) # stores names of predictors that are more that 90% correlated
ordered = melted[order(melted[,3], decreasing = TRUE),] # make orderd matrix for futher considerations
perc_0.9 = cor_0.9*100/(nrow(melted)-70)
perc_0.5 = cor_0.5*100/(nrow(melted)-70)
#perc_0.9 # 0.62%
#perc_0.5 # 16.77%

#ggplot(data = melted, aes(x=X1, y=X2, fill=value)) +  geom_tile() + ggtitle("Heatmap: All proteins expression levels") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

#heatmap visualization of correlating matrix for 10 genes
correlationMatrix2=(cor(cortex[,1:10]))
melted2 = melt(correlationMatrix2)
ggplot(data = melted2, aes(x=X1, y=X2, fill=value)) +  geom_tile() + ggtitle("Heatmap: 10 proteins expression levels") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


par(mfrow = c(1,2))

#trining the model
ridge.model <- glmnet(x.train,y.train, alpha=0, lambda=NULL, family = "binomial")
plot(ridge.model)


#Choosing optimal value for lambda
set.seed(1)
ridge.cv = cv.glmnet(x.train, y.train, alpha=0, family = "binomial", nfolds = length(train)) #lenght LOOCV
plot(ridge.cv)
best.lambda = ridge.cv$lambda.min

#Ridge model with best lambda selected by CV
model.ridge = glmnet(x.train, y.train, alpha=0, family = "binomial", lambda = best.lambda)
#Displaying the coefficietns
coef.ridge = as.matrix(coef(model.ridge))

#Displaying the coefficietns
coef.matrix = as.matrix(coef(model.ridge))
#getting non-zero coefficients
row_sub = apply(coef.matrix, 1, function(row) all(row !=0 ))
coef.matrix = sort(coef.matrix[row_sub,])
par(mfrow = c(1,1))
midpts = barplot(coef.matrix, names.arg = "", col = viridis(10))
title("Ridge coefficients for different proteins")
## Use grid to add the labels    
vps <- baseViewports()
pushViewport(vps$inner, vps$figure, vps$plot)

grid.text(names(coef.matrix),
          x = unit(midpts, "native"), y=unit(-1, "lines"),
          just="right", rot=50)

popViewport(3)

ridge.pred=predict(ridge.model ,s=best.lambda ,newx=x.test,type="response")
predictions=rep(0 ,length(y.test))
predictions[ridge.pred>0.5]= 1
table(y.test,predictions)
performanceRidge=length(which(predictions==y.test))/length(y.test)
performanceRidge #11+9)/24 = 83%

vals=predict(model.ridge,s= best.lambda,type="coefficients")
selected=colnames(x)[vals@i]



par(mfrow = c(1,2))

#trining the model
lasso.model <- glmnet(x.train,y.train, alpha=1, lambda=NULL, family = "binomial")
plot(lasso.model)


#Choosing optimal value for lambda
set.seed(1)
lasso.cv = cv.glmnet(x.train, y.train, alpha=1, family = "binomial", nfolds = length(train)) #lenght LOOCV
plot(lasso.cv)
best.lambda = lasso.cv$lambda.min

#Ridge model with best lambda selected by CV
model.lasso = glmnet(x.train, y.train, alpha=1, family = "binomial", lambda = best.lambda)


#Displaying the coefficietns
coef.matrix = as.matrix(coef(model.lasso))
#getting non-zero coefficients
row_sub = apply(coef.matrix, 1, function(row) all(row !=0 ))
coef.matrix = sort(coef.matrix[row_sub,])
par(mfrow = c(1,1))
midpts = barplot(coef.matrix, names.arg = "", col = viridis(10))
title("Lasso coefficients for different proteins")
## Use grid to add the labels    
vps <- baseViewports()
pushViewport(vps$inner, vps$figure, vps$plot)

grid.text(names(coef.matrix),
          x = unit(midpts, "native"), y=unit(-1, "lines"),
          just="right", rot=50)

popViewport(3)


lasso.pred=predict(lasso.model ,s=best.lambda ,newx=x.test,type="response")
predictions=rep(0 ,length(y.test))
predictions[lasso.pred>0.5]= 1
table(y.test,predictions)
performanceLasso=length(which(predictions==y.test))/length(y.test)
performanceLasso #14+9)/24 = 95,8%

vals=predict(model.lasso,s= best.lambda,type="coefficients")
selected=colnames(x)[vals@i]
selected

set.seed(1)
cortex$dummy = ifelse(cortex$Behavior == "C/S", 1, 0)
model.boost = gbm(cortex$dummy ~ ., data = cortex[train],  distribution = "bernoulli", 
                  n.trees = 3000, interaction.depth = 3)
summary(model.boost)


plot(model.boost, i = "DYRK1A_N")
plot(model.boost, i = "pS6_N")

boost.pred = predict(model.boost, cortex[test, ], n.trees = 3000)
predictions=rep(0 ,length(cortex[test, "dummy"]))
predictions[boost.pred>0.5]= 1
table(cortex[test, "dummy"],predictions)
performanceBoost=length(which(predictions==cortex[test, "dummy"]))/length(cortex[test, "dummy"])
performanceBoost # 100%

treatement.indexes = cortex$Treatment == "Memantine"
differentially_expressed = c()

for (i in 1:70){
  
  treated = cortex[treatement.indexes, i]
  untreated = cortex[!treatement.indexes, i]
  test = wilcox.test(treated, untreated)
  if(test$p.value < 0.05) {
    differentially_expressed = c(differentially_expressed, names(cortex)[i])
  }
}


par(mfrow = c(3,3))
for (i in differentially_expressed) {
  boxplot(cortex[, i] ~ cortex$Treatment, ylab = i, xlab = "treatement", col = viridis(4))
}


get_differentally_expressed <- function(df, category){
  
  genotype.indexes = cortex$Genotype == "Ts65Dn"
  behavior.indexes = cortex$Behavior == "C/S"
  differentially_expressed = c()
  treatement.indexes = cortex$Treatment == "Memantine"
  indexes= rep(1,70)
  for (i in 1:70){
    if (category == 1) {indexes[i] = ifelse(genotype.indexes[i] == T && behavior.indexes[i] == T, T, F) }
    if (category == 2) {indexes[i] = ifelse(genotype.indexes[i] == T && behavior.indexes[i] == F, T, F) }
    if (category == 3) {indexes[i] = ifelse(genotype.indexes[i] == F && behavior.indexes[i] == T, T, F) }
    if (category == 4) {indexes[i] = ifelse(genotype.indexes[i] == F && behavior.indexes[i] == F, T, F) }
  }
  
  data = df[indexes, ]
  
  
  for (i in 1:70){
    treated = data[treatement.indexes, i]
    untreated = data[!treatement.indexes, i]
    test = wilcox.test(treated, untreated, paired = FALSE)
    if(!is.na(test$p.value) && test$p.value < 0.05) {
      differentially_expressed = c(differentially_expressed, names(data)[i])
    }
  }
  differentially_expressed
}

cat1 = get_differentally_expressed(cortex, 1)
cat2 = get_differentally_expressed(cortex, 2)
cat3 = get_differentally_expressed(cortex, 3)
cat4 = get_differentally_expressed(cortex, 4)


get_differentally_expressed <- function(df, category){
  
  genotype.indexes = cortex$Genotype == "Ts65Dn"
  behavior.indexes = cortex$Behavior == "C/S"
  differentially_expressed = c()
  treatement.indexes = cortex$Treatment == "Memantine"
  
  if (category == 1) {indexes = genotype.indexes}
  if (category == 2) {indexes = !genotype.indexes}
  if (category == 3) {indexes = behavior.indexes}
  if (category == 4) {indexes = !behavior.indexes}
  
  
  data = df[indexes, ]
  
  
  for (i in 1:70){
    treated = data[treatement.indexes, i]
    untreated = data[!treatement.indexes, i]
    test = wilcox.test(treated, untreated, paired = FALSE)
    if(!is.null(test) && test$p.value < 0.05) {
      differentially_expressed = c(differentially_expressed, names(data)[i])
    }
  }
  differentially_expressed
}

cat1 = get_differentally_expressed(cortex, 1)
cat2 = get_differentally_expressed(cortex, 2)
cat3 = get_differentally_expressed(cortex, 3)
cat4 = get_differentally_expressed(cortex, 4)

