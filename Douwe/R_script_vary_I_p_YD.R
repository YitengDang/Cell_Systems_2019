# Here logistic models will be fitted to the dataset in which both the initial I and p were varied.

# Load the csv file produced in Matlab
setwd("K:\\bn\\hy\\Shared\\Douwe\\Recent Results\\Vary_p_and_I") # Change to folder containing csv file
d = read.csv("Vary_I_p.csv")

# Check if columns in data are of correct datatype:
str(d)
# p1,p2,I1,I2 should be numeric (num). TW should be a factor. path does not matter, will not be used. If nessasary,
# convert the data column. Likely, the only change that is required is changing the TW column
d$TW = as.factor(d$TW)
# Check if convertion worked:
str(d)
#-------------------------------------

# Fitting logistic models
M0 = glm(TW ~ 1, family = 'binomial', data = d) # null model, mean fraction wave
M0p1 = glm(TW ~ p1, family = 'binomial', data = d)
M0p1p2 = glm(TW ~ p1 + p2, family = 'binomial' ,data = d)
M1 = glm(TW ~ p1 + p2 + p1:p2, family = 'binomial', data = d) # Check if using p improves prediction. These terms
# should be significant, based on our results. Let's check:
summary(M1)
anova(M0,M0p1,M0p1p2,M1,test = 'Chisq')
# 
# Create models indluding I to compare with the M1
M2 = glm(TW ~ p1 + p2 + p1:p2 + I1 + I2, family = 'binomial',data = d)
summary(M2)
anova(M0,M1,M2,test = 'Chisq')
# 
M3 = glm(TW ~ p1 + p2 + p1:p2 + I1 + I2 + I1:I2, family = 'binomial',data = d)
summary(M3)
anova(M0,M1,M2,M3,test = 'Chisq')

M4 = glm(TW ~ p1 + p2 + p1:p2 + I1 + I2 + I1:I2 + I1:I2:p1, family = 'binomial',data = d)

M5 = glm(TW ~ p1 + p2 + p1:p2 + I1 + I2 + I1:I2 + I1:I2:p1:p2, family = 'binomial',data = d)
summary(M5)
# 
Mfull = glm(TW ~ p1*p2*I1*I2, family = 'binomial',data = d)
summary(Mfull)
anova(M3,M4,Mfull,test = 'Chisq') # Check if 3 and 4 way interactions are better
anova(M5,Mfull,test = "Chisq")

# Note about notation ':' = interaction, '*' = interaction and include the terms in model separately.
# For example, following two expressions equivalent:
# Y ~ X + Z + W + X:Z + X:W + Z:W + X:Z:W
# Y ~ X*Z*W

# Compare a few more models with each other
M0 = glm(TW ~ 1, family = 'binomial', data = d) # null model, mean fraction wave
M1I = glm(TW ~ I1*I2, family = 'binomial',data = d)
M1p = glm(TW ~ p1*p2, family = 'binomial',data = d)
M1pI = glm(TW ~ p1*p2*I1*I2, family = 'binomial',data = d)
anova(M0, M1p, M1I, M1pI, test = "Chisq")
#-----------------------------------
# Select the best model above to serve as classifier. The classifier predics whether a simulation will
# result in a TW based on the terms in de model. Here we will predict this for the trainingset. However,
# for rigorous testing the model should be tested on a independent test set. This could for instance be done
# by keeping a part of the original data seperate (dont use in model fitting, but only for testing).

TW_prob <- predict(Mfull, type = 'response') # Change model if necessary
Prediction <- rep("No wave",length(d$p1))
Prediction[TW_prob > 0.5] <- "wave" # if TW prob > 0.5, predict TW

confusion_matrix = table(Prediction,d$TW) # Rows are predictions. Columns the reality

# Calculating performance metrics
(confusion_matrix[1,1] + confusion_matrix[2,2]) / sum(confusion_matrix) # accuracy
(confusion_matrix[1,1])/(confusion_matrix[1,1] + confusion_matrix[2,1]) # specificity
(confusion_matrix[2,2])/(confusion_matrix[2,2] + confusion_matrix[1,2]) # recall
(confusion_matrix[2,2])/(confusion_matrix[2,2] + confusion_matrix[2,1]) #precision