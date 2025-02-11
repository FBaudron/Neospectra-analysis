#' ---
#' title: "PLS analysis Neospectra - soils Zimbabwe - Double Burden project"
#' author: "Frédéric Baudron"
#' date: "February 10th, 2025"
#' ---


# CLEAR ENVIRONMENT-------------------------------------------------------------

rm(list = ls())


# LOADING NECESSARY PACKAGES----------------------------------------------------

suppressPackageStartupMessages({
  require("pls")
  require("openxlsx")
  require("caret")
  require("chillR")
})


# SETTING UP THE DIRECTORY, LOADING THE DATA & PREPARING THE DATA---------------

setwd('D:\\Mes Donnees\\1. Cirad\\GAIA\\Neospectra\\')

data = read.xlsx("Data\\soils Zimbabwe - Double Burden project.xlsx")

# means every 5 rows (reps for the same samples)
data_mean = sapply(data[, c(4:36, 41:297)], function(x) colMeans(matrix(x, nrow = 5)))

data_mean = as.data.frame(data_mean)

# log-transform reflectance (to sort out issues of distribution; obtain something alike to absorbance?)
data_mean[, 34:290] = log(data_mean[, 34:290])


# create a list of 70% of the rows in the original dataset we can use for training
validation_index = createDataPartition(data_mean$pH, p = 0.7, list = FALSE)
# select 30% of the data for validation
valid = data_mean[-validation_index,]
# use the remaining 70% of data to training and testing the models
train = data_mean[validation_index,]


# MODEL PH----------------------------------------------------------------------

# make this example reproducible
set.seed(1)

# fit PCR model
model_ph = plsr(pH ~ ., ncomp = 100, data = train[, -c(2:33)], scale = TRUE, validation = "CV")

# view summary of model fitting
# summary(model_ph)

# choose the number of PLS components

# validationplot(model_ph, val.type = "MSEP")

msepcv.pls_ph = MSEP(model_ph, estimate = c("train","CV"))
plot(msepcv.pls_ph, lty = 1, type = "l", legendpos = "topright", main="")

ncomp.pls_ph = which.min(msepcv.pls_ph$val["CV",,])-1
ncomp.pls_ph

model_ph = plsr(pH ~ ., ncomp = ncomp.pls_ph, data = train[, -c(2:33)], scale = TRUE)

res_model_ph = residuals(model_ph)
plot(res_model_ph[,,ncomp.pls_ph], pch = 15, cex = .5 , ylab = "Residus", main = "")
abline(h = c(-2,0,2), lty = c(2,1,2))

# performance of the calibration
pred_ph_train = predict(model_ph, train, ncomp = ncomp.pls_ph)

# calculate R2 and RMSEP
R2_ph_train = R2(pred_ph_train, train$pH, na.rm = TRUE)
RMSEP_ph_train = RMSEP(pred_ph_train, train$pH, na.rm = TRUE)

png("Results\\Zimbabwe\\PLS\\1.1. pH calibration.png", units="in", width = 8, height = 6.5, res=1000)
plot(x = train$pH, y = pred_ph_train, xlab = "Observed", ylab = "Predicted", main = "pH calibration")
abline(0, 1, col = "red")
fit_train_ph = lm(pred_ph_train ~ train$pH)
sfit_train_ph = summary(fit_train_ph)
# fit_train_ph
abline(a = fit_train_ph$coefficients[1], b = fit_train_ph$coefficients[2], col = "black")
mtext(paste0("R2 = ", round(R2_ph_train, 2),"; RMSEP = ", round(RMSEP_ph_train, 3)), side = 3)
dev.off()


# use model to make predictions on a validation dataset

pred_ph_val = predict(model_ph, valid, ncomp = ncomp.pls_ph)

# calculate RMSEP and R2
RMSEP_ph_val = RMSEP(pred_ph_val, valid$pH, na.rm = TRUE)
R2_ph_val = R2(pred_ph_val, valid$pH, na.rm = TRUE)

png("Results\\Zimbabwe\\PLS\\1.2. pH validation.png", units="in", width = 8, height = 6.5, res=1000)
plot(x = valid$pH, y = pred_ph_val, xlab = "Observed", ylab = "Predicted", main = "pH validation")
abline(0, 1, col = "red")
fit_val_ph = lm(pred_ph_val ~ valid$pH)
sfit_val_ph = summary(fit_val_ph)
# fit_val_ph
abline(a = fit_val_ph$coefficients[1], b = fit_val_ph$coefficients[2], col = "black")
mtext(paste0("R2 = ", round(R2_ph_val, 2),"; RMSEP = ", round(RMSEP_ph_val, 3)), side = 3)
dev.off()



# MODEL OLSEN P-----------------------------------------------------------------

# make this example reproducible
set.seed(1)

# fit PCR model
model_po = plsr(OlsenP ~ ., ncomp = 100, data = train[, -c(1:3, 5:33)], scale = TRUE, validation = "CV")

# view summary of model fitting
# summary(model_po)

# choose the number of PLS components

# validationplot(model_po, val.type = "MSEP")

msepcv.pls_po = MSEP(model_po, estimate = c("train","CV"))
plot(msepcv.pls_po, lty = 1, type = "l", legendpos = "topright", main="")

ncomp.pls_po = which.min(msepcv.pls_po$val["CV",,])-1
ncomp.pls_po

model_po = plsr(OlsenP ~ ., ncomp = ncomp.pls_po, data = train[, -c(1:3, 5:33)], scale = TRUE)

res_model_po = residuals(model_po)
plot(res_model_po[,,ncomp.pls_po], pch = 15, cex = .5 , ylab = "Residus", main = "")
abline(h = c(-2,0,2), lty = c(2,1,2))

# performance of the calibration
pred_po_train = predict(model_po, train, ncomp = ncomp.pls_po)

# calculate R2 and RMSEP
R2_po_train = R2(pred_po_train, train$OlsenP, na.rm = TRUE)
RMSEP_po_train = RMSEP(pred_po_train, train$OlsenP, na.rm = TRUE)

png("Results\\Zimbabwe\\PLS\\2.1. P calibration.png", units="in", width = 8, height = 6.5, res=1000)
plot(x = train$OlsenP, y = pred_po_train, xlab = "Observed", ylab = "Predicted", main = "P Olsen calibration")
abline(0, 1, col = "red")
fit_train_po = lm(pred_po_train ~ train$OlsenP)
sfit_train_po = summary(fit_train_po)
# sfit_train_po
abline(a = sfit_train_po$coefficients[1], b = sfit_train_po$coefficients[2], col = "black")
mtext(paste0("R2 = ", round(R2_po_train, 2),"; RMSEP = ", round(RMSEP_po_train, 3)), side = 3)
dev.off()


# use model to make predictions on a validation dataset

pred_po_val = predict(model_po, valid, ncomp = ncomp.pls_po)

# calculate RMSEP and R2
RMSEP_po_val = RMSEP(pred_po_val, valid$OlsenP, na.rm = TRUE)
R2_po_val = R2(pred_po_val, valid$OlsenP, na.rm = TRUE)

png("Results\\Zimbabwe\\PLS\\2.1. P validationation.png", units="in", width = 8, height = 6.5, res=1000)
plot(x = valid$OlsenP, y = pred_po_val, xlab = "Observed", ylab = "Predicted", main = "P Olsen validation")
abline(0, 1, col = "red")
fit_val_po = lm(pred_po_val ~ valid$OlsenP)
sfit_val_po = summary(fit_val_po)
# sfit_val_po
abline(a = sfit_val_po$coefficients[1], b = sfit_val_po$coefficients[2], col = "black")
mtext(paste0("R2 = ", round(R2_po_val, 2),"; RMSEP = ", round(RMSEP_po_val, 3)), side = 3)
dev.off()



# MODEL Al----------------------------------------------------------------------

# make this example reproducible
set.seed(1)

# fit PCR model
model_al = plsr(Al ~ ., ncomp = 100, data = train[, -c(1:13, 15:33)], scale = TRUE, validation = "CV")

# view summary of model fitting
# summary(model_al)

# choose the number of PLS components

# validationplot(model_al, val.type = "MSEP")

msepcv.pls_al = MSEP(model_al, estimate = c("train","CV"))
plot(msepcv.pls_al, lty = 1, type = "l", legendpos = "topright", main="")

ncomp.pls_al = which.min(msepcv.pls_al$val["CV",,])-1
ncomp.pls_al

model_al = plsr(Al ~ ., ncomp = ncomp.pls_al, data = train[, -c(1:13, 15:33)], scale = TRUE)

res_model_al = residuals(model_al)
plot(res_model_al[,,ncomp.pls_al], pch = 15, cex = .5 , ylab = "Residus", main = "")
abline(h = c(-2,0,2), lty = c(2,1,2))

# performance of the calibration
pred_al_train = predict(model_al, train, ncomp = ncomp.pls_al)

# calculate R2 and RMSEP
R2_al_train = R2(pred_al_train, train$Al, na.rm = TRUE)
RMSEP_al_train = RMSEP(pred_al_train, train$Al, na.rm = TRUE)

png("Results\\Zimbabwe\\PLS\\3.1. Al calibration.png", units="in", width = 8, height = 6.5, res=1000)
plot(x = train$Al, y = pred_al_train, xlab = "Observed", ylab = "Predicted", main = "Al calibration")
abline(0, 1, col = "red")
fit_train_al = lm(pred_al_train ~ train$Al)
sfit_train_al = summary(fit_train_al)
# sfit_train_al
abline(a = sfit_train_al$coefficients[1], b = sfit_train_al$coefficients[2], col = "black")
mtext(paste0("R2 = ", round(R2_al_train, 2),"; RMSEP = ", round(RMSEP_al_train, 3)), side = 3)
dev.off()


# use model to make predictions on a validation dataset

pred_al_val = predict(model_al, valid, ncomp = ncomp.pls_al)

# calculate RMSEP and R2
RMSEP_al_val = RMSEP(pred_al_val, valid$Al, na.rm = TRUE)
R2_al_val = R2(pred_al_val, valid$Al, na.rm = TRUE)

png("Results\\Zimbabwe\\PLS\\3.1. Al validation.png", units="in", width = 8, height = 6.5, res=1000)
plot(x = valid$Al, y = pred_al_val, xlab = "Observed", ylab = "Predicted", main = "Al validation")
abline(0, 1, col = "red")
fit_val_al = lm(pred_al_val ~ valid$Al)
sfit_val_al = summary(fit_val_al)
# sfit_val_al
abline(a = sfit_val_al$coefficients[1], b = sfit_val_al$coefficients[2], col = "black")
mtext(paste0("R2 = ", round(R2_al_val, 2),"; RMSEP = ", round(RMSEP_al_val, 3)), side = 3)
dev.off()



# MODEL Ca----------------------------------------------------------------------

# make this example reproducible
set.seed(1)

# fit PCR model
model_ca = plsr(Ca.ex ~ ., ncomp = 100, data = train[, -c(1:4, 6:33)], scale = TRUE, validation = "CV")

# view summary of model fitting
# summary(model_ca)

# choose the number of PLS components

# validationplot(model_ca, val.type = "MSEP")

msepcv.pls_ca = MSEP(model_ca, estimate = c("train","CV"))
plot(msepcv.pls_ca, lty = 1, type = "l", legendpos = "topright", main="")

ncomp.pls_ca = which.min(msepcv.pls_ca$val["CV",,])-1
ncomp.pls_ca

model_ca = plsr(Ca.ex ~ ., ncomp = ncomp.pls_ca, data = train[, -c(1:4, 6:33)], scale = TRUE)

res_model_ca = residuals(model_ca)
plot(res_model_ca[,,ncomp.pls_ca], pch = 15, cex = .5 , ylab = "Residus", main = "")
abline(h = c(-2,0,2), lty = c(2,1,2))

# performance of the calibration
pred_ca_train = predict(model_ca, train, ncomp = ncomp.pls_ca)

# calculate R2 and RMSEP
R2_ca_train = R2(pred_ca_train, train$Ca.ex, na.rm = TRUE)
RMSEP_ca_train = RMSEP(pred_ca_train, train$Ca.ex, na.rm = TRUE)

png("Results\\Zimbabwe\\PLS\\4.1. Ca calibration.png", units="in", width = 8, height = 6.5, res=1000)
plot(x = train$Ca.ex, y = pred_ca_train, xlab = "Observed", ylab = "Predicted", main = "Ca calibration")
abline(0, 1, col = "red")
fit_train_ca = lm(pred_ca_train ~ train$Ca.ex)
sfit_train_ca = summary(fit_train_ca)
# sfit_train_ca
abline(a = sfit_train_ca$coefficients[1], b = sfit_train_ca$coefficients[2], col = "black")
mtext(paste0("R2 = ", round(R2_al_train, 2),"; RMSEP = ", round(RMSEP_al_train, 3)), side = 3)
dev.off()


# use model to make predictions on a validation dataset

pred_ca_val = predict(model_ca, valid, ncomp = ncomp.pls_ca)

# calculate RMSEP and R2
RMSEP_ca_val = RMSEP(pred_ca_val, valid$Ca.ex, na.rm = TRUE)
R2_ca_val = R2(pred_ca_val, valid$Ca.ex, na.rm = TRUE)

png("Results\\Zimbabwe\\PLS\\4.1. Ca validation.png", units="in", width = 8, height = 6.5, res=1000)
plot(x = valid$Ca.ex, y = pred_ca_val, xlab = "Observed", ylab = "Predicted", main = "Ca validation")
abline(0, 1, col = "red")
fit_val_ca = lm(pred_ca_val ~ valid$Ca.ex)
sfit_val_ca = summary(fit_val_ca)
# sfit_val_ca
abline(a = sfit_val_ca$coefficients[1], b = sfit_val_ca$coefficients[2], col = "black")
mtext(paste0("R2 = ", round(R2_ca_val, 2),"; RMSEP = ", round(RMSEP_ca_val, 3)), side = 3)
dev.off()



# MODEL Mg----------------------------------------------------------------------

# make this example reproducible
set.seed(1)

# fit PCR model
model_mg = plsr(Mg.ex ~ ., ncomp = 100, data = train[, -c(1:6, 8:33)], scale = TRUE, validation = "CV")

# view summary of model fitting
# summary(model_mg)

# choose the number of PLS components

# validationplot(model_mg, val.type = "MSEP")

msepcv.pls_mg = MSEP(model_mg, estimate = c("train","CV"))
plot(msepcv.pls_mg, lty = 1, type = "l", legendpos = "topright", main="")

ncomp.pls_mg = which.min(msepcv.pls_mg$val["CV",,])-1
ncomp.pls_mg

model_mg = plsr(Mg.ex ~ ., ncomp = ncomp.pls_mg, data = train[, -c(1:6, 8:33)], scale = TRUE)

res_model_mg = residuals(model_mg)
plot(res_model_mg[,,ncomp.pls_mg], pch = 15, cex = .5 , ylab = "Residus", main = "")
abline(h = c(-2,0,2), lty = c(2,1,2))

# performance of the calibration
pred_mg_train = predict(model_mg, train, ncomp = ncomp.pls_mg)

# calculate R2 and RMSEP
R2_mg_train = R2(pred_mg_train, train$Mg.ex, na.rm = TRUE)
RMSEP_mg_train = RMSEP(pred_mg_train, train$Mg.ex, na.rm = TRUE)

png("Results\\Zimbabwe\\PLS\\5.1. Mg calibration.png", units="in", width = 8, height = 6.5, res=1000)
plot(x = train$Mg.ex, y = pred_mg_train, xlab = "Observed", ylab = "Predicted", main = "Mg calibration")
abline(0, 1, col = "red")
fit_train_mg = lm(pred_mg_train ~ train$Mg.ex)
sfit_train_mg = summary(fit_train_mg)
# fit_train_mg
abline(a = fit_train_mg$coefficients[1], b = fit_train_mg$coefficients[2], col = "black")
mtext(paste0("R2 = ", round(R2_mg_train, 2),"; RMSEP = ", round(RMSEP_mg_train, 3)), side = 3)
dev.off()


# use model to make predictions on a validation dataset

pred_mg_val = predict(model_mg, valid, ncomp = ncomp.pls_mg)

# calculate RMSEP and R2
RMSEP_mg_val = RMSEP(pred_mg_val, valid$Mg.ex, na.rm = TRUE)
R2_mg_val = R2(pred_mg_val, valid$Mg.ex, na.rm = TRUE)

png("Results\\Zimbabwe\\PLS\\5.2. Mg validation.png", units="in", width = 8, height = 6.5, res=1000)
plot(x = valid$Mg.ex, y = pred_mg_val, xlab = "Observed", ylab = "Predicted", main = "Mg validation")
abline(0, 1, col = "red")
fit_val_mg = lm(pred_mg_val ~ valid$Mg.ex)
sfit_val_mg = summary(fit_val_mg)
# sfit_val_mg
abline(a = sfit_val_mg$coefficients[1], b = sfit_val_mg$coefficients[2], col = "black")
mtext(paste0("R2 = ", round(R2_mg_val, 2),"; RMSEP = ", round(RMSEP_mg_val, 3)), side = 3)
dev.off()



# MODEL CEC---------------------------------------------------------------------

# make this example reproducible
set.seed(1)

# fit PCR model
model_cec = plsr(eCEC.clean ~ ., ncomp = 100, data = train[, -c(1:8, 10:33)], scale = TRUE, validation = "CV")

# view summary of model fitting
# summary(model_cec)

# choose the number of PLS components

# validationplot(model_cec, val.type = "MSEP")

msepcv.pls_cec = MSEP(model_cec, estimate = c("train","CV"))
plot(msepcv.pls_cec, lty = 1, type = "l", legendpos = "topright", main="")

ncomp.pls_cec = which.min(msepcv.pls_cec$val["CV",,])-1
ncomp.pls_cec

model_cec = plsr(eCEC.clean ~ ., ncomp = ncomp.pls_cec, data = train[, -c(1:8, 10:33)], scale = TRUE)

res_model_cec = residuals(model_cec)
plot(res_model_cec[,,ncomp.pls_cec], pch = 15, cex = .5 , ylab = "Residus", main = "")
abline(h = c(-2,0,2), lty = c(2,1,2))

# performance of the calibration
pred_cec_train = predict(model_cec, train, ncomp = ncomp.pls_cec)

# calculate R2 and RMSEP
R2_cec_train = R2(pred_cec_train, train$eCEC.clean, na.rm = TRUE)
RMSEP_cec_train = RMSEP(pred_cec_train, train$eCEC.clean, na.rm = TRUE)

png("Results\\Zimbabwe\\PLS\\6.1. CEC calibration.png", units="in", width = 8, height = 6.5, res=1000)
plot(x = train$eCEC.clean, y = pred_cec_train, xlab = "Observed", ylab = "Predicted", main = "eCEC calibration")
abline(0, 1, col = "red")
fit_train_cec = lm(pred_cec_train ~ train$eCEC.clean)
sfit_train_cec = summary(fit_train_cec)
# sfit_train_cec
abline(a = sfit_train_cec$coefficients[1], b = sfit_train_cec$coefficients[2], col = "black")
mtext(paste0("R2 = ", round(R2_cec_train, 2),"; RMSEP = ", round(RMSEP_cec_train, 3)), side = 3)
dev.off()


# use model to make predictions on a validation dataset

pred_cec_val = predict(model_cec, valid, ncomp = ncomp.pls_cec)

# calculate RMSEP and R2
RMSEP_cec_val = RMSEP(pred_cec_val, valid$eCEC.clean, na.rm = TRUE)
R2_cec_val = R2(pred_cec_val, valid$eCEC.clean, na.rm = TRUE)

png("Results\\Zimbabwe\\PLS\\6.2. CEC validation.png", units="in", width = 8, height = 6.5, res=1000)
plot(x = valid$eCEC.clean, y = pred_cec_val, xlab = "Observed", ylab = "Predicted", main = "eCEC validation")
abline(0, 1, col = "red")
fit_val_cec = lm(pred_cec_val ~ valid$eCEC.clean)
sfit_val_cec = summary(fit_val_cec)
# sfit_val_cec
abline(a = sfit_val_cec$coefficients[1], b = sfit_val_cec$coefficients[2], col = "black")
mtext(paste0("R2 = ", round(R2_cec_val, 2),"; RMSEP = ", round(RMSEP_cec_val, 3)), side = 3)
dev.off()



# MODEL SOC---------------------------------------------------------------------

# make this example reproducible
set.seed(1)

# fit PCR model
model_soc = plsr(Total.soil.C ~ ., ncomp = 100, data = train[, -c(1:2, 4:33)], scale = TRUE, validation = "CV")

# view summary of model fitting
# summary(model_soc)

# choose the number of PLS components

# validationplot(model_soc, val.type = "MSEP")

msepcv.pls_soc = MSEP(model_soc, estimate = c("train","CV"))
plot(msepcv.pls_soc, lty = 1, type = "l", legendpos = "topright", main="")

ncomp.pls_soc = which.min(msepcv.pls_soc$val["CV",,])-1
ncomp.pls_soc

model_soc = plsr(Total.soil.C ~ ., ncomp = ncomp.pls_soc, data = train[, -c(1:2, 4:33)], scale = TRUE)

res_model_soc = residuals(model_soc)
plot(res_model_soc[,,ncomp.pls_soc], pch = 15, cex = .5 , ylab = "Residus", main = "")
abline(h = c(-2,0,2), lty = c(2,1,2))

# performance of the calibration
pred_soc_train = predict(model_soc, train, ncomp = ncomp.pls_soc)

# calculate R2 and RMSEP
R2_soc_train = R2(pred_soc_train, train$Total.soil.C, na.rm = TRUE)
RMSEP_soc_train = RMSEP(pred_soc_train, train$Total.soil.C, na.rm = TRUE)

png("Results\\Zimbabwe\\PLS\\7.1. SOC calibration.png", units="in", width = 8, height = 6.5, res=1000)
plot(x = train$Total.soil.C, y = pred_soc_train, xlab = "Observed", ylab = "Predicted", main = "SOC calibration")
abline(0, 1, col = "red")
fit_train_soc = lm(pred_soc_train ~ train$Total.soil.C)
sfit_train_soc = summary(fit_train_soc)
# sfit_train_soc
abline(a = sfit_train_soc$coefficients[1], b = sfit_train_soc$coefficients[2], col = "black")
mtext(paste0("R2 = ", round(R2_soc_train, 2),"; RMSEP = ", round(RMSEP_soc_train, 3)), side = 3)
dev.off()


# use model to make predictions on a validation dataset

pred_soc_val = predict(model_soc, valid, ncomp = ncomp.pls_soc)

# calculate RMSEP and R2
RMSEP_soc_val = RMSEP(pred_soc_val, valid$Total.soil.C, na.rm = TRUE)
R2_soc_val = R2(pred_soc_val, valid$Total.soil.C, na.rm = TRUE)

png("Results\\Zimbabwe\\PLS\\7.2. SOC validation.png", units="in", width = 8, height = 6.5, res=1000)
plot(x = valid$Total.soil.C, y = pred_soc_val, xlab = "Observed", ylab = "Predicted", main = "SOC validation")
abline(0, 1, col = "red")
fit_val_soc = lm(pred_soc_val ~ valid$Total.soil.C)
sfit_val_soc = summary(fit_val_soc)
# sfit_val_soc
abline(a = sfit_val_soc$coefficients[1], b = sfit_val_soc$coefficients[2], col = "black")
mtext(paste0("R2 = ", round(R2_soc_val, 2),"; RMSEP = ", round(RMSEP_soc_val, 3)), side = 3)
dev.off()





