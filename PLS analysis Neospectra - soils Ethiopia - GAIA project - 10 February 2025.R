#' ---
#' title: "PLS analysis Neospectra - soils Ethiopia - GAIA project"
#' author: "Frédéric Baudron"
#' date: "February 10th, 2025"
#' ---


# CLEAR ENVIRONMENT-------------------------------------------------------------

rm(list = ls())


# LOADING NECESSARY PACKAGES----------------------------------------------------

suppressPackageStartupMessages({
  require("openxlsx")
  require("tidyr")
  require("dplyr")
  require("pls")
  require("caret")
  require("chillR")
})


# SETTING UP THE DIRECTORY, LOADING THE DATA & PREPARING THE DATA---------------

setwd('D:\\Mes Donnees\\1. Cirad\\GAIA\\Neospectra\\')

jim_codes = read.xlsx("Data\\soils Ethiopia - GAIA project.xlsx", sheet = 1)
deb_codes = read.xlsx("Data\\soils Ethiopia - GAIA project.xlsx", sheet = 3)

jim_neo = read.xlsx("Data\\soils Ethiopia - GAIA project.xlsx", sheet = 2)
deb_neo = read.xlsx("Data\\soils Ethiopia - GAIA project.xlsx", sheet = 4)

jim_neo = read.xlsx("Data\\soils Ethiopia - GAIA project.xlsx", sheet = 2)
deb_neo = read.xlsx("Data\\soils Ethiopia - GAIA project.xlsx", sheet = 4)

icraf = read.xlsx("Data\\soils Ethiopia - GAIA project.xlsx", sheet = 5)

cropnut = read.xlsx("Data\\soils Ethiopia - GAIA project.xlsx", sheet = 6)


jim_codes = separate(jim_codes, col = NeoSpectra.scanning.codes, into = c("Neospectra", "Neospectra2"), sep = "-", remove = TRUE)
jim_codes = jim_codes[, -c(8)]
jim_codes$Neospectra = sub('.', '', jim_codes$Neospectra)
jim_codes$Neospectra = as.numeric(jim_codes$Neospectra)
jim_codes1 = jim_codes
jim_codes2 = jim_codes
jim_codes2$Neospectra = jim_codes2$Neospectra + 1
jim_codes3 = jim_codes
jim_codes3$Neospectra = jim_codes3$Neospectra + 2
jim_codes4 = jim_codes
jim_codes4$Neospectra = jim_codes4$Neospectra + 3
jim_codes5 = jim_codes
jim_codes5$Neospectra = jim_codes5$Neospectra + 4
jim_codes = rbind(jim_codes1, jim_codes2, jim_codes3, jim_codes4, jim_codes5)
jim_codes$Neospectra = paste0("J_", jim_codes$Neospectra)
jim_codes = jim_codes[, c(7, 9:10)]
names(jim_codes)[2:3] = c("ICRAF", "Spectral")

deb_codes = separate(deb_codes, col = Neo.Spectra.Scanned.Codes, into = c("Neospectra", "Neospectra2"), sep = "-", remove = TRUE)
deb_codes = deb_codes[, -c(6)]
deb_codes$Neospectra = gsub("^.{0,2}", "", deb_codes$Neospectra)
deb_codes$Neospectra = as.numeric(deb_codes$Neospectra)
deb_codes1 = deb_codes
deb_codes2 = deb_codes
deb_codes2$Neospectra = deb_codes2$Neospectra + 1
deb_codes3 = deb_codes
deb_codes3$Neospectra = deb_codes3$Neospectra + 2
deb_codes4 = deb_codes
deb_codes4$Neospectra = deb_codes4$Neospectra + 3
deb_codes5 = deb_codes
deb_codes5$Neospectra = deb_codes5$Neospectra + 4
deb_codes = rbind(deb_codes1, deb_codes2, deb_codes3, deb_codes4, deb_codes5)
deb_codes$Neospectra = paste0("DM_", deb_codes$Neospectra)
deb_codes = deb_codes[, c(5:7)]
names(deb_codes)[2:3] = c("ICRAF", "Spectral")


jim_neo = jim_neo[, -c(1, 3)]
names(jim_neo)[1] = "Neospectra"

deb_neo = deb_neo[, -c(2)]
names(deb_neo)[1] = "Neospectra"


icraf = icraf[, -c(1:3, 5:9)]
names(icraf)[1] = "ICRAF"


codes = rbind(jim_codes, deb_codes)

neo = rbind(jim_neo, deb_neo)

data = merge(codes, icraf, by = "ICRAF")
data = merge(data, neo, by = "Neospectra")


# create in ID to order Neospectra scan so the 5 reps can be averaged

data = separate(data, col = Neospectra, into = c("Neospectra1", "Neospectra2"), sep = "_", remove = FALSE)
data$Neospectra2 = sprintf("%03d", as.numeric(data$Neospectra2))
data$ID = paste0(data$Neospectra1, "_", data$Neospectra2)

data = arrange(data, ID)

data = data[, -c(2:3)]
data = data[, c(281, 1:280)]

data[, c(5:281)] = sapply(data[, c(5:281)], as.numeric)


# means every 5 rows (reps for the same samples)
data_mean = sapply(data[, c(5:281)], function(x) colMeans(matrix(x, nrow = 5)))

data_mean = as.data.frame(data_mean)

# log-transform reflectance (to sort out issues of distribution; obtain something alike to absorbance?)
data_mean[, 23:277] = log(data_mean[, 23:277])


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
model_ph = plsr(pH ~ ., ncomp = 77, data = train[, -c(2:20)], scale = TRUE, validation = "CV")

# view summary of model fitting
# summary(model_ph)

# choose the number of PLS components

# validationplot(model_ph, val.type = "MSEP")

msepcv.pls_ph = MSEP(model_ph, estimate = c("train","CV"))
plot(msepcv.pls_ph, lty = 1, type = "l", legendpos = "topright", main="")

ncomp.pls_ph = which.min(msepcv.pls_ph$val["CV",,])-1
ncomp.pls_ph

model_ph = plsr(pH ~ ., ncomp = ncomp.pls_ph, data = train[, -c(2:20)], scale = TRUE)

res_model_ph = residuals(model_ph)
plot(res_model_ph[,,ncomp.pls_ph], pch = 15, cex = .5 , ylab = "Residus", main = "")
abline(h = c(-2,0,2), lty = c(2,1,2))

# performance of the calibration
pred_ph_train = predict(model_ph, train, ncomp = ncomp.pls_ph)

# calculate R2 and RMSEP
R2_ph_train = R2(pred_ph_train, train$pH, na.rm = TRUE)
RMSEP_ph_train = RMSEP(pred_ph_train, train$pH, na.rm = TRUE)

png("Results\\Ethiopia\\PLS\\1.1. pH calibration.png", units="in", width = 8, height = 6.5, res=1000)
plot(x = train$pH, y = pred_ph_train, xlab = "Observed", ylab = "Predicted", main = "pH calibration")
abline(0, 1, col = "red")
fit_train_ph = lm(pred_ph_train ~ train$pH)
sfit_train_ph = summary(fit_train_ph)
# fit_train_ph
abline(a = sfit_train_ph$coefficients[1], b = sfit_train_ph$coefficients[2], col = "black")
mtext(paste0("R2 = ", round(R2_ph_train, 2),"; RMSEP = ", round(RMSEP_ph_train, 3)), side = 3)


# use model to make predictions on a validation dataset

pred_ph_val = predict(model_ph, valid, ncomp = ncomp.pls_ph)

# calculate RMSEP and R2
RMSEP_ph_val = RMSEP(pred_ph_val, valid$pH, na.rm = TRUE)
R2_ph_val = R2(pred_ph_val, valid$pH, na.rm = TRUE)

png("Results\\Ethiopia\\RF\\1.2. pH validation.png", units="in", width = 8, height = 6.5, res=1000)
plot(x = valid$pH, y = pred_ph_val, xlab = "Observed", ylab = "Predicted", main = "pH validation")
abline(0, 1, col = "red")
fit_val_ph = lm(pred_ph_val ~ valid$pH)
sfit_val_ph = summary(fit_val_ph)
# sfit_val_ph
abline(a = sfit_val_ph$coefficients[1], b = sfit_val_ph$coefficients[2], col = "black")
mtext(paste0("R2 = ", round(R2_ph_val, 2),"; RMSEP = ", round(RMSEP_ph_val, 3)), side = 3)



# MODEL PSI---------------------------------------------------------------------

# make this example reproducible
set.seed(1)

# fit PCR model
model_po = plsr(PSI ~ ., ncomp = 77, data = train[, -c(1:14, 16:20)], scale = TRUE, validation = "CV")

# view summary of model fitting
# summary(model_po)

# choose the number of PLS components

# validationplot(model_po, val.type = "MSEP")

msepcv.pls_po = MSEP(model_po, estimate = c("train","CV"))
plot(msepcv.pls_po, lty = 1, type = "l", legendpos = "topright", main="")

ncomp.pls_po = which.min(msepcv.pls_po$val["CV",,])-1
ncomp.pls_po

model_po = plsr(PSI ~ ., ncomp = ncomp.pls_po, data = train[, -c(1:14, 16:20)], scale = TRUE)

res_model_po = residuals(model_po)
plot(res_model_po[,,ncomp.pls_po], pch = 15, cex = .5 , ylab = "Residus", main = "")
abline(h = c(-2,0,2), lty = c(2,1,2))

# performance of the calibration
pred_po_train = predict(model_po, train, ncomp = ncomp.pls_po)

# calculate R2 and RMSEP
R2_po_train = R2(pred_po_train, train$PSI, na.rm = TRUE)
RMSEP_po_train = RMSEP(pred_po_train, train$PSI, na.rm = TRUE)

png("Results\\Ethiopia\\RF\\2.1. PSI calibration.png", units="in", width = 8, height = 6.5, res=1000)
plot(x = train$PSI, y = pred_po_train, xlab = "Observed", ylab = "Predicted", main = "PSI calibration")
abline(0, 1, col = "red")
fit_train_po = lm(pred_po_train ~ train$PSI)
sfit_train_po = summary(fit_train_po)
# sfit_train_po
abline(a = sfit_train_po$coefficients[1], b = sfit_train_po$coefficients[2], col = "black")
mtext(paste0("R2 = ", round(R2_po_train, 2),"; RMSEP = ", round(RMSEP_po_train, 3)), side = 3)


# use model to make predictions on a validation dataset

pred_po_val = predict(model_po, valid, ncomp = ncomp.pls_po)

# calculate RMSEP and R2
RMSEP_po_val = RMSEP(pred_po_val, valid$PSI, na.rm = TRUE)
R2_po_val = R2(pred_po_val, valid$PSI, na.rm = TRUE)

png("Results\\Ethiopia\\RF\\2.2. PSI validation.png", units="in", width = 8, height = 6.5, res=1000)
plot(x = valid$PSI, y = pred_po_val, xlab = "Observed", ylab = "Predicted", main = "PSI validation")
abline(0, 1, col = "red")
fit_val_po = lm(pred_po_val ~ valid$PSI)
sfit_val_po = summary(fit_val_po)
# sfit_val_po
abline(a = sfit_val_po$coefficients[1], b = sfit_val_po$coefficients[2], col = "black")
mtext(paste0("R2 = ", round(R2_po_val, 2),"; RMSEP = ", round(RMSEP_po_val, 3)), side = 3)



# MODEL ExAc--------------------------------------------------------------------

# make this example reproducible
set.seed(1)

# fit PCR model
model_exac = plsr(ExAc ~ ., ncomp = 77, data = train[, -c(1:12, 14:20)], scale = TRUE, validation = "CV")

# view summary of model fitting
# summary(model_exac)

# choose the number of PLS components

# validationplot(model_exac, val.type = "MSEP")

msepcv.pls_exac = MSEP(model_exac, estimate = c("train","CV"))
plot(msepcv.pls_exac, lty = 1, type = "l", legendpos = "topright", main="")

ncomp.pls_exac = which.min(msepcv.pls_exac$val["CV",,])-1
ncomp.pls_exac

model_exac = plsr(ExAc ~ ., ncomp = ncomp.pls_exac, data = train[, -c(1:12, 14:20)], scale = TRUE)

res_model_exac = residuals(model_exac)
plot(res_model_exac[,,ncomp.pls_exac], pch = 15, cex = .5 , ylab = "Residus", main = "")
abline(h = c(-2,0,2), lty = c(2,1,2))

# performance of the calibration
pred_exac_train = predict(model_exac, train, ncomp = ncomp.pls_exac)

# calculate R2 and RMSEP
R2_exac_train = R2(pred_exac_train, train$ExAc, na.rm = TRUE)
RMSEP_exac_train = RMSEP(pred_exac_train, train$ExAc, na.rm = TRUE)

png("Results\\Ethiopia\\RF\\3.1. ExAc calibration.png", units="in", width = 8, height = 6.5, res=1000)
plot(x = train$ExAc, y = pred_exac_train, xlab = "Observed", ylab = "Predicted", main = "ExAc calibration")
abline(0, 1, col = "red")
fit_train_exac = lm(pred_exac_train ~ train$ExAc)
sfit_train_exac = summary(fit_train_exac)
# sfit_train_exac
abline(a = sfit_train_exac$coefficients[1], b = sfit_train_exac$coefficients[2], col = "black")
mtext(paste0("R2 = ", round(R2_exac_train, 2),"; RMSEP = ", round(RMSEP_exac_train, 3)), side = 3)


# use model to make predictions on a validation dataset

pred_exac_val = predict(model_exac, valid, ncomp = ncomp.pls_exac)

# calculate RMSEP and R2
RMSEP_exac_val = RMSEP(pred_exac_val, valid$ExAc, na.rm = TRUE)
R2_exac_val = R2(pred_exac_val, valid$ExAc, na.rm = TRUE)

png("Results\\Ethiopia\\RF\\3.2. ExAc validation.png", units="in", width = 8, height = 6.5, res=1000)
plot(x = valid$ExAc, y = pred_exac_val, xlab = "Observed", ylab = "Predicted", main = "ExAc validation")
abline(0, 1, col = "red")
fit_val_exac = lm(pred_exac_val ~ valid$ExAc)
sfit_val_exac = summary(fit_val_exac)
# sfit_val_exac
abline(a = sfit_val_exac$coefficients[1], b = sfit_val_exac$coefficients[2], col = "black")
mtext(paste0("R2 = ", round(R2_exac_val, 2),"; RMSEP = ", round(RMSEP_exac_val, 3)), side = 3)



# MODEL Al----------------------------------------------------------------------

# make this example reproducible
set.seed(1)

# fit PCR model
model_al = plsr(m3.Al ~ ., ncomp = 77, data = train[, -c(1:4, 6:20)], scale = TRUE, validation = "CV")

# view summary of model fitting
# summary(model_al)

# choose the number of PLS components

# validationplot(model_al, val.type = "MSEP")

msepcv.pls_al = MSEP(model_al, estimate = c("train","CV"))
plot(msepcv.pls_al, lty = 1, type = "l", legendpos = "topright", main="")

ncomp.pls_al = which.min(msepcv.pls_al$val["CV",,])-1
ncomp.pls_al

model_al = plsr(m3.Al ~ ., ncomp = ncomp.pls_al, data = train[, -c(1:4, 6:20)], scale = TRUE)

res_model_al = residuals(model_al)
plot(res_model_al[,,ncomp.pls_al], pch = 15, cex = .5 , ylab = "Residus", main = "")
abline(h = c(-2,0,2), lty = c(2,1,2))

# performance of the calibration
pred_al_train = predict(model_al, train, ncomp = ncomp.pls_al)

# calculate R2 and RMSEP
R2_al_train = R2(pred_al_train, train$m3.Al, na.rm = TRUE)
RMSEP_al_train = RMSEP(pred_al_train, train$m3.Al, na.rm = TRUE)

png("Results\\Ethiopia\\RF\\4.1. Al calibration.png", units="in", width = 8, height = 6.5, res=1000)
plot(x = train$m3.Al, y = pred_al_train, xlab = "Observed", ylab = "Predicted", main = "Al calibration")
abline(0, 1, col = "red")
fit_train_al = lm(pred_al_train ~ train$m3.Al)
sfit_train_al = summary(fit_train_al)
# sfit_train_al
abline(a = sfit_train_al$coefficients[1], b = sfit_train_al$coefficients[2], col = "black")
mtext(paste0("R2 = ", round(R2_al_train, 2),"; RMSEP = ", round(RMSEP_al_train, 3)), side = 3)


# use model to make predictions on a validation dataset

pred_al_val = predict(model_al, valid, ncomp = ncomp.pls_al)

# calculate RMSEP and R2
RMSEP_al_val = RMSEP(pred_al_val, valid$m3.Al, na.rm = TRUE)
R2_al_val = R2(pred_al_val, valid$m3.Al, na.rm = TRUE)

png("Results\\Ethiopia\\RF\\4.2. Al validation.png", units="in", width = 8, height = 6.5, res=1000)
plot(x = valid$m3.Al, y = pred_al_val, xlab = "Observed", ylab = "Predicted", main = "Al validation")
abline(0, 1, col = "red")
fit_val_al = lm(pred_al_val ~ valid$m3.Al)
sfit_val_al = summary(fit_val_al)
# sfit_val_al
abline(a = sfit_val_al$coefficients[1], b = sfit_val_al$coefficients[2], col = "black")
mtext(paste0("R2 = ", round(R2_al_val, 2),"; RMSEP = ", round(RMSEP_al_val, 3)), side = 3)



# MODEL Ca----------------------------------------------------------------------

# make this example reproducible
set.seed(1)

# fit PCR model
model_ca = plsr(m3.Ca ~ ., ncomp = 77, data = train[, -c(1:6, 8:20)], scale = TRUE, validation = "CV")

# view summary of model fitting
# summary(model_ca)

# choose the number of PLS components

# validationplot(model_ca, val.type = "MSEP")

msepcv.pls_ca = MSEP(model_ca, estimate = c("train","CV"))
plot(msepcv.pls_ca, lty = 1, type = "l", legendpos = "topright", main="")

ncomp.pls_ca = which.min(msepcv.pls_ca$val["CV",,])-1
ncomp.pls_ca

model_ca = plsr(m3.Ca ~ ., ncomp = ncomp.pls_ca, data = train[, -c(1:6, 8:20)], scale = TRUE)

res_model_ca = residuals(model_ca)
plot(res_model_ca[,,ncomp.pls_ca], pch = 15, cex = .5 , ylab = "Residus", main = "")
abline(h = c(-2,0,2), lty = c(2,1,2))

# performance of the calibration
pred_ca_train = predict(model_ca, train, ncomp = ncomp.pls_ca)

# calculate R2 and RMSEP
R2_ca_train = R2(pred_ca_train, train$m3.Ca, na.rm = TRUE)
RMSEP_ca_train = RMSEP(pred_ca_train, train$m3.Ca, na.rm = TRUE)

png("Results\\Ethiopia\\RF\\5.1. Ca calibration.png", units="in", width = 8, height = 6.5, res=1000)
plot(x = train$m3.Ca, y = pred_ca_train, xlab = "Observed", ylab = "Predicted", main = "Ca calibration")
abline(0, 1, col = "red")
fit_train_ca = lm(pred_ca_train ~ train$m3.Ca)
sfit_train_ca = summary(fit_train_ca)
# sfit_train_ca
abline(a = sfit_train_ca$coefficients[1], b = sfit_train_ca$coefficients[2], col = "black")
mtext(paste0("R2 = ", round(R2_al_train, 2),"; RMSEP = ", round(RMSEP_al_train, 3)), side = 3)


# use model to make predictions on a validation dataset

pred_ca_val = predict(model_ca, valid, ncomp = ncomp.pls_ca)

# calculate RMSEP and R2
RMSEP_ca_val = RMSEP(pred_ca_val, valid$m3.Ca, na.rm = TRUE)
R2_ca_val = R2(pred_ca_val, valid$m3.Ca, na.rm = TRUE)

plot(x = valid$m3.Ca, y = pred_ca_val, xlab = "Observed", ylab = "Predicted", main = "Ca validation")
abline(0, 1, col = "red")
fit_val_ca = lm(pred_ca_val ~ valid$m3.Ca)
sfit_val_ca = summary(fit_val_ca)
# sfit_val_ca
abline(a = sfit_val_ca$coefficients[1], b = sfit_val_ca$coefficients[2], col = "black")
mtext(paste0("R2 = ", round(R2_ca_val, 2),"; RMSEP = ", round(RMSEP_ca_val, 3)), side = 3)



# MODEL Mg----------------------------------------------------------------------

# make this example reproducible
set.seed(1)

# fit PCR model
model_mg = plsr(m3.Mg ~ ., ncomp = 77, data = train[, -c(1:8, 10:20)], scale = TRUE, validation = "CV")

# view summary of model fitting
# summary(model_mg)

# choose the number of PLS components

# validationplot(model_mg, val.type = "MSEP")

msepcv.pls_mg = MSEP(model_mg, estimate = c("train","CV"))
plot(msepcv.pls_mg, lty = 1, type = "l", legendpos = "topright", main="")

ncomp.pls_mg = which.min(msepcv.pls_mg$val["CV",,])-1
ncomp.pls_mg

model_mg = plsr(m3.Mg ~ ., ncomp = ncomp.pls_mg, data = train[, -c(1:8, 10:20)], scale = TRUE)

res_model_mg = residuals(model_mg)
plot(res_model_mg[,,ncomp.pls_mg], pch = 15, cex = .5 , ylab = "Residus", main = "")
abline(h = c(-2,0,2), lty = c(2,1,2))

# performance of the calibration
pred_mg_train = predict(model_mg, train, ncomp = ncomp.pls_mg)

# calculate R2 and RMSEP
R2_mg_train = R2(pred_mg_train, train$m3.Mg, na.rm = TRUE)
RMSEP_mg_train = RMSEP(pred_mg_train, train$m3.Mg, na.rm = TRUE)

png("Results\\Ethiopia\\PLS\\5.1. Ca calibration.png", units="in", width = 8, height = 6.5, res=1000)
plot(x = train$m3.Mg, y = pred_mg_train, xlab = "Observed", ylab = "Predicted", main = "Mg calibration")
abline(0, 1, col = "red")
fit_train_mg = lm(pred_mg_train ~ train$m3.Mg)
sfit_train_mg = summary(fit_train_mg)
# sfit_train_mg
abline(a = sfit_train_mg$coefficients[1], b = sfit_train_mg$coefficients[2], col = "black")
mtext(paste0("R2 = ", round(R2_mg_train, 2),"; RMSEP = ", round(RMSEP_mg_train, 3)), side = 3)


# use model to make predictions on a validation dataset

pred_mg_val = predict(model_mg, valid, ncomp = ncomp.pls_mg)

# calculate RMSEP and R2
RMSEP_mg_val = RMSEP(pred_mg_val, valid$m3.Mg, na.rm = TRUE)
R2_mg_val = R2(pred_mg_val, valid$m3.Mg, na.rm = TRUE)

png("Results\\Ethiopia\\PLS\\5.2. Ca validation.png", units="in", width = 8, height = 6.5, res=1000)
plot(x = valid$m3.Mg, y = pred_mg_val, xlab = "Observed", ylab = "Predicted", main = "Mg validation")
abline(0, 1, col = "red")
fit_val_mg = lm(pred_mg_val ~ valid$m3.Mg)
sfit_val_mg = summary(fit_val_mg)
# sfit_val_mg
abline(a = sfit_val_mg$coefficients[1], b = sfit_val_mg$coefficients[2], col = "black")
mtext(paste0("R2 = ", round(R2_mg_val, 2),"; RMSEP = ", round(RMSEP_mg_val, 3)), side = 3)



# MODEL CEC---------------------------------------------------------------------

# make this example reproducible
set.seed(1)

# fit PCR model
model_cec = plsr(CEC ~ ., ncomp = 77, data = train[, -c(1:16, 18:20)], scale = TRUE, validation = "CV")

# view summary of model fitting
# summary(model_cec)

# choose the number of PLS components

# validationplot(model_cec, val.type = "MSEP")

msepcv.pls_cec = MSEP(model_cec, estimate = c("train","CV"))
plot(msepcv.pls_cec, lty = 1, type = "l", legendpos = "topright", main="")

ncomp.pls_cec = which.min(msepcv.pls_cec$val["CV",,])-1
ncomp.pls_cec

model_cec = plsr(CEC ~ ., ncomp = ncomp.pls_cec, data = train[, -c(1:16, 18:20)], scale = TRUE)

res_model_cec = residuals(model_cec)
plot(res_model_cec[,,ncomp.pls_cec], pch = 15, cex = .5 , ylab = "Residus", main = "")
abline(h = c(-2,0,2), lty = c(2,1,2))

# performance of the calibration
pred_cec_train = predict(model_cec, train, ncomp = ncomp.pls_cec)

# calculate R2 and RMSEP
R2_cec_train = R2(pred_cec_train, train$CEC, na.rm = TRUE)
RMSEP_cec_train = RMSEP(pred_cec_train, train$CEC, na.rm = TRUE)

png("Results\\Ethiopia\\PLS\\6.1. Mg calibration.png", units="in", width = 8, height = 6.5, res=1000)
plot(x = train$CEC, y = pred_cec_train, xlab = "Observed", ylab = "Predicted", main = "CEC calibration")
abline(0, 1, col = "red")
fit_train_cec = lm(pred_cec_train ~ train$CEC)
sfit_train_cec = summary(fit_train_cec)
# sfit_train_cec
abline(a = sfit_train_cec$coefficients[1], b = sfit_train_cec$coefficients[2], col = "black")
mtext(paste0("R2 = ", round(R2_cec_train, 2),"; RMSEP = ", round(RMSEP_cec_train, 3)), side = 3)


# use model to make predictions on a validation dataset

pred_cec_val = predict(model_cec, valid, ncomp = ncomp.pls_cec)

# calculate RMSEP and R2
RMSEP_cec_val = RMSEP(pred_cec_val, valid$CEC, na.rm = TRUE)
R2_cec_val = R2(pred_cec_val, valid$CEC, na.rm = TRUE)

png("Results\\Ethiopia\\PLS\\6.2. Mg validation.png", units="in", width = 8, height = 6.5, res=1000)
plot(x = valid$CEC, y = pred_cec_val, xlab = "Observed", ylab = "Predicted", main = "CEC validation")
abline(0, 1, col = "red")
fit_val_cec = lm(pred_cec_val ~ valid$CEC)
sfit_val_cec = summary(fit_val_cec)
# sfit_val_cec
abline(a = sfit_val_cec$coefficients[1], b = sfit_val_cec$coefficients[2], col = "black")
mtext(paste0("R2 = ", round(R2_cec_val, 2),"; RMSEP = ", round(RMSEP_cec_val, 3)), side = 3)



# MODEL SOC---------------------------------------------------------------------

# make this example reproducible
set.seed(1)

# fit PCR model
model_soc = plsr(SOC ~ ., ncomp = 77, data = train[, -c(1, 3:20)], scale = TRUE, validation = "CV")

# view summary of model fitting
# summary(model_soc)

# choose the number of PLS components

# validationplot(model_soc, val.type = "MSEP")

msepcv.pls_soc = MSEP(model_soc, estimate = c("train","CV"))
plot(msepcv.pls_soc, lty = 1, type = "l", legendpos = "topright", main="")

ncomp.pls_soc = which.min(msepcv.pls_soc$val["CV",,])-1
ncomp.pls_soc

model_soc = plsr(SOC ~ ., ncomp = ncomp.pls_soc, data = train[, -c(1, 3:20)], scale = TRUE)

res_model_soc = residuals(model_soc)
plot(res_model_soc[,,ncomp.pls_soc], pch = 15, cex = .5 , ylab = "Residus", main = "")
abline(h = c(-2,0,2), lty = c(2,1,2))

# performance of the calibration
pred_soc_train = predict(model_soc, train, ncomp = ncomp.pls_soc)

# calculate R2 and RMSEP
R2_soc_train = R2(pred_soc_train, train$SOC, na.rm = TRUE)
RMSEP_soc_train = RMSEP(pred_soc_train, train$SOC, na.rm = TRUE)

png("Results\\Ethiopia\\PLS\\7.1. CEC calibration.png", units="in", width = 8, height = 6.5, res=1000)
plot(x = train$SOC, y = pred_soc_train, xlab = "Observed", ylab = "Predicted", main = "SOC calibration")
abline(0, 1, col = "red")
fit_train_soc = lm(pred_soc_train ~ train$SOC)
sfit_train_soc = summary(fit_train_soc)
# sfit_train_soc
abline(a = sfit_train_soc$coefficients[1], b = sfit_train_soc$coefficients[2], col = "black")
mtext(paste0("R2 = ", round(R2_soc_train, 2),"; RMSEP = ", round(RMSEP_soc_train, 3)), side = 3)


# use model to make predictions on a validation dataset

pred_soc_val = predict(model_soc, valid, ncomp = ncomp.pls_soc)

# calculate RMSEP and R2
RMSEP_soc_val = RMSEP(pred_soc_val, valid$SOC, na.rm = TRUE)
R2_soc_val = R2(pred_soc_val, valid$SOC, na.rm = TRUE)

png("Results\\Ethiopia\\PLS\\7.2. CEC validation.png", units="in", width = 8, height = 6.5, res=1000)
plot(x = valid$SOC, y = pred_soc_val, xlab = "Observed", ylab = "Predicted", main = "SOC validation")
abline(0, 1, col = "red")
fit_val_soc = lm(pred_soc_val ~ valid$SOC)
sfit_val_soc = summary(fit_val_soc)
# sfit_val_soc
abline(a = sfit_val_soc$coefficients[1], b = sfit_val_soc$coefficients[2], col = "black")
mtext(paste0("R2 = ", round(R2_soc_val, 2),"; RMSEP = ", round(RMSEP_soc_val, 3)), side = 3)


