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
  require("h2o")
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

# convert data to h2o dataframes
train_h2o = as.h2o(train[, -c(2:20)])
valid_h2o = as.h2o(valid[, -c(2:20)])

# identify predictors and response 
response_col = "pH"
predictors = setdiff(names(data_mean[, -c(2:20)]), c(response_col))

# define hyperparameter grid
hyper_params = list(ntrees = seq(500, 1500, by = 50), 
                    mtries = seq(1:6), 
                    max_depth = seq(10, 25, by = 1), 
                    min_rows = seq(5, 100, by = 5),   
                    sample_rate = c(0.55, 0.75, by = 0.05))

search_criteria = list(strategy = "RandomDiscrete", 
                       max_models = 50, 
                       max_runtime_secs = 1000)

# grid search function
rf_grid = h2o.grid(algorithm = "randomForest", 
                   training_frame = train_h2o, 
                   x = predictors, 
                   y = response_col, 
                   hyper_params = hyper_params, 
                   nfolds = 10, 
                   seed = 1, 
                   search_criteria = search_criteria)

# sort models based on rmse
rf_gridperf = h2o.getGrid(grid_id = rf_grid@grid_id, sort_by = "rmse", decreasing = F)

# get the best model
best_model = h2o.getModel(rf_gridperf@model_ids[[1]])

# get the best parameters from the best model
best_params = best_model@parameters
print(best_params)

# assign best parameters to variables
ntrees_rf_optim = best_params$ntrees
mtries_rf_optim = best_params$mtries
max_depth_rf_optim = best_params$max_depth 
sample_rate_rf_optim = best_params$sample_rate
min_rows_rf_optim = best_params$min_rows

# retrain model with train set and validate in test set with the best hyper parameters
ph_rf = h2o.randomForest(x = predictors, 
                         y = response_col, 
                         ntrees = ntrees_rf_optim, 
                         mtries = mtries_rf_optim, 
                         max_depth = max_depth_rf_optim, 
                         min_rows = min_rows_rf_optim, 
                         sample_rate = sample_rate_rf_optim, 
                         training_frame = train_h2o, 
                         validation_frame = valid_h2o, 
                         seed = 1)

summary(ph_rf)

pred_ph_train = h2o.predict(object = ph_rf, newdata = train_h2o)

pred_ph_train_df = as.data.frame(pred_ph_train)


# calculate R2 and RMSEP
R2_ph_train = R2(pred_ph_train_df, train$pH, na.rm = TRUE)
RMSEP_ph_train = RMSEP(pred_ph_train_df, train$pH, na.rm = TRUE)

png("Results\\Ethiopia\\RF\\1.1. pH calibration.png", units="in", width = 8, height = 6.5, res=1000)
plot(x = train$pH, y = pred_ph_train_df$predict, xlab = "Observed", ylab = "Predicted", main = "pH calibration")
abline(0, 1, col = "red")
fit_train_ph = lm(pred_ph_train_df$predict ~ train$pH)
sfit_train_ph = summary(fit_train_ph)
# sfit_train_ph
abline(a = fit_train_ph$coefficients[1], b = fit_train_ph$coefficients[2], col = "black")
mtext(paste0("R2 = ", round(R2_ph_train, 2),"; RMSEP = ", round(RMSEP_ph_train, 3)), side = 3)
dev.off()


# use model to make predictions on a validation dataset

pred_ph_val = h2o.predict(object = ph_rf, newdata = valid_h2o)

pred_ph_val_df = as.data.frame(pred_ph_val)

# calculate RMSEP and R2
RMSEP_ph_val = RMSEP(pred_ph_val_df$predict, valid$pH, na.rm = TRUE)
R2_ph_val = R2(pred_ph_val_df$predict, valid$pH, na.rm = TRUE)

png("Results\\Ethiopia\\RF\\1.2. pH validation.png", units="in", width = 8, height = 6.5, res=1000)
plot(x = valid$pH, y = pred_ph_val_df$predict, xlab = "Observed", ylab = "Predicted", main = "pH validation")
abline(0, 1, col = "red")
fit_val_ph = lm(pred_ph_val_df$predict ~ valid$pH)
sfit_val_ph = summary(fit_val_ph)
# sfit_val_ph
abline(a = fit_val_ph$coefficients[1], b = fit_val_ph$coefficients[2], col = "black")
mtext(paste0("R2 = ", round(R2_ph_val, 2),"; RMSEP = ", round(RMSEP_ph_val, 3)), side = 3)
dev.off()


# MODEL PSI---------------------------------------------------------------------

# make this example reproducible
set.seed(1)

# convert data to h2o dataframes
train_h2o = as.h2o(train[, -c(1:14, 16:20)])
valid_h2o = as.h2o(valid[, -c(1:14, 16:20)])

# identify predictors and response 
response_col = "PSI"
predictors = setdiff(names(data_mean[, -c(1:14, 16:20)]), c(response_col))

# define hyperparameter grid
hyper_params = list(ntrees = seq(500, 1500, by = 50), 
                    mtries = seq(1:6), 
                    max_depth = seq(10, 25, by = 1), 
                    min_rows = seq(5, 100, by = 5),   
                    sample_rate = c(0.55, 0.75, by = 0.05))

search_criteria = list(strategy = "RandomDiscrete", 
                       max_models = 50, 
                       max_runtime_secs = 1000)

# grid search function
rf_grid = h2o.grid(algorithm = "randomForest", 
                   training_frame = train_h2o, 
                   x = predictors, 
                   y = response_col, 
                   hyper_params = hyper_params, 
                   nfolds = 10, 
                   seed = 1, 
                   search_criteria = search_criteria)

# sort models based on rmse
rf_gridperf = h2o.getGrid(grid_id = rf_grid@grid_id, sort_by = "rmse", decreasing = F)

# get the best model
best_model = h2o.getModel(rf_gridperf@model_ids[[1]])

# get the best parameters from the best model
best_params = best_model@parameters
print(best_params)

# assign best parameters to variables
ntrees_rf_optim = best_params$ntrees
mtries_rf_optim = best_params$mtries
max_depth_rf_optim = best_params$max_depth 
sample_rate_rf_optim = best_params$sample_rate
min_rows_rf_optim = best_params$min_rows

# retrain model with train set and validate in test set with the best hyper parameters
po_rf = h2o.randomForest(x = predictors, 
                         y = response_col, 
                         ntrees = ntrees_rf_optim, 
                         mtries = mtries_rf_optim, 
                         max_depth = max_depth_rf_optim, 
                         min_rows = min_rows_rf_optim, 
                         sample_rate = sample_rate_rf_optim, 
                         training_frame = train_h2o, 
                         validation_frame = valid_h2o, 
                         seed = 1)

summary(po_rf)

pred_po_train = h2o.predict(object = po_rf, newdata = train_h2o)

pred_po_train_df = as.data.frame(pred_po_train)


# calculate R2 and RMSEP
R2_po_train = R2(pred_po_train_df, train$PSI, na.rm = TRUE)
RMSEP_po_train = RMSEP(pred_po_train_df, train$PSI, na.rm = TRUE)

png("Results\\Ethiopia\\RF\\2.1. PSI calibration.png", units="in", width = 8, height = 6.5, res=1000)
plot(x = train$PSI, y = pred_po_train_df$predict, xlab = "Observed", ylab = "Predicted", main = "PSI calibration")
abline(0, 1, col = "red")
fit_train_po = lm(pred_po_train_df$predict ~ train$PSI)
sfit_train_po = summary(fit_train_po)
# sfit_train_po
abline(a = fit_train_po$coefficients[1], b = fit_train_po$coefficients[2], col = "black")
mtext(paste0("R2 = ", round(R2_po_train, 2),"; RMSEP = ", round(RMSEP_po_train, 3)), side = 3)
dev.off()


# use model to make predictions on a validation dataset

pred_po_val = h2o.predict(object = po_rf, newdata = valid_h2o)

pred_po_val_df = as.data.frame(pred_po_val)

# calculate RMSEP and R2
RMSEP_po_val = RMSEP(pred_po_val_df$predict, valid$PSI, na.rm = TRUE)
R2_po_val = R2(pred_po_val_df$predict, valid$PSI, na.rm = TRUE)

png("Results\\Ethiopia\\RF\\2.2. PSI validation.png", units="in", width = 8, height = 6.5, res=1000)
plot(x = valid$PSI, y = pred_po_val_df$predict, xlab = "Observed", ylab = "Predicted", main = "PSI validation")
abline(0, 1, col = "red")
fit_val_po = lm(pred_po_val_df$predict ~ valid$PSI)
sfit_val_po = summary(fit_val_po)
# sfit_val_po
abline(a = fit_val_po$coefficients[1], b = fit_val_po$coefficients[2], col = "black")
mtext(paste0("R2 = ", round(R2_po_val, 2),"; RMSEP = ", round(RMSEP_po_val, 3)), side = 3)
dev.off()


# MODEL ExAc--------------------------------------------------------------------

# make this example reproducible
set.seed(1)

# convert data to h2o dataframes
train_h2o = as.h2o(train[, -c(1:12, 14:20)])
valid_h2o = as.h2o(valid[, -c(1:12, 14:20)])

# identify predictors and response 
response_col = "ExAc"
predictors = setdiff(names(data_mean[, -c(1:12, 14:20)]), c(response_col))

# define hyperparameter grid
hyper_params = list(ntrees = seq(500, 1500, by = 50), 
                    mtries = seq(1:6), 
                    max_depth = seq(10, 25, by = 1), 
                    min_rows = seq(5, 100, by = 5),   
                    sample_rate = c(0.55, 0.75, by = 0.05))

search_criteria = list(strategy = "RandomDiscrete", 
                       max_models = 50, 
                       max_runtime_secs = 1000)

# grid search function
rf_grid = h2o.grid(algorithm = "randomForest", 
                   training_frame = train_h2o, 
                   x = predictors, 
                   y = response_col, 
                   hyper_params = hyper_params, 
                   nfolds = 10, 
                   seed = 1, 
                   search_criteria = search_criteria)

# sort models based on rmse
rf_gridperf = h2o.getGrid(grid_id = rf_grid@grid_id, sort_by = "rmse", decreasing = F)

# get the best model
best_model = h2o.getModel(rf_gridperf@model_ids[[1]])

# get the best parameters from the best model
best_params = best_model@parameters
print(best_params)

# assign best parameters to variables
ntrees_rf_optim = best_params$ntrees
mtries_rf_optim = best_params$mtries
max_depth_rf_optim = best_params$max_depth 
sample_rate_rf_optim = best_params$sample_rate
min_rows_rf_optim = best_params$min_rows

# retrain model with train set and validate in test set with the best hyper parameters
exac_rf = h2o.randomForest(x = predictors, 
                         y = response_col, 
                         ntrees = ntrees_rf_optim, 
                         mtries = mtries_rf_optim, 
                         max_depth = max_depth_rf_optim, 
                         min_rows = min_rows_rf_optim, 
                         sample_rate = sample_rate_rf_optim, 
                         training_frame = train_h2o, 
                         validation_frame = valid_h2o, 
                         seed = 1)

summary(exac_rf)

pred_exac_train = h2o.predict(object = exac_rf, newdata = train_h2o)

pred_exac_train_df = as.data.frame(pred_exac_train)


# calculate R2 and RMSEP
R2_exac_train = R2(pred_exac_train_df, train$ExAc, na.rm = TRUE)
RMSEP_exac_train = RMSEP(pred_exac_train_df, train$ExAc, na.rm = TRUE)

png("Results\\Ethiopia\\RF\\3.1. ExAc calibration.png", units="in", width = 8, height = 6.5, res=1000)
plot(x = train$ExAc, y = pred_exac_train_df$predict, xlab = "Observed", ylab = "Predicted", main = "ExAc calibration")
abline(0, 1, col = "red")
fit_train_exac = lm(pred_exac_train_df$predict ~ train$ExAc)
sfit_train_exac = summary(fit_train_exac)
# sfit_train_exac
abline(a = fit_train_exac$coefficients[1], b = fit_train_exac$coefficients[2], col = "black")
mtext(paste0("R2 = ", round(R2_exac_train, 2),"; RMSEP = ", round(RMSEP_exac_train, 3)), side = 3)
dev.off()


# use model to make predictions on a validation dataset

pred_exac_val = h2o.predict(object = exac_rf, newdata = valid_h2o)

pred_exac_val_df = as.data.frame(pred_exac_val)

# calculate RMSEP and R2
RMSEP_exac_val = RMSEP(pred_exac_val_df$predict, valid$ExAc, na.rm = TRUE)
R2_exac_val = R2(pred_exac_val_df$predict, valid$ExAc, na.rm = TRUE)

png("Results\\Ethiopia\\RF\\3.2. ExAc validation.png", units="in", width = 8, height = 6.5, res=1000)
plot(x = valid$ExAc, y = pred_exac_val_df$predict, xlab = "Observed", ylab = "Predicted", main = "ExAc validation")
abline(0, 1, col = "red")
fit_val_exac = lm(pred_exac_val_df$predict ~ valid$ExAc)
sfit_val_exac = summary(fit_val_exac)
# sfit_val_po
abline(a = fit_val_exac$coefficients[1], b = fit_val_exac$coefficients[2], col = "black")
mtext(paste0("R2 = ", round(R2_exac_val, 2),"; RMSEP = ", round(RMSEP_exac_val, 3)), side = 3)
dev.off()


# MODEL Al----------------------------------------------------------------------

# make this example reproducible
set.seed(1)

# convert data to h2o dataframes
train_h2o = as.h2o(train[, -c(1:4, 6:20)])
valid_h2o = as.h2o(valid[, -c(1:4, 6:20)])

# identify predictors and response 
response_col = "m3.Al"
predictors = setdiff(names(data_mean[, -c(1:4, 6:20)]), c(response_col))

# define hyperparameter grid
hyper_params = list(ntrees = seq(500, 1500, by = 50), 
                    mtries = seq(1:6), 
                    max_depth = seq(10, 25, by = 1), 
                    min_rows = seq(5, 100, by = 5),   
                    sample_rate = c(0.55, 0.75, by = 0.05))

search_criteria = list(strategy = "RandomDiscrete", 
                       max_models = 50, 
                       max_runtime_secs = 1000)

# grid search function
rf_grid = h2o.grid(algorithm = "randomForest", 
                   training_frame = train_h2o, 
                   x = predictors, 
                   y = response_col, 
                   hyper_params = hyper_params, 
                   nfolds = 10, 
                   seed = 1, 
                   search_criteria = search_criteria)

# sort models based on rmse
rf_gridperf = h2o.getGrid(grid_id = rf_grid@grid_id, sort_by = "rmse", decreasing = F)

# get the best model
best_model = h2o.getModel(rf_gridperf@model_ids[[1]])

# get the best parameters from the best model
best_params = best_model@parameters
print(best_params)

# assign best parameters to variables
ntrees_rf_optim = best_params$ntrees
mtries_rf_optim = best_params$mtries
max_depth_rf_optim = best_params$max_depth 
sample_rate_rf_optim = best_params$sample_rate
min_rows_rf_optim = best_params$min_rows

# retrain model with train set and validate in test set with the best hyper parameters
al_rf = h2o.randomForest(x = predictors, 
                         y = response_col, 
                         ntrees = ntrees_rf_optim, 
                         mtries = mtries_rf_optim, 
                         max_depth = max_depth_rf_optim, 
                         min_rows = min_rows_rf_optim, 
                         sample_rate = sample_rate_rf_optim, 
                         training_frame = train_h2o, 
                         validation_frame = valid_h2o, 
                         seed = 1)

summary(al_rf)

pred_al_train = h2o.predict(object = al_rf, newdata = train_h2o)

pred_al_train_df = as.data.frame(pred_al_train)


# calculate R2 and RMSEP
R2_al_train = R2(pred_al_train_df, train$m3.Al, na.rm = TRUE)
RMSEP_al_train = RMSEP(pred_al_train_df, train$m3.Al, na.rm = TRUE)

png("Results\\Ethiopia\\RF\\4.1. Al calibration.png", units="in", width = 8, height = 6.5, res=1000)
plot(x = train$m3.Al, y = pred_al_train_df$predict, xlab = "Observed", ylab = "Predicted", main = "Al calibration")
abline(0, 1, col = "red")
fit_train_al = lm(pred_al_train_df$predict ~ train$m3.Al)
sfit_train_al = summary(fit_train_al)
# sfit_train_al
abline(a = fit_train_al$coefficients[1], b = fit_train_al$coefficients[2], col = "black")
mtext(paste0("R2 = ", round(R2_al_train, 2),"; RMSEP = ", round(RMSEP_al_train, 3)), side = 3)
dev.off()


# use model to make predictions on a validation dataset

pred_al_val = h2o.predict(object = al_rf, newdata = valid_h2o)

pred_al_val_df = as.data.frame(pred_al_val)

# calculate RMSEP and R2
RMSEP_al_val = RMSEP(pred_al_val_df$predict, valid$m3.Al, na.rm = TRUE)
R2_al_val = R2(pred_al_val_df$predict, valid$m3.Al, na.rm = TRUE)

png("Results\\Ethiopia\\RF\\4.2. Al validation.png", units="in", width = 8, height = 6.5, res=1000)
plot(x = valid$m3.Al, y = pred_al_val_df$predict, xlab = "Observed", ylab = "Predicted", main = "Al validation")
abline(0, 1, col = "red")
fit_val_al = lm(pred_al_val_df$predict ~ valid$m3.Al)
sfit_val_al = summary(fit_val_al)
# sfit_val_po
abline(a = fit_val_al$coefficients[1], b = fit_val_al$coefficients[2], col = "black")
mtext(paste0("R2 = ", round(R2_al_val, 2),"; RMSEP = ", round(RMSEP_al_val, 3)), side = 3)
dev.off()


# MODEL Ca----------------------------------------------------------------------

# make this example reproducible
set.seed(1)

# convert data to h2o dataframes
train_h2o = as.h2o(train[, -c(1:6, 8:20)])
valid_h2o = as.h2o(valid[, -c(1:6, 8:20)])

# identify predictors and response 
response_col = "m3.Ca"
predictors = setdiff(names(data_mean[, -c(1:6, 8:20)]), c(response_col))

# define hyperparameter grid
hyper_params = list(ntrees = seq(500, 1500, by = 50), 
                    mtries = seq(1:6), 
                    max_depth = seq(10, 25, by = 1), 
                    min_rows = seq(5, 100, by = 5),   
                    sample_rate = c(0.55, 0.75, by = 0.05))

search_criteria = list(strategy = "RandomDiscrete", 
                       max_models = 50, 
                       max_runtime_secs = 1000)

# grid search function
rf_grid = h2o.grid(algorithm = "randomForest", 
                   training_frame = train_h2o, 
                   x = predictors, 
                   y = response_col, 
                   hyper_params = hyper_params, 
                   nfolds = 10, 
                   seed = 1, 
                   search_criteria = search_criteria)

# sort models based on rmse
rf_gridperf = h2o.getGrid(grid_id = rf_grid@grid_id, sort_by = "rmse", decreasing = F)

# get the best model
best_model = h2o.getModel(rf_gridperf@model_ids[[1]])

# get the best parameters from the best model
best_params = best_model@parameters
print(best_params)

# assign best parameters to variables
ntrees_rf_optim = best_params$ntrees
mtries_rf_optim = best_params$mtries
max_depth_rf_optim = best_params$max_depth 
sample_rate_rf_optim = best_params$sample_rate
min_rows_rf_optim = best_params$min_rows

# retrain model with train set and validate in test set with the best hyper parameters
ca_rf = h2o.randomForest(x = predictors, 
                         y = response_col, 
                         ntrees = ntrees_rf_optim, 
                         mtries = mtries_rf_optim, 
                         max_depth = max_depth_rf_optim, 
                         min_rows = min_rows_rf_optim, 
                         sample_rate = sample_rate_rf_optim, 
                         training_frame = train_h2o, 
                         validation_frame = valid_h2o, 
                         seed = 1)

summary(ca_rf)

pred_ca_train = h2o.predict(object = ca_rf, newdata = train_h2o)

pred_ca_train_df = as.data.frame(pred_ca_train)


# calculate R2 and RMSEP
R2_ca_train = R2(pred_ca_train_df, train$m3.Ca, na.rm = TRUE)
RMSEP_ca_train = RMSEP(pred_ca_train_df, train$m3.Ca, na.rm = TRUE)

png("Results\\Ethiopia\\RF\\5.1. Ca calibration.png", units="in", width = 8, height = 6.5, res=1000)
plot(x = train$m3.Ca, y = pred_ca_train_df$predict, xlab = "Observed", ylab = "Predicted", main = "ca calibration")
abline(0, 1, col = "red")
fit_train_ca = lm(pred_ca_train_df$predict ~ train$m3.Ca)
sfit_train_ca = summary(fit_train_ca)
# sfit_train_ca
abline(a = fit_train_ca$coefficients[1], b = fit_train_ca$coefficients[2], col = "black")
mtext(paste0("R2 = ", round(R2_ca_train, 2),"; RMSEP = ", round(RMSEP_ca_train, 3)), side = 3)
dev.off()


# use model to make predictions on a validation dataset

pred_ca_val = h2o.predict(object = ca_rf, newdata = valid_h2o)

pred_ca_val_df = as.data.frame(pred_ca_val)

# calculate RMSEP and R2
RMSEP_ca_val = RMSEP(pred_ca_val_df$predict, valid$m3.Ca, na.rm = TRUE)
R2_ca_val = R2(pred_ca_val_df$predict, valid$m3.Ca, na.rm = TRUE)

png("Results\\Ethiopia\\RF\\5.2. Ca validation.png", units="in", width = 8, height = 6.5, res=1000)
plot(x = valid$m3.Ca, y = pred_ca_val_df$predict, xlab = "Observed", ylab = "Predicted", main = "ca validation")
abline(0, 1, col = "red")
fit_val_ca = lm(pred_ca_val_df$predict ~ valid$m3.Ca)
sfit_val_al = summary(fit_val_ca)
# sfit_val_al
abline(a = fit_val_ca$coefficients[1], b = fit_val_ca$coefficients[2], col = "black")
mtext(paste0("R2 = ", round(R2_ca_val, 2),"; RMSEP = ", round(RMSEP_ca_val, 3)), side = 3)
dev.off()


# MODEL Mg----------------------------------------------------------------------

# make this example reproducible
set.seed(1)

# convert data to h2o dataframes
train_h2o = as.h2o(train[, -c(1:8, 10:20)])
valid_h2o = as.h2o(valid[, -c(1:8, 10:20)])

# identify predictors and response 
response_col = "m3.Mg"
predictors = setdiff(names(data_mean[, -c(1:8, 10:20)]), c(response_col))

# define hyperparameter grid
hyper_params = list(ntrees = seq(500, 1500, by = 50), 
                    mtries = seq(1:6), 
                    max_depth = seq(10, 25, by = 1), 
                    min_rows = seq(5, 100, by = 5),   
                    sample_rate = c(0.55, 0.75, by = 0.05))

search_criteria = list(strategy = "RandomDiscrete", 
                       max_models = 50, 
                       max_runtime_secs = 1000)

# grid search function
rf_grid = h2o.grid(algorithm = "randomForest", 
                   training_frame = train_h2o, 
                   x = predictors, 
                   y = response_col, 
                   hyper_params = hyper_params, 
                   nfolds = 10, 
                   seed = 1, 
                   search_criteria = search_criteria)

# sort models based on rmse
rf_gridperf = h2o.getGrid(grid_id = rf_grid@grid_id, sort_by = "rmse", decreasing = F)

# get the best model
best_model = h2o.getModel(rf_gridperf@model_ids[[1]])

# get the best parameters from the best model
best_params = best_model@parameters
print(best_params)

# assign best parameters to variables
ntrees_rf_optim = best_params$ntrees
mtries_rf_optim = best_params$mtries
max_depth_rf_optim = best_params$max_depth 
sample_rate_rf_optim = best_params$sample_rate
min_rows_rf_optim = best_params$min_rows

# retrain model with train set and validate in test set with the best hyper parameters
mg_rf = h2o.randomForest(x = predictors, 
                         y = response_col, 
                         ntrees = ntrees_rf_optim, 
                         mtries = mtries_rf_optim, 
                         max_depth = max_depth_rf_optim, 
                         min_rows = min_rows_rf_optim, 
                         sample_rate = sample_rate_rf_optim, 
                         training_frame = train_h2o, 
                         validation_frame = valid_h2o, 
                         seed = 1)

summary(mg_rf)

pred_mg_train = h2o.predict(object = mg_rf, newdata = train_h2o)

pred_mg_train_df = as.data.frame(pred_mg_train)


# calculate R2 and RMSEP
R2_mg_train = R2(pred_mg_train_df, train$m3.Mg, na.rm = TRUE)
RMSEP_mg_train = RMSEP(pred_mg_train_df, train$m3.Mg, na.rm = TRUE)

png("Results\\Ethiopia\\RF\\6.1. Mg calibration.png", units="in", width = 8, height = 6.5, res=1000)
plot(x = train$m3.Mg, y = pred_mg_train_df$predict, xlab = "Observed", ylab = "Predicted", main = "Mg calibration")
abline(0, 1, col = "red")
fit_train_mg = lm(pred_mg_train_df$predict ~ train$m3.Mg)
sfit_train_mg = summary(fit_train_mg)
# sfit_train_mg
abline(a = fit_train_mg$coefficients[1], b = fit_train_mg$coefficients[2], col = "black")
mtext(paste0("R2 = ", round(R2_mg_train, 2),"; RMSEP = ", round(RMSEP_mg_train, 3)), side = 3)
dev.off()


# use model to make predictions on a validation dataset

pred_mg_val = h2o.predict(object = mg_rf, newdata = valid_h2o)

pred_mg_val_df = as.data.frame(pred_mg_val)

# calculate RMSEP and R2
RMSEP_mg_val = RMSEP(pred_mg_val_df$predict, valid$m3.Mg, na.rm = TRUE)
R2_mg_val = R2(pred_mg_val_df$predict, valid$m3.Mg, na.rm = TRUE)

png("Results\\Ethiopia\\RF\\6.2. Mg validation.png", units="in", width = 8, height = 6.5, res=1000)
plot(x = valid$m3.Mg, y = pred_mg_val_df$predict, xlab = "Observed", ylab = "Predicted", main = "Mg validation")
abline(0, 1, col = "red")
fit_val_mg = lm(pred_mg_val_df$predict ~ valid$m3.Mg)
sfit_val_mg = summary(fit_val_mg)
# sfit_val_al
abline(a = fit_val_mg$coefficients[1], b = fit_val_mg$coefficients[2], col = "black")
mtext(paste0("R2 = ", round(R2_mg_val, 2),"; RMSEP = ", round(RMSEP_mg_val, 3)), side = 3)
dev.off()


# MODEL CEC---------------------------------------------------------------------

# make this example reproducible
set.seed(1)

# convert data to h2o dataframes
train_h2o = as.h2o(train[, -c(1:16, 18:20)])
valid_h2o = as.h2o(valid[, -c(1:16, 18:20)])

# identify predictors and response 
response_col = "CEC"
predictors = setdiff(names(data_mean[, -c(1:16, 18:20)]), c(response_col))

# define hyperparameter grid
hyper_params = list(ntrees = seq(500, 1500, by = 50), 
                    mtries = seq(1:6), 
                    max_depth = seq(10, 25, by = 1), 
                    min_rows = seq(5, 100, by = 5),   
                    sample_rate = c(0.55, 0.75, by = 0.05))

search_criteria = list(strategy = "RandomDiscrete", 
                       max_models = 50, 
                       max_runtime_secs = 1000)

# grid search function
rf_grid = h2o.grid(algorithm = "randomForest", 
                   training_frame = train_h2o, 
                   x = predictors, 
                   y = response_col, 
                   hyper_params = hyper_params, 
                   nfolds = 10, 
                   seed = 1, 
                   search_criteria = search_criteria)

# sort models based on rmse
rf_gridperf = h2o.getGrid(grid_id = rf_grid@grid_id, sort_by = "rmse", decreasing = F)

# get the best model
best_model = h2o.getModel(rf_gridperf@model_ids[[1]])

# get the best parameters from the best model
best_params = best_model@parameters
print(best_params)

# assign best parameters to variables
ntrees_rf_optim = best_params$ntrees
mtries_rf_optim = best_params$mtries
max_depth_rf_optim = best_params$max_depth 
sample_rate_rf_optim = best_params$sample_rate
min_rows_rf_optim = best_params$min_rows

# retrain model with train set and validate in test set with the best hyper parameters
cec_rf = h2o.randomForest(x = predictors, 
                          y = response_col, 
                          ntrees = ntrees_rf_optim, 
                          mtries = mtries_rf_optim, 
                          max_depth = max_depth_rf_optim, 
                          min_rows = min_rows_rf_optim, 
                          sample_rate = sample_rate_rf_optim, 
                          training_frame = train_h2o, 
                          validation_frame = valid_h2o, 
                          seed = 1)

summary(cec_rf)

pred_cec_train = h2o.predict(object = cec_rf, newdata = train_h2o)

pred_cec_train_df = as.data.frame(pred_cec_train)


# calculate R2 and RMSEP
R2_cec_train = R2(pred_cec_train_df, train$CEC, na.rm = TRUE)
RMSEP_cec_train = RMSEP(pred_cec_train_df, train$CEC, na.rm = TRUE)

png("Results\\Ethiopia\\RF\\7.1. CEC calibration.png", units="in", width = 8, height = 6.5, res=1000)
plot(x = train$CEC, y = pred_cec_train_df$predict, xlab = "Observed", ylab = "Predicted", main = "CEC calibration")
abline(0, 1, col = "red")
fit_train_cec = lm(pred_cec_train_df$predict ~ train$CEC)
sfit_train_cec = summary(fit_train_cec)
# sfit_train_cec
abline(a = fit_train_cec$coefficients[1], b = fit_train_cec$coefficients[2], col = "black")
mtext(paste0("R2 = ", round(R2_cec_train, 2),"; RMSEP = ", round(RMSEP_cec_train, 3)), side = 3)
dev.off()


# use model to make predictions on a validation dataset

pred_cec_val = h2o.predict(object = cec_rf, newdata = valid_h2o)

pred_cec_val_df = as.data.frame(pred_cec_val)

# calculate RMSEP and R2
RMSEP_cec_val = RMSEP(pred_cec_val_df$predict, valid$CEC, na.rm = TRUE)
R2_cec_val = R2(pred_cec_val_df$predict, valid$CEC, na.rm = TRUE)

png("Results\\Ethiopia\\RF\\7.2. CEC validation.png", units="in", width = 8, height = 6.5, res=1000)
plot(x = valid$CEC, y = pred_cec_val_df$predict, xlab = "Observed", ylab = "Predicted", main = "CEC validation")
abline(0, 1, col = "red")
fit_val_cec = lm(pred_cec_val_df$predict ~ valid$CEC)
sfit_val_cec = summary(fit_val_cec)
# sfit_val_cec
abline(a = fit_val_cec$coefficients[1], b = fit_val_cec$coefficients[2], col = "black")
mtext(paste0("R2 = ", round(R2_cec_val, 2),"; RMSEP = ", round(RMSEP_cec_val, 3)), side = 3)
dev.off()


# MODEL SOC---------------------------------------------------------------------

# make this example reproducible
set.seed(1)

# convert data to h2o dataframes
train_h2o = as.h2o(train[, -c(1, 3:20)])
valid_h2o = as.h2o(valid[, -c(1, 3:20)])

# identify predictors and response 
response_col = "SOC"
predictors = setdiff(names(data_mean[, -c(1, 3:20)]), c(response_col))

# define hyperparameter grid
hyper_params = list(ntrees = seq(500, 1500, by = 50), 
                    mtries = seq(1:6), 
                    max_depth = seq(10, 25, by = 1), 
                    min_rows = seq(5, 100, by = 5),   
                    sample_rate = c(0.55, 0.75, by = 0.05))

search_criteria = list(strategy = "RandomDiscrete", 
                       max_models = 50, 
                       max_runtime_secs = 1000)

# grid search function
rf_grid = h2o.grid(algorithm = "randomForest", 
                   training_frame = train_h2o, 
                   x = predictors, 
                   y = response_col, 
                   hyper_params = hyper_params, 
                   nfolds = 10, 
                   seed = 1, 
                   search_criteria = search_criteria)

# sort models based on rmse
rf_gridperf = h2o.getGrid(grid_id = rf_grid@grid_id, sort_by = "rmse", decreasing = F)

# get the best model
best_model = h2o.getModel(rf_gridperf@model_ids[[1]])

# get the best parameters from the best model
best_params = best_model@parameters
print(best_params)

# assign best parameters to variables
ntrees_rf_optim = best_params$ntrees
mtries_rf_optim = best_params$mtries
max_depth_rf_optim = best_params$max_depth 
sample_rate_rf_optim = best_params$sample_rate
min_rows_rf_optim = best_params$min_rows

# retrain model with train set and validate in test set with the best hyper parameters
soc_rf = h2o.randomForest(x = predictors, 
                          y = response_col, 
                          ntrees = ntrees_rf_optim, 
                          mtries = mtries_rf_optim, 
                          max_depth = max_depth_rf_optim, 
                          min_rows = min_rows_rf_optim, 
                          sample_rate = sample_rate_rf_optim, 
                          training_frame = train_h2o, 
                          validation_frame = valid_h2o, 
                          seed = 1)

summary(soc_rf)

pred_soc_train = h2o.predict(object = soc_rf, newdata = train_h2o)

pred_soc_train_df = as.data.frame(pred_soc_train)


# calculate R2 and RMSEP
R2_soc_train = R2(pred_soc_train_df, train$SOC, na.rm = TRUE)
RMSEP_soc_train = RMSEP(pred_soc_train_df, train$SOC, na.rm = TRUE)

png("Results\\Zimbabwe\\RF\\7.1. SOC calibration.png", units="in", width = 8, height = 6.5, res=1000)
plot(x = train$SOC, y = pred_soc_train_df$predict, xlab = "Observed", ylab = "Predicted", main = "SOC calibration")
abline(0, 1, col = "red")
fit_train_soc = lm(pred_soc_train_df$predict ~ train$SOC)
sfit_train_soc = summary(fit_train_soc)
# sfit_train_soc
abline(a = fit_train_soc$coefficients[1], b = fit_train_soc$coefficients[2], col = "black")
mtext(paste0("R2 = ", round(R2_soc_train, 2),"; RMSEP = ", round(RMSEP_soc_train, 3)), side = 3)
dev.off()

# use model to make predictions on a validation dataset

pred_soc_val = h2o.predict(object = soc_rf, newdata = valid_h2o)

pred_soc_val_df = as.data.frame(pred_soc_val)

# calculate RMSEP and R2
RMSEP_soc_val = RMSEP(pred_soc_val_df$predict, valid$SOC, na.rm = TRUE)
R2_soc_val = R2(pred_soc_val_df$predict, valid$SOC, na.rm = TRUE)

png("Results\\Zimbabwe\\RF\\7.2. SOC validation.png", units="in", width = 8, height = 6.5, res=1000)
plot(x = valid$SOC, y = pred_soc_val_df$predict, xlab = "Observed", ylab = "Predicted", main = "SOC validation")
abline(0, 1, col = "red")
fit_val_soc = lm(pred_soc_val_df$predict ~ valid$SOC)
sfit_val_soc = summary(fit_val_soc)
# sfit_val_soc
abline(a = fit_val_soc$coefficients[1], b = fit_val_soc$coefficients[2], col = "black")
mtext(paste0("R2 = ", round(R2_soc_val, 2),"; RMSEP = ", round(RMSEP_soc_val, 3)), side = 3)
dev.off()



