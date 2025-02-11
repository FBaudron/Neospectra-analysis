#' ---
#' title: "RF analysis Neospectra - soils Zimbabwe - Double Burden project"
#' author: "Frédéric Baudron"
#' date: "February 10th, 2025"
#' ---


# CLEAR ENVIRONMENT-------------------------------------------------------------

rm(list = ls())


# LOADING NECESSARY PACKAGES----------------------------------------------------

suppressPackageStartupMessages({
  require("h2o")
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


h2o.init()
h2o.clusterInfo()



# MODEL PH----------------------------------------------------------------------

# make this example reproducible
set.seed(1)

# convert data to h2o dataframes
train_h2o = as.h2o(train[, -c(2:33)])
valid_h2o = as.h2o(valid[, -c(2:33)])

# identify predictors and response 
response_col = "pH"
predictors = setdiff(names(data_mean[, -c(2:33)]), c(response_col))

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

png("Results\\Zimbabwe\\RF\\1.1. pH calibration.png", units="in", width = 8, height = 6.5, res=1000)
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

png("Results\\Zimbabwe\\RF\\1.2. pH validation.png", units="in", width = 8, height = 6.5, res=1000)
plot(x = valid$pH, y = pred_ph_val_df$predict, xlab = "Observed", ylab = "Predicted", main = "pH validation")
abline(0, 1, col = "red")
fit_val_ph = lm(pred_ph_val_df$predict ~ valid$pH)
sfit_val_ph = summary(fit_val_ph)
# sfit_val_ph
abline(a = fit_val_ph$coefficients[1], b = fit_val_ph$coefficients[2], col = "black")
mtext(paste0("R2 = ", round(R2_ph_val, 2),"; RMSEP = ", round(RMSEP_ph_val, 3)), side = 3)
dev.off()


# MODEL OLSEN P-----------------------------------------------------------------

# make this example reproducible
set.seed(1)

# convert data to h2o dataframes
train_h2o = as.h2o(train[, -c(1:3, 5:33)])
valid_h2o = as.h2o(valid[, -c(1:3, 5:33)])

# identify predictors and response 
response_col = "OlsenP"
predictors = setdiff(names(data_mean[, -c(1:3, 5:33)]), c(response_col))

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
R2_po_train = R2(pred_po_train_df, train$OlsenP, na.rm = TRUE)
RMSEP_po_train = RMSEP(pred_po_train_df, train$OlsenP, na.rm = TRUE)

png("Results\\Zimbabwe\\RF\\2.1. P calibration.png", units="in", width = 8, height = 6.5, res=1000)
plot(x = train$OlsenP, y = pred_po_train_df$predict, xlab = "Observed", ylab = "Predicted", main = "P calibration")
abline(0, 1, col = "red")
fit_train_po = lm(pred_po_train_df$predict ~ train$OlsenP)
sfit_train_po = summary(fit_train_po)
# sfit_train_po
abline(a = fit_train_po$coefficients[1], b = fit_train_po$coefficients[2], col = "black")
mtext(paste0("R2 = ", round(R2_po_train, 2),"; RMSEP = ", round(RMSEP_po_train, 3)), side = 3)
dev.off()

# use model to make predictions on a validation dataset

pred_po_val = h2o.predict(object = po_rf, newdata = valid_h2o)

pred_po_val_df = as.data.frame(pred_po_val)

# calculate RMSEP and R2
RMSEP_po_val = RMSEP(pred_po_val_df$predict, valid$OlsenP, na.rm = TRUE)
R2_po_val = R2(pred_po_val_df$predict, valid$OlsenP, na.rm = TRUE)

png("Results\\Zimbabwe\\RF\\2.1. P validationation.png", units="in", width = 8, height = 6.5, res=1000)
plot(x = valid$OlsenP, y = pred_po_val_df$predict, xlab = "Observed", ylab = "Predicted", main = "P validation")
abline(0, 1, col = "red")
fit_val_po = lm(pred_po_val_df$predict ~ valid$OlsenP)
sfit_val_po = summary(fit_val_po)
# sfit_val_po
abline(a = fit_val_po$coefficients[1], b = fit_val_po$coefficients[2], col = "black")
mtext(paste0("R2 = ", round(R2_po_val, 2),"; RMSEP = ", round(RMSEP_po_val, 3)), side = 3)
dev.off()


# MODEL Al----------------------------------------------------------------------

# make this example reproducible
set.seed(1)

# convert data to h2o dataframes
train_h2o = as.h2o(train[, -c(1:13, 15:33)])
valid_h2o = as.h2o(valid[, -c(1:13, 15:33)])

# identify predictors and response 
response_col = "Al"
predictors = setdiff(names(data_mean[, -c(1:13, 15:33)]), c(response_col))

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
R2_al_train = R2(pred_al_train_df, train$Al, na.rm = TRUE)
RMSEP_al_train = RMSEP(pred_al_train_df, train$Al, na.rm = TRUE)

png("Results\\Zimbabwe\\RF\\3.1. Al calibration.png", units="in", width = 8, height = 6.5, res=1000)
plot(x = train$Al, y = pred_al_train_df$predict, xlab = "Observed", ylab = "Predicted", main = "Al calibration")
abline(0, 1, col = "red")
fit_train_al = lm(pred_al_train_df$predict ~ train$Al)
sfit_train_al = summary(fit_train_al)
# sfit_train_al
abline(a = fit_train_al$coefficients[1], b = fit_train_al$coefficients[2], col = "black")
mtext(paste0("R2 = ", round(R2_al_train, 2),"; RMSEP = ", round(RMSEP_al_train, 3)), side = 3)
dev.off()

# use model to make predictions on a validation dataset

pred_al_val = h2o.predict(object = al_rf, newdata = valid_h2o)

pred_al_val_df = as.data.frame(pred_al_val)

# calculate RMSEP and R2
RMSEP_al_val = RMSEP(pred_al_val_df$predict, valid$Al, na.rm = TRUE)
R2_al_val = R2(pred_al_val_df$predict, valid$Al, na.rm = TRUE)

png("Results\\Zimbabwe\\RF\\3.1. Al validation.png", units="in", width = 8, height = 6.5, res=1000)
plot(x = valid$Al, y = pred_al_val_df$predict, xlab = "Observed", ylab = "Predicted", main = "Al validation")
abline(0, 1, col = "red")
fit_val_al = lm(pred_al_val_df$predict ~ valid$Al)
sfit_val_al = summary(fit_val_al)
# sfit_val_po
abline(a = fit_val_al$coefficients[1], b = fit_val_al$coefficients[2], col = "black")
mtext(paste0("R2 = ", round(R2_al_val, 2),"; RMSEP = ", round(RMSEP_al_val, 3)), side = 3)
dev.off()


# MODEL Ca----------------------------------------------------------------------

# make this example reproducible
set.seed(1)

# convert data to h2o dataframes
train_h2o = as.h2o(train[, -c(1:4, 6:33)])
valid_h2o = as.h2o(valid[, -c(1:4, 6:33)])

# identify predictors and response 
response_col = "Ca.ex"
predictors = setdiff(names(data_mean[, -c(1:4, 6:33)]), c(response_col))

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
R2_ca_train = R2(pred_ca_train_df, train$Ca.ex, na.rm = TRUE)
RMSEP_ca_train = RMSEP(pred_ca_train_df, train$Ca.ex, na.rm = TRUE)

png("Results\\Zimbabwe\\RF\\4.1. Ca calibration.png", units="in", width = 8, height = 6.5, res=1000)
plot(x = train$Ca.ex, y = pred_ca_train_df$predict, xlab = "Observed", ylab = "Predicted", main = "Ca calibration")
abline(0, 1, col = "red")
fit_train_ca = lm(pred_ca_train_df$predict ~ train$Ca.ex)
sfit_train_ca = summary(fit_train_ca)
# sfit_train_ca
abline(a = fit_train_ca$coefficients[1], b = fit_train_ca$coefficients[2], col = "black")
mtext(paste0("R2 = ", round(R2_ca_train, 2),"; RMSEP = ", round(RMSEP_ca_train, 3)), side = 3)
dev.off()

# use model to make predictions on a validation dataset

pred_ca_val = h2o.predict(object = ca_rf, newdata = valid_h2o)

pred_ca_val_df = as.data.frame(pred_ca_val)

# calculate RMSEP and R2
RMSEP_ca_val = RMSEP(pred_ca_val_df$predict, valid$Ca.ex, na.rm = TRUE)
R2_ca_val = R2(pred_ca_val_df$predict, valid$Ca.ex, na.rm = TRUE)

png("Results\\Zimbabwe\\RF\\4.1. Ca validation.png", units="in", width = 8, height = 6.5, res=1000)
plot(x = valid$Ca.ex, y = pred_ca_val_df$predict, xlab = "Observed", ylab = "Predicted", main = "Ca validation")
abline(0, 1, col = "red")
fit_val_ca = lm(pred_ca_val_df$predict ~ valid$Ca.ex)
sfit_val_al = summary(fit_val_ca)
# sfit_val_al
abline(a = fit_val_ca$coefficients[1], b = fit_val_ca$coefficients[2], col = "black")
mtext(paste0("R2 = ", round(R2_ca_val, 2),"; RMSEP = ", round(RMSEP_ca_val, 3)), side = 3)
dev.off()


# MODEL Mg----------------------------------------------------------------------

# make this example reproducible
set.seed(1)

# convert data to h2o dataframes
train_h2o = as.h2o(train[, -c(1:6, 8:33)])
valid_h2o = as.h2o(valid[, -c(1:6, 8:33)])

# identify predictors and response 
response_col = "Mg.ex"
predictors = setdiff(names(data_mean[, -c(1:6, 8:33)]), c(response_col))

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
R2_mg_train = R2(pred_mg_train_df, train$Mg.ex, na.rm = TRUE)
RMSEP_mg_train = RMSEP(pred_mg_train_df, train$Mg.ex, na.rm = TRUE)

png("Results\\Zimbabwe\\RF\\5.1. Mg calibration.png", units="in", width = 8, height = 6.5, res=1000)
plot(x = train$Mg.ex, y = pred_mg_train_df$predict, xlab = "Observed", ylab = "Predicted", main = "Mg calibration")
abline(0, 1, col = "red")
fit_train_mg = lm(pred_mg_train_df$predict ~ train$Mg.ex)
sfit_train_mg = summary(fit_train_mg)
# sfit_train_mg
abline(a = fit_train_mg$coefficients[1], b = fit_train_mg$coefficients[2], col = "black")
mtext(paste0("R2 = ", round(R2_mg_train, 2),"; RMSEP = ", round(RMSEP_mg_train, 3)), side = 3)
dev.off()

# use model to make predictions on a validation dataset

pred_mg_val = h2o.predict(object = mg_rf, newdata = valid_h2o)

pred_mg_val_df = as.data.frame(pred_mg_val)

# calculate RMSEP and R2
RMSEP_mg_val = RMSEP(pred_mg_val_df$predict, valid$Mg.ex, na.rm = TRUE)
R2_mg_val = R2(pred_mg_val_df$predict, valid$Mg.ex, na.rm = TRUE)

png("Results\\Zimbabwe\\RF\\5.2. Mg validation.png", units="in", width = 8, height = 6.5, res=1000)
plot(x = valid$Mg.ex, y = pred_mg_val_df$predict, xlab = "Observed", ylab = "Predicted", main = "Mg validation")
abline(0, 1, col = "red")
fit_val_mg = lm(pred_mg_val_df$predict ~ valid$Mg.ex)
sfit_val_mg = summary(fit_val_mg)
# sfit_val_al
abline(a = fit_val_mg$coefficients[1], b = fit_val_mg$coefficients[2], col = "black")
mtext(paste0("R2 = ", round(R2_mg_val, 2),"; RMSEP = ", round(RMSEP_mg_val, 3)), side = 3)
dev.off()


# MODEL CEC---------------------------------------------------------------------

# make this example reproducible
set.seed(1)

# convert data to h2o dataframes
train_h2o = as.h2o(train[, -c(1:8, 10:33)])
valid_h2o = as.h2o(valid[, -c(1:8, 10:33)])

# identify predictors and response 
response_col = "eCEC.clean"
predictors = setdiff(names(data_mean[, -c(1:8, 10:33)]), c(response_col))

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
R2_cec_train = R2(pred_cec_train_df, train$eCEC.clean, na.rm = TRUE)
RMSEP_cec_train = RMSEP(pred_cec_train_df, train$eCEC.clean, na.rm = TRUE)

png("Results\\Zimbabwe\\RF\\6.1. CEC calibration.png", units="in", width = 8, height = 6.5, res=1000)
plot(x = train$eCEC.clean, y = pred_cec_train_df$predict, xlab = "Observed", ylab = "Predicted", main = "CEC calibration")
abline(0, 1, col = "red")
fit_train_cec = lm(pred_cec_train_df$predict ~ train$eCEC.clean)
sfit_train_cec = summary(fit_train_cec)
# sfit_train_cec
abline(a = fit_train_cec$coefficients[1], b = fit_train_cec$coefficients[2], col = "black")
mtext(paste0("R2 = ", round(R2_cec_train, 2),"; RMSEP = ", round(RMSEP_cec_train, 3)), side = 3)
dev.off()

# use model to make predictions on a validation dataset

pred_cec_val = h2o.predict(object = cec_rf, newdata = valid_h2o)

pred_cec_val_df = as.data.frame(pred_cec_val)

# calculate RMSEP and R2
RMSEP_cec_val = RMSEP(pred_cec_val_df$predict, valid$eCEC.clean, na.rm = TRUE)
R2_cec_val = R2(pred_cec_val_df$predict, valid$eCEC.clean, na.rm = TRUE)

png("Results\\Zimbabwe\\RF\\6.2. CEC validation.png", units="in", width = 8, height = 6.5, res=1000)
plot(x = valid$eCEC.clean, y = pred_cec_val_df$predict, xlab = "Observed", ylab = "Predicted", main = "CEC validation")
abline(0, 1, col = "red")
fit_val_cec = lm(pred_cec_val_df$predict ~ valid$eCEC.clean)
sfit_val_cec = summary(fit_val_cec)
# sfit_val_cec
abline(a = fit_val_cec$coefficients[1], b = fit_val_cec$coefficients[2], col = "black")
mtext(paste0("R2 = ", round(R2_cec_val, 2),"; RMSEP = ", round(RMSEP_cec_val, 3)), side = 3)
dev.off()


# MODEL SOC---------------------------------------------------------------------

# make this example reproducible
set.seed(1)

# convert data to h2o dataframes
train_h2o = as.h2o(train[, -c(1:2, 4:33)])
valid_h2o = as.h2o(valid[, -c(1:2, 4:33)])

# identify predictors and response 
response_col = "Total.soil.C"
predictors = setdiff(names(data_mean[, -c(1:2, 4:33)]), c(response_col))

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
R2_soc_train = R2(pred_soc_train_df, train$Total.soil.C, na.rm = TRUE)
RMSEP_soc_train = RMSEP(pred_soc_train_df, train$Total.soil.C, na.rm = TRUE)

png("Results\\Zimbabwe\\RF\\7.1. SOC calibration.png", units="in", width = 8, height = 6.5, res=1000)
plot(x = train$Total.soil.C, y = pred_soc_train_df$predict, xlab = "Observed", ylab = "Predicted", main = "SOC calibration")
abline(0, 1, col = "red")
fit_train_soc = lm(pred_soc_train_df$predict ~ train$Total.soil.C)
sfit_train_soc = summary(fit_train_soc)
# sfit_train_soc
abline(a = fit_train_soc$coefficients[1], b = fit_train_soc$coefficients[2], col = "black")
mtext(paste0("R2 = ", round(R2_soc_train, 2),"; RMSEP = ", round(RMSEP_soc_train, 3)), side = 3)
dev.off()

# use model to make predictions on a validation dataset

pred_soc_val = h2o.predict(object = soc_rf, newdata = valid_h2o)

pred_soc_val_df = as.data.frame(pred_soc_val)

# calculate RMSEP and R2
RMSEP_soc_val = RMSEP(pred_soc_val_df$predict, valid$Total.soil.C, na.rm = TRUE)
R2_soc_val = R2(pred_soc_val_df$predict, valid$Total.soil.C, na.rm = TRUE)

png("Results\\Zimbabwe\\RF\\7.2. SOC validation.png", units="in", width = 8, height = 6.5, res=1000)
plot(x = valid$Total.soil.C, y = pred_soc_val_df$predict, xlab = "Observed", ylab = "Predicted", main = "SOC validation")
abline(0, 1, col = "red")
fit_val_soc = lm(pred_soc_val_df$predict ~ valid$Total.soil.C)
sfit_val_soc = summary(fit_val_soc)
# sfit_val_soc
abline(a = fit_val_soc$coefficients[1], b = fit_val_soc$coefficients[2], col = "black")
mtext(paste0("R2 = ", round(R2_soc_val, 2),"; RMSEP = ", round(RMSEP_soc_val, 3)), side = 3)
dev.off()


