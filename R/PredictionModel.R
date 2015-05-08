TrainModel <- function(train, method = "") {
  switch(method,
         gbm = gbm(y~., data = train, n.trees = 50, distribution = "gaussian",verbose = FALSE),
         linear = glm(y~., family = gaussian, data = train),
         xgboost = xgboost(data = train[,-ncol(train)], label = train[,"y"], nrounds = 1, objective = "reg:linear", max.depth = 2, verbose = 0))
}

TestPrediction <- function(model, test, method = "") {
  switch(method,
         gbm = predict.gbm(model, test, n.trees = 50),
         linear = predict.glm(model, test),
         xgboost = predict(object = model, as.matrix(test)))  
}