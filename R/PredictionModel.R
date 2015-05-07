TrainModel <- function(train, method = "") {
  switch(method,
         gbm = gbm(y~., data = train, n.trees = 50, distribution = "gaussian",verbose = FALSE),
         linear = glm(y~., family = gaussian, data = train))
}

TestPrediction <- function(model, test, method = "") {
  switch(method,
         gbm = predict.gbm(model, test, n.trees = 50),
         linear = predict.glm(model, test))  
}