rm(list = ls())
graphics.off()
library(caret)
library(mlbench)
data(Sonar)

set.seed(107)
inTrain <- createDataPartition(
  y = Sonar$Class,
  ## the outcome data are needed
  p = .75,
  ## The percentage of data in the
  ## training set
  list = FALSE
)
## The format of the results

## The output is a set of integers for the rows of Sonar
## that belong in the training set.
str(inTrain)
#>  int [1:157, 1] 1 2 3 4 5 7 10 11 12 13 ...
#>  - attr(*, "dimnames")=List of 2
#>   ..$ : NULL
#>   ..$ : chr "Resample1"

training <- Sonar[ inTrain,]
testing  <- Sonar[-inTrain,]

nrow(training)
#> [1] 157
nrow(testing)
#> [1] 51/
plsFit <- train(
  Class ~ .,
  data = training,
  method = "pls",
  ## Center and scale the predictors for the training
  ## set and all future samples.
  preProc = c("center", "scale")
)

plsFit <- train(
  Class ~ .,
  data = training,
  method = "pls",
  preProc = c("center", "scale"),
  ## added:
  tuneLength = 15
)
ctrl <- trainControl(method = "repeatedcv", repeats = 3)

plsFit <- train(
  Class ~ .,
  data = training,
  method = "pls",
  preProc = c("center", "scale"),
  tuneLength = 15,
  ## added:
  trControl = ctrl
)
ctrl <- trainControl(
  method = "repeatedcv", 
  repeats = 3,
  classProbs = TRUE, 
  summaryFunction = twoClassSummary
)

set.seed(123)
plsFit <- train(
  Class ~ .,
  data = training,
  method = "pls",
  preProc = c("center", "scale"),
  tuneLength = 15,
  trControl = ctrl,
  metric = "ROC"
)
plsFit
#> Partial Least Squares 
#> 
#> 157 samples
#>  60 predictor
#>   2 classes: 'M', 'R' 
#> 
#> Pre-processing: centered (60), scaled (60) 
#> Resampling: Cross-Validated (10 fold, repeated 3 times) 
#> Summary of sample sizes: 141, 141, 142, 142, 141, 142, ... 
#> Resampling results across tuning parameters:
#> 
#>   ncomp  ROC    Sens   Spec 
#>    1     0.805  0.726  0.690
#>    2     0.848  0.750  0.801
#>    3     0.849  0.764  0.748
#>    4     0.836  0.765  0.736
#>    5     0.812  0.748  0.755
#>    6     0.789  0.724  0.699
#>    7     0.794  0.744  0.689
#>    8     0.801  0.739  0.698
#>    9     0.793  0.758  0.677
#>   10     0.790  0.741  0.690
#>   11     0.787  0.742  0.710
#>   12     0.777  0.737  0.715
#>   13     0.772  0.738  0.700
#>   14     0.768  0.718  0.690
#>   15     0.768  0.715  0.690
#> 
#> ROC was used to select the optimal model using
#>  the largest value.
#> The final value used for the model was ncomp = 3.
ggplot(plsFit)
plsClasses <- predict(plsFit, newdata = testing)
str(plsClasses)
#>  Factor w/ 2 levels "M","R": 2 1 1 1 2 2 1 2 2 2 ...
plsProbs <- predict(plsFit, newdata = testing, type = "prob")
head(plsProbs)
#>        M     R
#> 6  0.288 0.712
#> 8  0.648 0.352
#> 9  0.659 0.341
#> 15 0.529 0.471
#> 26 0.430 0.570
#> 27 0.492 0.508
confusionMatrix(data = plsClasses, testing$Class)
