#GP Classification with kernlab

#Downloading the data

data = read.csv("https://github.com/STIMALiU/AdvMLCourse/raw/master/GaussianProcess/Code/banknoteFraud.csv", header=FALSE, sep=",")
names(data) = c("varWave","skewWave","kurtWave","entropyWave","fraud")
data[,5] = as.factor(data[,5])

#Dividing the data into training and testing
set.seed(111)
SelectTraining = sample(1:dim(data)[1],
                        size = 1000,
                        replace = FALSE)
trainingData = data[SelectTraining,]
testingData = data[-SelectTraining,]

GPfitFraud <- gausspr(fraud ~  varWave + skewWave, data=trainingData)
GPfitFraud

# predict on the training set
trainPredict = predict(GPfitFraud,trainingData[,1:2])
confusionMatrix = table(trainPredict, trainingData[,5]) 
trainAccuracy = sum(diag(confusionMatrix))/sum(confusionMatrix)
trainAccuracy

#Plotting class probabilities
# class probabilities 
probPreds <- predict(GPfitFraud, trainingData[,1:2], type="probabilities")
x1 <- seq(min(trainingData[,1]),max(trainingData[,1]),length=100)
x2 <- seq(min(trainingData[,2]),max(trainingData[,2]),length=100)
gridPoints <- meshgrid(x1, x2)
gridPoints <- cbind(c(gridPoints$x), c(gridPoints$y))

gridPoints <- data.frame(gridPoints)
names(gridPoints) <- names(trainingData)[1:2]
probPreds <- predict(GPfitFraud, gridPoints, type="probabilities")

# Plotting for Prob(setosa)
contour(x1,x2,matrix(probPreds[,1],100,byrow = TRUE), 20, xlab = "varWave", ylab = "skewWave", main = 'Prob(Fraud) - Fraud is blue')
points(trainingData[trainingData[,5]==1,1],trainingData[trainingData[,5]==1,1],col="blue")
points(trainingData[trainingData[,5]=='virginica',3],trainingData[trainingData[,5]=='virginica',4],col="blue")
points(trainingData[trainingData[,5]=='versicolor',3],trainingData[trainingData[,5]=='versicolor',4],col="green")

# Plotting for Prob(Versicolor)
contour(x1,x2,matrix(probPreds[,2],100,byrow = TRUE), 20, xlab = "Petal.Length", ylab = "Petal.Width", main = 'Prob(Versicolor) - Versicolor is green')
points(trainingData[trainingData[,5]=='setosa',3],trainingData[trainingData[,5]=='setosa',4],col="red")
points(trainingData[trainingData[,5]=='virginica',3],trainingData[trainingData[,5]=='virginica',4],col="blue")
points(trainingData[trainingData[,5]=='versicolor',3],trainingData[trainingData[,5]=='versicolor',4],col="green")












#Question2: calculating accuracy for the testing data###########################
#accuracy for the test data
testPredict = predict(GPfitFraud,testingData[,1:2])
confusionMatrix = table(testPredict, testingData[,5]) 
testAccuracy = sum(diag(confusionMatrix))/sum(confusionMatrix)
testAccuracy

#The accuracy for the training data is trainAccuracy and the test accuracy is testAccuracy




#Question3: training the model on the whole data ###############################

GPfitFraudAll <- gausspr(fraud ~  varWave + skewWave + kurtWave + entropyWave, data=trainingData)
GPfitFraudAll

testAllPredict = predict(GPfitFraudAll,testingData[,1:4])
confusionMatrixAll = table(testAllPredict, testingData[,5]) 
testAllAccuracy = sum(diag(confusionMatrixAll))/sum(confusionMatrixAll)
testAllAccuracy

#When fitting the model with all the data the accuracy increases up to 99 percent