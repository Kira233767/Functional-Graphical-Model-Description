# EEG Data, Full Application with functions for T_QDA


source('/mnt/home/lipeide/Classification/Project1/EEGApplication/T_QDA.R')



# Input Data
load('/mnt/home/lipeide/Classification/Project1/EEGApplication/Patient1.Rda')
load('/mnt/home/lipeide/Classification/Project1/EEGApplication/Patient2.Rda')
load('/mnt/home/lipeide/Classification/Project1/EEGApplication/Patient7.Rda')
load('/mnt/home/lipeide/Classification/Project1/EEGApplication/Patient4.Rda')
load('/mnt/home/lipeide/Classification/Project1/EEGApplication/Patient5.Rda')
load('/mnt/home/lipeide/Classification/Project1/EEGApplication/Patient6.Rda')



load('/mnt/home/lipeide/Classification/Project1/EEGApplication/Patient72.Rda')
load('/mnt/home/lipeide/Classification/Project1/EEGApplication/Patient67.Rda')
load('/mnt/home/lipeide/Classification/Project1/EEGApplication/Patient68.Rda')
load('/mnt/home/lipeide/Classification/Project1/EEGApplication/Patient69.Rda')
load('/mnt/home/lipeide/Classification/Project1/EEGApplication/Patient70.Rda')
load('/mnt/home/lipeide/Classification/Project1/EEGApplication/Patient71.Rda')
# load('Patient72.Rda')
# load('Patient73.Rda')
# load('Patient74.Rda')
# load('Patient75.Rda')
# load('Patient76.Rda')
# load('Patient77.Rda')



# Selecting 80 complete observations from each class

#60 different observations at the same point
Training.flag <- seq(1, 15360, length.out = 15360)  
Testing.flag <- seq(15361, length.out = 5120) 



# 360 * 2 training set and 20 testing set
X.Train <- rbind(Patient1[Training.flag, ], Patient2[Training.flag, ],
                 Patient7[Training.flag, ], Patient4[Training.flag, ],
                 Patient5[Training.flag, ], Patient6[Training.flag, ])
Y.Train <- rbind(Patient72[Training.flag, ], Patient67[Training.flag, ],
                 Patient68[Training.flag, ], Patient69[Training.flag, ],
                 Patient70[Training.flag, ], Patient71[Training.flag, ])



Testset <- rbind(Patient1[Testing.flag, ], Patient2[Testing.flag, ],
                 Patient7[Testing.flag, ], Patient4[Testing.flag, ],
                 Patient5[Testing.flag, ], Patient6[Testing.flag, ],
                 Patient72[Testing.flag, ], Patient67[Testing.flag, ],
                 Patient68[Testing.flag, ], Patient69[Testing.flag, ],
                 Patient70[Testing.flag, ], Patient71[Testing.flag, ])



# 
# # Ordinary procedure
X.Train$Time <- X.Train$Time * (1 / 256)
Y.Train$Time <- Y.Train$Time * (1 / 256)
Testset$Time <- Testset$Time * (1 / 256)
Prob.classX <- rep(0, nrow(Testset) / 256)

Kernel.Name <- c("Gaussian", "Epan", "Tri-Cube")
for(i in 1 : 3){
  for (j in 1 : length(Prob.classX)){
    Index.start <- (j - 1) * 256 + 1
    Index.end <- Index.start + 255
    Temp.patient <- Testset[c(Index.start : Index.end), ]
    Prob.classX[j] <- Time.QDA(Testing = Temp.patient, Training.X = X.Train,
                               Training.Y = Y.Train, Time.seq.new = Temp.patient$Time, flag.kernel = i)
  }
  
  
  TruePositive <- sum(Prob.classX[c(1 : 120)] > 0.5)
  TrueNegative <- sum(Prob.classX[-c(1 : 120)] <= 0.5)
  
  Sensitivity <- TurePostive / 120
  Specificity <- TrueNegative / (length(Prob.classX) - 120)
  Precision <- TruePositive / sum(Prob.classX > 0.5)
  Accuracy <- (TruePositive + TrueNegative) / length(Prob.classX)
  
  print(Kernel.Name[i])
  print(list('Sensitivity' = Sensitivity, 'Specificity' = Specificity, 'Precision' = Precision, 'Accuracy' = Accuracy))
}







