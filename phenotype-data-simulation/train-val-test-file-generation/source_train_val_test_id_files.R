#!/usr/bin/env Rscript

range <- 0:119999

#Split source into training / validation / test 
train = sample(range, 100000)
train_dummy = range[!(range %in% train)]
val = sample(train_dummy, 10000)
test = train_dummy[!(train_dummy %in% val)]

#Create text files for plink sample filtering 
train_data = matrix(nrow=100000, ncol=2)
for (i in 1:100000){
  train_data[i,1] = paste("id1_", train[i], sep="")
  train_data[i,2] = paste("id2_", train[i], sep="")
}

val_data = matrix(nrow=10000, ncol=2)
for (i in 1:10000){
  val_data[i,1] = paste("id1_", val[i], sep="")
  val_data[i,2] = paste("id2_", val[i], sep="")
}

test_data = matrix(nrow=10000, ncol=2)
for (i in 1:10000){
  test_data[i,1] = paste("id1_", test[i], sep="")
  test_data[i,2] = paste("id2_", test[i], sep="")
}

write.table(train_data, "train_samples.txt", sep=" ", quote = F, row.names = F)
write.table(val_data, "val_samples.txt", sep=" ", quote = F, row.names = F)
write.table(test_data, "test_samples.txt", sep=" ", quote = F, row.names = F)
