mad <- MDFS::madelon
n <- 100
mad_sample_X <- mad$data[1:n,1:3]
mad_sample_Y <- mad$decision[1:n]
mad_sample_XY <- data.frame(X = mad_sample_X, Y = mad_sample_Y)
write.csv(mad_sample_XY, file='mad_sample.csv', row.names=FALSE)
