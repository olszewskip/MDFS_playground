mad <- MDFS::madelon
n <- 300
mad_sample_X <- mad$data[1:n,301:500]
mad_sample_Y <- mad$decision[1:n]
mad_sample_XY <- data.frame(X = mad_sample_X, Y = mad_sample_Y)
write.csv(mad_sample_XY, file='mad_sample.csv', row.names=FALSE)

write.table(mad, file = "madelon.csv", row.names=FALSE, col.names=FALSE)

MDFS::MDFS(mad_sample_X, mad_sample_Y, dimension=2, divisions = 10, range=0, seed=123)
