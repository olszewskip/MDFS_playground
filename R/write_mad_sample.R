mad <- MDFS::madelon
n <- 2000
mad_sample_X <- mad$data[1:n, 451:500]
mad_sample_Y <- mad$decision[1:n]
mad_sample_XY <- data.frame(X = mad_sample_X, Y = mad_sample_Y)
write.csv(mad_sample_XY, file='mad_sample.csv', row.names=FALSE)

write.table(mad, file = "madelon.csv", row.names=FALSE, col.names=FALSE)

MDFS::ComputeMaxInfoGains(mad_sample_X, mad_sample_Y, dimensions=2, divisions = 10, range=0, pseudo.count=0.25, seed=123)

mad <- MDFS::madelon
result = MDFS::ComputeMaxInfoGains(mad$data, mad$decision, dimensions=2, divisions = 2, range=0, pseudo.count=0.25, seed=123)
order(-result$IG)[1:10]
result[order(-result$IG)[1:10],]

my_df = read.csv("my_df_1.csv", header=FALSE)
MDFS::ComputeMaxInfoGains(my_df[,1:2], my_df[,3], dimensions=3, divisions = 1, range=0, pseudo.count=1e-5, seed=123)

