mad <- read.csv('mad_sample.csv')
mad_X <- mad[,1:10]
mad_Y <- mad[,11]
divs <- 3
out <- MDFS::ComputeMaxInfoGains(mad_X,mad_Y, dimensions=2, divisions = divs, range=1, return.tuples = T)
write.csv(out, file='mad_sample_IG.csv', row.names=FALSE)
out
