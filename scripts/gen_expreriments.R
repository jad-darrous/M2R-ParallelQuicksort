# options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE)
out_file = args[1]
rm(args)

print (out_file)
library(DoE.base)
Design.1 <- fac.design(nfactors=2, replications=15, repeat.only=FALSE,
                       blocks=1, randomize=TRUE, seed=24625,
                       factor.names=list(size=(1:20)*50000, tlevel=c(2:6)))

# direct output to a file
sink(out_file, append=FALSE, split=FALSE)
df<-data.frame(Design.1)
print (df, row.names = FALSE)
# return output to the terminal
sink()
