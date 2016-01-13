# options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE)
expr1_out = args[1]
expr2_out = args[2]
rm(args)

library(DoE.base)


Design.1 <- fac.design(nfactors=2, replications=15, repeat.only=FALSE,
                       blocks=1, randomize=TRUE, seed=24625,
                       factor.names=list(size=(1:20)*50000, tlevel=c(2:6)))


Design.2 <- fac.design(nfactors=2, replications=5, repeat.only=FALSE,
                      blocks=1, randomize=TRUE, seed=24625,
                      factor.names=list(size=(1:10)*(10^6), tlevel=c(3:4)))


write_design_to_file <- function(design, fname) {
  # direct output to a file
  sink(fname, append=FALSE, split=FALSE)
  print (data.frame(design), row.names = FALSE)
  # return output to the terminal
  sink()
}

write_design_to_file(Design.1, expr1_out)
write_design_to_file(Design.2, expr2_out)
