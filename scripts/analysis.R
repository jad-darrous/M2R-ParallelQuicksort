
# setwd("/home/jad/dev/M2R-ParallelQuicksort/scripts/n2");  csv_file = "experiments_results.csv";
csv_file = "/home/jad/dev/M2R-ParallelQuicksort/data/jad-K52Jc_2016-01-10_13:44/experiments_results.csv";
nb_threads = 8

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop("usage: Rscript analysis.R <csv_file> <out_dir> <nb_threads>", call.=FALSE)
}
csv_file = args[1]
out_dir = args[2]
nb_threads = args[3]
rm(args)

# =================================================

allexp_df <- read.csv(csv_file, sep=" ", header=T, strip.white=TRUE)

# renaming for simplicity
names(allexp_df)[3:5] <- c("seq", "par", "blt")


# =================================================


# analysis the changes in execution time in relation to the size,
# for a specified number of thread in the parallel version
df_size <- allexp_df[allexp_df$threads == nb_threads,]
df_size <- df_size[df_size$size %% 200000 == 0,]


library(ggplot2)

p <- scale_colour_manual("Sorting", values = c("red", "blue", "purple"),
                         labels = c("Sequential", "Parallel", "Libc"))

o <- labs(title = "Mesearments", y = "Time", x = "Size")

# ggplot(df_size, aes(x=size)) +
#   geom_point(aes(y=seq), colour="red") +
#   geom_point(aes(y=par), colour="green") +
#   geom_point(aes(y=blt), colour="blue") +
#   labs(title = "Mesearments", y = "Time", x = "Size")

# Same as above but with a legend
# the Libc data can be removed to get higher resolution
ggplot(df_size, aes(x=size)) +
  geom_point(aes(y=seq, color="Sequential")) +
  geom_point(aes(y=par, color="Parallel")) +
  # geom_point(aes(y=blt, color="Libc")) +
  p + o

library(reshape2)

# the Libc data (col. 5) can be removed to get higher resolution when drawing
df_size_melt <- melt(df_size[,c(1,3,4,5)], id="size")

df_size_melt_par_seq <- df_size_melt[(df_size_melt$variable == 'par') | (df_size_melt$variable == 'seq'),]

# Same as above but in an elegant way
ggplot(df_size_melt, aes(size, value, colour=variable)) +
  geom_point() +
  p + o

ggsave(paste(out_dir, "size.png", sep="/"), width = 5, height = 5)


library(plyr)

compute_ci <- function(df_size_melt, sorting) {
  df <- df_size_melt[df_size_melt$variable == sorting,]
  df_ci <- ddply(df, "size", function(x) {
    mn <- mean(x$value)
    sd <- sd(x$value)
    se = 2*sd/sqrt(nrow(x))
    data.frame(mn, sd, se, variable = sorting)
  })
  # return(list("df" = df, "ci" = df_ci))
  return(df_ci)
}

df_size_melt_seq_ci <- compute_ci(df_size_melt, 'seq')
df_size_melt_par_ci <- compute_ci(df_size_melt, 'par')

# draw the CI
ggplot(df_size_melt_par_seq, aes(size, value, colour=variable)) +
  geom_point() +
  geom_errorbar(data=df_size_melt_par_ci, aes(y=mn,ymin=mn-se,ymax=mn+se)) +
  geom_errorbar(data=df_size_melt_seq_ci, aes(y=mn,ymin=mn-se,ymax=mn+se)) +
  p + o

ggsave(paste(out_dir, "size_ci.png", sep="/"), width = 5, height = 5)


# draw the regression line
ggplot(df_size_melt, aes(size, value, colour=variable)) +
  geom_point() +
  geom_smooth(method = "lm", formula = y ~ x+I(x^2), se = FALSE, size = 1) +
  p + o

ggsave(paste(out_dir, "size_regression.png", sep="/"), width = 5, height = 5)


reg_seq <- lm(seq ~ size+I(size^2), df_size)

# reg_par <- lm(par ~ size, df_size)
reg_par <- lm(par ~ size+I(size^2), df_size)
reg_par <- lm(par ~ size*I(log(size)), df_size)
reg_par <- lm(par ~ size+I(log(size)), df_size)
reg_par <- lm(par ~ size, df_size)
# summary(reg_par)

xv <- seq(100,10^6,10^3)
pred_y_par <- predict(reg_par, list(size=xv), se.fit = TRUE)
# sum(pred_y_par$se.fit)

pred_y_par <- predict(reg_par, list(size=xv))
pred_y_seq <- predict(reg_seq, list(size=xv))


# draw the regression line in a more complicated way
ggplot(data=df_size) + theme_bw() +
  labs(title = "Mesearments", y = "Time", x = "Size") +
  geom_hline(yintercept=0) + geom_vline(xintercept=0) +
  geom_point(aes(x=size,y=par), color="blue", size=2) +
  geom_point(aes(x=size,y=seq), color="red",  size=2) +
  geom_point(aes(x=size,y=blt), color="green",size=1) +
  geom_line(data=data.frame(x=xv,y=pred_y_par), aes(x=x,y=y), color="blue",size=1)  +
  geom_line(data=data.frame(x=xv,y=pred_y_seq), aes(x=x,y=y), color="red",size=1)



# =================================================

# analysis the parallel version w.r.t number of threads.
df_threads <- allexp_df[,c(1,2,4)]

# Simple plot
qplot (size, par, data=df_threads, color=factor(threads))

# same as above + regression line
ggplot(df_threads, aes(size, par, color=factor(threads))) +
  geom_point() +
  # formula = y ~ poly(x, 2)
  geom_smooth(method = "lm", size = 1, formula = y ~ x+I(x^2), se = FALSE) +
  scale_colour_manual("#Threads", values = c("red", "blue", "green", "purple", "yellow"))

ggsave(paste(out_dir, "threads.png", sep="/"), width = 5, height = 5)


# =================================================
# =================================================
# =================================================
# =================================================

#
# reg = reg_par
# summary(reg)
# plot(reg)
# par(mfrow=c(2,2));     plot(reg);     par(mfrow=c(1,1))
#
# wireframe(seq ~ size * threads, data=allexp_df)
# wireframe(par ~ size * threads, data=allexp_df)
#
# p <- wireframe(par ~ size * threads, data=allexp_df)
# npanel <- c(4, 2)
# rotx <- c(-50, -80)
# rotz <- seq(30, 300, length = npanel[1]+1)
# update(p[rep(1, prod(npanel))], layout = npanel,
#        panel = function(..., screen) {
#          panel.wireframe(..., screen = list(par = rotz[current.column()],
#                                             size = rotx[current.row()]))
#        })
#
#
# ggplot(df_size, aes(x=size)) + scale_x_log10() +
#   geom_line(aes(y=res), colour="red")  +
#   geom_point(aes(y=seq), colour="red") +
#   geom_point(aes(y=par), colour="green") +
#   geom_point(aes(y=blt), colour="blue") + ylab("Time") + xlab("size")
#
#
#
# ggplot(df_thread_M, aes(x=threads))  + scale_x_discrete(limits=c(1, 2, 4, 8, 16, 32)) +
#   geom_line(aes(y=seq), colour="red")  +
#   geom_line(aes(y=par), colour="green")  +
#   geom_line(aes(y=blt), colour="blue")  + ylab("Time") + xlab("size")


# library(plyr)
#
# aggr_df <- ddply(allexp_df, c("size", "threads"), summarize,
#            seq = mean(sequential_time), par = mean(parallel_time), blt = mean(builtin_time),
#            seq2 = seq^2, par2 = par^2, blt2 = blt^2);
