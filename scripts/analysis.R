
# setwd("/home/jad/dev/M2R-ParallelQuicksort/scripts/n2");  csv_file = "experiments_results.csv";
csv_file = "/home/jad/dev/M2R-ParallelQuicksort/data/jad-K52Jc_2016-01-13_07:46/experiments_results.csv";
nb_threads = 8

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) {
  stop("usage: Rscript analysis.R <csv_file> <csv_file_2> <out_dir> <nb_threads>", call.=FALSE)
}
csv_file = args[1]
csv_file_2 = args[2]
out_dir = args[3]
nb_threads = args[4]
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

o <- labs(title = "Measurements", y = "Time", x = "Size")

# ggplot(df_size, aes(x=size)) +
#   geom_point(aes(y=seq), colour="red") +
#   geom_point(aes(y=par), colour="green") +
#   geom_point(aes(y=blt), colour="blue") +
#   labs(title = "Measurements", y = "Time", x = "Size")

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

ggsave(paste(out_dir, "size_all.png", sep="/"), width = 5, height = 5)

# Same as above but without the LibC
ggplot(df_size_melt_par_seq, aes(size, value, colour=variable)) +
  geom_point() +
  p + o

ggsave(paste(out_dir, "size_par_seq.png", sep="/"), width = 5, height = 5)

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
  geom_errorbar(data=df_size_melt_par_ci, aes(y=mn, ymin=mn-se, ymax=mn+se)) +
  geom_errorbar(data=df_size_melt_seq_ci, aes(y=mn, ymin=mn-se, ymax=mn+se)) +
  p + o

ggsave(paste(out_dir, "size_ci_points.png", sep="/"), width = 5, height = 5)


df_size_melt_par_seq_mn <- ddply(df_size_melt_par_seq, c("size", "variable"), summarize,
           mn = mean(value));

# draw the CI with the mean
ggplot(df_size_melt_par_seq, aes(size, value, colour=variable)) +
  geom_line(data=df_size_melt_par_seq_mn, aes(y=mn, colour=variable)) +
  geom_errorbar(data=df_size_melt_par_ci, aes(y=mn, ymin=mn-se, ymax=mn+se)) +
  geom_errorbar(data=df_size_melt_seq_ci, aes(y=mn, ymin=mn-se, ymax=mn+se)) +
  p + o

ggsave(paste(out_dir, "size_ci_lines.png", sep="/"), width = 5, height = 5)


# draw the regression line
ggplot(df_size_melt, aes(size, value, colour=variable)) +
  geom_point() +
  geom_smooth(method = "lm", formula = y ~ x+I(x^2), se = FALSE, size = 1) +
  p + o

ggsave(paste(out_dir, "size_regression.png", sep="/"), width = 5, height = 5)


# Same as above but without the LibC
ggplot(df_size_melt_par_seq, aes(size, value, colour=variable)) +
  geom_point() +
  geom_smooth(method = "lm", formula = y ~ x+I(x^2), se = FALSE, size = 1) +
  p + o

ggsave(paste(out_dir, "size_regression_seq_par.png", sep="/"), width = 5, height = 5)


reg_seq <- lm(seq ~ size+I(size^2), df_size)

# reg_par <- lm(par ~ size, df_size)
reg_par <- lm(par ~ size+I(size^2), df_size)
reg_par <- lm(par ~ poly(size, 4), df_size)
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
  labs(title = "Measurements", y = "Time", x = "Size") +
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
  scale_colour_manual("#Threads", values = c("blue", "red", "purple", "green", "yellow"))

ggsave(paste(out_dir, "threads.png", sep="/"), width = 5, height = 5)

# same as above + regression line, just for 4 and 8
df_threads_4_8 <- df_threads[(df_threads$threads == 4) | (df_threads$threads == 8),]
ggplot(df_threads_4_8, aes(size, par, color=factor(threads))) +
  geom_point() +
  geom_smooth(method = "lm", size = 1, formula = y ~ poly(x, 2), se = FALSE) +
  # geom_smooth(method = "lm", size = 1, formula = y ~ x+I(x^2), se = FALSE) +
  scale_colour_manual("#Threads", values = c("red", "purple"))

ggsave(paste(out_dir, "threads_4_8_M.png", sep="/"), width = 5, height = 5)

# read the second expr.
df_tr48 <- read.csv(csv_file_2, sep=" ", header=T, strip.white=TRUE)

# renaming for simplicity
names(df_tr48)[3] <- c("par")

qplot (size, par, data=df_tr48, color=factor(threads))

ggplot(df_tr48, aes(size, par, color=factor(threads))) +
  geom_point() +
  geom_smooth(method = "lm", size = 1, formula = y ~ poly(x, 2), se = FALSE) +
  scale_colour_manual("#Threads", values = c("red", "purple"))

ggsave(paste(out_dir, "threads_4_8_10M.png", sep="/"), width = 5, height = 5)
