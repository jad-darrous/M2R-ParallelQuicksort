
out_dir = "/home/jad/dev/M2R-ParallelQuicksort/data/last";
# setwd(out_dir);  csv_file = "experiments_results.csv";

csv_file_medium = paste(out_dir, "experiments_medium_results.csv", sep="/");
csv_file_small = paste(out_dir, "experiments_small_results.csv", sep="/");
csv_file_big = paste(out_dir, "experiments_big_results.csv", sep="/");

nb_threads = 8

# args <- commandArgs(trailingOnly = TRUE)
# if (length(args) != 4) {
#   stop("usage: Rscript analysis.R <csv_file> <csv_file_2> <out_dir> <nb_threads>", call.=FALSE)
# }
# csv_file = args[1]
# csv_file_2 = args[2]
# out_dir = args[3]
# nb_threads = args[4]
# rm(args)

# =================================================

library(reshape2)
library(plyr)
library(ggplot2)



gg_color_hue <- function(n) {
  hues = seq(15, 375, length=n+1)
  hcl(h=hues, l=65, c=100)[1:n]
}
n = 3
p_colors = gg_color_hue(n)

# dev.new(width=4, height=4)
plot(1:n, pch=16, cex=2, col=p_colors)

# =================================================

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

get_error_layer <- function(df_melt_ci) {
  return (geom_errorbar(data=df_melt_ci, aes(y=mn, ymin=mn-se, ymax=mn+se)))
}


p_all <- scale_colour_manual("Sorting", values = p_colors,
                             labels = c("Sequential", "Parallel", "Libc"))
p_sq_pr <- scale_colour_manual("Sorting", values = p_colors[1:2],
                             labels = c("Sequential", "Parallel"))
p_sq_lc <- scale_colour_manual("Sorting", values = p_colors[c(1,3)],
                             labels = c("Sequential", "Libc"))

p_labs <- labs(title = "Measurements", y = "Time", x = "Size")

# =================================================


df_expr_medium <- read.csv(csv_file_medium, sep=" ", header=T, strip.white=TRUE)

# renaming for simplicity
names(df_expr_medium)[3:5] <- c("seq", "par", "blt")

# analysis the changes in execution time in relation to the size,
# for a specified number of thread in the parallel version
df_size <- df_expr_medium[df_expr_medium$threads == nb_threads,]
df_size <- df_size[df_size$size %% 200000 == 0,]


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
  geom_point(aes(y=blt, color="Libc")) +
  p_all + p_labs


# the Libc data (col. 5) can be removed to get higher resolution when drawing
df_size_melt <- melt(df_size[,c(1,3,4,5)], id="size")

df_size_melt_par_seq <- df_size_melt[df_size_melt$variable %in% c('seq', 'par'),]

# Same as above but in an elegant way
ggplot(df_size_melt, aes(size, value, colour=variable)) +
  geom_point() +
  p_all + p_labs

ggsave(paste(out_dir, "size_all.png", sep="/"), width = 5, height = 5)

# another plot
ggplot(df_size_melt, aes(size, value, fill=variable)) +
  geom_bar(stat="identity", position=position_dodge()) +
  p_all + p_labs


# Same as above but without the LibC
ggplot(df_size_melt_par_seq, aes(size, value, colour=variable)) +
  geom_point() +
  p_sq_pr + p_labs

ggsave(paste(out_dir, "size_par_seq.png", sep="/"), width = 5, height = 5)


df_size_melt_seq_ci <- compute_ci(df_size_melt, 'seq')
df_size_melt_par_ci <- compute_ci(df_size_melt, 'par')

# draw the CI
ggplot(df_size_melt_par_seq, aes(size, value, colour=variable)) +
  geom_point() +
  geom_errorbar(data=df_size_melt_par_ci, aes(y=mn, ymin=mn-se, ymax=mn+se)) +
  geom_errorbar(data=df_size_melt_seq_ci, aes(y=mn, ymin=mn-se, ymax=mn+se)) +
  p_sq_pr + p_labs

ggsave(paste(out_dir, "size_ci_points.png", sep="/"), width = 5, height = 5)


df_size_melt_par_seq_mn <- ddply(df_size_melt_par_seq, c("size", "variable"), summarize,
           mn = mean(value));

# draw the CI with the mean
ggplot(df_size_melt_par_seq, aes(size, value, colour=variable)) +
  geom_line(data=df_size_melt_par_seq_mn, aes(y=mn, colour=variable)) +
  geom_errorbar(data=df_size_melt_par_ci, aes(y=mn, ymin=mn-se, ymax=mn+se)) +
  geom_errorbar(data=df_size_melt_seq_ci, aes(y=mn, ymin=mn-se, ymax=mn+se)) +
  p_all + p_labs

ggsave(paste(out_dir, "size_ci_lines.png", sep="/"), width = 5, height = 5)


# draw the regression line
ggplot(df_size_melt, aes(size, value, colour=variable)) +
  geom_point() +
  geom_smooth(method = "lm", formula = y ~ x+I(x^2), se = FALSE, size = 1) +
  p_all + p_labs


ggsave(paste(out_dir, "size_regression.png", sep="/"), width = 5, height = 5)


# Same as above but without the LibC
ggplot(df_size_melt_par_seq, aes(size, value, colour=variable)) +
  geom_point() +
  geom_smooth(method = "lm", formula = y ~ x+I(x^2), se = FALSE, size = 1) +
  p_sq_pr + p_labs

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
df_threads <- df_expr_medium[,c(1,2,4)]

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
df_expr_big <- read.csv(csv_file_big, sep=" ", header=T, strip.white=TRUE)

# renaming for simplicity
names(df_expr_big)[3] <- c("par")

qplot (size, par, data=df_expr_big, color=factor(threads))

ggplot(df_expr_big, aes(size, par, color=factor(threads))) +
  geom_point() +
  geom_smooth(method = "lm", size = 1, formula = y ~ poly(x, 2), se = FALSE) +
  scale_colour_manual("#Threads", values = c("red", "purple"))

ggsave(paste(out_dir, "threads_4_8_10M.png", sep="/"), width = 5, height = 5)


# =================================================

df_expr_small <- read.csv(csv_file_small, sep=" ", header=T, strip.white=TRUE)

names(df_expr_small)[3:5] <- c("seq", "par", "blt")

df_expr_small <- df_expr_small[df_expr_small$threads == 2,]
df_expr_small <- df_expr_small[df_expr_small$size %% 60 == 0,]

df_lt1k_melt <- melt(df_expr_small[,c(1,3,4,5)], id="size")

# plot as bar plot
ggplot(df_lt1k_melt, aes(size, value, fill=variable)) +
  geom_bar(stat="identity", position=position_dodge()) +
  p_all + p_sq_lc + p_sq_pr + p_labs


df_small_melt_seq_ci <- compute_ci(df_lt1k_melt, 'seq')
df_small_melt_par_ci <- compute_ci(df_lt1k_melt, 'par')
df_small_melt_blt_ci <- compute_ci(df_lt1k_melt, 'blt')

# draw the CI

ggplot(df_lt1k_melt, aes(size, value, colour=variable)) +
  geom_point() +
  get_error_layer(df_small_melt_seq_ci) +
  get_error_layer(df_small_melt_blt_ci) +
  get_error_layer(df_small_melt_par_ci) +
  p_all + p_labs

df_lt1k_melt_wo_par <- df_lt1k_melt[df_lt1k_melt$variable %in% c('seq', 'blt'),]

ggplot(df_lt1k_melt_wo_par, aes(size, value, colour=variable)) +
  geom_point() +
  get_error_layer(df_small_melt_seq_ci) +
  get_error_layer(df_small_melt_blt_ci) +
  p_sq_lc + p_labs


#+++++++++++++++++++++++++
# Function to calculate the mean and the standard deviation
# for each group
#+++++++++++++++++++++++++
# data : a data frame
# varname : the name of a column containing the variable
#to be summariezed
# groupnames : vector of column names to be used as
# grouping variables
data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}

df3 <- data_summary(df_size_melt_par_seq, varname="value", 
                    groupnames=c("size", "variable"))

p <- ggplot(df3, aes(x=size, y=value, fill=variable)) + 
  geom_bar(stat="identity", position=position_dodge()) +
  geom_errorbar(aes(ymin=value-sd, ymax=value+sd), width=.4,
                position=position_dodge())


p + scale_fill_brewer(palette="Paired") + theme_minimal()
p + scale_fill_brewer(palette="Greens") + theme_minimal()
p + scale_fill_brewer(palette="Reds") + theme_minimal()
p + scale_fill_manual(values=c('red','blue','purple'))
