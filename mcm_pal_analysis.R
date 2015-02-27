setwd('C:/Users/Jeff/Documents/ducklow_lab/pal_mcm_paper')

#### prepare data ####

pal_compare <- read.table('pal_comparison_data.txt', sep = ',', header = T,
                          colClasses = c('character', 'numeric', 'numeric', 'numeric', 'numeric', 'numeric', 'numeric'))
mcm_compare <- read.table('mcm_comparison_data.txt', sep = ',', header = T,
                          colClasses = c('character', 'numeric', 'numeric', 'numeric', 'numeric', 'numeric', 'numeric'))

## remove negative values and nas

mcm_bp_no_na <- data.frame(mcm_compare[,1], mcm_compare[,2], mcm_compare[,9], mcm_compare[,3], mcm_compare[,4], mcm_compare[,5], mcm_compare[,6], mcm_compare[,7], mcm_compare[,8], stringsAsFactors = F)
pal_bp_no_na <- data.frame(pal_compare[,1], pal_compare[,2], pal_compare[,9], pal_compare[,3], pal_compare[,4], pal_compare[,5], pal_compare[,6], pal_compare[,7], pal_compare[,8], stringsAsFactors = F)
colnames(mcm_bp_no_na) <- c('run', 'depth', 'doy', 'td', 'leu', 'pp', 'abund', 'chl', 'doc')
colnames(pal_bp_no_na) <- c('run', 'depth', 'doy', 'td', 'leu', 'pp', 'abund', 'chl', 'doc')

mcm_bp_no_na[mcm_bp_no_na <= 0] <- NA
pal_bp_no_na[pal_bp_no_na <= 0] <- NA

mcm_bp_no_na <- na.omit(mcm_bp_no_na)
pal_bp_no_na <- na.omit(pal_bp_no_na)

mcm_leu_pp_no_na <- data.frame(mcm_compare[,1], mcm_compare[,2], mcm_compare[,9], mcm_compare[,4], mcm_compare[,5], mcm_compare[,6], mcm_compare[,7], mcm_compare[,8], stringsAsFactors = F)
pal_leu_pp_no_na <- data.frame(pal_compare[,1], pal_compare[,2], pal_compare[,9], pal_compare[,4], pal_compare[,5], pal_compare[,6], pal_compare[,7], pal_compare[,8], stringsAsFactors = F)
colnames(mcm_leu_pp_no_na) <- c('run', 'depth', 'doy', 'leu', 'pp', 'abund', 'chl', 'doc')
colnames(pal_leu_pp_no_na) <- c('run', 'depth', 'doy', 'leu', 'pp', 'abund', 'chl', 'doc')

mcm_leu_pp_no_na[mcm_leu_pp_no_na <= 0] <- NA
pal_leu_pp_no_na[pal_leu_pp_no_na <= 0] <- NA

mcm_leu_pp_no_na <- na.omit(mcm_leu_pp_no_na)
pal_leu_pp_no_na <- na.omit(pal_leu_pp_no_na)

#### evaluate specific parameters ####

## look at DOC distribution

pal_doc_hist <- hist(as.numeric(pal_leu_pp_no_na$doc), breaks = 100)
plot(pal_doc_hist$count, log="y", type='h')

mcm_doc_hist <- hist(as.numeric(mcm_leu_pp_no_na$doc), breaks = 100)

## plot Tdr ~ Leu for both sites

plot(pal_bp_no_na$td ~ pal_bp_no_na$leu,
     log = 'xy',
     type = 'n',
     xlab = 'Leu',
     ylab = 'Td',
     xlim = c(1e-5,3),
     ylim = c(1e-6,10.0))

points(mcm_bp_no_na$td ~ mcm_bp_no_na$leu,
       col = 'red')

points(pal_bp_no_na$td ~ pal_bp_no_na$leu,
       col = 'blue')

mcm_leu_td_lm <- lm(mcm_bp_no_na$td ~ mcm_bp_no_na$leu) 
pal_leu_td_lm <- lm(pal_bp_no_na$td ~ pal_bp_no_na$leu) 

abline(mcm_leu_td_lm,
       col = 'red',
       untf = T)

abline(pal_leu_td_lm,
       col = 'blue',
       untf = T)

legend('bottomleft',
       pch = 1,
       col = c('blue', 'red'),
       legend = c('PAL', 'MCM'))

## plot Leu ~ pp for both sites

plot(pal_leu_pp_no_na$leu ~ pal_leu_pp_no_na$pp,
     log = 'xy',
     type = 'n',
     xlab = 'PP',
     ylab = 'Leu',
     xlim = c(1e-5,1e4),
     ylim = c(1e-4,10))

points(pal_leu_pp_no_na$leu ~ pal_leu_pp_no_na$pp,
       col = 'blue')

points(mcm_leu_pp_no_na$leu ~ mcm_leu_pp_no_na$pp,
       col = 'red')

pal_leu_pp_lm <- lm(pal_leu_pp_no_na$leu ~ pal_leu_pp_no_na$pp)
mcm_leu_pp_lm <- lm(mcm_leu_pp_no_na$leu ~ mcm_leu_pp_no_na$pp)

abline(pal_leu_pp_lm,
       col = 'blue',
       untf = T)

abline(mcm_leu_pp_lm,
       col = 'red',
       untf = T)

legend('bottomleft',
       pch = 1,
       col = c('blue', 'red'),
       legend = c('PAL', 'MCM'))

## plot Leu ~ doc

plot(mcm_leu_pp_no_na$leu ~ mcm_leu_pp_no_na$doc,
     log = 'xy',
     type = 'n',
     xlab = 'doc',
     ylab = 'Leu',
     #xlim = c(1e-5,1550),
     ylim = c(1e-3,10)
     )

#CHANGE!!
points(pal_leu_pp_no_na$leu ~ c(as.numeric(pal_leu_pp_no_na$doc) / 1000),
       col = 'blue')

points(mcm_leu_pp_no_na$leu ~ mcm_leu_pp_no_na$doc,
       col = 'red')

#CHANGE!!
pal_leu_doc_lm <- lm(pal_leu_pp_no_na$leu ~ pal_leu_pp_no_na$doc/1000)
mcm_leu_doc_lm <- lm(mcm_leu_pp_no_na$leu ~ mcm_leu_pp_no_na$doc)

abline(pal_leu_doc_lm,
       col = 'blue',
       untf = T)

abline(mcm_leu_doc_lm,
       col = 'red',
       untf = T)

legend('bottomleft',
       pch = 1,
       col = c('blue', 'red'),
       legend = c('PAL', 'MCM'))

#### identify overlaps  - central tendency method ####
## in leu ~ pp space which pal points co-locate with mcm mean and sd?

library(plotrix)

## mcm mean and sd

mcm_leu_mean_log <- mean(log(mcm_leu_pp_no_na$leu))
mcm_pp_mean_log <- mean(log(mcm_leu_pp_no_na$pp))

mcm_leu_sd_log <- sd(log(mcm_leu_pp_no_na$leu))
mcm_pp_sd_log <- sd(log(mcm_leu_pp_no_na$pp))

## pal mean and sd

pal_leu_mean_log <- mean(log(pal_leu_pp_no_na$leu))
pal_pp_mean_log <- mean(log(pal_leu_pp_no_na$pp))

pal_leu_sd_log <- sd(log(pal_leu_pp_no_na$leu))
pal_pp_sd_log <- sd(log(pal_leu_pp_no_na$pp))

## get the pal points that look like mcm points

pal_in_mcm <- pal_leu_pp_no_na[intersect(which(log(pal_leu_pp_no_na$leu) < mcm_leu_mean_log + mcm_leu_sd_log & log(pal_leu_pp_no_na$leu) > mcm_leu_mean_log - mcm_leu_sd_log),
                                         which(log(pal_leu_pp_no_na$pp) < mcm_pp_mean_log + mcm_pp_sd_log & log(pal_leu_pp_no_na$pp) > mcm_pp_mean_log - mcm_pp_sd_log)),]

## get the mcm points that look like pal points

mcm_in_pal <- mcm_leu_pp_no_na[intersect(which(log(mcm_leu_pp_no_na$leu) < pal_leu_mean_log + pal_leu_sd_log & log(mcm_leu_pp_no_na$leu) > pal_leu_mean_log - pal_leu_sd_log),
                                         which(log(mcm_leu_pp_no_na$pp) < pal_pp_mean_log + pal_pp_sd_log & log(mcm_leu_pp_no_na$pp) > pal_pp_mean_log - pal_pp_sd_log)),]

plot(log(pal_leu_pp_no_na$leu) ~ log(pal_leu_pp_no_na$pp),
     #log = 'xy',
     type = 'n',
     xlab = 'log(PP)',
     ylab = 'log(leu)',
     xlim = c(-10,10),
     ylim = c(-10,4))

abline(lm(log(mcm_leu_pp_no_na$leu) ~ log(mcm_leu_pp_no_na$pp)),
       untf = T,
       col = 'red')

abline(lm(log(pal_leu_pp_no_na$leu) ~ log(pal_leu_pp_no_na$pp)),
       untf = T,
       col = 'blue')

points(log(pal_leu_pp_no_na$leu) ~ log(pal_leu_pp_no_na$pp),
       col = 'blue')

points(log(mcm_leu_pp_no_na$leu) ~ log(mcm_leu_pp_no_na$pp),
       col = 'red')

points(mcm_pp_mean_log, mcm_leu_mean_log, pch = 19, col = 'orange')

draw.ellipse(mcm_pp_mean_log, mcm_leu_mean_log,
             mcm_pp_sd_log, mcm_leu_sd_log,
             border = 'orange',
             lwd = 2)

points(pal_pp_mean_log, pal_leu_mean_log, pch = 19, col = 'orange')

draw.ellipse(pal_pp_mean_log, pal_leu_mean_log,
             pal_pp_sd_log, pal_leu_sd_log,
             border = 'orange',
             lwd = 2)

legend('bottomleft',
       pch = c(1,1,19),
       col = c('blue', 'red', 'orange'),
       legend = c('PAL', 'MCM', 'mean'))

#### identify overlaps - slope method ####
## in leu ~ pp space which pal points fall along mcm linear model?

## get the pal points that look like mcm points

mcm_leu_pp_log_lm <- lm(log(mcm_leu_pp_no_na$leu) ~ log(mcm_leu_pp_no_na$pp))
mcm_leu_pp_log_m <- as.numeric(mcm_leu_pp_log_lm$coefficients[2])
mcm_leu_pp_log_b <- as.numeric(mcm_leu_pp_log_lm$coefficients[1])

pal_in_mcm <- matrix(ncol = length(pal_leu_pp_no_na[1,]), nrow = 0)
colnames(pal_in_mcm) <- colnames(pal_leu_pp_no_na)

for(i in 1:length(pal_leu_pp_no_na$leu)){
  leu <- log(pal_leu_pp_no_na$leu[i])
  pp <- log(pal_leu_pp_no_na$pp[i])
  predict <- mcm_leu_pp_log_m * pp + mcm_leu_pp_log_b
  if(abs(leu - predict) < abs(0.2 * predict)){
    pal_in_mcm <- rbind(pal_in_mcm, pal_leu_pp_no_na[i,])
  }
}

## get the mcm points that look like pal points

pal_leu_pp_log_lm <- lm(log(pal_leu_pp_no_na$leu) ~ log(pal_leu_pp_no_na$pp))
pal_leu_pp_log_m <- as.numeric(pal_leu_pp_log_lm$coefficients[2])
pal_leu_pp_log_b <- as.numeric(pal_leu_pp_log_lm$coefficients[1])

mcm_in_pal <- matrix(ncol = length(mcm_leu_pp_no_na[1,]), nrow = 0)
colnames(mcm_in_pal) <- colnames(mcm_leu_pp_no_na)

for(i in 1:length(mcm_leu_pp_no_na$leu)){
  leu <- log(mcm_leu_pp_no_na$leu[i])
  pp <- log(mcm_leu_pp_no_na$pp[i])
  predict <- pal_leu_pp_log_m * pp + pal_leu_pp_log_b
  if(abs(leu - predict) < abs(0.2 * predict)){
    mcm_in_pal <- rbind(mcm_in_pal, mcm_leu_pp_no_na[i,])
  }
}

## make plots

plot(log(pal_leu_pp_no_na$leu) ~ log(pal_leu_pp_no_na$pp),
     #log = 'xy',
     type = 'n',
     xlab = 'log(PP)',
     ylab = 'log(leu)',
     xlim = c(-10,4),
     ylim = c(-10,4))

abline(lm(log(mcm_leu_pp_no_na$leu) ~ log(mcm_leu_pp_no_na$pp)),
       untf = T,
       col = 'red')

abline(lm(log(pal_leu_pp_no_na$leu) ~ log(pal_leu_pp_no_na$pp)),
       untf = T,
       col = 'blue')

points(log(pal_leu_pp_no_na$leu) ~ log(pal_leu_pp_no_na$pp),
       col = 'blue')

points(log(mcm_leu_pp_no_na$leu) ~ log(mcm_leu_pp_no_na$pp),
       col = 'red')

points(log(mcm_in_pal$leu) ~ log(mcm_in_pal$pp),
       pch = 19,
       col = 'green')

points(log(pal_in_mcm$leu) ~ log(pal_in_mcm$pp),
       pch = 19,
       col = 'orange')

legend('topleft',
       pch = c(1,1,19,19),
       col = c('blue', 'red', 'orange', 'green'),
       legend = c('PAL', 'MCM', 'PAL w MCM trend', 'MCM w PAL trend'))

## what mcm dates, locations are pal-like?

select_limno_runs <- unique(mcm_bp[which(mcm_bp$LIMNO_RUN %in% mcm_in_pal$run),c(2,4,7)])

mcm_in_pal$loc <- vector(length = length(mcm_in_pal[,1]))
mcm_in_pal$date <- vector(length = length(mcm_in_pal[,1]))

for(i in seq(1:length(mcm_in_pal$run))){
  run <- mcm_in_pal$run[i]
  mcm_in_pal$loc[i] <- as.character(select_limno_runs[which(select_limno_runs[,1] == run), 2])
  mcm_in_pal$date[i] <- as.character(select_limno_runs[which(select_limno_runs[,1] == run), 3])
}

## what pal dates, locations are mcm-like?

pal_in_mcm$run <- as.character(pal_in_mcm$run)

select_pal_runs <- unique(pal[which(pal$Run %in% pal_in_mcm$run),c(24,4,5)])

pal_in_mcm$loc <- vector(length = length(pal_in_mcm[,1]))
pal_in_mcm$date <- vector(length = length(pal_in_mcm[,1]))

for(i in seq(1:length(pal_in_mcm$run))){
  run <- as.character(pal_in_mcm$run[i])
  pal_in_mcm$loc[i] <- as.character(select_pal_runs[which(select_pal_runs[,1] == run), 3])
  pal_in_mcm$date[i] <- as.character(select_pal_runs[which(select_pal_runs[,1] == run), 2])
}

********************************************************************************

## any coherence to date/depth for pal?

old_pal_dates <- as.data.frame(strsplit(pal_in_mcm$date, c('[/ -]')))

doys <- NULL

for(i in seq(1,length(old_pal_dates[1,]))){
  if(old_pal_dates[2,i] == 'Jan'){
    doy <- old_pal_dates[1,i]}else(doy <- old_pal_dates[2,i])
  doys <- append(doys, as.numeric(as.character(doy)))
}

pal_date_depth <- cbind(pal_in_mcm$depth, doys, paste(pal_in_mcm$depth, doys, sep = '-'))

pal_date_depth_tally <- as.data.frame(table(pal_date_depth[,3]))

strsplit(as.character(pal_date_depth_tally[,1]), '-')

******************* YOU ARE HERE, THE ABOVE IS SHIT ******************************

#### per cell comparisons ####

pal_leu_pp_no_na$spec_leu <- pal_leu_pp_no_na$leu/pal_leu_pp_no_na$abund
mcm_leu_pp_no_na$spec_leu <- mcm_leu_pp_no_na$leu/mcm_leu_pp_no_na$abund

plot(pal_leu_pp_no_na$spec_leu ~ pal_leu_pp_no_na$pp,
     log = 'xy',
     type = 'n',
     xlab = 'PP',
     ylab = 'Specific Leu',
     xlim = c(1e-5,2000),
     ylim = c(1e-10,0.00015))

points(pal_leu_pp_no_na$spec_leu ~ pal_leu_pp_no_na$pp,
       col = 'blue')

points(mcm_leu_pp_no_na$spec_leu ~ mcm_leu_pp_no_na$pp,
       col = 'red')

legend('bottomleft',
       pch = 1,
       col = c('blue', 'red'),
       legend = c('PAL', 'MCM'))

## leu ~ doc

#Change!!
plot(pal_leu_pp_no_na$spec_leu ~ c(as.numeric(pal_leu_pp_no_na$doc)/1000),
     log = 'xy',
     type = 'n',
     xlab = 'DOC',
     ylab = 'Specific Leu')

#Change!!
points(pal_leu_pp_no_na$spec_leu ~ c(as.numeric(pal_leu_pp_no_na$doc)/1000),
       col = 'blue')

points(mcm_leu_pp_no_na$spec_leu ~ mcm_leu_pp_no_na$doc,
       col = 'red')

legend('bottomleft',
       pch = 1,
       col = c('blue', 'red'),
       legend = c('PAL', 'MCM'))

lm_pal_spec_leu_doc <- lm(as.numeric(pal_leu_pp_no_na$spec_leu) ~ as.numeric(pal_leu_pp_no_na$doc))
abline(lm_pal_spec_leu_doc)

#### identify cell specific overlaps: central tendency method ####
## in leu ~ pp space which pal points co-locate with mcm mean and sd?

library(plotrix)

## mcm mean and sd

mcm_specific_leu_mean_log <- mean(log(mcm_leu_pp_no_na$spec_leu))
mcm_specific_pp_mean_log <- mean(log(mcm_leu_pp_no_na$pp))

mcm_specific_leu_sd_log <- sd(log(mcm_leu_pp_no_na$spec_leu))
mcm_specific_pp_sd_log <- sd(log(mcm_leu_pp_no_na$pp))

## pal mean and sd

pal_specific_leu_mean_log <- mean(log(pal_leu_pp_no_na$spec_leu))
pal_specific_pp_mean_log <- mean(log(pal_leu_pp_no_na$pp))

pal_specific_leu_sd_log <- sd(log(pal_leu_pp_no_na$spec_leu))
pal_specific_pp_sd_log <- sd(log(pal_leu_pp_no_na$pp))

## get the pal points that look like mcm points

pal_in_mcm_specific <- pal_leu_pp_no_na[intersect(which(log(pal_leu_pp_no_na$spec_leu) < mcm_specific_leu_mean_log + mcm_specific_leu_sd_log & log(pal_leu_pp_no_na$spec_leu) > mcm_specific_leu_mean_log - mcm_specific_leu_sd_log),
                                         which(log(pal_leu_pp_no_na$pp) < mcm_specific_pp_mean_log + mcm_specific_pp_sd_log & log(pal_leu_pp_no_na$pp) > mcm_specific_pp_mean_log - mcm_specific_pp_sd_log)),]

## get the mcm points that look like pal points

mcm_in_pal_specific <- mcm_leu_pp_no_na[intersect(which(log(mcm_leu_pp_no_na$spec_leu) < pal_specific_leu_mean_log + pal_specific_leu_sd_log & log(mcm_leu_pp_no_na$spec_leu) > pal_specific_leu_mean_log - pal_specific_leu_sd_log),
                                         which(log(mcm_leu_pp_no_na$pp) < pal_specific_pp_mean_log + pal_specific_pp_sd_log & log(mcm_leu_pp_no_na$pp) > pal_specific_pp_mean_log - pal_specific_pp_sd_log)),]

plot(log(pal_leu_pp_no_na$spec_leu) ~ log(pal_leu_pp_no_na$pp),
     #log = 'xy',
     type = 'n',
     xlab = 'log(PP)',
     ylab = 'log(leu)',
     xlim = c(-8,10),
     ylim = c(-30,0)
     )

abline(lm(log(mcm_leu_pp_no_na$spec_leu) ~ log(mcm_leu_pp_no_na$pp)),
       untf = T,
       col = 'red')

abline(lm(log(pal_leu_pp_no_na$spec_leu) ~ log(pal_leu_pp_no_na$pp)),
       untf = T,
       col = 'blue')

points(log(pal_leu_pp_no_na$spec_leu) ~ log(pal_leu_pp_no_na$pp),
       col = 'blue')

points(log(mcm_leu_pp_no_na$spec_leu) ~ log(mcm_leu_pp_no_na$pp),
       col = 'red')

points(mcm_specific_pp_mean_log, mcm_specific_leu_mean_log, pch = 19, col = 'orange')

draw.ellipse(mcm_specific_pp_mean_log, mcm_specific_leu_mean_log,
             mcm_specific_pp_sd_log, mcm_specific_leu_sd_log,
             border = 'orange',
             lwd = 2)

points(pal_specific_pp_mean_log, pal_specific_leu_mean_log, pch = 19, col = 'orange')

draw.ellipse(pal_specific_pp_mean_log, pal_specific_leu_mean_log,
             pal_specific_pp_sd_log, pal_specific_leu_sd_log,
             border = 'orange',
             lwd = 2)

legend('bottomleft',
       pch = c(1,1,19),
       col = c('blue', 'red', 'orange'),
       legend = c('PAL', 'MCM', 'mean'))

hist(pal_in_mcm_specific$depth)

#### bp climatology ####

plot(pal_compare$depth ~ pal_compare$leu,
     type = 'n',
     log = 'x',
     ylim = c(500, 0))

for(run in unique(pal_compare$run)){
  print(run)
  points(pal_compare$depth[which(pal_compare$run == run)] ~ pal_compare$leu[which(pal_compare$run == run)],
         pch = 19,
         cex = 0.5)
}

#### pal heatmaps ####

heat_color <- colorRampPalette(c('white', 'blue', 'red'))(100)

library(gplots)

### pp ###

pal_int_pp <- as.matrix(read.table('pal_int_pp.txt', sep = ',', header = T))

pal_int_pp[which(pal_int_pp < 0)] <- NA

heatmap.2(log(pal_int_pp),
          col = heat_color,
          Colv = NA,
          Rowv = NA,
          na.color = 'grey',
          symkey = F,
          margins= c(5, 10),
          sepcolor = 'white',
          trace = 'none',
          main = 'pal depth-int pp',
          keysize = 1,
          density.info = c('none'),
          cex.main = 0.5)

### leu ###

pal_int_leu <- as.matrix(read.table('pal_int_leu.txt', sep = ',', header = T))

pal_int_leu[which(pal_int_leu <= 0)] <- NA

heatmap.2(pal_int_leu,
          col = heat_color,
          Colv = NA,
          Rowv = NA,
          na.color = 'grey',
          symkey = F,
          margins= c(5, 10),
          sepcolor = 'white',
          trace = 'none',
          main = 'pal depth-int leu',
          keysize = 1,
          density.info = c('none'),
          cex.main = 0.5)

### chl ###

pal_int_chl <- as.matrix(read.table('pal_int_chl.txt', sep = ',', header = T))

pal_int_chl[which(pal_int_chl <= 0)] <- NA

heatmap.2(log(pal_int_chl),
          col = heat_color,
          Colv = NA,
          Rowv = NA,
          na.color = 'grey',
          symkey = F,
          margins= c(5, 10),
          sepcolor = 'white',
          trace = 'none',
          main = 'pal log(depth-int chl)',
          keysize = 1,
          density.info = c('none'),
          cex.main = 0.5)

#### mcm heatmaps ####

library(gplots)

mcm_int_tdr <- as.matrix(read.table('mcm_int_tdr.txt', sep = ',', header = T))
mcm_int_leu <- as.matrix(read.table('mcm_int_leu.txt', sep = ',', header = T))
mcm_int_pp <- as.matrix(read.table('mcm_int_pp.txt', sep = ',', header = T))
mcm_int_chl <- as.matrix(read.table('mcm_int_chl.txt', sep = ',', header = T))

heat_color <- colorRampPalette(c('white', 'blue', 'red'))(100)

## for each site need to integrate across season or take mean.  If integrate will need to come up
## with mechanism for labeling by day of year

### td ###

## label rows by year

#mcm_dates <- t(as.data.frame(strsplit(row.names(mcm_int_tdr), '/')))
#row.names(mcm_dates) <- NULL
#row.names(mcm_int_tdr) <- mcm_dates[,3]

## populate matrix with integrated values by site and year

mcm_mean_tdr <- matrix(nrow = length(colnames(mcm_int_tdr)), ncol = length(unique(row.names(mcm_int_tdr))))
colnames(mcm_mean_tdr) <- unique(row.names(mcm_int_tdr))
row.names(mcm_mean_tdr) <- colnames(mcm_int_tdr)

for(site in colnames(mcm_int_tdr)){
  for(season in unique(row.names(mcm_int_tdr))){
    select <- mcm_int_tdr[which(row.names(mcm_int_tdr) == season), which(colnames(mcm_int_tdr) == site)]
    if(length(select[is.na(select) == F] > 0)){
      temp <- mean(select[is.na(select) == F])
      mcm_mean_tdr[which(row.names(mcm_mean_tdr) == site), which(colnames(mcm_mean_tdr) == season)] <- temp
      print(paste(season, site, temp))
    }
  }
}

heatmap.2(log(mcm_mean_tdr),
          col = heat_color,
          Colv = NA,
          Rowv = NA,
          na.color = 'grey',
          symkey = F,
          margins= c(5, 10),
          sepcolor = 'white',
          trace = 'none',
          main = 'MCM depth-int mean log(tdr)',
          keysize = 1,
          density.info = c('none'),
          cex.main = 0.5,
          labCol = seq(1993,2011) + 1) # plus one to get in agreement with PAL seasons

### leu ###

## label rows by year

mcm_dates <- t(as.data.frame(strsplit(row.names(mcm_int_leu), '/')))
row.names(mcm_dates) <- NULL
row.names(mcm_int_leu) <- mcm_dates[,3]

## populate matrix with integrated values by site and year

mcm_mean_leu <- matrix(nrow = length(colnames(mcm_int_leu)), ncol = length(unique(row.names(mcm_int_leu))))
colnames(mcm_mean_leu) <- unique(row.names(mcm_int_leu))
row.names(mcm_mean_leu) <- colnames(mcm_int_leu)

for(site in colnames(mcm_int_leu)){
  for(year in unique(row.names(mcm_int_leu))){
    select <- mcm_int_leu[which(row.names(mcm_int_leu) == year), which(colnames(mcm_int_leu) == site)]
    if(length(select[is.na(select) == F] > 0)){
      temp <- mean(select[is.na(select) == F])
      mcm_mean_leu[which(row.names(mcm_mean_leu) == site), which(colnames(mcm_mean_leu) == year)] <- temp
      print(paste(year, site, temp))
    }
  }
}

heatmap.2(log(mcm_mean_leu),
          col = heat_color,
          Colv = NA,
          Rowv = NA,
          na.color = 'grey',
          symkey = F,
          margins= c(5, 10),
          sepcolor = 'white',
          trace = 'none',
          main = 'MCM depth-int mean log(leu)',
          keysize = 1,
          density.info = c('none'),
          cex.main = 0.5,
          labCol = seq(1993,2011) + 1)

### pp ###

## label rows by year

# mcm_dates <- t(as.data.frame(strsplit(row.names(mcm_int_pp), '/')))
# row.names(mcm_dates) <- NULL
# row.names(mcm_int_pp) <- mcm_dates[,3]

## populate matrix with integrated values by site and year

mcm_mean_pp <- matrix(nrow = length(colnames(mcm_int_pp)), ncol = length(unique(row.names(mcm_int_pp))))
colnames(mcm_mean_pp) <- unique(row.names(mcm_int_pp))
row.names(mcm_mean_pp) <- colnames(mcm_int_pp)

for(site in colnames(mcm_int_pp)){
  for(year in unique(row.names(mcm_int_pp))){
    select <- mcm_int_pp[which(row.names(mcm_int_pp) == year), which(colnames(mcm_int_pp) == site)]
    if(length(select[is.na(select) == F] > 0)){
      temp <- mean(select[is.na(select) == F])
      mcm_mean_pp[which(row.names(mcm_mean_pp) == site), which(colnames(mcm_mean_pp) == year)] <- temp
      print(paste(year, site, temp))
    }
  }
}

heatmap.2(mcm_mean_pp,,
          col = heat_color,
          Colv = NA,
          Rowv = NA,
          na.color = 'grey',
          symkey = F,
          margins= c(5, 10),
          sepcolor = 'white',
          trace = 'none',
          main = 'MCM depth-int mean pp',
          keysize = 1,
          density.info = c('none'),
          cex.main = 0.5,
          labCol = seq(1993,2011) + 1)

### chl ###

mcm_mean_chl <- matrix(nrow = length(colnames(mcm_int_chl)), ncol = length(unique(row.names(mcm_int_chl))))
colnames(mcm_mean_chl) <- unique(row.names(mcm_int_chl))
row.names(mcm_mean_chl) <- colnames(mcm_int_chl)

for(site in colnames(mcm_int_chl)){
  for(year in unique(row.names(mcm_int_chl))){
    select <- mcm_int_chl[which(row.names(mcm_int_chl) == year), which(colnames(mcm_int_chl) == site)]
    if(length(select[is.na(select) == F] > 0)){
      temp <- mean(select[is.na(select) == F])
      mcm_mean_chl[which(row.names(mcm_mean_chl) == site), which(colnames(mcm_mean_chl) == year)] <- temp
      print(paste(year, site, temp))
    }
  }
}

heatmap.2(mcm_mean_chl,,
          col = heat_color,
          Colv = NA,
          Rowv = NA,
          na.color = 'grey',
          symkey = F,
          margins= c(5, 10),
          sepcolor = 'white',
          trace = 'none',
          main = 'MCM depth-int mean chl',
          keysize = 1,
          density.info = c('none'),
          cex.main = 0.5,
          labCol = seq(1993,2011) + 1)

## check for correlation between pal and mcm years.  no correlation.

plot(apply(mcm_mean_pp, 2, mean)[10:19] ~ apply(pal_int_pp, 2, mean, na.rm = T))

#### NMDS analysis ####

## not currently working, need better normalization

library(vegan)

combined_leu_pp_no_na <- as.matrix(rbind(pal_leu_pp_no_na, mcm_leu_pp_no_na))

row.names(combined_leu_pp_no_na) <- NULL

combined_leu_pp_no_na <- combined_leu_pp_no_na[,-c(1:3)]

combined_leu_pp_no_na <- apply(combined_leu_pp_no_na, 2, as.numeric)

combined_leu_pp_no_na_std <- combined_leu_pp_no_na
combined_leu_pp_no_na_std[,1] <- log(combined_leu_pp_no_na_std[,1])
combined_leu_pp_no_na_std[,2] <- log(combined_leu_pp_no_na_std[,2])
combined_leu_pp_no_na_std[,3] <- sqrt(combined_leu_pp_no_na_std[,3])

combined_leu_pp_no_na_std <- decostand(combined_leu_pp_no_na_std, method = 'standardize', margin = 2)
#combined_leu_pp_no_na_std <- decostand(combined_leu_pp_no_na_std, method = 'range', margin = 2)

small = sample(c(1:length(combined_leu_pp_no_na_std[,1])), 200)

combined_leu_pp_no_na_mds <- metaMDS(combined_leu_pp_no_na_std, distance = 'euclidean', k = 2, autotransform = F)

plot(combined_leu_pp_no_na_mds$points[,1],
     combined_leu_pp_no_na_mds$points[,2],
     xlim = c(-0.1,0.1))


#### PCA ####

combined_leu_pp_no_na_pca <- prcomp(combined_leu_pp_no_na_std)

plot(combined_leu_pp_no_na_pca$x,
     type = 'n',
     ylim = c(-10,10),
     xlim = c(-10,10))

points(combined_leu_pp_no_na_pca$x[1:length(pal_leu_pp_no_na[,1]),1],
       combined_leu_pp_no_na_pca$x[1:length(pal_leu_pp_no_na[,1]),2],
       col = 'blue')

points(combined_leu_pp_no_na_pca$x[length(pal_leu_pp_no_na[,1]):length(combined_leu_pp_no_na[,1]),1],
       combined_leu_pp_no_na_pca$x[length(pal_leu_pp_no_na[,1]):length(combined_leu_pp_no_na[,1]),2],
       col = 'red')

loadings <- combined_leu_pp_no_na_pca$rotation

added_loadings <- abs(loadings[,1]) + abs(loadings[,2])

ordered_loadings <- loadings[order(added_loadings, decreasing = T),]

arrows(0,
       0,
       loadings[,1] * 10,
       loadings[,2] * 10
)

text(loadings[,1] * 11,
     loadings[,2] * 11,
     labels = c(row.names(loadings)))



