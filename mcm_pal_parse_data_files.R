setwd('C:/Users/Jeff/Documents/ducklow_lab/pal_mcm_paper')

## significant errors in source files, had to eliminate all comments in "clean" files using excel
## also added SEASON column to all, and DOY column to some.  DOY is not used in current analysis.

mcm_bp <- read.table('mcm_data/mcm_bp_clean.txt', sep = '\t', header = T)
mcm_pp <- read.table('mcm_data/mcm_pp_clean.txt', sep = '\t', header = T)
mcm_enum <- read.table('mcm_data/98_08_bact_enum_clean.txt', sep = '\t', header = T)
mcm_chl <- read.table('mcm_data/mcm_chl_clean.txt', sep = '\t', header = T)
mcm_doc <- read.table('mcm_data/mcm_doc_clean.txt', sep = '\t', header = T)

mcm_bp_recent <- read.table('mcm_data/mcm_bp_recent_clean.txt', sep = '\t', header = T)
mcm_pp_recent <- read.table('mcm_data/mcm_pp_recent_clean.txt', sep = '\t', header = T)
mcm_chl_recent <- read.table('mcm_data/mcm_chl_recent_clean.txt', sep = '\t', header = T)
mcm_doc_recent <- read.table('mcm_data/mcm_doc_recent_clean.txt', sep = '\t', header = T)

## enum is not available for recent MCM data
mcm_bp <- rbind(mcm_bp, mcm_bp_recent)
mcm_pp <- rbind(mcm_pp, mcm_pp_recent)
mcm_chl <- rbind(mcm_chl, mcm_chl_recent)
mcm_doc <- rbind(mcm_doc, mcm_doc_recent)

#### MCM interpolate across cast ####
## goal is interpolated BP/PP/Bact abund for comparison

library(pracma)

mcm_analyze <- seq(4, 70) # depths to be mcm_compared across parameters after interpolation
mcm_compare <- matrix(nrow = length(unique(mcm_bp$LIMNO_RUN)) * length(mcm_analyze), ncol = 9)

## populate matrix with limno run, depth, doy

mcm_doy <- vector(length = length(mcm_bp$DATE_TIME))

i <- 0
for(d in mcm_bp$DATE_TIME){
  i <- i + 1
  temp_date <- as.Date(d, "%m/%d/%y")
  mcm_doy[i] <- as.numeric(strftime(temp_date, format = "%j"))
}

mcm_bp$real_DOY <- mcm_doy

i <- 0
for(id in unique(as.character(mcm_bp$LIMNO_RUN))){
  doy <- mcm_bp$real_DOY[which(mcm_bp$LIMNO_RUN == id)[1]]
  for(a in mcm_analyze){
    i <- i + 1
    mcm_compare[i,1] <- id
    mcm_compare[i,2] <- a
    mcm_compare[i,9] <- doy
  }
}

## MCM tdr method switches to centrifuge in 2006, line 1401 on spreadsheet.  combine tdr columns
mcm_bp$TDR..nM.TDR.day.[1400:length(mcm_bp$TDR..nM.TDR.day.)] <- mcm_bp$TDR.CENT..nM.TDR.day.[1400:length(mcm_bp$TDR.CENT..nM.TDR.day.)]

# for testing run <- 'EB9394L1'

## MCM thymidine incorporation: mcm_compare[,3] - units are nM L-1 day-1

mcm_int_tdr <- matrix(ncol = length(unique(mcm_bp$LOCATION.CODE)), nrow = length(unique(mcm_bp$DATE_TIME)))
row.names(mcm_int_tdr) <- unique(mcm_bp$DATE_TIME)
colnames(mcm_int_tdr) <- unique(mcm_bp$LOCATION.CODE)

for(run in unique(mcm_bp$LIMNO_RUN)){
  select <- which(mcm_bp$LIMNO_RUN == run)
  station <- mcm_bp$LOCATION.CODE[select[1]]
  date <- mcm_bp$DATE_TIME[select[1]]
  if(length(select) > 1){
    temp_depth <- mcm_bp$DEPTH..m.[select]
    temp_bp <- mcm_bp$TDR..nM.TDR.day.[select]
    xi <- seq(min(temp_depth), max(temp_depth), length = length(select) * 10) # 10 to convert L to M
    temp <- interp1(x = temp_depth, y = temp_bp, xi = xi, method = 'spline')
    for(a in mcm_analyze){
      keep <- NA
      if(a > max(temp_depth)){
        keep <- NA
      }else{
        if(a < min(temp_depth)){
          keep <- NA
        }else{
          try(
            if(temp[which(abs(xi - a) == min(abs(xi - a)))][1] < 0){
              keep <- 0
            }else{
              keep <- temp[which(abs(xi - a) == min(abs(xi - a)))][1]
            },
            silent = T
          )
        }
      }
      if(is.na(keep) == F){
        print(paste('tdr', run, a, keep))
      }
      mcm_compare[intersect(which(mcm_compare[,1] == run), which(mcm_compare[,2] == a)),3] <- keep
    }
    sum_int <- sum(temp) 
    mcm_int_tdr[which(row.names(mcm_int_tdr) == date), which(colnames(mcm_int_tdr) == station)] <- sum_int
  }
}

## relabel row names with season

for(name in row.names(mcm_int_tdr)){
  row.names(mcm_int_tdr)[which(row.names(mcm_int_tdr) == name)] <- mcm_bp$SEASON[which(mcm_bp$DATE_TIME == name)][1]
}

## MCM leucine incorporation: mcm_compare[,4] - units are nM L-1 day-1

mcm_int_leu <- matrix(ncol = length(unique(mcm_bp$LOCATION.CODE)), nrow = length(unique(mcm_bp$DATE_TIME)))
row.names(mcm_int_leu) <- unique(mcm_bp$DATE_TIME)
colnames(mcm_int_leu) <- unique(mcm_bp$LOCATION.CODE)

for(run in unique(mcm_bp$LIMNO_RUN)){
  select <- intersect(which(mcm_bp$LIMNO_RUN == run), which(is.na(mcm_bp$LEU.CENT..nM.Leu.day.) == F))
  date <- mcm_bp$DATE_TIME[select[1]]
  station <- mcm_bp$LOCATION.CODE[select[1]]
  if(length(select) > 1){
    temp_depth <- mcm_bp$DEPTH..m.[select]
    temp_leu <- mcm_bp$LEU.CENT..nM.Leu.day.[select]
    xi <- seq(min(temp_depth), max(temp_depth), length = length(select) * 10) # 10 to convert L to M
    temp <- interp1(x = temp_depth, y = temp_leu, xi = xi, method = 'spline')
    for(a in mcm_analyze){
      if(a > max(temp_depth)){
        keep <- NA 
      }else{
        if(a < min(temp_depth)){
          keep <- NA
        }else{
          try(
            if(temp[which(abs(xi - a) == min(abs(xi - a)))][1] < 0){
              keep <- 0
            }else{
              keep <- temp[which(abs(xi - a) == min(abs(xi - a)))][1]
            },
            silent = T
          )
        }
      }
      if(is.na(keep) == F){
        print(paste('leu', run, a, keep))
      }
      mcm_compare[intersect(which(mcm_compare[,1] == run), which(mcm_compare[,2] == a)),4] <- keep
    }
    sum_int <- sum(temp) 
    mcm_int_leu[which(row.names(mcm_int_leu) == date), which(colnames(mcm_int_leu) == station)] <- sum_int
  }
}

## relabel row names with season

for(name in row.names(mcm_int_leu)){
  row.names(mcm_int_leu)[which(row.names(mcm_int_leu) == name)] <- mcm_bp$SEASON[which(mcm_bp$DATE_TIME == name)][1]
}

## MCM primary production: mcm_compare[,5] - units are ug C L-1 day-1

mcm_int_pp <- matrix(ncol = length(unique(mcm_pp$LOCATION.CODE)), nrow = length(unique(mcm_pp$DATE_TIME)))
row.names(mcm_int_pp) <- unique(mcm_pp$DATE_TIME)
colnames(mcm_int_pp) <- unique(mcm_pp$LOCATION.CODE)

for(run in unique(mcm_pp$LIMNO_RUN)){
  select <- which(mcm_pp$LIMNO_RUN == run)
  date <- mcm_pp$DATE_TIME[select[1]]
  station <- mcm_pp$LOCATION.CODE[select[1]]
  if(length(select) > 1){
    temp_depth <- mcm_pp$DEPTH..m.[select]
    temp_pp <- mcm_pp$PPR..ug.C.L.d.[select]
    xi <- seq(min(temp_depth), max(temp_depth), length = length(select) * 10) # 10 to convert L to M
    temp <- interp1(x = temp_depth, y = temp_pp, xi = xi, method = 'spline')
    for(a in mcm_analyze){
      if(a > max(temp_depth)){
        keep <- NA 
      }else{
        if(a < min(temp_depth)){
          keep <- NA
        }else{
          try(
            if(temp[which(abs(xi - a) == min(abs(xi - a)))][1] < 0){
              keep <- 0
            }else{
              keep <- temp[which(abs(xi - a) == min(abs(xi - a)))][1]
            },
            silent = T
          )
        }
      }
      if(is.na(keep) == F){
        print(paste('pp', run, a, keep))
      }
      mcm_compare[intersect(which(mcm_compare[,1] == run), which(mcm_compare[,2] == a)),5] <- keep
    }
    
    sum_int <- sum(temp) 
    mcm_int_pp[which(row.names(mcm_int_pp) == date), which(colnames(mcm_int_pp) == station)] <- sum_int
  }
}

## relabel row names with season

for(name in row.names(mcm_int_pp)){
  row.names(mcm_int_pp)[which(row.names(mcm_int_pp) == name)] <- mcm_pp$SEASON[which(mcm_pp$DATE_TIME == name)][1]
}

## MCM bacterial enumeration: mcm_compare[,6] - units are cells ml-1

# for testing run <- "HOR0001L3"

mcm_int_enum <- matrix(ncol = length(unique(mcm_enum$LOCATION.CODE)), nrow = length(unique(mcm_enum$DATE_TIME)))
row.names(mcm_int_enum) <- unique(mcm_enum$DATE_TIME)
colnames(mcm_int_enum) <- unique(mcm_enum$LOCATION.CODE)

for(run in unique(mcm_enum$LIMNO_RUN)){
  select <- intersect(which(mcm_enum$LIMNO_RUN == run), which(is.na(mcm_enum$TOTAL_BACT..cells.ml.) == F))
  date <- mcm_enum$DATE_TIME[select[1]]
  station <- mcm_enum$LOCATION.CODE[select[1]]
  if(length(select) > 1){
    temp_depth <- mcm_enum$DEPTH..m.[select]
    temp_enum <- mcm_enum$TOTAL_BACT..cells.ml.[select]
    xi <- seq(min(temp_depth), max(temp_depth), length = length(select) * 100) # 100 to convert ml to M
    temp <- interp1(x = temp_depth, y = temp_enum, xi = xi, method = 'spline')
    for(a in mcm_analyze){
      if(a > max(temp_depth)){
        keep <- NA 
      }else{
        if(a < min(temp_depth)){
          keep <- NA
        }else{
          try(
            if(temp[which(abs(xi - a) == min(abs(xi - a)))][1] < 0){
              keep <- NA
            }else{
              keep <- temp[which(abs(xi - a) == min(abs(xi - a)))][1]
            },
            silent = T
          )
        }
      }
      if(is.na(keep) == F){
        print(paste('enum', run, a, keep))
      }
      mcm_compare[intersect(which(mcm_compare[,1] == run), which(mcm_compare[,2] == a)),6] <- keep
    }
    sum_int <- sum(temp) # this is incorrect integration
    mcm_int_enum[which(row.names(mcm_int_enum) == date), which(colnames(mcm_int_enum) == station)] <- sum_int
  }
}

## relabel row names with season

for(name in row.names(mcm_int_enum)){
  row.names(mcm_int_enum)[which(row.names(mcm_int_enum) == name)] <- mcm_enum$SEASON[which(mcm_enum$DATE_TIME == name)][1]
}

mcm_compare[mcm_compare[,6] == 0,6] <- NA # can't have 0 enumeration values

## MCM chl: mcm_compare[,7] - units are ug C L-1

mcm_int_chl <- matrix(ncol = length(unique(mcm_chl$LOCATION.CODE)), nrow = length(unique(mcm_chl$DATE_TIME)))
row.names(mcm_int_chl) <- unique(mcm_chl$DATE_TIME)
colnames(mcm_int_chl) <- unique(mcm_chl$LOCATION.CODE)

for(run in unique(mcm_chl$LIMNO_RUN)){
  select <- which(mcm_chl$LIMNO_RUN == run)
  date <- mcm_chl$DATE_TIME[select[1]]
  station <- mcm_chl$LOCATION.CODE[select[1]]
  if(length(select) > 1){
    temp_depth <- mcm_chl$DEPTH..m.[select]
    temp_chl <- mcm_chl$CHL..ug.chl.a.L.[select]
    xi <- seq(min(temp_depth), max(temp_depth), length = length(select) * 10) # 10 to convert L to M
    temp <- interp1(x = temp_depth, y = temp_chl, xi = xi, method = 'spline')
    for(a in mcm_analyze){
      if(a > max(temp_depth)){
        keep <- NA 
      }else{
        if(a < min(temp_depth)){
          keep <- NA
        }else{
          try(
            if(temp[which(abs(xi - a) == min(abs(xi - a)))][1] < 0){
              keep <- 0
            }else{
              keep <- temp[which(abs(xi - a) == min(abs(xi - a)))][1]
            },
            silent = T
          )
        }
      }
      if(is.na(keep) == F){
        print(paste('chl', run, a, keep))
      }
      mcm_compare[intersect(which(mcm_compare[,1] == run), which(mcm_compare[,2] == a)),7] <- keep
    }
    
    sum_int <- sum(temp) 
    mcm_int_chl[which(row.names(mcm_int_chl) == date), which(colnames(mcm_int_chl) == station)] <- sum_int
  }
}

## relabel row names with season

for(name in row.names(mcm_int_chl)){
  row.names(mcm_int_chl)[which(row.names(mcm_int_chl) == name)] <- mcm_chl$SEASON[which(mcm_chl$DATE_TIME == name)][1]
}

## MCM doc: mcm_compare[,7] - units are mM C L-1

mcm_int_doc <- matrix(ncol = length(unique(mcm_doc$LOCATION.CODE)), nrow = length(unique(mcm_doc$DATE_TIME)))
row.names(mcm_int_doc) <- unique(mcm_doc$DATE_TIME)
colnames(mcm_int_doc) <- unique(mcm_doc$LOCATION.CODE)

for(run in unique(mcm_doc$LIMNO_RUN)){
  select <- which(mcm_doc$LIMNO_RUN == run)
  date <- mcm_doc$DATE_TIME[select[1]]
  station <- mcm_doc$LOCATION.CODE[select[1]]
  if(length(select) > 1){
    temp_depth <- mcm_doc$DEPTH..m.[select]
    temp_doc <- mcm_doc$DOC..mM.[select]
    xi <- seq(min(temp_depth), max(temp_depth), length = length(select) * 10) # 10 to convert L to M
    temp <- interp1(x = temp_depth, y = temp_doc, xi = xi, method = 'spline')
    for(a in mcm_analyze){
      if(a > max(temp_depth)){
        keep <- NA 
      }else{
        if(a < min(temp_depth)){
          keep <- NA
        }else{
          try(
            if(temp[which(abs(xi - a) == min(abs(xi - a)))][1] < 0){
              keep <- 0
            }else{
              keep <- temp[which(abs(xi - a) == min(abs(xi - a)))][1]
            },
            silent = T
          )
        }
      }
      if(is.na(keep) == F){
        print(paste('doc', run, a, keep))
      }
      mcm_compare[intersect(which(mcm_compare[,1] == run), which(mcm_compare[,2] == a)),8] <- keep
    }
    
    sum_int <- sum(temp) 
    mcm_int_doc[which(row.names(mcm_int_doc) == date), which(colnames(mcm_int_doc) == station)] <- sum_int
  }
}

## relabel row names with season

for(name in row.names(mcm_int_doc)){
  row.names(mcm_int_doc)[which(row.names(mcm_int_doc) == name)] <- mcm_doc$SEASON[which(mcm_doc$DATE_TIME == name)][1]
}

##########################################################################################################
############################### PAL interpolate across cast ##############################################

pal <- read.table('pal_data/pal_master_data.csv', sep = ',', header = T)
pal_doc <- read.table('pal_data/Dissolved Organic Carbon.csv', sep = ',', header = T)
pal_doc <- pal_doc[is.na(pal_doc$Depth..m.) == F,]

## generate doy column
pal_doy <- vector(length = length(pal$Date))

i <- 0
for(d in as.character(pal$Date)){
  temp_date <- NULL
  i <- i + 1
  if(length(grep('-', d)) == 1){
    temp_date <- as.Date(d, "%d-%b-%y")
  }else{
    if(length(grep('/', d)) == 1){
      temp_date <- as.Date(d, "%m/%d/%y")
      }else{
      temp_date <- as.Date(d, "%b %d %Y")
      }
    }

  pal_doy[i] <- as.numeric(strftime(temp_date, format = "%j"))
  print(pal_doy[i])
}

pal$doy <- pal_doy

## pal abundance data prior to year 2012 us in L, change to mL - other conversions take place later
## but you need to index by year for this
pal$Abundance[which(pal$Cruise != 2012)] <- pal$Abundance[which(pal$Cruise != 2012)] / 1000

## generate unique identifies for each cast
pal$Run <- paste0(pal$Cruise, '.', pal$Event)

## do the same for the doc dataset
#year_fun <- function(s) strsplit(s, "-")[[1]][1]
#pal_doc$year <- sapply(as.character(pal_doc$Datetime.GMT), year_fun)
pal_doc$Run <- paste0(pal_doc$Year, '.', pal_doc$Event)

#### NEED TO FIX #### pal_doc$Station <- pal_doc$Grid.Line.Intended + 0.01 * pal_doc$Grid.Station.Intended

## NA depths will break script downstream, and are useless anyway.  remove.
pal <- pal[which(is.na(pal$DepSM) == F),]
pal_analyze <- seq(0, 500, 10)
pal_compare <- matrix(nrow = (length(pal$Depth) * length(pal_analyze)), ncol = 9)

## populate with event,depth,doy

i <- 0
for(id in unique(as.character(pal$Run))){
  doy <- pal$doy[which(pal$Run == id)[1]]
  for(a in pal_analyze){
    i <- i + 1
    pal_compare[i,1] <- id
    pal_compare[i,2] <- a
    pal_compare[i,9] <- doy
  }
}

## Tdr: pal_compare[,3], pmol l-1 h-1

pal_int_tdr <- matrix(ncol = length(unique(pal$Cruise)), nrow = length(unique(pal$Station)))
colnames(pal_int_tdr) <- seq(min(as.numeric(unique(pal$Cruise))), max(as.numeric(unique(pal$Cruise))))
row.names(pal_int_tdr) <- unique(pal$Station)

for(run in unique(pal$Run)){
  select <- intersect(which(pal$Run == run), which(is.na(pal$Thymidine) == F))
  year <- pal$Cruise[select[1]]
  station <- pal$Station[select[1]]
  if(length(select) > 1){
    temp_depth <- pal$DepSM[select]
    temp_bp <- pal$Thymidine[select]
    xi <- seq(min(temp_depth), max(temp_depth), length = (max(temp_depth) - min(temp_depth)) * 10) # 10 to convert L to M
    temp <- interp1(x = temp_depth, y = temp_bp, xi = xi, method = 'spline')
    for(a in pal_analyze){
      keep <- NA
      if(a > max(temp_depth)){
        keep <- NA 
      }else{
        if(a < min(temp_depth)){
          keep <- NA
        }else{
          try(
            if(temp[which(abs(xi - a) == min(abs(xi - a)))][1] < 0){
              keep <- 0
            }else{
              keep <- temp[which(abs(xi - a) == min(abs(xi - a)))][1]
            },
            silent = T
          )
        }
      }
      if(is.na(keep) == F){
        print(paste('tdr', run, a, keep))
      }
      pal_compare[intersect(which(pal_compare[,1] == run), which(pal_compare[,2] == a)),3] <- keep
    }
    sum_int <- sum(temp) 
    pal_int_tdr[which(row.names(pal_int_tdr) == station), which(colnames(pal_int_tdr) == year)] <- sum_int
  }
}

pal_int_tdr <- (pal_int_tdr / 1000) * 24 # convert so consitent with MCM units

## Leu: pal_compare[,4], pmol l-1 h-1

pal_int_leu <- matrix(ncol = length(unique(pal$Cruise)), nrow = length(unique(pal$Station)))
colnames(pal_int_leu) <- seq(min(as.numeric(unique(pal$Cruise))), max(as.numeric(unique(pal$Cruise))))
row.names(pal_int_leu) <- unique(pal$Station)


for(run in unique(pal$Run)){
  select <- intersect(which(pal$Run == run), which(is.na(pal$Leucine) == F))
  year <- pal$Cruise[select[1]]
  station <- pal$Station[select[1]]
  #  temp_depth <- pal$DepSM[select]
  #  temp_bp <- pal$Thymidine[select]
  if(length(select) > 1){
    temp_depth <- pal$DepSM[select]
    temp_bp <- pal$Leucine[select]
    xi <- seq(min(temp_depth), max(temp_depth), length = (max(temp_depth) - min(temp_depth)) * 10) # 10 to convert L to M
    temp <- interp1(x = temp_depth, y = temp_bp, xi = xi, method = 'spline')
    for(a in pal_analyze){
      if(a > max(temp_depth)){
        keep <- NA 
      }else{
        if(a < min(temp_depth)){
          keep <- NA
        }else{
          try(
            if(temp[which(abs(xi - a) == min(abs(xi - a)))][1] < 0){
              keep <- 0
            }else{
              keep <- temp[which(abs(xi - a) == min(abs(xi - a)))][1]
            },
            silent = T
          )
        }
      }
      if(is.na(keep) == F){
        print(paste('leu', run, a, keep))
      }
      pal_compare[intersect(which(pal_compare[,1] == run), which(pal_compare[,2] == a)),4] <- keep
    }
    
    sum_int <- sum(temp) 
    pal_int_leu[which(row.names(pal_int_leu) == station), which(colnames(pal_int_leu) == year)] <- sum_int
  }
}

pal_int_leu <- (pal_int_leu / 1000) * 24 # convert so consistent with MCM units

## pp: pal_compare[,5], mg m-3 day-1

pal_int_pp <- matrix(ncol = length(unique(pal$Cruise)), nrow = length(unique(pal$Station)))
colnames(pal_int_pp) <- seq(min(as.numeric(unique(pal$Cruise))), max(as.numeric(unique(pal$Cruise))))
row.names(pal_int_pp) <- unique(pal$Station)

for(run in unique(pal$Run)){
  select <- intersect(which(pal$Run == run), which(is.na(pal$PP) == F))
  year <- pal$Cruise[select[1]]
  station <- pal$Station[select[1]]
  if(length(select) > 2){
    temp_depth <- pal$DepSM[select]
    temp_pp <- pal$PP[select]
    xi <- seq(min(temp_depth), max(temp_depth), length = (max(temp_depth) - min(temp_depth))) # already in m space
    temp <- interp1(x = temp_depth, y = temp_pp, xi = xi, method = 'spline')
    for(a in pal_analyze){
      if(a > max(temp_depth)){
        keep <- NA 
      }else{
        if(a < min(temp_depth)){
          keep <- NA
        }else{
          try(
            if(temp[which(abs(xi - a) == min(abs(xi - a)))][1] < 0){
              keep <- 0
            }else{
              keep <- temp[which(abs(xi - a) == min(abs(xi - a)))][1]
            },
            silent = T
          )
        }
      }
      if(is.na(keep) == F){
        print(paste('pp', run, a, keep))
      }
      pal_compare[intersect(which(pal_compare[,1] == run), which(pal_compare[,2] == a)),5] <- keep
    }
    
    sum_int <- sum(temp) 
    pal_int_pp[which(row.names(pal_int_pp) == station), which(colnames(pal_int_pp) == year)] <- sum_int
  }
}

## enum: pal_compare[,6], bact ml-1

pal_int_enum <- matrix(ncol = length(unique(pal$Cruise)), nrow = length(unique(pal$Station)))
colnames(pal_int_enum) <- seq(min(as.numeric(unique(pal$Cruise))), max(as.numeric(unique(pal$Cruise))))
row.names(pal_int_enum) <- unique(pal$Station)

for(run in unique(pal$Run)){
  select <- intersect(which(pal$Run == run), which(is.na(pal$Abundance) == F))
  year <- pal$Cruise[select[1]]
  station <- pal$Station[select[1]]
  if(length(select) > 3){
    temp_depth <- pal$DepSM[select]
    temp_bp <- pal$Abundance[select]
    xi <- seq(min(temp_depth), max(temp_depth), length = (max(temp_depth) - min(temp_depth)) * 100) # 100 to convert ml to M
    temp <- interp1(x = temp_depth, y = temp_bp, xi = xi, method = 'spline')
    for(a in pal_analyze){
      if(a > max(temp_depth)){
        keep <- NA 
      }else{
        if(a < min(temp_depth)){
          keep <- NA
        }else{
          try(
            if(temp[which(abs(xi - a) == min(abs(xi - a)))][1] < 0){
              keep <- NA
            }else{
              keep <- temp[which(abs(xi - a) == min(abs(xi - a)))][1]
            },
            silent = T
          )
        }
      }
      if(is.na(keep) == F){
        print(paste('enum', run, a, keep))
      }
      pal_compare[intersect(which(pal_compare[,1] == run), which(pal_compare[,2] == a)),6] <- keep
    }
    
    sum_int <- sum(temp) 
    pal_int_enum[which(row.names(pal_int_enum) == station), which(colnames(pal_int_enum) == year)] <- sum_int
  }
}
pal_compare[pal_compare[,6] == 0,6] <- NA # can't have 0 enumeration values

## chl: pal_compare[,7], mg m^-3

pal_int_chl <- matrix(ncol = length(unique(pal$Cruise)), nrow = length(unique(pal$Station)))
colnames(pal_int_chl) <- seq(min(as.numeric(unique(pal$Cruise))), max(as.numeric(unique(pal$Cruise))))
row.names(pal_int_chl) <- unique(pal$Station)

for(run in unique(pal$Run)){
  select <- intersect(which(pal$Run == run), which(is.na(pal$Chl) == F))
  year <- pal$Cruise[select[1]]
  station <- pal$Station[select[1]]
  if(length(select) > 3){
    temp_depth <- pal$DepSM[select]
    temp_bp <- pal$Chl[select]
    xi <- seq(min(temp_depth), max(temp_depth), length = (max(temp_depth) - min(temp_depth)) * 100) # 100 to convert ml to M
    temp <- interp1(x = temp_depth, y = temp_bp, xi = xi, method = 'spline')
    for(a in pal_analyze){
      if(a > max(temp_depth)){
        keep <- NA 
      }else{
        if(a < min(temp_depth)){
          keep <- NA
        }else{
          try(
            if(temp[which(abs(xi - a) == min(abs(xi - a)))][1] < 0){
              keep <- NA
            }else{
              keep <- temp[which(abs(xi - a) == min(abs(xi - a)))][1]
            },
            silent = T
          )
        }
      }
      if(is.na(keep) == F){
        print(paste('chl', run, a, keep))
      }
      pal_compare[intersect(which(pal_compare[,1] == run), which(pal_compare[,2] == a)),7] <- keep
    }
    
    sum_int <- sum(temp) 
    pal_int_chl[which(row.names(pal_int_chl) == station), which(colnames(pal_int_chl) == year)] <- sum_int
  }
}

## DOC: pal_compare[,8], umol L-1

for(run in unique(pal_doc$Run)){
  select <- intersect(which(pal_doc$Run == run), which(is.na(pal_doc$DOC..µmol.L.) == F))
  year <- pal_doc$year[select[1]]
  #  station <- pal_doc$Station[select[1]]
  if(length(select) > 3){
    temp_depth <- pal_doc$Depth..m.[select]
    temp_doc <- pal_doc$DOC..µmol.L.[select]
    xi <- seq(min(temp_depth), max(temp_depth), length = (max(temp_depth) - min(temp_depth)) * 10) # 10 to convert L to M
    temp <- interp1(x = temp_depth, y = temp_doc, xi = xi, method = 'spline')
    for(a in pal_analyze){
      if(a > max(temp_depth)){
        keep <- NA 
      }else{
        if(a < min(temp_depth)){
          keep <- NA
        }else{
          try(
            if(temp[which(abs(xi - a) == min(abs(xi - a)))][1] < 0){
              keep <- NA
            }else{
              keep <- temp[which(abs(xi - a) == min(abs(xi - a)))][1]
            },
            silent = T
          )
        }
      }
      if(is.na(keep) == F){
        print(paste('doc', run, a, keep))
      }
      pal_compare[intersect(which(pal_compare[,1] == run), which(pal_compare[,2] == a)),8] <- keep
    }
    
    sum_int <- sum(temp) 
  }
}

#### unit conversions - get PAL data into same units as MCM data ####

## convert pmol l-1 h-1 to nmol l-1 d-1 for PAL bp data
## divide by 1000, multiply by 24

pal_compare[,3] <- (as.numeric(pal_compare[,3]) / 1000) * 24
pal_compare[,4] <- (as.numeric(pal_compare[,4]) / 1000) * 24

## convery mg m-3 day-1 ug L-1 day-1 for PAL pp
## no conversion necessary

pal_compare[,5] <- (as.numeric(pal_compare[,5])) * 1

## convery umol L-1 DOC to mmol L-1 DOC
pal_compare[,8] <- (as.numeric(pal_compare[,8])) / 1000

#### write tables ####

write.table(mcm_int_tdr, 'mcm_int_tdr.txt', quote = F, sep = ',')
write.table(mcm_int_leu, 'mcm_int_leu.txt', quote = F, sep = ',')
write.table(mcm_int_pp, 'mcm_int_pp.txt', quote = F, sep = ',')
write.table(mcm_int_chl, 'mcm_int_chl.txt', quote = F, sep = ',')

write.table(pal_int_tdr, 'pal_int_tdr.txt', quote = F, sep = ',')
write.table(pal_int_leu, 'pal_int_leu.txt', quote = F, sep = ',')
write.table(pal_int_pp, 'pal_int_pp.txt', quote = F, sep = ',')
write.table(pal_int_enum, 'pal_int_enum.txt', quote = F, sep = ',')
write.table(pal_int_chl, 'pal_int_chl.txt', quote = F, sep = ',')

colnames(pal_compare) <- c('run', 'depth', 'thymidine', 'leucine', 'pp', 'abund', 'chl', 'doc', 'doy')
write.table(pal_compare, file = 'pal_comparison_data.txt', sep = ',', quote = F, row.names = F)

colnames(mcm_compare) <- c('run', 'depth', 'thymidine', 'leucine', 'pp', 'abund', 'chl', 'doc', 'doy')
write.table(mcm_compare, file = 'mcm_comparison_data.txt', sep = ',', quote = F, row.names = F)