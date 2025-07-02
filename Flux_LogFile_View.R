# R script to import, format, and graph CO2 fluxes from LI-840A
# Use this script to visualise all the records from flux log files. 

# all files must follow the naming convention: 
# site_samplenumber_postwaterh_fluxrep.txt
# eg: H_T1_0h.txt or H_T1_0h_1.txt for Holloman sample T1 0h after water addition and the first flux file (no flux rep number) and second flux file (rep: _1)

# All txt log files must be in the same folder from the sampling date yyyymmdd

# load libraries
library(tidyr)
library(dplyr)
library(ggplot2)
library(lubridate)
library(stringr)
library(data.table)
library(plyr)
library(cowplot)

# list all .txt flux files logged with LI-840 and list files in date folder
# include file path in name
flux.files <- list.files("./Data/CFluxes/20250701/", full.names=TRUE)

# create column names that are easier to type
# Original:
# [1] "Date.Y.M.D."         "Time.H.M.S."         "CO2.ppm."           
# [4] "H2O.ppt."            "H2O.C."              "Cell_Temperature.C."
# [7] "Cell_Pressure.kPa."  "CO2_Absorption"      "H2O_Absorption"   
flux.file.colnames <- c("date","time","co2","h2o.ppt","h2o.c","cell_temp","cell_press","co2_absorption","h2o_absorption")

read_txt_colname <- function(colname){
  ret <- read.table(colname, sep=" ", skip=2, header=FALSE, stringsAsFactors=FALSE,col.names = flux.file.colnames)
  ret$ind_count <- seq(from=1,to=nrow(ret),by=1) # create a count that runs the length of each measurement
  obj_name <- tools::file_path_sans_ext(basename(colname))
  ret$sampleID <- obj_name #EDIT
  ret
}

data <- ldply(flux.files, read_txt_colname)

# format timestamp of data
data <- data %>%
  mutate(datetime = ymd_hms(paste(date,time, sep=" ")))

# add sample size
data <- data %>%
  mutate(samplesize="10g")

# use the sample ID to create a site label,
# sample number (replicate of jars),
# post_water (hours since water addition),
# flux rep (redo of a flux measurement)
data <- data %>%
  separate_wider_delim(sampleID, delim = "_",
                       names = c("site","sample_num","post_water","flux_rep"),
                       too_few="align_start", too_many="merge", cols_remove=FALSE)


# graph the data
# All CO2 fluxes by sample ID
data %>%
  ggplot(.,aes(ind_count,co2))+
  geom_point(size=0.5)+
  geom_line()+
  facet_wrap(sampleID~.)

# All CO2 fluxes by sample ID
# all in same graph but with colors
data %>%
  ggplot(.,aes(ind_count,co2,color=sampleID))+
  geom_point(size=0.5)+
  geom_line()

# Almost all samples have linear CO2 increase section between time 50-225
# exclude potting soil
data %>%
  filter(!(str_detect(sampleID, "pottingsoil")))%>%
  ggplot(.,aes(ind_count,co2,color=sampleID))+
  geom_point(size=0.5)+
  geom_line()+
  xlim(c(50,225))

# exclude potting soil
# grup by sample ID, hour, rep
data %>%
  filter(!(str_detect(sampleID, "pottingsoil")))%>%
  ggplot(.,aes(ind_count,co2,color=sample_num, shape=factor(flux_rep)))+
  geom_point(size=0.5)+
  geom_line()+
 # xlim(c(50,250))+
  facet_grid(site~post_water)

# exclude potting soil
# grup by sample ID, hour, rep
# show only for 50-230 time points (linear section except P_T1)
data %>%
  filter(!(str_detect(sampleID, "pottingsoil")))%>%
  ggplot(.,aes(ind_count,co2,color=sample_num, shape=factor(flux_rep)))+
  geom_point(size=0.5)+
  geom_line()+
  xlim(c(50,225))+
  facet_grid(site~post_water)


# Function to calculate slopes
# calculate keeling slopes and intercepts
flux_calc = function(Z) 
{
  Flux.slope=coef(summary(lm(co2~ind_count, Z)))[2];
  Flux.r2=summary(lm(co2~ind_count, Z))$r.squared;
  #Time=mean(Z$count);
  #Plot=unique(Z$Plot)
  return (cbind(Flux.slope, Flux.r2))
}

# subset only data and time period for flux calculation
data.flux.sub <- data %>%
  filter(!(str_detect(sampleID, "pottingsoil")))%>%
  filter(ind_count> 49 & ind_count < 250)

fluxes=ddply(data.flux.sub, .(sampleID, site, sample_num,post_water,flux_rep),  flux_calc)

# graph slopes by sample ID
ggplot(fluxes, aes(sampleID,Flux.slope))+
  geom_point()

# graph slopes by post water, sample number, and site
ggplot(fluxes, aes(post_water,Flux.slope,color=sample_num))+
  geom_point()+
  facet_grid(site~.)

# graph r2 by post water, sample number, and site
ggplot(fluxes, aes(post_water,Flux.r2,color=sample_num))+
  geom_point()+
  geom_hline(yintercept=0.8)+
  facet_grid(site~.)
