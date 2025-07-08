# R script to import, format, and graph CO2 fluxes from LI-840A
# Use this script to visualise all the records from flux log files. 

# all files must follow the naming convention: 
# site_samplenumber_postwaterh_fluxrep.txt
# eg: H_T1_0h.txt or H_T1_0h_1.txt for Holloman sample T1 0h after water addition and the first flux file (no flux rep number) and second flux file (rep: _1)

# All txt log files must be in the same folder from the sampling date yyyymmdd

# Folders needed: 
# Data/CFluxes: raw data from licor organised by date
# CalcuatedFluxes: save the calculated flux slopes and R2 here, Can append and add to file or save by date or by batch
# FluxFitFigures: save a pdf of the flux fits graph by sampleID with equations and r2

### To Do: 
# flux_calc function calculates slopes but need to refine time interval
# some files need custom time intervals for flux calculation

# load libraries
library(tidyr)
library(dplyr)
library(ggplot2)
library(lubridate)
library(stringr)
library(plyr)
library(ggpmisc)

# list all .txt flux files logged with LI-840 and list files in date folder
# include file path in name
flux.files1 <- list.files("./Data/CFluxes/20250707/", full.names=TRUE)
flux.files2<- list.files("./Data/CFluxes/20250708/", full.names=TRUE)
#flux.files3 <- list.files("./Data/CFluxes/20250703/", full.names=TRUE)

# create column names that are easier to type
# Original:
# [1] "Date.Y.M.D."         "Time.H.M.S."         "CO2.ppm."           
# [4] "H2O.ppt."            "H2O.C."              "Cell_Temperature.C."
# [7] "Cell_Pressure.kPa."  "CO2_Absorption"      "H2O_Absorption"   
flux.file.colnames <- c("date","time","co2","h2o.ppt","h2o.c","cell_temp","cell_press","co2_absorption","h2o_absorption")

# function that reads files, uses the filename to create sample ID and adds and index count to rows for each file
read_txt_colname <- function(colname){
  ret <- read.table(colname, sep=" ", skip=2, header=FALSE, stringsAsFactors=FALSE,col.names = flux.file.colnames)
  ret$ind_count <- seq(from=1,to=nrow(ret),by=1) # create a count that runs the length of each measurement
  obj_name <- tools::file_path_sans_ext(basename(colname))
  ret$sampleID <- obj_name #EDIT
  ret
}

# use the read_txt_colnames function to read and combine all files from the flux.files lists
data1 <- ldply(flux.files1, read_txt_colname)
data2 <- ldply(flux.files2, read_txt_colname)
#data3 <- ldply(flux.files3, read_txt_colname)

# combine data from multiple day folders (list as many days as you want to combine)
data <- rbind(data1, data2)

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
                       names = c("site","sample_num","post_water_h","flux_rep"),
                       too_few="align_start", too_many="merge", cols_remove=FALSE)

# for post_water, create a numeric variable
data <- data %>%
  separate_wider_delim(post_water_h, delim = "h",
                       names = "post_water",too_many="drop", cols_remove=FALSE) %>%
  mutate(post_water = as.numeric(post_water))

# make the flux_rep 0 when NA
data <- data %>%
  mutate(flux_rep = case_when (is.na(flux_rep) ~ 0L,
                                       .default = TRUE))

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
  geom_line() +
  theme(legend.position="none")

# Almost all test samples have linear CO2 increase section between time 50-225
# for first round of 10g measurements timing was shorter, select 50-190
# exclude potting soil
data %>%
  filter(!(str_detect(sampleID, "pottingsoil")))%>%
  ggplot(.,aes(ind_count,co2,color=sampleID))+
  geom_point(size=0.5)+
  geom_line()+
  xlim(c(50,190))

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


# visualise potential line fits: 
# grup by sample ID, hour, rep
# show only for 50-230 time points (linear section except P_T1)
# use stat_poly from ggpmisc::
# https://cran.r-project.org/web/packages/ggpmisc/vignettes/model-based-annotations.html
formula <- y ~ x

data %>%
  filter(!(str_detect(sampleID, "pottingsoil")))%>%
  filter(ind_count> 49 & ind_count < 191) %>% 
  #filter(post_water==8 | post_water==24) %>%
  ggplot(.,aes(ind_count,co2,color=sample_num, shape=factor(flux_rep)))+
  geom_point(size=0.7)+
  #stat_poly_line(formula = formula, linewidth=0.3) +
  #stat_poly_eq(formula = formula)+
  #xlim(c(50,225))+
  facet_grid(site~post_water,scales="free_y")

# graph line fits for each sample ID in a separate facet to see more clearly
graph.flux.fits <- data %>%
  filter(!(str_detect(sampleID, "pottingsoil")))%>%
  filter(ind_count> 49 & ind_count < 191) %>% 
  filter(post_water==8 | post_water==24) %>%
  ggplot(.,aes(ind_count,co2,color=sample_num, shape=factor(flux_rep)))+
  geom_point(size=0.7)+
  stat_poly_line(formula = formula, linewidth=0.3) +
  stat_poly_eq(formula = formula, size=2,label.x = "right", label.y = "bottom") +
  stat_poly_eq(mapping = use_label("eq"),formula = formula, size=2,label.x = "left", label.y = "top")+
  facet_wrap(sampleID~.)

graph.flux.fits

# option to save the flux fits by sample ID graph for reference
ggsave("GraphFluxFits_10gBatch1.pdf",
       plot=graph.flux.fits,
       path="./FluxFitFigures",device="pdf",
       width=12, height=7, units="in")

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
  filter(ind_count> 49 & ind_count < 191)

fluxes=ddply(data.flux.sub, .(date, samplesize,sampleID, site, sample_num,post_water,flux_rep),  flux_calc)


# graph slopes by sample ID
ggplot(fluxes, aes(sampleID,Flux.slope))+
  geom_point()

# graph r2 by post water, sample number, and site
ggplot(fluxes, aes(post_water,Flux.r2,color=sample_num))+
  geom_point()+
  geom_hline(yintercept=0.8)+
  facet_grid(site~.)

# graph slopes by post water, sample number, and site
ggplot(fluxes, aes(post_water,Flux.slope,color=sample_num))+
  geom_point()+
  geom_hline(yintercept=0)+
  facet_grid(site~.)


##### TO SAVE CALCULATED FLUXES ######

# Add code to export fluxes as csv file
# each time this is run, it should add to the existing file
write.table(fluxes, file="./CalculatedFluxes/FluxSlopes.csv",append=TRUE,sep=",",dec=".",row.names=FALSE)

# to save a new file by date or by group of samples,
# chance the EDIT for FluxSlopes to either:
# the date: yyyymmdd
# the batch: eg: 10gBatch1
write.table(fluxes, file="./CalculatedFluxes/CalculatedFluxes_EDIT.csv",append=TRUE,sep=",",dec=".",row.names=FALSE)

