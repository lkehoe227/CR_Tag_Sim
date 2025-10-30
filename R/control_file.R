rm(list=ls())
options(warn=1)

library('zoo')
library('RMark')
library('snowfall')
library('parallel')
library('PBSmodelling')
library('snow')
library('foreach')
library('doSNOW')
library('tcltk2')
library('matrixStats')

#set working directory
setwd("~/GitHub/CR_Tag_Sim")
dir.base = getwd()
dir.out = paste0(dir.base,'/output')

#assign a label for each scenario 
runname = 'test' 


####################################################################################################
# inputs
####################################################################################################
#Switches
#turn on simulation  (1=on, 0=off)
run_simulation = 1
#turn on age structure  (1=on, 0=off) {switches in init_pop}
age_on = 1
#turn on Time variability  (1=on, 0=off) {switches in run_sim}
time_on = 1
#Turn on estimations (1=on, 0=off)
est_on = 1
#estimation type # 'pre-post' "july" "high-low" "may-dec" "june-dec"
est_type = 'pre-post' #july is the only one for 2 and 3 year , july_yr2 includes the first july in high effort season
#type of effort and cpue data input
eff_data = 'srfs' #can be "mrip" or "srfs"
tag_release = 'monthly' #monthly (even split in tag season); fixed, set number of tags by month
prev_tag_seas = c(0) #when tags were previously deployed
prev_tag = 0 #total tags deployed prior to this study, not included in nmarks

#setup model dimensions-----------------------------------------------------------------------------
nyrs = 1
ntsteps = nyrs*12 

maxage = 20
ages.m  = round(seq(1,maxage+11/12,1/12),2)


source('R/init_population.R') #sources data and BAM inputs
source('R/run_sim.R') 
source('R/run_est.R')
source('R/tag_functions.R')



#tagging parameters---------------------------------------------------------------------------------
handle.mort = .234 #handling mortality (release mortality)
Nmarks=c(seq(500,3000,500)) #number of tags
rr=.95#reporting rate (decimal)
p.trips = 0.02 #Percentage of participating trips in ibgt,  
season = c(3,4) #month(s) of tag releases
wd = dir.base
timestep = 'monthly' #writes data in either annual or monthly timestep for estimation 
tag_type = 'conv'#type of tagging (Conventional(conv) or Genetic(gene) ) 
tag_after_cap = "clip" #this can be clip for a tag removed after capture or "keep" for the tag stays in the fish


#number of simulation runs per number of tags 
Ntrials = 10


if (run_simulation == 1) {
#makes run wd
dir.run=paste0(dir.out,"/",runname)

if(!dir.exists(dir.run)) dir.create(dir.run)
#creates run list
if(tag_type == "conv"){
  run.list = data.frame(marks = rep(Nmarks,each = Ntrials))
  run.list$trial = rep(1:Ntrials, length(Nmarks))
  run.list$ntimes = ntsteps
  run.list$run.num <- 1:nrow(run.list)
  run.list$fold = paste0("Nmarks_",run.list$marks,'/run_',run.list$trial)
  run.list$dir = paste0(dir.run,"/Nmarks_",run.list$marks)
  run.list$full.dir = paste0(run.list$dir,'/run_',run.list$trial)
}
if(tag_type == "gene"){
  run.list = data.frame(marks = rep(p.trips,each = Ntrials))
  run.list$trial = rep(1:Ntrials, length(p.trips))
  run.list$ntimes = ntsteps
  run.list$run.num <- 1:nrow(run.list)
  run.list$fold = paste0("p.trips_",run.list$marks,'/run_',run.list$trial)
  run.list$dir = paste0(dir.run,"/p.trips_",run.list$marks)
  run.list$full.dir = paste0(run.list$dir,'/run_',run.list$trial)
}



#clears and creates folders 
unlink(list.dirs(dir.run,recursive=F),recursive=T) #clears folders and deletes contents!!!!!
lapply(run.list$dir,function(x) if(!dir.exists(x)) dir.create(x))
lapply(paste0(run.list$dir,'/run_',run.list$trial),function(x) if(!dir.exists(x)) dir.create(x))
#make a csv run.list
write.csv(run.list,paste0(dir.run,'/',runname,'-run_list.csv'),row.names = F)

##################################################
#Simulation 
##################################################
#setup parallel clusters
cl <- makeSOCKcluster(detectCores()-1)
clusterExport(cl, c("wd","est_on","run.list","run.sim","run_estimation",'time_on','age_on')) #LK,'time_on','age_on'
registerDoSNOW(cl)

pbar <- winProgressBar("Running Simulations",label=paste0("Simulation 0 of ",nrow(run.list)),max=100)
prog <- function(n) setWinProgressBar(pbar,(n/nrow(run.list)*100),label=paste("Simulation Run", n,"of", nrow(run.list),"Completed"))
opts <- list(progress=prog)

#first do the simulations-----------------------------------------------------------
est_on <<-est_on
time_on <<- time_on 
age_on <<- age_on 
source('R/run_sim.R')
source('R/run_est.R')

mod <- foreach(d=1:length(run.list$marks),.errorhandling='pass',.packages=c('snowfall','parallel','tcltk2','RMark')
) %dopar% {
  run.sim(run.list,d) 
}
close(pbar)

#check for missing runs - checks if pdf file for sim plots is written
done.sims = paste0(basename(dirname(dirname(list.files(dir.run, recursive=TRUE, full.names=TRUE,pattern="\\.pdf$")))),'/',basename((dirname(list.files(dir.run, recursive=TRUE, full.names=TRUE,pattern="\\.pdf$")))))
miss.sims = run.list[which(!paste0(run.list$fold)%in%done.sims),]
print(paste(length(miss.sims$marks),"simulations did not write to file"))

pbar <- winProgressBar("Running Simulations",label=paste0("Simulation 0 of ",nrow(run.list)),max=100)
prog <- function(n) setWinProgressBar(pbar,(n/nrow(run.list)*100),label=paste("Simulation Run", n,"of", nrow(run.list),"Completed"))
opts <- list(progress=prog)


#redo missing simulation runs, repeat until no more
if(length(miss.sims$marks)>0){
  repeat{
    miss.runs = miss.sims
    est_on <<-est_on
    time_on <<- time_on #LK
    age_on <<- age_on #LK
    mod<- foreach(d=1:length(miss.runs$marks),.errorhandling='pass',.packages=c('snowfall','parallel','tcltk2','RMark')
    ) %dopar% {
      run.sim(miss.runs,d) 
    }
      
    #check for missing runs - checks if pdf file for sim plots is written
    done.sims = paste0(basename(dirname(dirname(list.files(dir.run, recursive=TRUE, full.names=TRUE,pattern="\\.pdf$")))),'/',basename((dirname(list.files(dir.run, recursive=TRUE, full.names=TRUE,pattern="\\.pdf$")))))
    miss.sims = run.list[which(!paste0(run.list$fold)%in%done.sims),]
    print(paste(length(miss.sims$marks),"simulations did not write to file"))
    
    #continue or stop?
    if(length(miss.sims$marks)==0){
      break
    }
  }
}
close(pbar)


##################################################
#Estimation 
##################################################
est_on <<-est_on
time_on <<- time_on 
age_on <<- age_on 
if(est_on == 1){
mod <- foreach(d=1:length(run.list$marks),.errorhandling='pass',.packages=c('snowfall','parallel','tcltk2','RMark')
) %dopar% {
  run_estimation(run.list,d)
}



#check for missing runs - checks if .vcv file for est is written
done.ests = paste0(basename(dirname(dirname(dirname(list.files(dir.run, recursive=TRUE, full.names=TRUE,pattern="\\.vcv$"))))),'/',basename(dirname(dirname(list.files(dir.run, recursive=TRUE, full.names=TRUE,pattern="\\.vcv$")))))
miss.ests = run.list[which(!paste0(run.list$fold)%in%done.ests),]
print(paste(length(miss.ests$marks),"estimations did not write to file"))

#redo missing simulation runs, repeat until no more
if(length(miss.ests$marks)>0){
  repeat{
    miss.estruns = miss.ests
    est_on <<-est_on
    time_on <<- time_on 
    age_on <<- age_on 
    mod<- foreach(d=1:length(miss.estruns$marks),.errorhandling='pass',.packages=c('snowfall','parallel','tcltk2','RMark')
    ) %dopar% {
      run_estimation(miss.estruns,d)
    }
    
    #check for missing runs - checks if .vcv file for est is written
    done.ests = paste0(basename(dirname(dirname(dirname(list.files(dir.run, recursive=TRUE, full.names=TRUE,pattern="\\.vcv$"))))),'/',basename(dirname(dirname(list.files(dir.run, recursive=TRUE, full.names=TRUE,pattern="\\.vcv$")))))
    miss.ests = run.list[which(!paste0(run.list$fold)%in%done.ests),]
    print(paste(length(miss.ests$marks),"estimations did not write to file"))
    
    #continue or stop?
    if(length(miss.ests$marks)==0){
      break
    }
  }
}
}
stopCluster(cl);  closeAllConnections();  output3<-gc()

}
if (run_simulation == 0){
  print('Simulation Off<>< <>< <>< <>< ')
}
if(est_on == 0){
  print('Estimation Off<>< <>< <>< <>< ')
}


