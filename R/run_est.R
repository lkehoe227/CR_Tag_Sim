###################################################################################################################################
#### Run Estimation ####
###################################################################################################################################


run_estimation = function(run.list,d){
  setwd(dir.base)
  if(est_on == 1){

    setwd(dir.base)
Tdat2 = read.csv(paste0(run.list$full.dir[d],'/barker_dat.csv'))


col.idx = which(colnames(Tdat2)=="L1")

Tdat2$ch = apply(Tdat2[,col.idx:ncol(Tdat2)],1,function(x)paste(x,collapse=""))  #create capture history string



#-----------------------------------
#annual estimation 
#-----------------------------------
if (timestep == 'annual'){
Tdat2$age.marked.intyr = as.factor(Tdat2$age.marked.intyr)
#make Barker data file-----------------------------------------------------------------------------
barker.dat = process.data(data=Tdat2[,c(1:c(col.idx-1),ncol(Tdat2))], model='Barker',groups=c("age.marked.intyr"),
                          age.var=1, initial.ages=as.numeric(levels(Tdat2$age.marked.intyr))) 
barker.des1 = make.design.data(data=barker.dat)
barker.des1 = add.design.data(data=barker.dat, barker.des1, parameter="S",type="age",bins=c(0:6,25),name='agebin')#bins=c(1:6,25)
barker.des1 = add.design.data(data=barker.dat, barker.des1, parameter="r",type="age",bins=c(0:6,25),name='agebin')
barker.des1 = add.design.data(data=barker.dat, barker.des1, parameter="R",type="age",bins=c(0:6,25),name='agebin')
barker.des1 = add.design.data(data=barker.dat, barker.des1, parameter="p",type="age",bins=c(0:6,25),name='agebin') 
#age.bins are [0,1] (1,2]  (2,3]  (3,4]  (4,5]  (5,6] (6,25] 
names(barker.des1)
barker.des1$S$fix=NULL
barker.des1$p$fix=NULL
barker.des1$r$fix=NULL
barker.des1$R$fix=NULL
barker.des1$RPrime$fix=NULL
barker.des1$F$fix=NULL
barker.des1$FPrime$fix=NULL

#model 1: S.,p.-------------------------------------------------------------------------------------
dir.est = paste0(run.list$full.dir[d],"/Barker1_Sdotpdot")


if(!dir.exists(dir.est)) dir.create(dir.est);
unlink(dir.est);setwd(dir.est)

barker.mod1 = make.mark.model(data=barker.dat, ddl=barker.des1,title=basename(dir.run),
                              parameters=list(S=list(formula=~agebin+time), 
                                              p=list(formula=~time), 
                                              r=list(formula=~agebin+time), 
                                              R=list(formula=~1),
                                              RPrime=list(formula=~1,fixed=0),
                                              F=list(formula=~1,fixed=1),
                                              FPrime=list(formula=~1,fixed=0))) 
barker.results1 = run.mark.model(model=barker.mod1) 
est_out = barker.results1$results$real

write.csv(est_out,paste0(run.list$full.dir[d],"/est_results.csv"))

}
  
  
#-----------------------------------
#monthly estimation 
#-----------------------------
  if (timestep == 'monthly'){
    if(age_on == 1&time_on == 1&nyrs==1){
#make Barker data file-----------------------------------------------------------------------------
barker.dat = process.data(data=Tdat2[,c(1:c(col.idx-1),ncol(Tdat2))], model='Barker',groups=c("age.marked.mo"), #can add an age marked covariate 
                          age.var=1, initial.ages=as.numeric(unique(Tdat2$age.marked.mo)))# time.intervals = rep(1/12,12) , begin.time = 4
barker.des1 = make.design.data(data=barker.dat) #makes all possible combinations of covariates (age, year, cohort, etc)

barker.des1 = add.design.data(data=barker.dat, barker.des1, parameter="S",type="age",bins=c(12,24,72,25*12),name='agebin') 
barker.des1 = add.design.data(data=barker.dat, barker.des1, parameter="r",type="age",bins=c(12,24,72,25*12),name='agebin')
barker.des1 = add.design.data(data=barker.dat, barker.des1, parameter="R",type="age",bins=c(12,24,72,25*12),name='agebin')
barker.des1 = add.design.data(data=barker.dat, barker.des1, parameter="p",type="age",bins=c(12,24,72,25*12),name='agebin')
#age bins in months   [12,24]  (24,36]  (36,48]  (48,60]  (60,72] (72,300]

if(est_type == "pre-post"){
barker.des1$S$season=ifelse(barker.des1$S$time%in%c(1:4),0,ifelse(barker.des1$S$time%in%c(5:8),1,2)) #0 = pre, 1=during, 2=after fishing
barker.des1$r$season=ifelse(barker.des1$r$time%in%c(1:4),0,ifelse(barker.des1$r$time%in%c(5:8),1,2)) #0 = pre, 1=during, 2=after fishing
barker.des1$R$season=ifelse(barker.des1$R$time%in%c(1:4),0,ifelse(barker.des1$R$time%in%c(5:8),1,2)) #0 = pre, 1=during, 2=after fishing

barker.des1$S$season = as.factor(barker.des1$S$season)
barker.des1$r$season = as.factor(barker.des1$r$season)
barker.des1$R$season = as.factor(barker.des1$R$season)

barker.des1$r$fix[barker.des1$r$season=='0']=0.0001 
barker.des1$r$fix[barker.des1$r$season=='2']=0.0001 
}
if(est_type == "july"){
#to make july season
barker.des1$S$season=ifelse(barker.des1$S$time%in%c(1:4,10:12),0,ifelse(barker.des1$S$time%in%c(7),2,1)) #0 = low pressure, 1 = high pressure, 2= july og = (1:4,9:12)
barker.des1$r$season=ifelse(barker.des1$r$time%in%c(1:4,10:12),0,ifelse(barker.des1$r$time%in%c(7),2,1)) #0 = low pressure, 1 = high pressure, 2= july
barker.des1$R$season=ifelse(barker.des1$R$time%in%c(1:4,10:12),0,ifelse(barker.des1$R$time%in%c(7),2,1)) #0 = low pressure, 1 = high pressure, 2= july

barker.des1$S$season = as.factor(barker.des1$S$season)
barker.des1$r$season = as.factor(barker.des1$r$season)
barker.des1$R$season = as.factor(barker.des1$R$season)

barker.des1$r$fix[barker.des1$r$season=='0']=0.0001 

}
if(est_type == "high-low"){
  barker.des1$S$season=ifelse(barker.des1$S$time%in%c(1:4,10:12),0,1) #0 = low pressure, 1 = high pressure
  barker.des1$r$season=ifelse(barker.des1$r$time%in%c(1:4,10:12),0,1) #0 = low pressure, 1 = high pressure
  barker.des1$R$season=ifelse(barker.des1$R$time%in%c(1:4,10:12),0,1) #0 = low pressure, 1 = high pressure
  
  barker.des1$S$season = as.factor(barker.des1$S$season)
  barker.des1$r$season = as.factor(barker.des1$r$season)
  barker.des1$R$season = as.factor(barker.des1$R$season) 
  barker.des1$r$fix[barker.des1$r$season=='0']=0.0001 
  
}
if(est_type == "may-dec"){
  #to make may 1st tagging start season
  barker.des1$S$season=ifelse(barker.des1$S$time%in%c(1:5),0,ifelse(barker.des1$S$time%in%c(6:9),1,2)) #0 = pre, 1=during, 2=after fishing
  barker.des1$r$season=ifelse(barker.des1$r$time%in%c(1:5),0,ifelse(barker.des1$r$time%in%c(6:9),1,2)) #0 = pre, 1=during, 2=after fishing
  barker.des1$R$season=ifelse(barker.des1$R$time%in%c(1:5),0,ifelse(barker.des1$R$time%in%c(6:9),1,2)) #0 = pre, 1=during, 2=after fishing
  barker.des1$S$season = as.factor(barker.des1$S$season)
  barker.des1$r$season = as.factor(barker.des1$r$season)
  barker.des1$R$season = as.factor(barker.des1$R$season)
  
  barker.des1$r$fix[barker.des1$r$season=='0']=0.0001 
  barker.des1$r$fix[barker.des1$r$season=='2']=0.0001 
  
}
if(est_type == "june-dec"){
  #to make may 1st tagging start season
  barker.des1$S$season=ifelse(barker.des1$S$time%in%c(1:6),0,ifelse(barker.des1$S$time%in%c(7:9),1,2)) #0 = pre, 1=during, 2=after fishing
  barker.des1$r$season=ifelse(barker.des1$r$time%in%c(1:6),0,ifelse(barker.des1$r$time%in%c(7:9),1,2)) #0 = pre, 1=during, 2=after fishing
  barker.des1$R$season=ifelse(barker.des1$R$time%in%c(1:6),0,ifelse(barker.des1$R$time%in%c(7:9),1,2)) #0 = pre, 1=during, 2=after fishing
  barker.des1$S$season = as.factor(barker.des1$S$season)
  barker.des1$r$season = as.factor(barker.des1$r$season)
  barker.des1$R$season = as.factor(barker.des1$R$season)
  
  barker.des1$r$fix[barker.des1$r$season=='0']=0.0001 
  barker.des1$r$fix[barker.des1$r$season=='2']=0.0001 
  
}

barker.des1$S$fix=NULL
barker.des1$p$fix=NULL
barker.des1$r$fix[barker.des1$r$agebin=='[12,24]']=0.0001  #fixes age 0-1 at a very small number 

barker.des1$R$fix=NULL
barker.des1$RPrime$fix=NULL
barker.des1$F$fix=NULL
barker.des1$FPrime$fix=NULL

#model 1: S.,p.-------------------------------------------------------------------------------------
dir.est = paste0(run.list$full.dir[d],"/Barker1_Sdotpdot")

if(!dir.exists(dir.est)) dir.create(dir.est);
unlink(dir.est);setwd(dir.est)

barker.mod1 = make.mark.model(data=barker.dat, ddl=barker.des1,title=basename(dir.run),
                              parameters=list(S=list(formula=~agebin+season),
                                              p=list(formula=~1,fixed=0),
                                              r=list(formula=~agebin+season),
                                              R=list(formula=~agebin+season),
                                              RPrime=list(formula=~1,fixed=0),
                                              F=list(formula=~1,fixed=1),
                                              FPrime=list(formula=~1,fixed=1))) 
barker.results1 = run.mark.model(model=barker.mod1) 
est_out = barker.results1$results$real

#two year
}
    if(age_on == 1&time_on == 1&nyrs==2){
      #make Barker data file-----------------------------------------------------------------------------
      barker.dat = process.data(data=Tdat2[,c(1:c(col.idx-1),ncol(Tdat2))], model='Barker',groups=c("age.marked.mo"),
                                age.var=1, initial.ages=as.numeric(unique(Tdat2$age.marked.mo)))# time.intervals = rep(1/12,12) , begin.time = 4
      barker.des1 = make.design.data(data=barker.dat) #makes all possible combinations of covariates (age, year, cohort, etc)
      barker.des1 = add.design.data(data=barker.dat, barker.des1, parameter="S",type="age",bins=c(12,24,72,25*12),name='agebin') 
      barker.des1 = add.design.data(data=barker.dat, barker.des1, parameter="r",type="age",bins=c(12,24,72,25*12),name='agebin')
      barker.des1 = add.design.data(data=barker.dat, barker.des1, parameter="R",type="age",bins=c(12,24,72,25*12),name='agebin')
      #age bins in months   [12,24]  (24,72] (72,300]
      
      if(est_type == "july"){ 
        barker.des1$S$season=ifelse(barker.des1$S$time%in%c(1:4,9:12, 13:16,21:24),0,ifelse(barker.des1$S$time%in%c(7,19),2,1)) #0 = low pressure, 1 = high pressure, 2= july
        barker.des1$r$season=ifelse(barker.des1$r$time%in%c(1:4,9:12, 13:16,21:24),0,ifelse(barker.des1$r$time%in%c(7,19),2,1)) #0 = low pressure, 1 = high pressure, 2= july
        barker.des1$R$season=ifelse(barker.des1$R$time%in%c(1:4,9:12, 13:16,21:24),0,ifelse(barker.des1$R$time%in%c(7,19),2,1)) #0 = low pressure, 1 = high pressure, 2= july
        
        barker.des1$S$season = as.factor(barker.des1$S$season)
        barker.des1$r$season = as.factor(barker.des1$r$season)
        barker.des1$R$season = as.factor(barker.des1$R$season)
        names(barker.des1)
        barker.des1$S$fix=NULL
        barker.des1$p$fix=NULL
        barker.des1$r$fix[barker.des1$r$agebin=='[12,24]']=0.0001  #fixes age 0-1 at a very small number 
        barker.des1$r$fix[barker.des1$r$season=='0']=0.0001
        barker.des1$r$fix[barker.des1$r$season=='1']=0.0001
      }
      if(est_type == "high-low_yr2"){ 
        barker.des1$S$season=ifelse(barker.des1$S$time%in%c(1:16),0,ifelse(barker.des1$S$time%in%c(17:20),2,1)) #0 = yr1-2 tagging season , 1 = yr2 fishing season, 2= yr2 post season
        barker.des1$r$season=ifelse(barker.des1$r$time%in%c(1:16),0,ifelse(barker.des1$r$time%in%c(17:20),2,1)) #0 = low pressure, 1 = high pressure, 2= july
        barker.des1$R$season=ifelse(barker.des1$R$time%in%c(1:16),0,ifelse(barker.des1$R$time%in%c(17:20),2,1)) #0 = low pressure, 1 = high pressure, 2= july
        
        barker.des1$S$season = as.factor(barker.des1$S$season)
        barker.des1$r$season = as.factor(barker.des1$r$season)
        barker.des1$R$season = as.factor(barker.des1$R$season)
        names(barker.des1)
        barker.des1$S$fix=NULL
        barker.des1$p$fix=NULL
        barker.des1$r$fix[barker.des1$r$agebin=='[12,24]']=0.0001  #fixes age 0-1 at a very small number (0 didn't seem to work right)
        barker.des1$r$fix[barker.des1$r$season=='2']=0.0001
      }
      barker.des1$R$fix=NULL
      barker.des1$RPrime$fix=NULL
      barker.des1$F$fix=NULL
      barker.des1$FPrime$fix=NULL
      
      #model 1: S.,p.-------------------------------------------------------------------------------------
      dir.est = paste0(run.list$full.dir[d],"/Barker1_Sdotpdot")
      
      if(!dir.exists(dir.est)) dir.create(dir.est);
      unlink(dir.est);setwd(dir.est)
      
      barker.mod1 = make.mark.model(data=barker.dat, ddl=barker.des1,title=basename(dir.run),
                                    parameters=list(S=list(formula=~agebin+season),
                                                    p=list(formula=~1,fixed=0),
                                                    r=list(formula=~agebin+season),
                                                    R=list(formula=~agebin+season),
                                                    RPrime=list(formula=~1,fixed=0),
                                                    F=list(formula=~1,fixed=1),
                                                    FPrime=list(formula=~1,fixed=1))) 
      barker.results1 = run.mark.model(model=barker.mod1) 
      est_out = barker.results1$results$real
      
      
    }
    #three year 
if(age_on == 1&time_on == 1&nyrs==3){
  #make Barker data file-----------------------------------------------------------------------------
  barker.dat = process.data(data=Tdat2[,c(1:c(col.idx-1),ncol(Tdat2))], model='Barker',groups=c("age.marked.mo"), 
                            age.var=1, initial.ages=as.numeric(unique(Tdat2$age.marked.mo)))
  barker.des1 = make.design.data(data=barker.dat) #makes all possible combinations of covariates (age, year, cohort, etc)

  barker.des1 = add.design.data(data=barker.dat, barker.des1, parameter="S",type="age",bins=c(12,24,72,25*12),name='agebin') 
  barker.des1 = add.design.data(data=barker.dat, barker.des1, parameter="r",type="age",bins=c(12,24,72,25*12),name='agebin')
  barker.des1 = add.design.data(data=barker.dat, barker.des1, parameter="R",type="age",bins=c(12,24,72,25*12),name='agebin')

  if(est_type == "july"){ 
  barker.des1$S$season=ifelse(barker.des1$S$time%in%c(1:4,9:12, 13:16,21:24, 25:28,33:36),0,ifelse(barker.des1$S$time%in%c(7,19,31),2,1)) #0 = low pressure, 1 = high pressure, 2= july
  barker.des1$r$season=ifelse(barker.des1$r$time%in%c(1:4,9:12, 13:16,21:24, 25:28,33:36),0,ifelse(barker.des1$r$time%in%c(7,19,31),2,1)) #0 = low pressure, 1 = high pressure, 2= july
  barker.des1$R$season=ifelse(barker.des1$R$time%in%c(1:4,9:12, 13:16,21:24, 25:28,33:36),0,ifelse(barker.des1$R$time%in%c(7,19,31),2,1)) #0 = low pressure, 1 = high pressure, 2= july
  
  barker.des1$S$season = as.factor(barker.des1$S$season)
  barker.des1$r$season = as.factor(barker.des1$r$season)
  barker.des1$R$season = as.factor(barker.des1$R$season)
  names(barker.des1)
  barker.des1$S$fix=NULL
  barker.des1$p$fix=NULL
  barker.des1$r$fix[barker.des1$r$agebin=='[12,24]']=0.0001  #fixes age 0-1 at a very small number (0 didn't seem to work right)
  barker.des1$r$fix[barker.des1$r$season=='0']=0.0001
  barker.des1$r$fix[barker.des1$r$season=='2']=0.0001
  }
  
  barker.des1$R$fix=NULL
  barker.des1$RPrime$fix=NULL
  barker.des1$F$fix=NULL
  barker.des1$FPrime$fix=NULL
  
  #model 1: S.,p.-------------------------------------------------------------------------------------
  dir.est = paste0(run.list$full.dir[d],"/Barker1_Sdotpdot")
  
  if(!dir.exists(dir.est)) dir.create(dir.est);
  unlink(dir.est);setwd(dir.est)
  
  barker.mod1 = make.mark.model(data=barker.dat, ddl=barker.des1,title=basename(dir.run),
                                parameters=list(S=list(formula=~agebin+season),
                                                p=list(formula=~1,fixed=0),
                                                r=list(formula=~agebin+season),
                                                R=list(formula=~agebin+season),
                                                RPrime=list(formula=~1,fixed=0),
                                                F=list(formula=~1,fixed=1),
                                                FPrime=list(formula=~1,fixed=1))) 
  barker.results1 = run.mark.model(model=barker.mod1) 
  est_out = barker.results1$results$real
  
  
    }
    if(age_on == 1&time_on==0){
      #make Barker data file-----------------------------------------------------------------------------
      barker.dat = process.data(data=Tdat2[,c(1:c(col.idx-1),ncol(Tdat2))], model='Barker',groups=c("age.marked.mo"), 
                                age.var=1, initial.ages=as.numeric(unique(Tdat2$age.marked.mo)))
      barker.des1 = make.design.data(data=barker.dat) #makes all possible combinations of covariates (age, year, cohort, etc)
      barker.des1 = add.design.data(data=barker.dat, barker.des1, parameter="S",type="age",bins=c(1:6*12,25*12),name='agebin') 
      barker.des1 = add.design.data(data=barker.dat, barker.des1, parameter="r",type="age",bins=c(1:6*12,25*12),name='agebin')
      barker.des1 = add.design.data(data=barker.dat, barker.des1, parameter="R",type="age",bins=c(1:6*12,25*12),name='agebin')
      barker.des1 = add.design.data(data=barker.dat, barker.des1, parameter="p",type="age",bins=c(1:6*12,25*12),name='agebin')
      #age bins in months  [0,12] (12,24]  (24,36]  (36,48]  (48,60]  (60,72] (72,300]
      names(barker.des1)
      barker.des1$S$fix=NULL
      barker.des1$p$fix=NULL
      barker.des1$r$fix[barker.des1$r$agebin=='[0,12]']=0  #fixes age 0-1 at 0
      barker.des1$R$fix=NULL
      barker.des1$RPrime$fix=NULL
      barker.des1$F$fix=NULL
      barker.des1$FPrime$fix=NULL
      
      #model 1: S.,p.-------------------------------------------------------------------------------------
      dir.est = paste0(run.list$full.dir[d],"/Barker1_Sdotpdot")
      
      if(!dir.exists(dir.est)) dir.create(dir.est);
      unlink(dir.est);setwd(dir.est)
      
      barker.mod1 = make.mark.model(data=barker.dat, ddl=barker.des1,title=basename(dir.run),
                                    parameters=list(S=list(formula=~agebin),
                                                    p=list(formula=~1,fixed=0),
                                                    r=list(formula=~agebin),
                                                    R=list(formula=~agebin),
                                                    RPrime=list(formula=~1,fixed=0),
                                                    F=list(formula=~1,fixed=1),
                                                    FPrime=list(formula=~1,fixed=1))) 
      barker.results1 = run.mark.model(model=barker.mod1) 
      est_out = barker.results1$results$real
    }
    if(age_on == 0){
      #make Barker data file-----------------------------------------------------------------------------
      barker.dat = process.data(data=Tdat2[,c(1:c(col.idx-1),ncol(Tdat2))], model='Barker')
      barker.des1 = make.design.data(data=barker.dat) #makes all possible combinations of covariates (age, year, cohort, etc)
      names(barker.des1)
      barker.des1$S$fix=NULL
      barker.des1$p$fix=NULL
      barker.des1$r$fix=NULL
      barker.des1$R$fix=NULL
      barker.des1$RPrime$fix=NULL
      barker.des1$F$fix=NULL
      barker.des1$FPrime$fix=NULL
      
      #model 1: S.,p.-------------------------------------------------------------------------------------
      dir.est = paste0(run.list$full.dir[d],"/Barker1_Sdotpdot")
      
      if(!dir.exists(dir.est)) dir.create(dir.est);
      unlink(dir.est);setwd(dir.est)
      
      barker.mod1 = make.mark.model(data=barker.dat, ddl=barker.des1,title=basename(dir.run),
                                    parameters=list(S=list(formula=~1),
                                                    p=list(formula=~1,fixed=0),
                                                    r=list(formula=~1),
                                                    R=list(formula=~1),
                                                    RPrime=list(formula=~1,fixed=0),
                                                    F=list(formula=~1,fixed=1),
                                                    FPrime=list(formula=~1,fixed=1))) 
      barker.results1 = run.mark.model(model=barker.mod1) 
      est_out = barker.results1$results$real
    }
    
write.csv(est_out,paste0(run.list$full.dir[d],"/est_results.csv"))
  }
  print("Done Estimation")
  }
  if (est_on ==0)
    print('Estimation Off <>< <>< <>< <>< <>< <><')
}

