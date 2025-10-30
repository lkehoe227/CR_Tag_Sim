#Data and BAM processing script

#from BAM-------------------------------------------------------------------------------------------
#bam = dget("data/s73.rdat") # can extract directly from bam 
#or read in the supplied file and not run this script
load("data/sim_data.RData")

if(exists("bam")){
#ages
ages    = 1:maxage
nages   = length(ages)

#numbers at age
Na.bam  = colMeans(bam$N.age[67:71,])  #numbers at age, start of year, avg last five years, 
#need monthly age (use spline from mid point of year)
Na.spline = smooth.spline((ages+0.5)*12,(Na.bam/12))
Na.calc = predict(Na.spline,1:(nages*12))

#size at age
La    = bam$a.series$length /10 
La.m  = (bam$parms$Linf*(1-exp(-bam$parms$K*(ages.m-bam$parms$t0))))/10  
Wa.m  = 1000*bam$parms$wgt.a*(La.m*10)^bam$parms$wgt.b  #in grams   


#mortality
Ma.bam    = bam$a.series$M
Fa.bam    = colMeans(bam$F.age[66:70,])  #total F-at-age, average over last five years
Za.bam    = Ma.bam+Fa.bam
dmort.bam = bam$parms$D.mort.GR3*.9

#natural mortality - calculated
Ma.m1 = (3.69*Wa.m^-.305)/12
Mtarg = mean(Ma.bam)/12
Ma.m  = Mtarg*(Ma.m1/mean(Ma.m1))

#selex at age
selex.bam = bam$sel.age$sel.v.wgted.tot #combined fleets
smspl = smooth.spline((ages-0.5)*12,selex.bam)
selex.m = predict(smspl,1:(nages*12))$y 
selex.m = selex.m/max(selex.m)


#recruitment
Ro        = median(bam$t.series$recruits)     #median recruits for projections=439823 
Rsigma    = bam$parms$R.sigma.par   #standard deviation of recruitment residuals in log space, set to 0 for constant, fixed recruitment
Rbiascor  = bam$parms$Mean.biascorr  

#monthly recruitment proportion
spawn = dnorm(1:12,mean=5,sd=1.5)
spawn = spawn/sum(spawn)

#Initialize-----------------------------------------------------------------------------------------
#F-at-age
F30 = 0.21  
Fcur = mean(tail(bam$t.series$F.full,6),na.rm=T)  #F full from BAM
Fcur.m = Fcur/12
Fa.m = Fcur.m*selex.m
#new formulation numbers approach (average 2017-2019 live releases and total catch (2017-2019))
tot_releases = (mean(bam$t.series$total.D.knum[68:70])/0.234)*1000
tot_harvest = (mean(bam$t.series$total.L.knum[68:70]))*1000
Ucur = (tot_releases+tot_harvest)/sum(colMeans(bam$N.age[68:70,])) #0.9283913
#est exploitation approach - 
#release_e = (mean(bam$t.series$E.D.num[68:70])*1000)/0.234
#harv_e = mean(bam$t.series$E.L.num[68:70])*1000
#Ucur = (release_e+harv_e)/ 1000 #0.9295519

#Z-at-age 
Za.m = Fa.m+Ma.m

#survivorship
surv0 = survF = rep(NA,nages*12)
surv0[1] = survF[1] = 1
for(a in 2:length(surv0)){
  surv0[a] = surv0[a-1] * exp(-Ma.m[a-1])
  survF[a] = survF[a-1] * exp(-Za.m[a-1])
}
  
#monthly fishing effort and cpue inputs-------------------------------------------------------------
if(eff_data == 'mrip'){
eff           = read.csv('data/mrip effort and cpue.csv')
dangtrips     = eff$RS.dangtrips.pvt #angler trips that caught or targeted red snapper NEFL
dtrips        = eff$RS.dtrips.pvt #trips that caught or targeted red snapper NEFL
cpue.dangtrip = eff$med.cpue.dangtrips.pvt #catch per angler trip NE FL
cpue.dtrip    = eff$med.cpue.dtrip.pvt #catch per trip NE FL
dprop         = eff$dprop.pvt #proportion of catch that is discarded, NEFL (includes out of seasons proportions because of avg years when harvest was outside july)
}
if(eff_data == 'srfs'){ #default
eff           = read.csv('data/srfs_effort_dat.csv')
dtrips        = eff$dtrips #trips that caught or targeted red snapper NEFL
cpue.dtrip    = eff$cpue #catch per trip  FL
dprop = c(rep(1,6),eff$dprop[7], rep(1,5)) #makes only july having harvest

}

#Effort and model dimensions 
Etype = 'dtrips'  #can only use dangtrips for mrip

#time dynamic forcing vectors-----------------------------------------------------------------------
#spawning seasonality
spawn.m = rep(spawn,nyrs)

#set effort and cpue type
if(Etype=='dtrips'){effort = dtrips; cpue = cpue.dtrip}

#scale effort so when multiplied by cpue you get catch from BAM
bam.totC = bam$L.age.pred.knum+bam$D.age.pred.knum/dmort.bam
Cbam = (sum(colMeans(bam.totC[66:70,]))*1000)
Chat = sum(effort*cpue)
q = Cbam/Chat
effort = ceiling(effort*q)

#repeat effort and cpue vectors for each year
effort.m = rep(effort,nyrs)
cpue.m = rep(cpue,nyrs)
dprop.m = rep(dprop,nyrs)

#catchability, fishing mortality, exploitation rate
qF.base = Fcur/sum(effort)
qU.base = Ucur/sum(effort)

qF.m = qF.base*(cpue/mean(cpue))
qF.m = qF.m*Fcur/(sum(qF.m*effort))
qU.m = qU.base*(cpue/mean(cpue))
qU.m = qU.m*Ucur/(sum(qU.m*effort))
qF.m = rep(qF.m,nyrs)
qU.m = rep(qU.m,nyrs)
Ft = effort.m*qF.m
Ut = effort.m*qU.m

}else{
  print('Populaiton Intialized')
}
