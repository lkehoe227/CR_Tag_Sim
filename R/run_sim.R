#### Run Simulation 

####################################################################################################
# Simulation function 
####################################################################################################
  run.sim <- function(run.list,d){
    
      
  setwd(dir.base) 
    
    #tagging parameters---------------------------------------------------------------------------------
    if(tag_type == "gene"){
      #number of effort units reporting each year
      Edot = p.trips[d]*sum(effort)
      Edot.m = round(Edot*(effort/sum(effort)))
      Edot.m = rep(Edot.m,nyrs)
      Edotprop.m = Edot.m/effort.m
    }
    if(tag_type == "conv"){
     ntags = run.list$marks[d]
    
      #number of effort units reporting each year ~~~~~~~~~~~~~~~~~~~~~~REPORTING RATE~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~@@@@@@@@@@@@@@
      reporting_rate <- rr 
      rep_yr <-rep(reporting_rate,nyrs) #annual reporting rate 
      rep_mo <-rep(reporting_rate,length(effort)*nyrs) #monthly reporting rate 
    
      #tagging distribution 
      tag_seas <- season #indexing what months tags go out (4 = April, 1:12 = all months)
      if(tag_release == 'monthly'){
      tag_dist <- rep(0,length(effort)) 
      tag_dist[tag_seas] <- ntags/length(season) 
      tag_dist <- rep(tag_dist,nyrs)
      }
      if(tag_release == 'fixed'){
        tag_dist <- rep(0,ntsteps)
        tag_dist[tag_seas] <- ntags/length(season)
        tag_dist[prev_tag_seas] <- prev_tag/length(prev_tag_seas)
      }
      Edot.m = rep_mo #this is just a place holder for conv
    }

    ####################################################################################################
    #  RUN TAGGING SIMULATION
    ####################################################################################################

    for(t in 1:ntsteps){
      #initialize (t==1)==================================================================================
      if(t==1){
        gc();
        time.start <-Sys.time()
        #setup storage objects------------------------------------------------------------------------------
        #fishing mortality, natural mortality, numbers N, Prop N tagged, exploitation, total mortality, catch, discards, retained yield,tagged N
        Fta <- Mta <- Nta <- Pta <- Uta <- Zta <- Cta <- Dta <- Tta <- Rta <- Lta <- Kta <- matrix(0,nrow = ntsteps,ncol = length(ages.m))
        Hit = Sit = Lit = Kit = Wit = Tit = Dit = Vit = Pit= Rit = Xit = (matrix(nrow=0,ncol=ntsteps))  
        colnames(Hit) = colnames(Sit) =colnames(Lit) =colnames(Kit) =colnames(Wit) = colnames(Tit) =c(paste('t',formatC(1:ntsteps,digits=2,flag='0')))
        Rt <-rep(0,ntsteps)
        Tdat <-(matrix(0,nrow=0,ncol=3)) ; colnames(Tdat) = c('tstep.marked', 'age.marked','age.t') 
        
        #initialize population------------------------------------------------------------------------------
        Nta[t,] <- Na.calc$y #population numbers (based on BAM)
        Fta[t,] <- qF.m[t]*effort.m[t]*selex.m #fishing mortality 
        Uta[t,] <- qU.m[t]*effort.m[t]*selex.m #exploitation rate
        Zta[t,] <- Fta[t,] + Ma.m #total mortality
        Cta[t,] <- ceiling(Nta[t,]*Uta[t,]) #catch
        
        #first marking event--------------------------------------------------------------------------------
        if(tag_type == "gene"){
          #participant caps at age
          caps = ceiling(qU.m[t]*Edot.m[t]*sum(selex.m*Nta[t,])) #total captures for timestep by participants, given U, effort, selex, and Na
          recaps = sum(Pit[,t],na.rm=T)  #number of recaps in timestep
          newcaps = max(0,caps-recaps) #number of new marked fish released by participants (total caps minus recaps)
          Rt[t] <- newcaps  #total number of new marks
        }
        if(tag_type == "conv"){
          Rt[t] <-tag_dist[t]  #total number of new marks
        }
        print('--determining fates of new marked fish');flush.console()
        if(Rt[t]>0){
          Rta[t,] <-rmultinom(n=1, size=Rt[t], Cta[t,]/sum(Cta[t,])) #new caps at age
        } else{
          Rta[t,] <-0
        }
        
        #newly marked individuals
        T1 <- data.frame(tstep = rep(t,sum(Rta[t,])),age = rep(ages.m,times=Rta[t,]))
        T1$a.t <- ages.m[match(T1$age,ages.m)]
        Tdat <- rbind(Tdat,as.matrix(T1)) 
        
        #expand matrices to store fates
        Sit <- as.matrix(rbind(Sit,(matrix(nrow=sum(Rta[t,]),ncol=ntsteps))))
        Tit <- as.matrix(rbind(Tit,(matrix(nrow=sum(Rta[t,]),ncol=ntsteps))))
        Lit <- as.matrix(rbind(Lit,(matrix(nrow=sum(Rta[t,]),ncol=ntsteps))))
        Kit <- as.matrix(rbind(Kit,(matrix(nrow=sum(Rta[t,]),ncol=ntsteps))))
        Wit <- as.matrix(rbind(Wit,(matrix(nrow=sum(Rta[t,]),ncol=ntsteps))))
        Dit <- as.matrix(rbind(Dit,(matrix(nrow=sum(Rta[t,]),ncol=ntsteps))))
        Vit <- as.matrix(rbind(Vit,(matrix(nrow=sum(Rta[t,]),ncol=ntsteps))))
        Pit <- as.matrix(rbind(Pit,(matrix(nrow=sum(Rta[t,]),ncol=ntsteps))))
        Rit <- as.matrix(rbind(Rit,(matrix(nrow=sum(Rta[t,]),ncol=ntsteps))))
        Xit <- as.matrix(rbind(Xit,(matrix(nrow=sum(Rta[t,]),ncol=ntsteps))))
        
        #fates of new releases: Survival to next timestep is based on release mortality.
        i.idx <- numeric(0)
        if(Rt[t]>0) i.idx <-(nrow(Tdat)-nrow(T1)+1):nrow(Tdat)
        Tit[i.idx,t] <-1  #set initial tag state to 1

        
        if(length(i.idx)>0){
          #run fates function for new marks
          if (tag_type == 'conv'){
            sapply(i.idx, function(i){
            a.t <<- Tdat[i,'age.t']
            a.idx <<- match(a.t,ages.m)
            
            #first timestep 
            if(t==1){
              Pit[i,t] <<- 1  #capture probability = 1 as this is a new mark   
              Xit[i,t] <<- rbinom(1,1,handle.mort*Pit[i,t]) #handling mortality, 1=yes
              Sit[i,t] <<- (1-Xit[i,t])*rbinom(1,1,exp(-(Zta[t,a.idx])))    #survive timestep, 1=yes
              Tit[i,t] <<- Sit[i,t]  #update tag state
            } 
            if(t>1){
              #new marks in subsequent timesteps
              if(is.na(Tit[i,t-1])){
                Pit[i,t] <<- 1 #capture probability = 1 as this is a new mark   
                Xit[i,t] <<- rbinom(1,1,handle.mort*Pit[i,t]) #handling mortality, 1=yes
                Sit[i,t] <<- (1-Xit[i,t])*rbinom(1,1,exp(-(Zta[t,a.idx])))       #survive current timestep, 1=yes
                Tit[i,t] <<- Sit[i,t]  #update tag state
                #existing marks
              } else {
                Pit[i,t] <<- Tit[i,t-1]*rbinom(1,1,Pta[t,a.idx])*rbinom(1,1,Rta[t,a.idx]/Nta[t,a.idx]) #recapture probability, 1=yes
                Xit[i,t] <<- rbinom(1,1,handle.mort*Pit[i,t]) #handling mortality, 1=yes
                Sit[i,t] <<- (1-Xit[i,t])*Tit[i,t-1]*rbinom(1,1,exp(-(Zta[t,a.idx])))      #survive current timestep, 1=yes
                Tit[i,t] <<- Tit[i,t-1]*Sit[i,t]  #update tag state
              }
            }
            #if survived
            Kit[i,t] <<- Sit[i,t]*rbinom(1,1,1-exp(-Uta[t,a.idx]))*rbinom(1,1,dprop.m[t]) #*rbinom(1,1,Pta[t,a.idx])#if survived, was it captured and released alive, 1=yes #added *rbinom(1,1,Pta[t,a.idx])
            #if died
            Lit[i,t] <<- (1-Xit[i,t])*(1-Sit[i,t])*rbinom(1,1,Fta[t,a.idx]/Zta[t,a.idx]) #if died, was it was due to fishing? F/Z
            Dit[i,t] <<- Lit[i,t]*rbinom(1,1,min(c(1,dmort.bam*Uta[t,a.idx]/Fta[t,a.idx])))*rbinom(1,1,dprop.m[t]) #if died due to fishing, was it discarding? Fdiscard/Ftotal, Fdiscard=U*Dmort
            Wit[i,t] <<- Lit[i,t]*(1-Dit[i,t]) #if not discarded, then it was harvested
            Rit[i,t] <<- (1-Xit[i,t])*(1-Sit[i,t])*rbinom(1,1,1-exp(-Uta[t,a.idx]))*rbinom(1,1,dprop.m[t]) #probability resighted alive before death
            Vit[i,t] <<- ifelse(Kit[i,t]==1|Lit[i,t]==1|Rit[i,t]==1,rbinom(1,1,rep_mo[t]),0) #was it reported if resighted or harvested
          })
          }
          if (tag_type == 'gene'){
            #run fates function for new marks
            sapply(i.idx,function(i){
              a.t <<- Tdat[i,'age.t']
              a.idx <<- match(a.t,ages.m)
              
              #first timestep
              if(t==1){
                Pit[i,t] <<- 1  #capture probability = 1 as this is a new mark   
                Xit[i,t] <<- rbinom(1,1,handle.mort*Pit[i,t]) #handling mortality, 1=yes
                Sit[i,t] <<- (1-Xit[i,t])*rbinom(1,1,exp(-(Zta[t,a.idx])))    #survive timestep, 1=yes
                Tit[i,t] <<- Sit[i,t]  #update tag state
              } 
              if(t>1){
                #new marks in subsequent timesteps
                if(is.na(Tit[i,t-1])){
                  Pit[i,t] <<- 1 #capture probability = 1 as this is a new mark   
                  Xit[i,t] <<- rbinom(1,1,handle.mort*Pit[i,t]) #handling mortality, 1=yes
                  Sit[i,t] <<- (1-Xit[i,t])*rbinom(1,1,exp(-(Zta[t,a.idx])))       #survive current timestep, 1=yes
                  Tit[i,t] <<- Sit[i,t]  #update tag state
                  #existing marks
                } else {
                  Pit[i,t] <<- Tit[i,t-1]*rbinom(1,1,Pta[t,a.idx])*rbinom(1,1,Rta[t,a.idx]/Nta[t,a.idx]) #recapture probability, 1=yes
                  Xit[i,t] <<- rbinom(1,1,handle.mort*Pit[i,t]) #handling mortality, 1=yes
                  Sit[i,t] <<- (1-Xit[i,t])*Tit[i,t-1]*rbinom(1,1,exp(-(Zta[t,a.idx])))      #survive current timestep, 1=yes
                  Tit[i,t] <<- Tit[i,t-1]*Sit[i,t]  #update tag state
                }
              }
              #if survived
              Kit[i,t] <<- Sit[i,t]*rbinom(1,1,1-exp(-Uta[t,a.idx]))*rbinom(1,1,dprop.m[t]) #if survived, was it captured and released alive, 1=yes
              #if died
              Lit[i,t] <<- (1-Xit[i,t])*(1-Sit[i,t])*rbinom(1,1,Fta[t,a.idx]/Zta[t,a.idx]) #if died, was it was due to fishing? F/Z
              Dit[i,t] <<- Lit[i,t]*rbinom(1,1,min(c(1,dmort.bam*Uta[t,a.idx]/Fta[t,a.idx])))*rbinom(1,1,dprop.m[t]) #if died due to fishing, was it discarding? Fdiscard/Ftotal, Fdiscard=U*Dmort
              Wit[i,t] <<- Lit[i,t]*(1-Dit[i,t]) #if not discarded, then it was harvested
              Rit[i,t] <<- (1-Xit[i,t])*(1-Sit[i,t])*rbinom(1,1,1-exp(-Uta[t,a.idx]))*rbinom(1,1,dprop.m[t]) #probability resighted alive before death
              Vit[i,t] <<- ifelse(Kit[i,t]==1|Lit[i,t]==1|Rit[i,t]==1,rbinom(1,1,Edotprop.m[t]),0) #was it reported if resighted or harvested
            })
            }
        }
        #update tagged pop and proportion with those that survived this timestep----------------------------
        Tta[t,match(unique(Tdat[,'age.t']),ages.m)] <- as.numeric(by(Tit[,t],Tdat[,'age.t'],sum,na.rm=T))
        Pta[t,] <- Tta[t,]/Nta[t,]
        Pta[Pta>1] <- 1
      }
      
      #remaining time steps (t>1)=========================================================================
      if(t>1){
        print(paste('running time step',t,' ----------'));flush.console()
        Nta[t,1] <- exp(rnorm(1,mean=log(Ro),sd=Rsigma))*spawn.m[t] #new age-1 recruits for time t
        Fta[t,] <- qF.m[t]*effort.m[t]*selex.m #fishing mortality at age
        Uta[t,] <- qU.m[t]*effort.m[t]*selex.m #exploitation rate at age
        Zta[t,] <- Fta[t,] + Ma.m #total mortality at age
        
        #population age structure accounting----------------------------------------------------------------
        for(a in 2:length(ages.m)){
          Nta[t,a] <- ceiling(Nta[t-1,a-1] * exp(-Zta[t-1,a-1])) #total pop
        } #end age loop
        Cta[t,] <- ceiling(Nta[t,]*Uta[t,]) #catch (includes discards)
        
        #individual tag accounting for existing tags at large-----------------------------------------------
        print('--determining fates of existing marked fish');flush.console()
        #identify fish that are still at large
        doT.idx <- which(Tit[,t-1]==1)
        #for all others, set values to zero
        Tit[-doT.idx,t] <- Pit[-doT.idx,t] <- Xit[-doT.idx,t] <- Sit[-doT.idx,t] <- Vit[-doT.idx,t] <- Wit[-doT.idx,t] <- Rit[-doT.idx,t] <- Dit[-doT.idx,t] <- Kit[-doT.idx,t] <- Lit[-doT.idx,t] <- Hit[-doT.idx,t] <- 0                 #expired tags not captured, so this is 0 always
        
        #advance age for fish at large - this could be changed, but working for now 
        Tdat[,'age.t'] <- ages.m[ifelse((match(Tdat[,'age.marked'],ages.m)+(t-1))>=length(ages.m),length(ages.m),
                                       match(Tdat[,'age.marked'],ages.m)+(t-1))]
        
        if (tag_type == 'conv'){
          sapply(doT.idx, function(i){
          a.t <<- Tdat[i,'age.t']
          a.idx <<- match(a.t,ages.m)
          
          #first timestep 
          if(t==1){
            Pit[i,t] <<- 1  #capture probability = 1 as this is a new mark   
            Xit[i,t] <<- rbinom(1,1,handle.mort*Pit[i,t]) #handling mortality, 1=yes
            Sit[i,t] <<- (1-Xit[i,t])*rbinom(1,1,exp(-(Zta[t,a.idx])))    #survive timestep, 1=yes
            Tit[i,t] <<- Sit[i,t]  #update tag state
          } 
          if(t>1){
            #new marks in subsequent timesteps
            if(is.na(Tit[i,t-1])){
              Pit[i,t] <<- 1 #capture probability = 1 as this is a new mark   
              Xit[i,t] <<- rbinom(1,1,handle.mort*Pit[i,t]) #handling mortality, 1=yes
              Sit[i,t] <<- (1-Xit[i,t])*rbinom(1,1,exp(-(Zta[t,a.idx])))       #survive current timestep, 1=yes
              Tit[i,t] <<- Sit[i,t]  #update tag state
              #existing marks
            } else {
              Pit[i,t] <<- Tit[i,t-1]*rbinom(1,1,Pta[t,a.idx])*rbinom(1,1,Rta[t,a.idx]/Nta[t,a.idx]) #recapture probability, 1=yes
              Xit[i,t] <<- rbinom(1,1,handle.mort*Pit[i,t]) #handling mortality, 1=yes
              Sit[i,t] <<- (1-Xit[i,t])*Tit[i,t-1]*rbinom(1,1,exp(-(Zta[t,a.idx])))      #survive current timestep, 1=yes
              Tit[i,t] <<- Tit[i,t-1]*Sit[i,t]  #update tag state
            }
          }
          #if survived
          Kit[i,t] <<- Sit[i,t]*rbinom(1,1,1-exp(-Uta[t,a.idx]))*rbinom(1,1,dprop.m[t]) 
          #if died
          Lit[i,t] <<- (1-Xit[i,t])*(1-Sit[i,t])*rbinom(1,1,Fta[t,a.idx]/Zta[t,a.idx]) #if died, was it was due to fishing? F/Z
          Dit[i,t] <<- Lit[i,t]*rbinom(1,1,min(c(1,dmort.bam*Uta[t,a.idx]/Fta[t,a.idx])))*rbinom(1,1,dprop.m[t]) #if died due to fishing, was it discarding? Fdiscard/Ftotal, Fdiscard=U*Dmort
          Wit[i,t] <<- Lit[i,t]*(1-Dit[i,t]) #if not discarded, then it was harvested
          Rit[i,t] <<- (1-Xit[i,t])*(1-Sit[i,t])*rbinom(1,1,1-exp(-Uta[t,a.idx]))*rbinom(1,1,dprop.m[t]) #probability resighted alive before death
          Vit[i,t] <<- ifelse(Kit[i,t]==1|Lit[i,t]==1|Rit[i,t]==1,rbinom(1,1,rep_mo[t]),0) #was it reported if resighted or harvested
          })
        }
        if (tag_type == 'gene'){
          #run fates function for each tag at large
          sapply(doT.idx,function(i){
          a.t <<- Tdat[i,'age.t']
          a.idx <<- match(a.t,ages.m)
          
          #first timestep
          if(t==1){
            Pit[i,t] <<- 1  #capture probability = 1 as this is a new mark   
            Xit[i,t] <<- rbinom(1,1,handle.mort*Pit[i,t]) #handling mortality, 1=yes
            Sit[i,t] <<- (1-Xit[i,t])*rbinom(1,1,exp(-(Zta[t,a.idx])))    #survive timestep, 1=yes
            Tit[i,t] <<- Sit[i,t]  #update tag state
          } 
          if(t>1){
            #new marks in subsequent timesteps
            if(is.na(Tit[i,t-1])){
              Pit[i,t] <<- 1 #capture probability = 1 as this is a new mark   
              Xit[i,t] <<- rbinom(1,1,handle.mort*Pit[i,t]) #handling mortality, 1=yes
              Sit[i,t] <<- (1-Xit[i,t])*rbinom(1,1,exp(-(Zta[t,a.idx])))       #survive current timestep, 1=yes
              Tit[i,t] <<- Sit[i,t]  #update tag state
              #existing marks
            } else {
              Pit[i,t] <<- Tit[i,t-1]*rbinom(1,1,Pta[t,a.idx])*rbinom(1,1,Rta[t,a.idx]/Nta[t,a.idx]) #recapture probability, 1=yes
              Xit[i,t] <<- rbinom(1,1,handle.mort*Pit[i,t]) #handling mortality, 1=yes
              Sit[i,t] <<- (1-Xit[i,t])*Tit[i,t-1]*rbinom(1,1,exp(-(Zta[t,a.idx])))      #survive current timestep, 1=yes
              Tit[i,t] <<- Tit[i,t-1]*Sit[i,t]  #update tag state
            }
          }
          #if survived
          Kit[i,t] <<- Sit[i,t]*rbinom(1,1,1-exp(-Uta[t,a.idx]))*rbinom(1,1,dprop.m[t]) #if survived, was it captured and released alive, 1=yes
          #if died
          Lit[i,t] <<- (1-Xit[i,t])*(1-Sit[i,t])*rbinom(1,1,Fta[t,a.idx]/Zta[t,a.idx]) #if died, was it was due to fishing? F/Z
          Dit[i,t] <<- Lit[i,t]*rbinom(1,1,min(c(1,dmort.bam*Uta[t,a.idx]/Fta[t,a.idx])))*rbinom(1,1,dprop.m[t]) #if died due to fishing, was it discarding? Fdiscard/Ftotal, Fdiscard=U*Dmort
          Wit[i,t] <<- Lit[i,t]*(1-Dit[i,t]) #if not discarded, then it was harvested
          Rit[i,t] <<- (1-Xit[i,t])*(1-Sit[i,t])*rbinom(1,1,1-exp(-Uta[t,a.idx]))*rbinom(1,1,dprop.m[t]) #probability resighted alive before death
          Vit[i,t] <<- ifelse(Kit[i,t]==1|Lit[i,t]==1|Rit[i,t]==1,rbinom(1,1,Edotprop.m[t]),0) #was it reported if resighted or harvested
          })
          
        }
        
        #new captures---------------------------------------------------------------------------------------
        print('--determining fates of new marked fish');flush.console()
        if(tag_type == "conv"){
        Rt[t] <- tag_dist[t] 
        }
        if(tag_type == "gene"){
          #participant recaps at age
          caps = ceiling(qU.m[t]*Edot.m[t]*sum(selex.m*Nta[t,])) #total captures for timestep by participants, given U, effort, selex, and Na
          recaps = sum(Pit[,t],na.rm=T)  #number of recaps in timestep
          newcaps = max(0,caps-recaps) #number of new marked fish released by participants (total caps minus recaps)
          Rt[t] = newcaps 
        }
        if(Rt[t]>0){
          Rta[t,] <- rmultinom(n=1, size=Rt[t], Cta[t,]/sum(Cta[t,])) #new caps at age
        } else{
          Rta[t,] <- 0
        }
        
        #newly marked individuals
        T1 <- data.frame(tstep = rep(t,sum(Rta[t,])),age = rep(ages.m,times=Rta[t,]))
        T1$a.t <- ages.m[match(T1$age,ages.m)]
        Tdat <- rbind(Tdat,as.matrix(T1))
        
        #expand matrices to store fates
        Sit <- as.matrix(rbind(Sit,(matrix(nrow=sum(Rta[t,]),ncol=ntsteps))))
        Tit <- as.matrix(rbind(Tit,(matrix(nrow=sum(Rta[t,]),ncol=ntsteps))))
        Lit <- as.matrix(rbind(Lit,(matrix(nrow=sum(Rta[t,]),ncol=ntsteps))))
        Kit <- as.matrix(rbind(Kit,(matrix(nrow=sum(Rta[t,]),ncol=ntsteps))))
        Wit <- as.matrix(rbind(Wit,(matrix(nrow=sum(Rta[t,]),ncol=ntsteps))))
        Dit <- as.matrix(rbind(Dit,(matrix(nrow=sum(Rta[t,]),ncol=ntsteps))))
        Vit <- as.matrix(rbind(Vit,(matrix(nrow=sum(Rta[t,]),ncol=ntsteps))))
        Pit <- as.matrix(rbind(Pit,(matrix(nrow=sum(Rta[t,]),ncol=ntsteps))))
        Rit <- as.matrix(rbind(Rit,(matrix(nrow=sum(Rta[t,]),ncol=ntsteps))))
        Xit <- as.matrix(rbind(Xit,(matrix(nrow=sum(Rta[t,]),ncol=ntsteps))))
        
        #fates of new releases: Survival to next timestep is based on release mortality.
        #needs to be written to fucntion "fn.newcaps"
        i.idx <- numeric(0)
        if(Rt[t]>0) i.idx <- (nrow(Tdat)-nrow(T1)+1):nrow(Tdat)
        Tit[i.idx,t] <- 1  #set initial tag state to 1
            
        #run fates function for new marks
        if(length(i.idx)>0){
         
          if (tag_type == 'conv'){
             sapply(i.idx, function(i){
               a.t <<- Tdat[i,'age.t']
            a.idx <<- match(a.t,ages.m)
            
            #first timestep 
            if(t==1){
              Pit[i,t] <<- 1  #capture probability = 1 as this is a new mark   
              Xit[i,t] <<- rbinom(1,1,handle.mort*Pit[i,t]) #handling mortality, 1=yes
              Sit[i,t] <<- (1-Xit[i,t])*rbinom(1,1,exp(-(Zta[t,a.idx])))    #survive timestep, 1=yes
              Tit[i,t] <<- Sit[i,t]  #update tag state
            } 
            if(t>1){
              #new marks in subsequent timesteps
              if(is.na(Tit[i,t-1])){
                Pit[i,t] <<- 1 #capture probability = 1 as this is a new mark   
                Xit[i,t] <<- rbinom(1,1,handle.mort*Pit[i,t]) #handling mortality, 1=yes
                Sit[i,t] <<- (1-Xit[i,t])*rbinom(1,1,exp(-(Zta[t,a.idx])))       #survive current timestep, 1=yes
                Tit[i,t] <<- Sit[i,t]  #update tag state
                #existing marks
              } else {
                Pit[i,t] <<- Tit[i,t-1]*rbinom(1,1,Pta[t,a.idx])*rbinom(1,1,Rta[t,a.idx]/Nta[t,a.idx]) #recapture probability, 1=yes
                Xit[i,t] <<- rbinom(1,1,handle.mort*Pit[i,t]) #handling mortality, 1=yes
                Sit[i,t] <<- (1-Xit[i,t])*Tit[i,t-1]*rbinom(1,1,exp(-(Zta[t,a.idx])))      #survive current timestep, 1=yes
                Tit[i,t] <<- Tit[i,t-1]*Sit[i,t]  #update tag state
              }
            }
            #if survived
            Kit[i,t] <<- Sit[i,t]*rbinom(1,1,1-exp(-Uta[t,a.idx]))*rbinom(1,1,dprop.m[t]) #*rbinom(1,1,Pta[t,a.idx]) #if survived, was it captured and released alive, 1=yes  #LK added *rbinom(1,1,Pta[t,a.idx]) new fate func
            #if died
            Lit[i,t] <<- (1-Xit[i,t])*(1-Sit[i,t])*rbinom(1,1,Fta[t,a.idx]/Zta[t,a.idx]) #if died, was it was due to fishing? F/Z
            Dit[i,t] <<- Lit[i,t]*rbinom(1,1,min(c(1,dmort.bam*Uta[t,a.idx]/Fta[t,a.idx])))*rbinom(1,1,dprop.m[t]) #if died due to fishing, was it discarding? Fdiscard/Ftotal, Fdiscard=U*Dmort
            Wit[i,t] <<- Lit[i,t]*(1-Dit[i,t]) #if not discarded, then it was harvested
            Rit[i,t] <<- (1-Xit[i,t])*(1-Sit[i,t])*rbinom(1,1,1-exp(-Uta[t,a.idx]))*rbinom(1,1,dprop.m[t]) #probability resighted alive before death
            Vit[i,t] <<- ifelse(Kit[i,t]==1|Lit[i,t]==1|Rit[i,t]==1,rbinom(1,1,rep_mo[t]),0) #was it reported if resighted or harvested
             }) 
            }
          if (tag_type == 'gene'){
            #run fates function for each tag at large
            sapply(i.idx,function(i){
              a.t <<- Tdat[i,'age.t']
              a.idx <<- match(a.t,ages.m)
              
              #first timestep
              if(t==1){
                Pit[i,t] <<- 1  #capture probability = 1 as this is a new mark   
                Xit[i,t] <<- rbinom(1,1,handle.mort*Pit[i,t]) #handling mortality, 1=yes
                Sit[i,t] <<- (1-Xit[i,t])*rbinom(1,1,exp(-(Zta[t,a.idx])))    #survive timestep, 1=yes
                Tit[i,t] <<- Sit[i,t]  #update tag state
              } 
              if(t>1){
                #new marks in subsequent timesteps
                if(is.na(Tit[i,t-1])){
                  Pit[i,t] <<- 1 #capture probability = 1 as this is a new mark   
                  Xit[i,t] <<- rbinom(1,1,handle.mort*Pit[i,t]) #handling mortality, 1=yes
                  Sit[i,t] <<- (1-Xit[i,t])*rbinom(1,1,exp(-(Zta[t,a.idx])))       #survive current timestep, 1=yes
                  Tit[i,t] <<- Sit[i,t]  #update tag state
                  #existing marks
                } else {
                  Pit[i,t] <<- Tit[i,t-1]*rbinom(1,1,Pta[t,a.idx])*rbinom(1,1,Rta[t,a.idx]/Nta[t,a.idx]) #recapture probability, 1=yes
                  Xit[i,t] <<- rbinom(1,1,handle.mort*Pit[i,t]) #handling mortality, 1=yes
                  Sit[i,t] <<- (1-Xit[i,t])*Tit[i,t-1]*rbinom(1,1,exp(-(Zta[t,a.idx])))      #survive current timestep, 1=yes
                  Tit[i,t] <<- Tit[i,t-1]*Sit[i,t]  #update tag state
                }
              }
              #if survived
              Kit[i,t] <<- Sit[i,t]*rbinom(1,1,1-exp(-Uta[t,a.idx]))*rbinom(1,1,dprop.m[t]) #if survived, was it captured and released alive, 1=yes
              #if died
              Lit[i,t] <<- (1-Xit[i,t])*(1-Sit[i,t])*rbinom(1,1,Fta[t,a.idx]/Zta[t,a.idx]) #if died, was it was due to fishing? F/Z
              Dit[i,t] <<- Lit[i,t]*rbinom(1,1,min(c(1,dmort.bam*Uta[t,a.idx]/Fta[t,a.idx])))*rbinom(1,1,dprop.m[t]) #if died due to fishing, was it discarding? Fdiscard/Ftotal, Fdiscard=U*Dmort
              Wit[i,t] <<- Lit[i,t]*(1-Dit[i,t]) #if not discarded, then it was harvested
              Rit[i,t] <<- (1-Xit[i,t])*(1-Sit[i,t])*rbinom(1,1,1-exp(-Uta[t,a.idx]))*rbinom(1,1,dprop.m[t]) #probability resighted alive before death
              Vit[i,t] <<- ifelse(Kit[i,t]==1|Lit[i,t]==1|Rit[i,t]==1,rbinom(1,1,Edotprop.m[t]),0) #was it reported if resighted or harvested
            })
          }
          
   
        }
        #update tagged pop and proportion with those that survived this timestep----------------------------
        Tta[t,match(unique(Tdat[,'age.t']),ages.m)] <- as.numeric(by(Tit[,t],Tdat[,'age.t'],sum,na.rm=T))
        Pta[t,] <- Tta[t,]/Nta[t,]
        Pta[Pta>1] <- 1
        Kta[t,] <- ceiling(Cta[t,]*Pta[t,]) 
        
        #runtime
        if(t==ntsteps){
          time.end <- Sys.time()
          print(time.end-time.start);flush.console()
        }
        
      } #end if t>1
    } #end t loop 
    
  
    #capture histories==================================================================================
    Hit <- Tit; Hit[] = NA 
    sapply(1:ncol(Hit),function(t){ 
    Hit[which(Tit[,t]==1 & Sit[,t]==1 & Xit[,t]==0 & Kit[,t]==0 & Lit[,t]==0 & Dit[,t]==0 & Wit[,t]==0 & Rit[,t]==0 & Vit[,t]==0),t] <<-   0  #survive, not encountered
    Hit[which(Tit[,t]==1 & Sit[,t]==1 & Xit[,t]==0 & Kit[,t]==1 & Lit[,t]==0 & Dit[,t]==0 & Wit[,t]==0 & Rit[,t]==0 & Vit[,t]==1),t] <<-   1  #survive, resighted, reported
    Hit[which(Tit[,t]==1 & Sit[,t]==1 & Xit[,t]==0 & Kit[,t]==1 & Lit[,t]==0 & Dit[,t]==0 & Wit[,t]==0 & Rit[,t]==0 & Vit[,t]==0),t] <<-  -1 #survive, resighted, not reported
  
    Hit[which(Tit[,t]==0 & Sit[,t]==0 & Xit[,t]==0 & Kit[,t]==0 & Lit[,t]==1 & Dit[,t]==1 & Wit[,t]==0 & Rit[,t]==0 & Vit[,t]==1),t] <<-   2 #dead, discard, reported
    Hit[which(Tit[,t]==0 & Sit[,t]==0 & Xit[,t]==0 & Kit[,t]==0 & Lit[,t]==1 & Dit[,t]==1 & Wit[,t]==0 & Rit[,t]==0 & Vit[,t]==0),t] <<-  -2 #dead, discard, not reported
    Hit[which(Tit[,t]==0 & Sit[,t]==0 & Xit[,t]==0 & Kit[,t]==0 & Lit[,t]==1 & Dit[,t]==1 & Wit[,t]==0 & Rit[,t]==1 & Vit[,t]==1),t] <<-  21 #dead, discard, resighted prior, reported 
    Hit[which(Tit[,t]==0 & Sit[,t]==0 & Xit[,t]==0 & Kit[,t]==0 & Lit[,t]==1 & Dit[,t]==1 & Wit[,t]==0 & Rit[,t]==1 & Vit[,t]==0),t] <<- -21 #dead, discard, resighted prior, not reported 
  
    Hit[which(Tit[,t]==0 & Sit[,t]==0 & Xit[,t]==0 & Kit[,t]==0 & Lit[,t]==1 & Dit[,t]==0 & Wit[,t]==1 & Rit[,t]==0 & Vit[,t]==1),t] <<-   3 #dead, harvest, reported
    Hit[which(Tit[,t]==0 & Sit[,t]==0 & Xit[,t]==0 & Kit[,t]==0 & Lit[,t]==1 & Dit[,t]==0 & Wit[,t]==1 & Rit[,t]==0 & Vit[,t]==0),t] <<-  -3 #dead, harvest, not reported
    Hit[which(Tit[,t]==0 & Sit[,t]==0 & Xit[,t]==0 & Kit[,t]==0 & Lit[,t]==1 & Dit[,t]==0 & Wit[,t]==1 & Rit[,t]==1 & Vit[,t]==1),t] <<-  31 #dead, harvest, resighted prior, reported
    Hit[which(Tit[,t]==0 & Sit[,t]==0 & Xit[,t]==0 & Kit[,t]==0 & Lit[,t]==1 & Dit[,t]==0 & Wit[,t]==1 & Rit[,t]==1 & Vit[,t]==0),t] <<- -31 #dead, harvest, resighted prior, not reported
  
    Hit[which(Tit[,t]==0 & Sit[,t]==0 & Xit[,t]==0 & Kit[,t]==0 & Lit[,t]==0 & Dit[,t]==0 & Wit[,t]==0 & Rit[,t]==0 & Vit[,t]==0),t] <<-   4 #natural mortality, not resighted
    Hit[which(Tit[,t]==0 & Sit[,t]==0 & Xit[,t]==0 & Kit[,t]==0 & Lit[,t]==0 & Dit[,t]==0 & Wit[,t]==0 & Rit[,t]==1 & Vit[,t]==1),t] <<-  41 #natural mortality, resighted prior, reported
    Hit[which(Tit[,t]==0 & Sit[,t]==0 & Xit[,t]==0 & Kit[,t]==0 & Lit[,t]==0 & Dit[,t]==0 & Wit[,t]==0 & Rit[,t]==1 & Vit[,t]==0),t] <<- -41 #natural mortality, resighted prior, not reported
    Hit[which(Tit[,t]==0 & Sit[,t]==0 & Xit[,t]==1 & Kit[,t]==0 & Lit[,t]==0 & Dit[,t]==0 & Wit[,t]==0 & Rit[,t]==0 & Vit[,t]==0),t] <<-   5 #handling mortality
    if(t>1){
      Hit[which(Tit[,t-1]==0),t] <<- 0 
    }
    })
    #} #end t loop 
    

    
    
    caphist = data.frame(tstep.marked = Tdat[,'tstep.marked'],age.marked=Tdat[,'age.marked'],
                         tstep.last=NA, age.dead=NA, Tit.sum=NA, Sit.sum = NA, Kit.sum=NA, Lit.sum=NA, Dit.sum=NA, Wit.sum=NA, Rit.sum=NA, 
                         Vit.sum=NA, Pit.sum=NA, Xit.sum=NA, Tit.last=NA, fate=NA, fate2=NA)
    for(i in 1:nrow(caphist)){
      caphist$tstep.last[i] = ifelse(length(which(Tit[i,]==0))==0,ntsteps,min(which(Tit[i,]==0)))
      caphist$age.dead[i] = ages.m[min(match(caphist$age.marked[i],ages.m) + (caphist$tstep.last[i]-caphist$tstep.marked[i]),length(ages.m))]
      caphist$Tit.sum[i] = sum(Tit[i,],na.rm=T)
      caphist$Sit.sum[i] = sum(Sit[i,],na.rm=T)
      caphist$Kit.sum[i] = sum(Kit[i,],na.rm=T)
      caphist$Lit.sum[i] = sum(Lit[i,],na.rm=T)
      caphist$Dit.sum[i] = sum(Dit[i,],na.rm=T)
      caphist$Wit.sum[i] = sum(Wit[i,],na.rm=T)
      caphist$Rit.sum[i] = sum(Rit[i,],na.rm=T)
      caphist$Vit.sum[i] = sum(Vit[i,],na.rm=T)
      caphist$Pit.sum[i] = sum(Pit[i,],na.rm=T)
      caphist$Xit.sum[i] = sum(Xit[i,],na.rm=T)
      caphist$Tit.last[i] = Tit[i,ncol(Tit)]
      caphist$fate[i] = Hit[i,caphist$tstep.last[i]]  
    }
    
    caphist$fate2 = ifelse(caphist$fate %in% c(3,-3,31,-31),'harvest',
                           ifelse(caphist$fate %in% c(2,-2,21,-21),'dmort',#dead discard
                                  ifelse(caphist$fate %in% c(4,-4,41,-41),'natmort',
                                         ifelse(caphist$Tit.last==1 & caphist$Kit.sum>=1,'surv.resighted', #tit.last ==1 means tag is still at large in last timestep
                                                ifelse(caphist$Tit.last==1 & caphist$Kit.sum==0,'surv.notresighted',
                                                       ifelse(caphist$fate==5, 'handmort','survived'))))))
    
    caphist$reporting = ifelse(caphist$fate %in% c(-1,-2,-3,-21,-31,-41),'cap_notrep', #re-sighted - not reported (dead or alive)
                               ifelse(caphist$fate %in% c(1,2,3,21,31,41),'cap_rep', #re-sighted - reported (dead or alive)
                                      ifelse(caphist$fate %in% c(0),'not_recap', #survived - not encountered (survived)
                                             ifelse(caphist$fate%in% c(4,5),'not_recap','err')))) #this needs to be updated, 3 is tagged and killed in the same time step, similar to handling mortality 
    
    
    #write.csv(caphist,paste0("~/GitHub/tag_sim/output/update-srfs-tagging-95rep","/capture_history_update_95rep.csv"),row.names = F)
    #write.csv(caphist,paste0("~/GitHub/tag_sim/output/update_gene_full","/capture_history_update_ptrip0.01.csv"),row.names = F)
    
    #timeseries=========================================================================================
    rm(t.series)
    t.series = data.frame(
      tstep = 1:ntsteps,
      Nt = rowSums(Nta),
      NtVul = rowSums(sweep(Nta,2,selex.m,"*")),
      Ct = rowSums(Cta),
      Ut = Ut,
      Ft = Ft,
      E = effort.m,
      Edot = Edot.m,
      recaps.efp = NA,
      recaps.tot  = apply(Hit,2,FUN=function(x)length(which(abs(x)%in%c(11,12,13)))),
      P.recap.efp   = NA,
      P.recap.tot   = NA,
      PTit = round(colSums(Tit,na.rm=T)/rowSums(Nta),4),
      PTitVul = round(colSums(Tit,na.rm=T)/rowSums(sweep(Nta,2,selex.m,"*")),4),
      cummarks      = apply(Tit,2,FUN=function(x)length(x[!is.na(x)])),
      Sit.sum = colSums(Sit,na.rm=T),
      Lit.sum = colSums(Lit,na.rm=T),
      Wit.sum = colSums(Wit,na.rm=T),
      Kit.sum = colSums(Kit,na.rm=T),
      Dit.sum = colSums(Dit,na.rm=T),
      Vit.sum = colSums(Vit,na.rm=T),
      Xit.sum = colSums(Xit,na.rm=T),
      Pit.sum = colSums(Pit,na.rm=T),
      Rit.sum = colSums(Rit,na.rm=T),
      Tit.sum = colSums(Tit,na.rm=T)
    )

    t.series$P.recap.tot = t.series$recaps.tot/t.series$Ct
    t.series$recaps.efp[1] = t.series$recaps.tot[1] = t.series$P.recap.efp[1] = t.series$P.recap.tot[1] = 0
    
    #output=============================================================================================

    write.csv(t.series,paste0(run.list$full.dir[d],"/tseries.csv"),row.names = F)
    write.csv(Hit,paste0(run.list$full.dir[d],"/capture_history.csv"),row.names = F)
    
    
    ####################################################################################################
    ##                          AGGREGATE FOR ANNUAL INTERVAL
    ####################################################################################################
    
    if(timestep == 'annual'){ #needs reworking to be functional with gene and not 3 years
      
      dim(Hit)
      tsteps.sub = 1:36
      Tdat2 = as.data.frame(Tdat)
      Tdat2 = Tdat2[Tdat2$tstep.marked%in%tsteps.sub,]  #subset time steps
      names(Tdat2)[2] = 'age.marked.yr'
      Tdat2$age.t = NULL
      Tdat2$age.marked.mo = as.factor(round(Tdat2$age.marked.yr*12,0))
      Tdat2$age.marked.intyr = as.factor(floor(Tdat2$age.marked.yr))

      #build LDLD capture histories for an annual interval model----------------------------------
      col.idx = ncol(Tdat2)+1
      Tdat2 = cbind(Tdat2,matrix(NA,nrow=nrow(Tdat2),ncol=nyrs*2))
      names(Tdat2)[((ncol(Tdat2)-nyrs*2)+1):ncol(Tdat2)] = paste0(rep(c('L','D'),nyrs),rep(1:nyrs,each=2))
      
      Lmat = Dmat = matrix(NA,nrow=nrow(Tdat2),ncol=nyrs)
      
      #these two line recode Hit fates to L and D
      for(y in 1:nyrs){
        y.idx = ((y-1)*12+1):(y*12)
        Pit2 = Pit[,y.idx]
        Hit2 = Hit[,y.idx]
        
        Lmat[,y] = apply(Pit2,1,function(x) ifelse(sum(x,na.rm=T)>=1,1,0))          #captured or recaptured during tagging, L=1
        table(Lmat[,y])
        

        #in reality, this is how it would be received:  All reported discards are considered live resightings,
        #and only reported harvest gets a one in the capture history 
        Dmat[,y] = t(apply(Hit2,1,function(x) ifelse(length(which(x%in%c(3,31)))>0,1,  #harvest, D=1
                                                     ifelse(length(which(x%in%c(1,2,21,41)))>0,2,0)))) #live resighting, D=2
      }
      
      
      Lmat[is.na(Lmat)] = 0
      Dmat[is.na(Dmat)] = 0
      dim(Lmat)
      Tdat2[,seq(col.idx,ncol(Tdat2),2)] = Lmat
      Tdat2[,seq(col.idx+1,ncol(Tdat2),2)] = Dmat
      
    } 
    
    ####################################################################################################
    ##                          AGGREGATE FOR MONTHLY INTERVAL
    ####################################################################################################
    
    if (timestep == "monthly") {
      tsteps.sub = 1:ntsteps
      Tdat2 = as.data.frame(Tdat)
      Tdat2 = Tdat2[Tdat2$tstep.marked%in%tsteps.sub,]  #subset time steps
      names(Tdat2)[2] = 'age.marked.yr'
      Tdat2$age.t = NULL
      Tdat2$age.marked.mo = as.factor(round(Tdat2$age.marked.yr*12,0))
      Tdat2$age.marked.intyr = as.factor(floor(Tdat2$age.marked.yr))
      Hit2 = Hit[,tsteps.sub]
      Hit2 = Hit2[rowSums(is.na(Hit2)) != ncol(Hit2),]
      Pit2 = Pit[1:nrow(Hit2),1:ncol(Hit2)]
      
      #build LDLD capture histories for a 1-year, monthly interval model----------------------------------
      col.idx = ncol(Tdat2)+1
      Tdat2 = cbind(Tdat2,matrix(NA,nrow=nrow(Tdat2),ncol=ntsteps*2))
      names(Tdat2)[(ncol(Tdat2)-(ntsteps*2)+1):ncol(Tdat2)] = paste0(rep(c('L','D'),ntsteps),rep(1:ntsteps,each=2))
      #these two line recode Hit fates to L and D
      Lmat = t(apply(Pit2,1,function(x) ifelse(x==1,1,0)))          #captured, L=1
      
      #in reality, this is how it would be received:
      Dmat = t(apply(Hit2,1,function(x) ifelse(x%in%c(3,31),1,  #harvest, D=1
                                               ifelse(x%in%c(1,2,21,41),2,0)))) #live resighting, D=2
      
      Lmat[is.na(Lmat)] = 0
      Dmat[is.na(Dmat)] = 0
      
      Tdat2[,seq(col.idx,ncol(Tdat2),2)] = Lmat
      Tdat2[,seq(col.idx+1,ncol(Tdat2),2)] = Dmat
    }
    
    # this removes all recaptures prior to run
    if(tag_after_cap == "clip"){
      Tdat2 = clipped_tag(Tdat2)
    }
    write.csv(Tdat2,paste0(run.list$full.dir[d],'/barker_dat.csv'),row.names = F)
    
    ################ PLOTTING ################
        pdf(file = paste0(run.list$full.dir[d],'/figures_marks-',run.list$marks[d],'run_',run.list$trial[d],'.pdf'))
    
    par(mfrow=c(3,2))
    plot(NtVul~tstep,t.series,type='l',ylim=c(.8*min(t.series$NtVul,na.rm=T),1.2*max(t.series$NtVul,na.rm=T)),main='A: Vulnerable Numbers',lwd=2, ylab = '# Vulnerable')
    plot(t.series$tstep,round(Nta[,1],0),type='l',main='B: Age-1.0 Recruits',lwd=2, ylab = '# Recruits', xlab = 'tstep')
    plot(Ct~tstep,t.series,type='l',main='C: Total Fishery Catch',lwd=2,ylab = 'Catch' )
    plot(Tit.sum~tstep,t.series,type='l',main='E: Marks at Large',ylab='N Marked Fish')
    
    plot(PTitVul~tstep,t.series,type='l',main='F: Proportion Marked',ylab='Proportion',
         ylim=c(0,max(t.series$PTitVul,t.series$P.recap.efp,t.series$P.recap.tot,na.rm=T)))
    lines(P.recap.tot~tstep,t.series,col='gray')
    #fate barplot---------------------------------------------------------------------------------------
    fate.barplot.dat = prop.table(table(caphist$fate2,useNA='ifany'))
    par(mfrow=c(1,1),mar = c(7.5,4,3,2))
    barplot(fate.barplot.dat,ylim=c(0,max(fate.barplot.dat)*1.1),main='Known Fates', las = 2)#,xlab='Fate Code'
    #reporting bar plot --------------------------------------------
    reporting.barplot.dat = prop.table(table(caphist$reporting,useNA='ifany'))
    par(mfrow=c(1,1))
    barplot(reporting.barplot.dat,ylim=c(0,max(reporting.barplot.dat)*1.1),main='Reporting')#,xlab='Fate Code')
    
    dev.off()
    
  } 









