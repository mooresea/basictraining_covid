buildNetworkFun = function(i){

  # cocoon size 30-65
  cocoon.size = P[i,12]
  numcocoons = ceiling(numrecruits/cocoon.size)
  
  # company size up to 240, shared common areas
  company.size = 120
  numcompany = ceiling(numcocoons/round(company.size/cocoon.size))
  
  # assign recruits to cocoons and companies
  recruit.cocoon = sample(rep(1:numcocoons,each = cocoon.size))[1:numrecruits]
  
  # assign recruits to companies
  # set up so that companies are made up of cocoons
  platoon.company = sample(rep(1:numcompany,each = ceiling(numcocoons/numcompany)))[1:numcocoons] 
  recruit.company = numeric()
  for (i in 1:numcompany){
    j = which(platoon.company == i)
    k = which(recruit.cocoon %in% j)
    recruit.company[k] = i
  }
  

  # build contact network
  # simulate contacts among cocoons
  network.cocoon = matrix(0,numrecruits,numrecruits)
  for(ii in 1:numcocoons){
    which.in.cocoon = expand.grid(
      which(recruit.cocoon==ii),which(recruit.cocoon==ii))
    network.cocoon[as.matrix(which.in.cocoon)] =
      rbinom(nrow(which.in.cocoon),1,P[i,13])
  }
  contacts.cocoon = apply(network.cocoon,1,function(x)which(x>0)); rm(network.cocoon)
  
  # simulate contacts among companies
  network.company = matrix(0,numrecruits,numrecruits)
  for(ii in 1:numcompany){
    which.in.company = expand.grid(
      which(recruit.company==ii),which(recruit.company==ii))
    network.company[as.matrix(which.in.company)] =
      rbinom(nrow(which.in.company),1,P[i,14])
  }
  contacts.company = apply(network.company,1,function(x)which(x>0)); rm(network.company)
  
  # simulate random contacts
  network.random = matrix(
    rbinom(numrecruits^2,1,P[i,15]),
    numrecruits,numrecruits)
  network.random = network.random + t(network.random)
  contacts.random = apply(network.random,1,function(x)which(x>0)); rm(network.random)
  
  
  # tally the total contacts of each recruit
  total.contacts =
    unlist(lapply(contacts.random,length)) +
    unlist(lapply(contacts.cocoon,length)) +
    unlist(lapply(contacts.company,length)) 
  
  result = list(contacts.cocoon,contacts.company,contacts.random,total.contacts)
  return(result)

}


###################################################################


simOutbreakFun = function(i) {

  network = buildNetworkFun(i)
  contacts.cocoon = network[[1]]
  contacts.company = network[[2]]
  contacts.random = network[[3]]
  total.contacts = network[[4]]
  
  # parameter values for infection processes
  R0.mag = P[i,1]
  R0.tim = P[i,2]
  symp.mag = P[i,3]
  symp.tim = P[i,4]
  symp.mean = 10
  asymp.adjust = P[i,5]
  initial.infected = P[i,6] # Fort Benning outbreak - 4/640 recruits positive upon arrival
  initial.immune = P[i,6]
  
  R0.prob = R0.mag / mean(total.contacts) / (symp.mag + (1 - symp.mag) * asymp.adjust)
  
  # isolation/quarantine and contact tracing
  BOOL_clinicalrelease = 1
  isolate.max = 14
  isolate.min = 10
  isolate.nosymp = 1
  isolate.effectiveness = P[i,7]
  quarantine.max = 14
  quarantine.effectiveness = P[i,8]
  quarantine.contacts = P[i,9]
  
  # masks
  mask.protection = P[i,10]
  mask.compliance = P[i,11]
  compliance.avg = mask.protection * mask.compliance
  compliance.spread = 1e1
  compliance = rbeta(
    numrecruits,
    compliance.avg * compliance.spread,
    (1 - compliance.avg) * compliance.spread
  )
  
  # specify secondary infection distribution
  R0 = R0.prob * dpois(1:21, R0.tim) / sum(dpois(1:21, R0.tim))
  
  # specify delay between infection and symptom onset
  prop.inf.symp = symp.mag * dpois(1:28, symp.tim) / sum(dpois(1:28, symp.tim))
  
  # previously infected and immune
  immune = rep(0, numrecruits)
  immune[sample(numrecruits, rbinom(1, numrecruits, initial.immune))] = 1
  
  # infected status
  infected = rep(0, numrecruits)
  infected[sample(numrecruits, rbinom(1, sum(immune), initial.infected))] = 1
  init.infect = sum(infected)
  
  # simulate day of infection of initial infecteds
  R0_init = R0.prob * dpois(1:21, R0.tim/2) / sum(dpois(1:21, R0.tim/2)) # specify secondary infection distribution for initial infecteds
  date.infected = rep(NA, numrecruits)
  for (ii in which(infected > 0)) {
    date.infected[ii] = arrivebase[ii] - sample(1:length(R0_init), 1, prob = R0_init)
  }
  
  # add immunes to infected and record date in far past
  infected[which(immune > 0)] = 1
  for (ii in which(immune > 0)) {
    date.infected[ii] = -1e3
  }
  
  # simulate date of symptom onset and duration
  date.symptoms = rep(NA, numrecruits)
  symptom.duration = rep(NA, numrecruits)
  for (ii in which(infected > 0)) {
    if (rbinom(1, 1, symp.mag) == 1) {
      date.symptoms[ii] =
        date.infected[ii] +
        sample(1:length(prop.inf.symp), 1, prob = prop.inf.symp)
      symptom.duration[ii] = rpois(1,symp.mean)
    }
  }
  
  # individual status storage
  testpos = rep(NA, numrecruits)
  tested.date = rep(NA, numrecruits)
  isolate.date = rep(NA, numrecruits)
  quarantine.date = rep(NA, numrecruits)
  init.pos = 0
  
  # timespan of simulation
  time = 1:84
  
  # aggregate numbers over time
  numIsolated = rep(0, length(time))
  numQuarantine = rep(0, length(time))
  numTested = rep(0, length(time))
  numInfected = rep(0, length(time))
  numSymptoms = rep(0, length(time))
  
  # keep track of who infected whom
  edges = cbind(rep(0, sum(infected)), which(infected > 0))
  
  # loop over each day of basic training
  for (tt in time) {
    
    # test on arrival
    if (tt %in% unique(arrivebase) & BOOL_testAfterArrival == 1) {
      needtesting = which(arrivebase == tt)
      trueneg = needtesting[is.na(date.infected[needtesting])]
      truepos = setdiff(needtesting, trueneg)
      testpos = rep(NA, length(needtesting))
      if (length(trueneg) > 0) {
        testpos[which(needtesting %in% trueneg)] =
          rbinom(length(trueneg), 1, 1 - pcr.spec)
      }
      if (length(truepos) > 0) {
        testpos[which(needtesting %in% truepos)] =
          rbinom(length(truepos), 1, pcr.sens(tt - date.infected[truepos] +
                                                1))
      }
      tested.date[needtesting] = tt
      needtracing = needtesting[which(testpos == 1)]
      isolate.date[needtracing] = tt + testDelayArrival + testReturn
      numTested[tt] = numTested[tt] + length(needtesting)
      init.pos = init.pos + length(needtracing)
    }
    
    
    # test on specified days
    if (tt %in% testdates) {
      needtesting = 1:numrecruits
      trueneg = needtesting[is.na(date.infected[needtesting])]
      truepos = setdiff(needtesting, trueneg)
      testpos = rep(NA, length(needtesting))
      if (length(trueneg) > 0) {
        testpos[which(needtesting %in% trueneg)] =
          rbinom(length(trueneg), 1, 1 - pcr.spec)
      }
      if (length(truepos) > 0) {
        testpos[which(needtesting %in% truepos)] =
          rbinom(length(truepos), 1, pcr.sens(tt - date.infected[truepos] +
                                                1))
      }
      tested.date[needtesting] = tt
      needtracing = needtesting[which(testpos == 1)]
      isolate.date[needtracing] = tt + testReturn
      numTested[tt] = numTested[tt] + length(needtesting)
    }
    
    # loop through those who are capable of infecting others
    infectious = which(infected > 0 &
                         date.infected < tt &
                         date.infected > (tt - length(R0)) &
                         arrivebase <= tt)
    for (ii in infectious) {
      # look up probability of infection today
      infect.today =
        R0[tt - date.infected[ii] + 1] *
        ifelse(is.na(date.symptoms[ii]), asymp.adjust, 1) *
        ifelse(
          is.na(isolate.date[ii]),
          1,
          1 - isolate.effectiveness
        ) *
        ifelse(
          is.na(quarantine.date[ii]),
          1,
          1 - quarantine.effectiveness
        ) *
        (1 - compliance[ii])
      
      # look up this person's contacts today
      if (tt <= 14){
        contacts.today = unique(contacts.random[[ii]],contacts.cocoon[[ii]])
      } else {
        contacts.today = unique(c(contacts.random[[ii]],contacts.cocoon[[ii]],contacts.company[[ii]]))
      }
      
      # determine who becomes newly infected
      infect.who = rbinom(length(contacts.today),
                          1,
                          infect.today * (1 - compliance[contacts.today]))
      if (sum(infect.who) > 0) {
        infect.who = contacts.today[which(infect.who == 1)]
        
        # update their status if they're not already infected
        if (sum(infected[infect.who] == 0) > 0) {
          infect.new = infect.who[which(infected[infect.who] == 0)]
          infected[infect.new] = 1
          date.infected[infect.new] = tt
          date.symptoms[infect.new] =
            ifelse(
              rbinom(length(infect.new), 1, symp.mag) == 1,
              date.infected[infect.new] +
                sample(
                  1:length(prop.inf.symp),
                  length(infect.new),
                  prob = prop.inf.symp,
                  replace = T
                ),
              NA
            )
          symptom.duration[infect.new] = rpois(1,symp.mean)
          edges = rbind(edges, cbind(ii, infect.new))
        }
      }
    } # end infectious loop
    
    # test symptomatics, and perform isolation and quarantine accordingly
    needtesting = which(date.symptoms == tt)
    needtesting = needtesting[which(arrivebase[needtesting] < tt)]
    needtesting = needtesting[which(is.na(isolate.date[needtesting]))]
    needtesting = needtesting[which(is.na(quarantine.date[needtesting]))]
    numTested[tt] = numTested[tt] + length(needtesting)
    if (length(needtesting) > 0) {
      trueneg = needtesting[is.na(date.infected[needtesting])]
      truepos = setdiff(needtesting, trueneg)
      testpos = rep(NA, length(needtesting))
      if (length(trueneg) > 0) {
        testpos[which(needtesting %in% trueneg)] =
          rbinom(length(trueneg), 1, 1 - pcr.spec)
      }
      if (length(truepos) > 0) {
        testpos[which(needtesting %in% truepos)] =
          rbinom(length(truepos), 1, pcr.sens(tt - date.infected[truepos] +
                                                1))
      }
      tested.date[needtesting] = tt
      needtracing = needtesting[which(testpos == 1)]
      isolate.date[needtracing] = tt + testReturn
      for (ii in needtracing) {
        if (tt <= 14){
          contacts.all = unique(contacts.random[[ii]],contacts.cocoon[[ii]])
        } else {
          contacts.all = unique(c(contacts.random[[ii]],contacts.cocoon[[ii]],contacts.company[[ii]]))
        }
        contacts.all = setdiff(contacts.all, which(!is.na(quarantine.date)))
        quarantine.date[sample(contacts.all, min(
          rpois(1, quarantine.contacts),
          length(contacts.all)
        ), replace = F)] =
          tt + 1
      }
    }
    
    # test those who were recently quarantined and release if negative
    if (tt > testDelayQuarantine) {
      needtesting = which(quarantine.date == (tt - testDelayQuarantine))
      numTested[tt - testDelayQuarantine] =
        numTested[tt - testDelayQuarantine] + length(needtesting)
      if (length(needtesting) > 0) {
        trueneg = needtesting[is.na(tt - date.infected[needtesting] + 1)]
        truepos = setdiff(needtesting, trueneg)
        testpos = rep(NA, length(needtesting))
        if (length(trueneg) > 0) {
          testpos[which(needtesting %in% trueneg)] =
            rbinom(length(trueneg), 1, 1 - pcr.spec)
        }
        if (length(truepos) > 0) {
          testpos[which(needtesting %in% truepos)] =
            rbinom(length(truepos),
                   1,
                   pcr.sens(tt - date.infected[truepos] + 1))
        }
        tested.date[needtesting] = tt
        release = which(testpos == 0)
        if (length(release) > 0) {
          quarantine.date[needtesting[release]] = NA
        }
        needtracing = needtesting[which(testpos == 1)]
        for (ii in needtracing) {
          if ( tt <= 14){
            contacts.all = unique(contacts.random[[ii]],contacts.cocoon[[ii]])
          } else {
            contacts.all = unique(c(contacts.random[[ii]],contacts.cocoon[[ii]],contacts.company[[ii]]))
          }
          contacts.all = setdiff(contacts.all, which(!is.na(quarantine.date)))
          quarantine.date[sample(contacts.all, min(
            rpois(1, quarantine.contacts),
            length(contacts.all)
          ), replace = F)] =
            tt + 1
        }
      }
    }
    
    # record numbers for the day
    numIsolated[tt] = sum(!is.na(isolate.date))
    numQuarantine[tt] = sum(!is.na(quarantine.date))
    numInfected[tt] = sum(date.infected == tt, na.rm = T)
    numSymptoms[tt] = sum(date.symptoms == tt, na.rm = T)
    
    # release from isolation
    if (BOOL_clinicalrelease){ # release based on clinical criteria
      release = c(which(max(date.symptoms+isolate.min, # 10 days from symptom onset 
                            date.symptoms+symptom.duration+isolate.nosymp) == tt), # delay from last day of fever
                  which(tested.date[which(is.na(date.symptoms))]+isolate.min == tt)) # 10 days from test date for asymptomatics
    } else { # release based on testing
      needtesting = which(isolate.date <= (tt - isolate.max - testDelayQuarantine))
      numTested[tt - testDelayQuarantine] = numTested[tt - testDelayQuarantine] + length(needtesting)
      if (length(needtesting) > 0) {
        trueneg = needtesting[is.na(tt - date.infected[needtesting] + 1)]
        truepos = setdiff(needtesting, trueneg)
        testpos = rep(NA, length(needtesting))
        if (length(trueneg) > 0) {
          testpos[which(needtesting %in% trueneg)] =
            rbinom(length(trueneg), 1, 1 - pcr.spec)
        }
        if (length(truepos) > 0) {
          testpos[which(needtesting %in% truepos)] =
            rbinom(length(truepos),
                   1,
                   pcr.sens(tt - date.infected[truepos] + 1))
        }
        tested.date[needtesting] = tt
        release = which(testpos == 0)
      } else {
        release = numeric()
      }
    }
    
    if (length(release) > 0) {
      isolate.date[release] = NA
    }
    
    # release from quarantine if time up
    release = which(quarantine.date == (tt - quarantine.max))
    if (length(release) > 0) {
      quarantine.date[release] = NA
    }
    
  } # end time loop
  
  # outputs
  summary_out = c(init.infect,
                  init.pos,
                  sum(numInfected), 
                  sum(numSymptoms),
                  max(numIsolated + numQuarantine),
                  sum(numTested))
  #timeseries_out = data.frame(cbind(numInfected,numSymptoms,numIsolated,numQuarantine,numTested))
  #indivudal_out = data.frame(cbind(date.symptoms,isolate.date,quarantine.date))
  
  #return(list(summary_out,timeseries_out,indivudal_out))
  return(summary_out)
  
} # end function definition