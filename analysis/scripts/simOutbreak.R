## simOutbreak.R
## 
## Simulates SARS-CoV-2 transmission and interventions in a basic training setting, tracks a single 
## cohort of recruits through a 12-week basic training session
## 
## INPUTS: contact network built by network_BT.R, parameter values for infection processes, 
## parameter values for control processes(testing, isolation/quarantine/contact tracing, masks), 
## number of replicates to simulate
##
## OUTPUTS: output contains summary and time series data on number of infections, number of symptomatic 
## infections, number in isolation and quarantine, number of tests given, who infected who, and individual 
## level dates of infection/symptoms/test
##
## FEATURES:
##
## staggered arrival dates
## quarantine with coccoon from day 4 to day 14
## random contacts allowed throughout 
## test symptomatics
## test return delay
## isolate positive cases
## contact tracing and quarantine test-positive cases
## test quarantine after 3 days and release if test-negative
## individual compliance
## recent pre-arrival infection
## clinical and test-based criteria to leave isolation
## https://www.who.int/news-room/commentaries/detail/criteria-for-releasing-covid-19-patients-from-isolation 
## https://www.journalofinfection.com/article/S0163-4453(20)30119-5/fulltext#seccesectitle0014
## https://www.cdc.gov/coronavirus/2019-ncov/hcp/duration-isolation.html
## Contact with drill sergeants/staff

## TO DO LIST:
##
## [x] Test trainers
## [x] Implement TOST and incubation period
## [x] Implement a more thorough sensitivity

simOutbreak = function(){
    ## specify secondary infection distribution    
    R0 = R0.prob * genint
        
    ## specify delay between infection and symptom onset
    ## Replacing this line with the incub function
    ##prop.inf.symp = symp.mag * dpois(1:28, symp.tim) / sum(dpois(1:28, symp.tim))
    ## INCUBATION PERIOD
    incub = pgamma(1:21,shape=k_inc,scale=0.948) - pgamma(0:20,shape=k_inc,scale=0.948)
    prop.inf.symp = incub / sum(incub)
    
    ## previously infected and immune
    initial.immune = initial.immune + initial.infected
    immune = rep(0, numrecruits+numdrillserg+numStaff)
    immune[sample(numrecruits+numdrillserg+numStaff, rbinom(1, numrecruits+numdrillserg+numStaff, initial.immune), replace = F)] = 1
    
    ## infected status
    ## which(immune > 0, replace = F)
    infected = rep(0, numrecruits+numdrillserg+numStaff)
    if(length(which(immune > 0)) > 0){
        infected[sample(which(immune > 0), rbinom(1, sum(immune), initial.infected/initial.immune), replace = F)] = 1
    }
    init.infect = sum(infected)
    
    ## simulate day of infection of initial infecteds
    ## We don't need the foor loop
    ## arrivebase[something] - sample(1:21, init.infect)
    ## R0_init = R0.prob * dpois(1:21, R0.tim/2) / sum(dpois(1:21, R0.tim/2)) # specify secondary infection distribution for initial infecteds to force some recent infections
    
    date.infected = rep(NA, numrecruits+numdrillserg+numStaff)
    secondary.infected = rep(NA, numrecruits+numdrillserg+numStaff)
    incub.period.inf = rep(NA, numrecruits+numdrillserg+numStaff)
    
    ## currently infected
    date.infected[which(infected > 0)] = arrivebase[which(infected > 0)] - sample(1:39, init.infect, replace = T)

    ## already immune
    date.infected[immune > 0 & infected == 0] = -1e3       
    infected[immune > 0] = 1

    
    ## simulate date of symptom onset and duration
    date.symptoms = rep(NA, numrecruits+numdrillserg+numStaff)
    symptom.duration = rep(NA, numrecruits+numdrillserg+numStaff)
    for(ii in which(infected > 0)){
        if(rbinom(1, 1, symp.mag) == 1){
            date.symptoms[ii] =
                date.infected[ii] +
                sample(1:length(prop.inf.symp), 1, prob = prop.inf.symp)
            symptom.duration[ii] = rpois(1,symp.mean)
            incub.period.inf[ii] = date.symptoms[ii]
        }else{
            incub.period.inf[ii] = date.infected[ii] +
                sample(1:length(prop.inf.symp),1,replace = T,prob = prop.inf.symp)
        }
    }

    
    ## individual status storage
    tested.date = rep(NA, numrecruits+numdrillserg+numStaff)
    isolate.date = rep(NA, numrecruits+numdrillserg+numStaff)
    quarantine.date = rep(NA, numrecruits+numdrillserg+numStaff)    
    testneg.date = rep(NA, numrecruits+numdrillserg+numStaff)
    init.pos = 0
    retest.pos = 0
        
    ## aggregate numbers over time
    numIsolated = rep(0, length(time))
    numQuarantine = rep(0, length(time))
    numTested = rep(0, length(time))
    numInfected = rep(0, length(time))
    numSymptoms = rep(0, length(time))
    numImported = rep(0, length(time))
    numTestPositive = rep(0, length(time))
    
    ## keep track of who infected whom
    edges = cbind(rep(0, sum(infected)), which(infected > 0))
    
    ## loop over each day of basic training
    for(tt in time){
        compliance.mult = 1
        if(BOOL_compliance_time == 1){
            compliance.mult = ((compliance.final.avg - compliance.avg)/(length(time) - 1) * (tt-1) + compliance.avg) / compliance.avg
        }
        ## test on arrival
        if(tt %in% unique(arrivebase) & BOOL_testOnArrival == 1){
            needtesting = which(arrivebase == tt)
            trueneg = needtesting[is.na(date.infected[needtesting])]
            truepos = setdiff(needtesting, trueneg)
            testpos = rep(NA, length(needtesting))
            if(length(trueneg) > 0){ # false pos
                testpos[which(needtesting %in% trueneg)] =
                    rbinom(length(trueneg), 1, 1 - pcr.spec.commercial)
            }
            ## Upon arrival, always test with PCR
            if(length(truepos) > 0){ # true pos
                testpos[which(needtesting %in% truepos)] =
                    rbinom(length(truepos), 1,
                           pcr.sens.commercial(tt - date.infected[truepos] + 1))
            }
            tested.date[needtesting] = tt
            needtracing = needtesting[which(testpos == 1)] # all pos
            isolate.date[needtracing] = tt + testReturn
            testneg.date[setdiff(needtesting,needtracing)] = tt
            numTested[tt] = numTested[tt] + length(needtesting)
            init.pos = init.pos + length(needtracing)
        }

        ## Only if want to test everyday, for bookkeeping
        if(BOOL_testDaily == TRUE){
            needtesting = 1:numrecruits
            truepos = needtesting[!is.na(date.infected[needtesting])]
            trueneg = numrecruits - length(truepos)
            testpos = 0
            if(trueneg > 0){
                testpos = sum(rbinom(trueneg, 1, 1 - pcr.spec.commercial))
            }
            if(length(truepos) > 0){
                testpos = testpos + sum(
                    rbinom(length(truepos), 1,
                           pcr.sens.commercial(tt - date.infected[truepos] + 1)
                           ))
            }
            numTestPositive[tt] = testpos
        }
        
        ## test on specified days
        if(tt %in% testdates){
            needtesting = 1:numrecruits
            trueneg = needtesting[is.na(date.infected[needtesting])]
            truepos = setdiff(needtesting, trueneg)
            testpos = rep(NA, length(needtesting))

            ind_testdates = which(testdates == tt)
            tmp_test_type = 'pcr'
            if(ind_testdates <= length(testdates_type)){
                tmp_test_type = testdates_type[ind_testdates]
            }
            ## can be antigen or PCR, if not antigen, assume pcr
            if(tmp_test_type == "antigen"){
                if(length(trueneg) > 0){
                    testpos[which(needtesting %in% trueneg)] =
                        rbinom(length(trueneg), 1, 1 - pcr.spec.screen)
                }
                if(length(truepos) > 0){
                    testpos[which(needtesting %in% truepos)] =
                        rbinom(length(truepos), 1, pcr.sens.screen(tt - date.infected[truepos] +
                                                                   1))
                }
            }else{
                if(length(trueneg) > 0){
                    testpos[which(needtesting %in% trueneg)] =
                        rbinom(length(trueneg), 1, 1 - pcr.spec.commercial)
                }
                if(length(truepos) > 0){
                    testpos[which(needtesting %in% truepos)] =
                        rbinom(length(truepos), 1, pcr.sens.commercial(tt - date.infected[truepos] +
                                                                   1))
                }
            }
            tested.date[needtesting] = tt
            needtracing = needtesting[which(testpos == 1)]
            isolate.date[needtracing] = tt + testReturn
            testneg.date[setdiff(needtesting,needtracing)] = tt
            numTested[tt] = numTested[tt] + length(needtesting)
            retest.pos = retest.pos + length(needtracing)
            
            for(ii in needtracing){
                if(ii <= numrecruits){ # recruits                    
                    if(tt <= 14){
                        contacts.all = unique(c(contacts.random[[ii]],contacts.cocoon[[ii]]))
                    } else {
                        contacts.all = unique(c(contacts.random[[ii]],contacts.company[[ii]]))
                    } 
                } else if(ii > numrecruits + numdrillserg){ # staff
                    contacts.all = contacts.staff.recruit[[ii - numrecruits - numdrillserg]]
                }else { # drill sergeants
                    contacts.all = contacts.drillserg.recruit[[ii - numrecruits]]
                }

                contacts.all = contacts.all[contacts.all != ii]
                contacts.all = setdiff(contacts.all, which(!is.na(quarantine.date)))

                ## CONTACT TRACING!! quarantine.contacts and set date to isolate if positive
                contacts.needtest = sample(contacts.all,
                       min(
                           rpois(1, quarantine.contacts),
                           length(contacts.all)
                       ), replace = F)
                
                trueneg = contacts.needtest[is.na(date.infected[contacts.needtest])]
                truepos = setdiff(contacts.needtest, trueneg)
                testpos = rep(NA, length(contacts.needtest))

                ## can be antigen or PCR, if not antigen, assume pcr
                if(tmp_test_type == "antigen"){
                    if(length(trueneg) > 0){
                        testpos[which(contacts.needtest %in% trueneg)] =
                            rbinom(length(trueneg), 1, 1 - pcr.spec.screen)
                    }
                    if(length(truepos) > 0){
                        testpos[which(contacts.needtest %in% truepos)] =
                            rbinom(length(truepos), 1, pcr.sens.screen(tt - date.infected[truepos] +
                                                                       1))
                    }
                }else{
                    if(length(trueneg) > 0){
                        testpos[which(contacts.needtest %in% trueneg)] =
                        rbinom(length(trueneg), 1, 1 - pcr.spec.commercial)
                    }
                    if(length(truepos) > 0){
                        testpos[which(contacts.needtest %in% truepos)] =
                            rbinom(length(truepos), 1, pcr.sens.commercial(tt - date.infected[truepos] +
                                                                           1))
                    }
                }
                tested.date[contacts.needtest] = tt                
                quarantine.date[contacts.needtest[testpos == 1]] = tt + testReturn # Can also be isolate.date
                testneg.date[contacts.needtest[testpos != 1]] = tt
                numTested[tt] = numTested[tt] + length(contacts.needtest)
                retest.pos = retest.pos + length(which(testpos == 1))            
            }
        }


        ## test staff on specified frequencies
        if((testStaffFreq > 0) && ((tt - 1) %% testStaffFreq == 0)){
            needtesting = (numrecruits + 1):length(infected) # Is this for staff or drill sergants too?
            trueneg = needtesting[is.na(date.infected[needtesting])]
            truepos = setdiff(needtesting, trueneg)
            testpos = rep(NA, length(needtesting))
           
            ## can be antigen or PCR, if not antigen, assume pcr
            if(testStaffType == "antigen"){
                if(length(trueneg) > 0){
                    testpos[which(needtesting %in% trueneg)] =
                        rbinom(length(trueneg), 1, 1 - pcr.spec.screen)
                }
                if(length(truepos) > 0){
                    testpos[which(needtesting %in% truepos)] =
                        rbinom(length(truepos), 1, pcr.sens.screen(tt - date.infected[truepos] +
                                                                   1))
                }
            }else{
                if(length(trueneg) > 0){
                    testpos[which(needtesting %in% trueneg)] =
                        rbinom(length(trueneg), 1, 1 - pcr.spec.commercial)
                }
                if(length(truepos) > 0){
                    testpos[which(needtesting %in% truepos)] =
                        rbinom(length(truepos), 1, pcr.sens.commercial(tt - date.infected[truepos] +
                                                                       1))
                }
            }
            tested.date[needtesting] = tt
            
            needtracing = needtesting[which(testpos == 1)]
            isolate.date[needtracing] = tt + testReturn
            testneg.date[setdiff(needtesting,needtracing)] = tt
            numTested[tt] = numTested[tt] + length(needtesting)
            retest.pos = retest.pos + length(needtracing)                        
        }
        
        ## importation from staff        
        infected.offcampus = which(rbinom(numrecruits+numdrillserg+numStaff, 1,(1-compliance*compliance.mult)*importation) == 1)
        numImported[tt] = length(infected.offcampus)
        if(length(infected.offcampus) > 0){
            if(sum(is.na(date.infected[infected.offcampus])) > 0){
                edges = rbind(
                    edges,
                    cbind(0,infected.offcampus[which(is.na(date.infected[infected.offcampus]))]))
            }
            infected[infected.offcampus] = 1
            date.infected[infected.offcampus] = ifelse(
                is.na(date.infected[infected.offcampus]),
                tt,
                date.infected[infected.offcampus])
            date.symptoms[infected.offcampus] =
                ifelse(
                    is.na(date.symptoms[infected.offcampus]),
                ifelse(
                    rbinom(length(infected.offcampus), 1, symp.mag) == 1,
                    date.infected[infected.offcampus] +
                    sample(
                        1:length(prop.inf.symp),
                        length(infected.offcampus),
                        prob = prop.inf.symp,
                        replace = T),
                    NA),
                date.symptoms[infected.offcampus])
            incub.period.inf[infected.offcampus] = date.symptoms[infected.offcampus]

            
            incub.period.inf[infected.offcampus[is.na(date.symptoms[infected.offcampus])]] = date.infected[infected.offcampus[is.na(date.symptoms[infected.offcampus])]] +
                sample(1:length(prop.inf.symp),
                       length(infected.offcampus[is.na(date.symptoms[infected.offcampus])]),
                       replace = T, prob = prop.inf.symp)
        }
        
        ## loop through those who are capable of infecting others
        max_infectious_period = 21
        infectious = which(infected > 0 &
                           date.infected < tt &
                           date.infected > (tt - max_infectious_period) &
                           arrivebase <= tt)


        for(ii in infectious){
            infect.today =
                R0[(tt - date.infected[ii]) + 1] *
                ifelse(is.na(date.symptoms[ii]), asymp.adjust, 1) * (1 - compliance[ii]*compliance.mult)
            
            ## We need to prevent immune people from becoming infected again
            ## look up this person's contacts today
            if(!is.na(isolate.date[ii]) | !is.na(quarantine.date[ii])){
                if(ii <= numrecruits){
                    which.in.isolation = unique(c(which(!is.na(isolate.date)), which(!is.na(quarantine.date))))
                    contacts.all = which.in.isolation
                }else{
                    contacts.all = numeric()
                }
            } else {
                if(ii <= numrecruits){ # recruits
                    if(tt <= 14){
                        contacts.all = unique(c(contacts.random[[ii]],contacts.cocoon[[ii]]))
                    } else {
                        contacts.all = unique(c(contacts.random[[ii]],contacts.company[[ii]]))
                    }
                } else if(ii > numrecruits + numdrillserg){ # staff
                    contacts.all = unique(c(contacts.staff.recruit[[ii - numrecruits - numdrillserg]]))
                }else { # drill sergeants
                    contacts.all = unique(c(contacts.drillserg.recruit[[ii - numrecruits]]))
                }
            }
            contacts.all = contacts.all[contacts.all != ii]

            ## determine who becomes newly infected
            infect.who = rbinom(length(contacts.all),
                                1,
                                infect.today *(1 - compliance[contacts.all]*compliance.mult))
            if(sum(is.na(infect.who)) > 0){
                print(infect.today)
                print(date.infected[ii])
                print(incub.period.inf[ii])
                print((tt - date.infected[ii]) + 1)
                stop()
            }
            if(sum(infect.who) > 0){
                infect.who = contacts.all[which(infect.who == 1)]
                
                ## update their status if they're not already infected
                if(sum(infected[infect.who] == 0) > 0){
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
                    incub.period.inf[infect.new] = date.symptoms[infect.new]                    
                    incub.period.inf[infect.new[is.na(date.symptoms[infect.new])]] = date.infected[infect.new[is.na(date.symptoms[infect.new])]] +
                        sample(1:length(prop.inf.symp),
                               length(infect.new[is.na(date.symptoms[infect.new])]),
                               replace = T,
                               prob = prop.inf.symp)
                    symptom.duration[infect.new] = rpois(1,symp.mean)
                    edges = rbind(edges, cbind(ii, infect.new))
                }
            }
        } # end infectious loop
        
        ## test symptomatics, and perform isolation and quarantine accordingly
        needtesting = which(date.symptoms == tt)
        needtesting = needtesting[which(arrivebase[needtesting] < tt)]
        needtesting = needtesting[which(is.na(isolate.date[needtesting]))]
        needtesting = needtesting[which(is.na(quarantine.date[needtesting]))]
        needtesting = needtesting[which(tested.date[needtesting] != tt)]
        numTested[tt] = numTested[tt] + length(needtesting)
        if(length(needtesting) > 0){
            trueneg = needtesting[is.na(date.infected[needtesting])]
            truepos = setdiff(needtesting, trueneg)
            testpos = rep(NA, length(needtesting))
            if(length(trueneg) > 0){
                testpos[which(needtesting %in% trueneg)] =
                    rbinom(length(trueneg), 1, 1 - pcr.spec.commercial)
            }
            if(length(truepos) > 0){
                testpos[which(needtesting %in% truepos)] =
                    rbinom(length(truepos), 1, pcr.sens.commercial(tt - date.infected[truepos] +
                                                                   1))
            }
            tested.date[needtesting] = tt
            needtracing = needtesting[which(testpos == 1)]
            isolate.date[needtracing] = tt # isolate on test date
            testneg.date[setdiff(needtesting,needtracing)] = tt
            for(ii in needtracing){
                if(ii <= numrecruits){ # recruits
                    if(tt <= 14){
                        contacts.all = unique(c(contacts.random[[ii]],contacts.cocoon[[ii]]))
                    } else {
                        contacts.all = unique(c(contacts.random[[ii]],contacts.company[[ii]]))
                    } 
                } else if(ii > numrecruits + numdrillserg){ # staff
                    contacts.all = contacts.staff.recruit[[ii - numrecruits - numdrillserg]]
                }else { # drill sergeants
                    contacts.all = contacts.drillserg.recruit[[ii - numrecruits]]
                }
                contacts.all = contacts.all[contacts.all != ii]
                contacts.all = setdiff(contacts.all, which(!is.na(quarantine.date)))
                quarantine.date[sample(contacts.all, min(
                                                         rpois(1, quarantine.contacts),
                                                         length(contacts.all)
                                                     ), replace = F)] =
                    tt + 1
            }
        }
        
        ## test those who were recently quarantined and release if negative        
        if(tt > testDelayQuarantine){
            needtesting = which(quarantine.date ==(tt - testDelayQuarantine))
            needtesting = needtesting[which(tested.date[needtesting] != tt)]
            numTested[tt] = numTested[tt] + length(needtesting)
            if(length(needtesting) > 0){
                trueneg = needtesting[is.na(date.infected[needtesting] + 1)]
                truepos = setdiff(needtesting, trueneg)
                testpos = rep(NA, length(needtesting))
                if(length(trueneg) > 0){
                    testpos[which(needtesting %in% trueneg)] =
                        rbinom(length(trueneg), 1, 1 - pcr.spec.commercial)
                }
                if(length(truepos) > 0){
                    testpos[which(needtesting %in% truepos)] =
                        rbinom(length(truepos),
                               1,
                               pcr.sens.commercial(tt - date.infected[truepos] + 1))
                }
                tested.date[needtesting] = tt
                release = which(testpos == 0)
                if(length(release) > 0){
                    quarantine.date[needtesting[release]] = NA
                }
                needtracing = needtesting[which(testpos == 1)]
                for(ii in needtracing){
                    if(ii <= numrecruits){ # recruits
                        if( tt <= 14){
                            contacts.all = unique(c(contacts.random[[ii]],contacts.cocoon[[ii]]))
                        } else {
                            contacts.all = unique(c(contacts.random[[ii]],contacts.company[[ii]]))
                        } 
                    } else if(ii > numrecruits + numdrillserg){ # staff
                        contacts.all = contacts.staff.recruit[[ii - numrecruits - numdrillserg]]
                    }else { # drill sergeants
                        contacts.all = contacts.drillserg.recruit[[ii - numrecruits]]
                    }
                    contacts.all = contacts.all[contacts.all != ii]
                    contacts.all = setdiff(contacts.all, which(!is.na(quarantine.date)))
                    quarantine.date[sample(contacts.all, min(
                                                             rpois(1, quarantine.contacts),
                                                             length(contacts.all)
                                                         ), replace = F)] =
                        tt + 1
                }
            }
        }
        
        ## record numbers for the day
        numIsolated[tt] = sum(!is.na(isolate.date))
        numQuarantine[tt] = sum(!is.na(quarantine.date))
        numInfected[tt] = sum(date.infected == tt, na.rm = T)
        numSymptoms[tt] = sum(date.symptoms == tt, na.rm = T)
        
        ## release from isolation
        if(BOOL_clinicalrelease){ # release based on clinical criteria 
            release = c(which(max(date.symptoms+isolate.length, # 10 days from symptom onset 
                                  date.symptoms+symptom.duration+isolate.nosymp) == tt), # delay from last day of fever
                        which(tested.date[which(is.na(date.symptoms))]+isolate.length == tt), # 10 days from test date for asymptomatics
                        which(testneg.date == tt - testReturn)) # negative tests returned 
        } else { # release based on testing
            needtesting = which(isolate.date <=(tt - isolate.length))
            numTested[tt] = numTested[tt] + length(needtesting)
            if(length(needtesting) > 0){
                trueneg = needtesting[is.na(tt - date.infected[needtesting] + 1)]
                truepos = setdiff(needtesting, trueneg)
                testpos = rep(NA, length(needtesting))
                if(length(trueneg) > 0){
                    testpos[which(needtesting %in% trueneg)] =
                        rbinom(length(trueneg), 1, 1 - pcr.spec.commercial)
                }
                if(length(truepos) > 0){
                    testpos[which(needtesting %in% truepos)] =
                        rbinom(length(truepos),
                               1,
                               pcr.sens.commercial(tt - date.infected[truepos] + 1))
                }
                tested.date[needtesting] = tt
                release = which(testpos == 0)
            } else {
                release = numeric()
            }
        }
        
        if(length(release) > 0){
            isolate.date[release] = NA
        }
        
        ## release from quarantine if time up
        release = which(quarantine.date ==(tt - quarantine.max))
        if(length(release) > 0){
            quarantine.date[release] = NA
        }
        
    } # end time loop
    
    ## outputs
    summary_out = c(init.infect,
                    init.pos,
                    retest.pos,
                    sum(numInfected), 
                    sum(numSymptoms),
                    max(numIsolated + numQuarantine),
                    sum(numTested))
    timeseries_out = data.frame(cbind(numInfected,numSymptoms,numIsolated,numQuarantine,numTested, numImported, numTestPositive))
    indivudal_out = data.frame(cbind(date.symptoms,isolate.date,quarantine.date))
    return(list(summary_out,timeseries_out,indivudal_out,edges))
    
} # end function definition

