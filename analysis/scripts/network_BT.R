## network_BCT.R
## 
## Builds a contact network for a single cohort of basic training recruits
## 
## INPUTS: numrecruits, cocoon size, company size, numStaff, drillsergPerCocoon
##
## OUTPUTS: contact lists for each recruit, saved as a .RData file to be loaded by simOutbreak.R
##
## FEATURES:
##
## mixing among cocoons
## mixing among companies
## random mixing
## drill sergeant contact with recruits
## staff contact with recruits
##=============================================#
## Setup network variables------------
##=============================================#
## number of recruits in basic training
training_site = 'full'
numrecruits_list = list('full' = 1200, 'FB' = 640, 'FLW' = 500)

numrecruits = numrecruits_list[[training_site]]

## number of staff interacting with recruits
numStaff = 0.05*numrecruits ## assume numStaff = 5% of numrecruits
RecruitContactsPerStaff = 20

## cocoon size 30-65
cocoon.size = 60
numcocoons = ceiling(numrecruits/cocoon.size)

## company size up to 240, shared common areas
company.size = 240
numcompany = ceiling(numcocoons/round(company.size/cocoon.size))

## number of drill sergeants
drillsergPerCocoon = 2
numdrillserg = numcocoons*drillsergPerCocoon

##=============================================#
## assign recruits ---------
##=============================================#
## assign recruits to cocoons
recruit.cocoon = sample(rep(1:numcocoons,each = cocoon.size))[1:numrecruits]

## companies are made up of cocoons
cocoon.company = sample(rep(1:numcompany,each = ceiling(numcocoons/numcompany)))[1:numcocoons] 
recruit.company = numeric()
for (i in 1:numcompany){
  j = which(cocoon.company == i)
  k = which(recruit.cocoon %in% j)
  recruit.company[k] = i
}

## assign drill sergeants to cocoons
drillserg.cocoon = sample(rep(1:numcocoons,each = drillsergPerCocoon))[1:numdrillserg]

##=============================================#
## build contact network---------------
##=============================================#
## simulate contacts among cocoons
network.cocoon = matrix(0,numrecruits,numrecruits)
for(ii in 1:numcocoons){
  which.in.cocoon = expand.grid(
    which(recruit.cocoon==ii),which(recruit.cocoon==ii))
  network.cocoon[as.matrix(which.in.cocoon)] =
    rbinom(nrow(which.in.cocoon),1,0.2)
}
contacts.cocoon = apply(network.cocoon,1,function(x)which(x>0)); rm(network.cocoon)

## add drill sergeant contacts
for (ii in 1:numrecruits){
  contacts.cocoon[[ii]] = c(contacts.cocoon[[ii]],
                      numrecruits + which(drillserg.cocoon == recruit.cocoon[[ii]]))
}

contacts.drillserg.recruit = list()
for(ii in 1:numdrillserg){
  contacts.drillserg.recruit[[ii]] = which(recruit.cocoon == drillserg.cocoon[[ii]])
}

## simulate contacts among companies
network.company = matrix(0,numrecruits,numrecruits)
for(ii in 1:numcompany){
  which.in.company = expand.grid(
    which(recruit.company==ii),which(recruit.company==ii))
  network.company[as.matrix(which.in.company)] =
    rbinom(nrow(which.in.company),1,0.1)
}
contacts.company = apply(network.company,1,function(x)which(x>0)); rm(network.company)

## simulate random contacts
## 0% of random contacts
## network.random = matrix(
##   rbinom(numrecruits^2,1,0.1),
##   numrecruits,numrecruits)
## network.random = network.random + t(network.random)
## contacts.random = apply(network.random,1,function(x)which(x>10000)); rm(network.random)
contacts.random = lapply(1:numrecruits, function(x){c()})



## add staff contacts
staffContact = rbeta(numStaff,1,1)
staffRecruitContacts = rmultinom(1,numStaff*RecruitContactsPerStaff,staffContact)

mat.staffrecruit = matrix(0,numStaff,numrecruits)
for(ii in 1:numStaff){
  mat.staffrecruit[ii,sample(numrecruits,staffRecruitContacts[ii],replace=F)] = 1
}
contacts.staff.recruit = list()
for(ii in 1:numStaff){
  contacts.staff.recruit[[ii]] = which(mat.staffrecruit[ii,]>0)
}

## No random contacts right now
## for(ii in 1:numrecruits){
##   contacts.random[[ii]] = c(
##     contacts.random[[ii]],
##     numrecruits + numdrillserg +
##       which(mat.staffrecruit[,ii]>0))
## }

## tally the total contacts of each recruit
total.contacts =
  unlist(lapply(contacts.random,length)) +
  unlist(lapply(contacts.cocoon,length)) +
  unlist(lapply(contacts.company,length)) 


## save results
if(training_site == 'full'){
    save(list=ls(),file='../output/network.RData')
}else{
    save(list=ls(),file=sprintf('../output/network_%s.RData', training_site))
}
    

