# script to plot outputs of a single simOutbreak replicate

# load libraries
library(igraph)

pdf('output/outbreak.pdf',width=6.5,height=8.5)

layout(matrix(1:6,3,2,byrow=T))

# incident infections
eval(parse(text=paste("
barplot(
  numInfected,las=1,space=0,
  xlab='Days since first recruit on base',
  ylab='Daily number of new infections',
  main='Total infections: ",sum(numInfected),"')",sep='')))
#axis(1)

# incident cases
eval(parse(text=paste("
barplot(
  numSymptoms,las=1,space=0,
  xlab='Days since first recruit on base',
  ylab='Daily number of new symptomatic infections',
  main='Total symptomatic infections: ",sum(numSymptoms),"')",sep='')))
#axis(1)

# daily number in isolation
eval(parse(text=paste("
  barplot(
  numIsolated+numQuarantine,las=1,space=0,
  xlab='Days since first recruit on base',
  ylab='Daily number of beds filled',
  main='Peak beds: ",max(numIsolated+numQuarantine),"')",sep='')))
#axis(1)

# daily number in isolation
eval(parse(text=paste("
barplot(
  numIsolated,las=1,space=0,
  xlab='Days since first recruit on base',
  ylab='Daily number in isolation',
  main='Peak in isolation: ",max(numIsolated),"')",sep='')))
#axis(1)

# daily number in quarantine
eval(parse(text=paste("
barplot(
  numQuarantine,las=1,space=0,
  xlab='Days since first recruit on base',
  ylab='Daily number in quarantine',
  main='Peak in quarantine: ",max(numQuarantine),"')",sep='')))
#axis(1)

# daily number of tests
eval(parse(text=paste("
barplot(
  numTested,las=1,space=0,
  xlab='Days since first recruit on base',
  ylab='Daily number of tests performed',
  main='Total tests needed: ",sum(numTested),"')",sep='')))
#axis(1)

# plot network of infections on campus
# campusCases = edges[which(edges[,1]!=0),]
# symptomatic = !is.na(date.symptoms[campusCases])
# isolated = !is.na(isolate.date[campusCases])
# quarantined = !is.na(quarantine.date[campusCases])
# contactGraph = graph_from_data_frame(campusCases)
# plot(contactGraph, edge.arrow.size = 0, edge.color = 'black',
#      vertex.size = 2, # ifelse(symptomatic,2,1),
#      vertex.label=NA, edge.width = 0.3, vertex.frame.color = NA,
#      vertex.color = 'DarkBlue', # ifelse(symptomatic,rgb(1,0,0,0.5),rgb(0,0,0,0.5)),
#      layout = layout_as_tree)

dev.off()
