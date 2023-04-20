

#####################################################################
# figure: histograms of outbreak size for testing days scenarios
dev.new()
par(mfrow = c(4,1))

load('./output/results13Aug/testS1.RData')
hist(summary_out[,4], breaks = 10, main = 'Scenario 1', xlab = 'total number of infections',cex.lab = 1.2,col = 'lightblue',ylim = c(1,350))
text(400,250,eval(paste('probability of outbreak =',length(which(summary_out[,4]>10))/500*100,'%' )),cex = 1.4)
text(400,170,eval(paste('median outbreak size =',median(summary_out[which(summary_out[,4]>10),4]))),cex = 1.4)

load('./output/results13Aug/testS2.RData')
hist(summary_out[,4], breaks = 10, main = 'Scenario 2', xlab = 'total number of infections',cex.lab = 1.2,col = 'lightblue',ylim = c(1,350))
text(400,250,eval(paste('probability of outbreak =',length(which(summary_out[,4]>10))/500*100,'%' )),cex = 1.4)
text(400,170,eval(paste('median outbreak size =',median(summary_out[which(summary_out[,4]>10),4]))),cex = 1.4)

load('./output/results13Aug/testS3.RData')
hist(summary_out[,4], breaks = 10, main = 'Scenario 3', xlab = 'total number of infections',cex.lab = 1.2,col = 'lightblue',ylim = c(1,350))
text(400,250,eval(paste('probability of outbreak =',length(which(summary_out[,4]>10))/500*100,'%' )),cex = 1.4)
text(400,170,eval(paste('median outbreak size =',median(summary_out[which(summary_out[,4]>10),4]))),cex = 1.4)

load('./output/results13Aug/testS4.RData')
hist(summary_out[,4], breaks = 10, main = 'Scenario 4', xlab = 'total number of infections',cex.lab = 1.2,col = 'lightblue',ylim = c(1,350))
text(400,250,eval(paste('probability of outbreak =',length(which(summary_out[,4]>10))/500*100,'%' )),cex = 1.4)
text(400,170,eval(paste('median outbreak size =',median(summary_out[which(summary_out[,4]>10),4]))),cex = 1.4)

dev.off()

#####################################################################
# figure: boxplot of outbreak size and tests given for test return time scenarios
# dev.new()
# 
# load('./output/outbreak_testreturnS1.RData')
# numinf = summary_out[,3]
# numtests = summary_out[,6]
# 
# load('./output/outbreak_testreturnS2.RData')
# numinf = c(numinf,summary_out[,3])
# numtests = c(numtests,summary_out[,6])
# 
# load('./output/outbreak_testreturnS3.RData')
# numinf = c(numinf,summary_out[,3])
# numtests = c(numtests,summary_out[,6])
# 
# par(mfrow = c(2,1))
# scenario = rep(1:3,each = 250)
# id = which(numinf>0)
# boxplot(numinf[id]~scenario[id], ylab = 'number of infections',xlab = 'scenario',col = 'lightblue')
# boxplot(numtests[id]~scenario[id], ylab = 'number of tests')
# 
# dev.off()

#####################################################################
# figure: boxplot of outbreak size and tests given for cocoon size scenarios

# dev.new()
# 
# load('./output/outbreak_cocoonS1.RData')
# numinf = summary_out[,3]
# numtests = summary_out[,6]
# 
# load('./output/outbreak_cocoonS2.RData')
# numinf = c(numinf,summary_out[,3])
# numtests = c(numtests,summary_out[,6])
# 
# load('./output/outbreak_cocoonS3.RData')
# numinf = c(numinf,summary_out[,3])
# numtests = c(numtests,summary_out[,6])
# 
# par(mfrow = c(2,1))
# scenario = rep(1:3,each = 250)
# id = which(numinf>0)
# boxplot(numinf[id]~scenario[id], ylab = 'number of infections',xlab = 'scenario',col = 'lightblue')
# boxplot(numtests[id]~scenario[id], ylab = 'number of tests')
# 
# dev.off()

#####################################################################
# figure: histograms of outbreak size for mask compliance scenarios

dev.new()
par(mfrow = c(4,1))

load('./output/results13Aug/maskS1.RData')
numinf1 = summary_out[,4]
numtests1 = summary_out[,7]
hist(summary_out[,4], breaks = 10, main = 'Scenario 1', xlab = 'total number of infections',cex.lab = 1.2,col = 'lightblue',ylim = c(1,350))
text(400,250,eval(paste('probability of outbreak =',length(which(summary_out[,4]>10))/500*100,'%' )),cex = 1.4)

load('./output/results13Aug/maskS2.RData')
numinf1 = c(numinf1,summary_out[,4])
numtests1 = c(numtests1,summary_out[,7])
hist(summary_out[,4], breaks = 10, main = 'Scenario 2', xlab = 'total number of infections',cex.lab = 1.2,col = 'lightblue',ylim = c(1,350))
text(400,250,eval(paste('probability of outbreak =',length(which(summary_out[,4]>10))/500*100,'%' )),cex = 1.4)

load('./output/results13Aug/maskS3.RData')
numinf1 = c(numinf1,summary_out[,4])
numtests1 = c(numtests1,summary_out[,7])
hist(summary_out[,4], breaks = 10, main = 'Scenario 3', xlab = 'total number of infections',cex.lab = 1.2,col = 'lightblue',ylim = c(1,350))
text(400,250,eval(paste('probability of outbreak =',length(which(summary_out[,4]>10))/500*100,'%' )),cex = 1.4)

load('./output/results13Aug/maskS4.RData')
numinf1 = c(numinf1,summary_out[,4])
numtests1 = c(numtests1,summary_out[,7])
hist(summary_out[,4], breaks = 10, main = 'Scenario 4', xlab = 'total number of infections',cex.lab = 1.2,col = 'lightblue',ylim = c(1,350))
text(400,250,eval(paste('probability of outbreak =',length(which(summary_out[,4]>10))/500*100,'%' )),cex = 1.4)

dev.off()

#####################################################################
# figure: histograms of outbreak size for contact tracing scenarios

dev.new()
par(mfrow = c(3,1))

load('./output/results13Aug/cont1.RData')
numinf1 = summary_out[,4]
numtests1 = summary_out[,7]
hist(summary_out[,4], breaks = 10, main = 'Scenario 1', xlab = 'total number of infections',cex.lab = 1.2,col = 'lightblue',ylim = c(1,350))
text(400,250,eval(paste('probability of outbreak =',length(which(summary_out[,4]>10))/500*100,'%' )),cex = 1.4)
text(400,200,eval(paste('median outbreak size =',median(summary_out[which(summary_out[,4]>10),4]))),cex = 1.4)

load('./output/results13Aug/cont2.RData')
numinf1 = c(numinf1,summary_out[,4])
numtests1 = c(numtests1,summary_out[,7])
hist(summary_out[,4], breaks = 10, main = 'Scenario 2', xlab = 'total number of infections',cex.lab = 1.2,col = 'lightblue',ylim = c(1,350))
text(400,250,eval(paste('probability of outbreak =',length(which(summary_out[,4]>10))/500*100,'%' )),cex = 1.4)
text(400,200,eval(paste('median outbreak size =',median(summary_out[which(summary_out[,4]>10),4]))),cex = 1.4)

load('./output/results13Aug/cont3.RData')
numinf1 = c(numinf1,summary_out[,4])
numtests1 = c(numtests1,summary_out[,7])
hist(summary_out[,4], breaks = 10, main = 'Scenario 3', xlab = 'total number of infections',cex.lab = 1.2,col = 'lightblue',ylim = c(1,350))
text(400,250,eval(paste('probability of outbreak =',length(which(summary_out[,4]>10))/500*100,'%' )),cex = 1.4)
text(400,200,eval(paste('median outbreak size =',median(summary_out[which(summary_out[,4]>10),4]))),cex = 1.4)

dev.off()

#####################################################################

