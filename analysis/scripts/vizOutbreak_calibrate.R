par(mfrow = c(2,3))

######################################################################
# Fort Benning model
load('FB_fit.RData')

hist(summary_out[,2], breaks = 10, main = 'positive tests upon arrival \n(obs = 4)', xlab = 'number of positive tests',col = "lightblue", cex.lab = 1.5, cex.main =  1.5, cex.axis = 1.5)
#abline(v = 4, col = "red", lwd = 3)

hist(summary_out[,3], breaks = 10, main = 'positive tests on day 22 \n(obs = 142)', xlab = 'number of positive tests',col = "lightblue", cex.lab = 1.5, cex.main =  1.5, cex.axis = 1.5)
#abline(v = 142, col = "red", lwd = 3)

hist(summary_out[,4], breaks = 10, main = 'total infections over \n12-week period', xlab = 'number of infections',col = "lightblue", cex.lab = 1.5, cex.main =  1.5, cex.axis = 1.5)

outbreakprob = length(which(summary_out[,4]>10))/length(summary_out[,4])
print(outbreakprob)


######################################################################
# Fort Leonard Wood model
load('FLW_fit.RData')

hist(summary_out[,2], breaks = 10, main = 'positive tests upon arrival \n(obs = 0)', xlab = 'number of positive tests',col = "lightblue", cex.lab = 1.5, cex.main =  1.5, cex.axis = 1.5)
#abline(v = 0, col = "red", lwd = 3)

hist(summary_out[,3], breaks = 10, main = 'positive tests on day 18 \n(obs = 70)', xlab = 'number of positive tests',col = "lightblue", cex.lab = 1.5, cex.main =  1.5, cex.axis = 1.5)
#abline(v = 70, col = "red", lwd = 3)

hist(summary_out[,4], breaks = 10, main = 'total infections over \n12-week period', xlab = 'number of infections',col = "lightblue", cex.lab = 1.5, cex.main =  1.5, cex.axis = 1.5)

outbreakprob = length(which(summary_out[,4]>10))/length(summary_out[,4])
print(outbreakprob)