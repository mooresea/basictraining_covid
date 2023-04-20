library(igraph)

load('FB_fit.RData')

############################################################################
# plot network graph
x = 102
edges1 = edges[[x]]
date.symptoms1 = date.symptoms[x,]
isolate.date1 = isolate.date[x,]
quarantine.date1 = quarantine.date[x,]

campusCases = edges1[which(edges1[,1]!=0),]
symptomatic = !is.na(date.symptoms1[campusCases])
isolated = !is.na(isolate.date1[campusCases])
quarantined = !is.na(quarantine.date1[campusCases])

samecocoon = recruit.cocoon[campusCases[,1]]==recruit.cocoon[campusCases[,2]]
samecompany = recruit.company[campusCases[,1]]==recruit.company[campusCases[,2]]
staffinf = which(campusCases[,1]>numrecruits)


edgecol = rep('gray',length(campusCases[,1]))
edgecol[which(samecompany)] = 'red'
edgecol[which(samecocoon)] = 'blue'
edgecol[staffinf] = 'green'


contactGraph = graph_from_data_frame(campusCases)
plot(contactGraph, edge.arrow.size = 0, edge.color = edgecol,
     vertex.size = 5, 
     vertex.label=NA, edge.width = 4, vertex.frame.color = NA,
     vertex.color = 'black', 
     layout = layout_as_tree)


recruit.cocoon[campusCases[,1]]
recruit.cocoon[campusCases[,2]]
length(which(recruit.cocoon[campusCases[,1]]==recruit.cocoon[campusCases[,2]]))/(length(campusCases)/2)
length(which(recruit.company[campusCases[,1]]==recruit.company[campusCases[,2]]))/(length(campusCases)/2)

legend('topright', c('cocoon','company','staff','other'), col = c('blue','red','green','gray'), lty = 1, lwd = 4)

############################################################################
# stacked bar chart for contact breakdown

samecocoon = vector()
samecompany = vector()
staffinf = vector()
other = vector()

for (x in 1:500){
        #print(x)
        
        edges1 = edges[[x]]
        
        campusCases = edges1[which(edges1[,1]!=0),]
        
        if (length(campusCases)/2 > 1){
                samecocoon[x] = length(which(
                        recruit.cocoon[campusCases[,1]]==recruit.cocoon[campusCases[,2]]))/(length(campusCases)/2)
                samecompany[x] = length(which(
                        recruit.company[campusCases[,1]]==recruit.company[campusCases[,2]]))/(length(campusCases)/2)
                staffinf[x] = length(which(campusCases[,1]>numrecruits))/(length(campusCases)/2)
                other[x] = 1-samecocoon[x]-samecompany[x]-staffinf[x]
        } else if (length(campusCases)/2== 1){
                samecocoon[x] = length(which(recruit.cocoon[campusCases[1]]==recruit.cocoon[campusCases[2]]))/(length(campusCases)/2)
                samecompany[x] = length(which(recruit.company[campusCases[1]]==recruit.company[campusCases[2]]))/(length(campusCases)/2)
                staffinf[x] = length(which(campusCases[1]>numrecruits))/(length(campusCases)/2)
                other[x] = 1-samecocoon[x]-samecompany[x]-staffinf[x]
        }
}

idx = c(11,52,100,202,450) # selection of simulations to plot
dat = matrix(c(samecocoon[idx],samecompany[idx],staffinf[idx],other[idx]),ncol = 5,byrow = T)
colnames(dat) = c('sim1','sim2','sim3','sim4','sim5')
par(mar=c(5, 4, 4, 10), xpd=TRUE)
barplot(dat, col = c('blue','red','green','gray'))

legend('topright', inset=c(-0.3,0), c('cocoon','company','staff','other'), 
       col = c('blue','red','green','gray'),lty = 1, lwd = 4)

