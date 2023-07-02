library(stats)
library(babynames)
library(dplyr)
library(naniar)
library(tidyverse)
require(boot)
require(matrixStats)
require(lubridate)
require(dplyr)
require(reshape2)
require(plotly)
require(ggrepel)
library(stringr)
library(ggplot2)

#paths
pathresults = "../results/uqrMaster2102_starchDynamic/"
lateDry = "10to18to25to25dry10871/" 
earlyDry = "10to11to18to25dry10870/" 
baseline = "10to18to25to25wet10869/" 

getVal <-function(myPath, valName, multFactor)#"transrate.txt"
{
  #setwd(myPath)
  x <- readLines(valName, n=-1)
  x = strsplit(x, ",");
  ncol =length(x[[length(x)]])
  trans = read.table(pathresults+myPath+valName, sep=",",header = FALSE,fill=TRUE, 
                     col.names = paste0("V",seq_len(ncol)))
  transum =rowSums(trans, na.rm = TRUE)
  transcumsum = cumsum(transum * multFactor) 
  trans2 = as.data.frame(transum)
  time =t(read.table(pathresults+myPath+"time.txt", sep=",",header = FALSE))
  trans2$time = time[1,]
  trans2$cumsum = transcumsum
  return(as.data.frame(trans2))
}
getValAll <-function(valName,fileName, multFactor)
{
  transearlyDry = getVal(earlyDry,fileName, multFactor )
  transearlyDry$scenario = "drySpell"
  transearlyDry$growth = "a) 11-18d"
  translateDry = getVal(lateDry,fileName, multFactor )
  translateDry$scenario = "drySpell"
  translateDry$growth = "b) 18-25d"
  transebaseline = getVal(baseline,fileName , multFactor)
  transebaseline$scenario = "baseline"
  transebaseline$growth = "b) 18-25d"
  
  trans = rbind(transearlyDry,translateDry,transebaseline)
  trans$var = valName
  return(trans)
}


############
#
# transpiration and assimilation (rate and cumsum)
#
############
trans = getValAll("Transpiration","transrate.txt", 1)
Ag = getValAll("Ag","AgPhl.txt", 1/24)

transAg = rbind(trans,Ag)

###############
#
# Graphic
#
###############
#water =18.01528 g/mol
transtemp = transAg
transtemp$cumsumtemp = NA
##old conversion:
#transtemp[transtemp$var == "Transpiration",]$cumsumtemp = 
# transtemp[transtemp$var == "Transpiration",]$cumsum *1*(1/0.01802)/100 #*10
#cumsum_cm3 * (g/cm3) * (mol/g) * (mmol/mol) * ratiotrans(for visualisation)
ratioTrans = 250 
transtemp[transtemp$var == "Transpiration",]$cumsum = 
  transtemp[transtemp$var == "Transpiration",]$cumsum *1000*(1/18.01528)/ratioTrans
#cumsum_cm3 * (mg/cm3) * (mmol/mg)  * scalingFactor
#transtemp[transtemp$var == "Ag",]$cumsum = transtemp[transtemp$var == "Ag",]$cumsum *5
transtemp$growth = ifelse(transtemp$growth==unique(transtemp$growth)[1],
                          "a) early dry spell",
                          "b) late dry spell")


var.lab3 = c("Transpiration" = "cumulative transpiration", 
             "Ag" = bquote("cumulative"~A["g"]))

vlabeller <- function (variable, value) {
  return(var.lab3[value])
}


both=ggplot(data=transtemp,#[transtemp$scenario!="baseline",], 
            aes(x=time, y=cumsum,col =var,linetype =scenario,
                                size=scenario,alpha=scenario))+
  geom_line()+
  scale_size_manual(breaks=c("drySpell","baseline"),name="",
                    values = c(1.5,0.5), drop = FALSE)+ 
  scale_alpha_manual(breaks=c("drySpell","baseline"),name="",
                     values = c(0.3,1), drop = FALSE)+ 
  scale_x_continuous(breaks=seq(from =6,to=28,by=2))+
  scale_y_continuous(sec.axis = sec_axis(~.*ratioTrans, name = "cumulative transpiration (mmol)"))+
  xlab('day of growth (d)')+
  facet_wrap(~growth,scales = "fixed", ncol = 1)+
  ylab(bquote(A["g"]~" (mmol)"))+#'cumulative Ag (mmol)')+ 
  geom_vline(data=filter(transtemp, growth==unique(transtemp$growth)[1]), aes(xintercept=11), colour="black") + 
  geom_vline(data=filter(transtemp, growth==unique(transtemp$growth)[1]), aes(xintercept=18), colour="black") + 
  geom_vline(data=filter(transtemp, growth==unique(transtemp$growth)[2]), aes(xintercept=18), colour="black") + 
  geom_vline(data=filter(transtemp, growth==unique(transtemp$growth)[2]), aes(xintercept=25), colour="black") + 
  scale_color_manual(breaks=c("Transpiration","Ag"),name="",
                     labels=vlabeller,#c("cumulative transpiration","cumulative Ag"),
                     values = c("#8DA0CB","#FC8D62"),
                     guide = "none")+  
  labs(colour = "", linetype ="")+
  theme(legend.text = element_text( size=15),
        legend.title = element_text( size=15),
        legend.position=c(0.2,0.3),
        #legend.direction = "horizontal",
        legend.box = "vertical",
        panel.grid.minor = element_line(colour="white"),
        panel.grid.major = element_line(colour="white"),
        panel.background = element_rect(fill="white", color="black"),
        strip.text.y =element_text( size=15) ,
        axis.text = element_text( size=17),
        axis.title.x = element_text( size=21),
        axis.title.y.left = element_text( size=21,color = "#FC8D62"),
        axis.title.y.right = element_text( size=21,color = "#8DA0CB"),
        axis.line.y.left = element_line(color = "#FC8D62"),
        axis.ticks.y.left = element_line(color = "#FC8D62"),
        axis.text.y.left = element_text(color = "#FC8D62"),
        axis.line.y.right = element_line(color = "#8DA0CB"),
        axis.ticks.y.right = element_line(color = "#8DA0CB"),
        axis.text.y.right = element_text(color = "#8DA0CB"),
        strip.background = element_rect(fill="white"),
        strip.text. = element_text(size = 21, hjust = 0,face ="bold"),
        strip.placement = "inside",
        legend.margin=margin()) ;both
ggsave('Agdryvswetwhite_all.png',plot=both,  width = 8, height = 10, dpi = 500)
