# Nom               : Article 1
# Type              : Propre
# Objet             : Community analysis through Bray Curtis,CTA, PTA, Costatis
# Input             : REPHY + RHLN datasets, Environmental data (see article for references) - treated datasets
# Output            : Given in article "Decadal trajectories of phytoplankton communities in contrasted estuarine systems in an epicontinental Sea"
# Auteur            : ANGELINE LEFRAN
# R version         : 3.5.1, (ggplot2 : V3.1.0)
# Encoding          : ISO 8859-1 (default)
# Date de création  : 15 MARCH 2020
# Date de modification : 28 MAY 2021
#______________________________________________________________________________________________

#### ANALYSIS START HERE ####

#Macros used :
source("Macros/Sp99percent.R")



#Fonction imported
library(ggplot2) #graph
library(ggmap)
library(lubridate) #date
library(vegan)
library(ade4)
library(vegclust)
library(smacof)
library(adegraphics) #costatis representation
library(indicspecies)
library(factoextra)
library(NbClust)
library(viridis) #color
library(RColorBrewer)
library(tidyr) #data-play
library(dplyr)
library(plyr)
library(Kendall) #trends
library(fume) #si autocor
library(grid)

# Before hand :

Sites<-data.frame(Num=c("010-P-001","010-P-109","014-P-023","018-P-054","006-P-001","022-P-002"),
                  Nom=c("Antifer","Cabourg", "Géfosse", "Donville","At so", "StCast_Hebihens"),
                  Code=c("ANT","CAB","GEF","DONV","ATSO","STCA"))

Dataset<- read.table(file = "Derived_datasets/PhytoData_complete_final-05-2021.csv",header=T,sep= ";", dec=".")
Dataset$Passage...Date<-as.Date(Dataset$Passage...Date,"%d/%m/%Y") #Warning for NAs ! Go back to csv file to change the column type to date
Dataset$Résultat...Valeur.de.la.mesure<-as.numeric(as.character(Dataset$Résultat...Valeur.de.la.mesure)) #as.numeric isn't enough to convert a factor class
#if NAs be sure the csv only contains numeric and not 1E6 that is character
Dataset<-Dataset[Dataset$Passage...Année >=2008, ]
Dataset<-Dataset[Dataset$Lieu.de.surveillance...Mnémonique %in% Sites$Num[1:6],]

for (i in (1:length(Dataset$Passage...Date))){
  Dataset$Station[i]<-as.character(Sites$Code[Sites$Num == as.character(Dataset$Lieu.de.surveillance...Mnémonique[i])])
} #Creating a column with sites' code


#Data selection according to sampling strategy -> will be called again at the beginning of each analysis
Tot<-Dataset[Dataset$Résultat...Code.paramètre == "FLORTOT",] #


#### ****************
###### I) Anomalies from environement variable ####
#### ****************

# Preparation of the dataset
DataEnv<-read.csv("Derived_datasets/DataEnv_Costatis-Nut-09-2020.csv",header = T,sep= ";", dec=".") #DataEnv_Costatis-Nut-04-2020
DataEnv<-DataEnv[DataEnv$Lieu.de.surveillance...Libellé %in% Sites$Nom,]

DataEnv<-DataEnv[DataEnv$Passage...Année >=2008, ]

#***:Changing the winter's reference year for example : Winter 2007 = Dec 2007 + Jan et Feb 2008 and not Dec 2007 + Jan et Feb 2007 ***
DataEnv$Passage...Année[DataEnv$Passage...Mois %in% c(1,2)]<-DataEnv$Passage...Année[DataEnv$Passage...Mois %in% c(1,2)]-1
#*********

for (i in (1:dim(DataEnv)[1])){
  DataEnv$Station[i]<-as.character(Sites$Code[Sites$Nom == as.character(DataEnv$Lieu.de.surveillance...Libellé[i])])
} #Creating a column with sites' code

{
  # all parameters available : c("NH4","NO3+NO2","PO4","SALI","SIOH","TEMP","Précipitations..mm.","Durée.d.insolation...min.Vitesse..ms.","Débit..m3s.")
  DataAdj<-DataEnv[DataEnv$Résultat...Code.paramètre %in% c("NH4","NO3+NO2","PO4"),]
  DataNsP<-ddply(DataAdj,.(Saison, Station, Passage...Année, Lieu.de.surveillance...Libellé),function(x)data.frame(
    Mean = sum(x$Mean[x$Résultat...Code.paramètre %in% c("NH4","NO3+NO2")]) / x$Mean[x$Résultat...Code.paramètre %in% c("PO4")],
    Med = sum(x$Med[x$Résultat...Code.paramètre %in% c("NH4","NO3+NO2")]) / x$Med[x$Résultat...Code.paramètre %in% c("PO4")],
    Qone = sum(x$Qone[x$Résultat...Code.paramètre %in% c("NH4","NO3+NO2")]) / x$Qone[x$Résultat...Code.paramètre %in% c("PO4")],
    Qthree = sum(x$Qthree[x$Résultat...Code.paramètre %in% c("NH4","NO3+NO2")]) / x$Qthree[x$Résultat...Code.paramètre %in% c("PO4")] ))
  
  
  DataNsP<-tibble::add_column(DataNsP, Résultat...Code.paramètre = rep("N/P",dim(DataNsP)[1]), .after = "Passage...Année")  
  
  DataEnv<-rbind(DataEnv,DataNsP)
}# Creation of the N/P parameter and selection of the abiotic variables  

{ 
  Opti_SW_Wind <- 270
  DataEnv$Mean[DataEnv$Résultat...Code.paramètre == "Direction..degres."]<- 
    abs(100 - abs(DataEnv$Mean[DataEnv$Résultat...Code.paramètre == "Direction..degres."] - Opti_SW_Wind) / 180 *100)
  
  DataEnv$Med[DataEnv$Résultat...Code.paramètre == "Direction..degres."]<- 
    abs(100 - abs(DataEnv$Med[DataEnv$Résultat...Code.paramètre == "Direction..degres."] - Opti_SW_Wind) / 180 *100)
  
  DataEnv$Qone[DataEnv$Résultat...Code.paramètre == "Direction..degres."]<- 
    abs(100 - abs(DataEnv$Qone[DataEnv$Résultat...Code.paramètre == "Direction..degres."] - Opti_SW_Wind) / 180 *100)
  
  DataEnv$Qthree[DataEnv$Résultat...Code.paramètre == "Direction..degres."]<- 
    abs(100 - abs(DataEnv$Qthree[DataEnv$Résultat...Code.paramètre == "Direction..degres."] - Opti_SW_Wind) / 180 *100)
  
}#transformation of the wind direction values so that it is a similarity to an optimal wind direction


ChosenVar<- c("N/P","SALI","TEMP","Précipitations..mm.","Durée.d.insolation...min.", "Débit..m3s.","Vitesse..ms.","Direction..degres.","OXYGENE","TURB")

DataEnv<-DataEnv[DataEnv$Résultat...Code.paramètre %in% ChosenVar,]

Missing_value_replacement<-mean(DataEnv$Mean[DataEnv$Saison=="Winter" & DataEnv$Station == "STCA" & DataEnv$Résultat...Code.paramètre == "N/P" ])


#Changing the parameters' name
DataEnv$Résultat...Code.paramètre<-as.character(DataEnv$Résultat...Code.paramètre)
DataEnv$Résultat...Code.paramètre[DataEnv$Résultat...Code.paramètre %in% "Débit..m3s."] <- "Flow (m^3/s)"
DataEnv$Résultat...Code.paramètre[DataEnv$Résultat...Code.paramètre %in% "Durée.d.insolation...min."] <- "Daylight (min)"
DataEnv$Résultat...Code.paramètre[DataEnv$Résultat...Code.paramètre %in% "Précipitations..mm."] <- "Rainfall (mm)"
DataEnv$Résultat...Code.paramètre[DataEnv$Résultat...Code.paramètre %in% "Vitesse..ms."] <- "Wind speed (m/s)"
DataEnv$Résultat...Code.paramètre[DataEnv$Résultat...Code.paramètre %in% "Direction..degres."] <- "West Wind (%)"
DataEnv$Résultat...Code.paramètre[DataEnv$Résultat...Code.paramètre %in% "OXYGENE"] <- "Oxygen (mg/l)"
DataEnv$Résultat...Code.paramètre[DataEnv$Résultat...Code.paramètre %in% "SALI"] <- "Salinity (PSU)" 
DataEnv$Résultat...Code.paramètre[DataEnv$Résultat...Code.paramètre %in% "TEMP"] <- "Water temperature (°C)"
DataEnv$Résultat...Code.paramètre[DataEnv$Résultat...Code.paramètre %in% "TURB"] <- "Turbidity (NTU)"


##Facet 9 parameters anomaly
{
  DataIn<-DataEnv
  library(lubridate)
  library(ggplot2)
  library(plyr)
  
  Table<-data.frame(date=DataIn$Saison,
                    year=DataIn$Passage...Année,
                    Station=DataIn$Station,
                    Val=DataIn$Mean,
                    Para=DataIn$Résultat...Code.paramètre)
  
  
  #per site and season we create the anomaly value
  for (p in (1:length(unique(Table$Para)))) {
    Parameter <- unique(Table$Para)[p]
    
    for (s in (unique(Table$Station))) {
      means<-c()
      n<-0
      for (i in unique(Table$date)) {
        n<-n+1
        means<-c(means,mean(na.omit(Table$Val[Table$date==i & Table$Para==Parameter & Table$Station == s ])))
        Table$Val[Table$date==i & Table$Para==Parameter & Table$Station == s]<-Table$Val[Table$date==i & 
                                                                                           Table$Para==Parameter & 
                                                                                           Table$Station == s] - means[n]
        #exception of the flow, too different between station so had to be calculated in percentage
        if (Parameter == "Flow (m^3/s)") {
          Table$Val[Table$date==i & Table$Para=="Flow (m^3/s)" & Table$Station == s]<- (Table$Val[Table$date==i & 
                                                                                                    Table$Para=="Flow (m^3/s)" & 
                                                                                                    Table$Station == s] / means[n]) * 100
        }
      }
    }
  }
  
  
  Table$date<-factor(Table$date, levels = c("Winter", "Spring", "Summer","Autumn"))
  
  Table$date<-as.character(Table$date)
  Table$datenum[Table$date == "Winter"] <- 4
  Table$datenum[Table$date == "Spring"] <- 1 
  Table$datenum[Table$date == "Summer"] <- 2
  Table$datenum[Table$date == "Autumn"] <- 3
  
  Table$period<-paste(Table$year,Table$datenum,sep="-")
  Table$period2<-paste(Table$year,Table$date,sep="-")
  
  
  Table$col<-ifelse(Table$Val>=0,"red","blue")
  
  Table<-Table[order(Table$period),]
  
  Compl<- ddply (Table, .(date,year,period,period2,Para), 	function(x)data.frame(Mean = mean(x$Val, na.rm = T),
                                                                                 Med = median(x$Val,na.rm = T),
                                                                                 Qone = quantile(x$Val,na.rm = T) [2],
                                                                                 Qthree= quantile(x$Val,na.rm = T) [4]))
  
  Compl<-Compl[order(Compl$period),] 
  Compl$col<-ifelse(Compl$Med>=0,"blue","red")
  
  Table<-Table[-1,]
  Table$period<-as.factor(Table$period)
  
  pp<-ggplot(data = Compl,aes(x=period,y=Med, group=1))+ theme_bw() +
    geom_line(lwd=1)+
    #scale_color_continuous(low = "yellow", high = "brown")+
    #scale_color_manual(values = c("blue","red"))+
    geom_ribbon(aes(x=period, ymin = Qone,ymax = Med),fill="grey", alpha =0.5)+
    geom_ribbon(aes(x= period, ymin = Med ,ymax = Qthree),fill="grey", alpha =0.5) +
    facet_grid(rows = vars(Para),scales = "free")+
    
    labs(x = "Year-Season (1=spring, 2=summer, 3=autumn, 4=winter)", y = paste("Anomalies of parameters")) +
    geom_vline(xintercept=which(levels(Table$period) %in% c("2008-1","2009-1","2010-1","2011-1","2012-1","2013-1","2014-1","2015-1","2016-1",
                                                            "2017-1","2018-1","2019-1")),linetype="dashed",col="darkgrey") +
    geom_hline(yintercept =0,linetype="dashed",col="salmon")+
    geom_point(data = Table,mapping = aes(x=period,y=Val, col=Station),size=2 )+  #colour="grey", pch=21,,
    scale_color_manual(values=c(brewer.pal(6,"Dark2"))) +
    theme(axis.text.x=element_text(angle=60, hjust=1))
  
}    

#pdf(paste("Figures/anomalies_facet3_legend",format(Sys.Date(), "%d%m%Y"),".pdf") , width = 10 , height = 15)
png(paste("Figures/anomalies_facet4",format(Sys.Date(), "%d%m%Y"),".png") , width = 300 , height = 450, units ="mm",res =300,family ="Arial",pointsize =10)
pp
dev.off()

#gridExtra::grid.arrange(pp, bottom="Figure 2: Seasonal anomalies of environmental variables over the period of 2008-2019 \n(black line: median value between all stations; grey ribbon: Q1 and Q3 limits). ")


###### II) Dissimilarity : Bray Curtis ####

Tabsite<-Tot[order(as.Date(Tot$Passage...Date, format="%d/%m/%Y")),]
Tabsite<-Tabsite[Tabsite$Passage...Année >=2008, ] #study period
Tabsite<-Tabsite[Tabsite$Passage...Année <=2019, ] 
Tabsite<-Tabsite[Tabsite$Résultat...Valeur.de.la.mesure >100, ] #limit of detection

Tabtouse<-Sp99(Tabsite,Sites)

Tabtouse<-Tabtouse[-which(Tabtouse$Résultat...Nom.du.taxon %in% c("Pennées","Centriques","Peridiniales")) ,]

Plotlist<-vector('list',6)

for (o in 1:6) {
  Lieu<-Sites$Num[o]
  Tabsite<-Tabtouse[Tabtouse$Lieu.de.surveillance...Mnémonique == as.character(Lieu),]
  Tabsite$Résultat...Valeur.de.la.mesure<-round(log10(Tabsite$Résultat...Valeur.de.la.mesure),digits=2)
  L<-length(unique(Tabsite$Passage...Date))
  S<-length(unique(Tabsite$Résultat...Nom.du.taxon))
  
  Tabsite<-Tabsite[colnames(Tabsite) %in% c("Passage...Date","Résultat...Nom.du.taxon","Résultat...Valeur.de.la.mesure")]
  matrixfauna <-   tidyr::pivot_wider(data=Tabsite, id_cols = c(Passage...Date), names_from=Résultat...Nom.du.taxon, values_from=Résultat...Valeur.de.la.mesure,values_fill = list(Résultat...Valeur.de.la.mesure=0)) 
  #Saving the first column for dates
  Temps<-unique(matrixfauna$Passage...Date)
  matrixfauna<-matrixfauna[,-(1)]
  
  #####Bray-Curtis index####
  matrixfauna.diss <- vegdist(matrixfauna,"bray") #givve distancy list
  matrixfauna.diss2<- as.matrix(vegdist(matrixfauna,"bray")) #create the dissimilarity matrix
  matrixfauna.sim2<-1-matrixfauna.diss2
  #hist(matrixfauna.sim2)
  
  
  ####Building the parallel matrix for time intervals####
  matrixfauna.tmp <-matrix(0,ncol=L, nrow=L)
  
  for (x in 1:L) {
    
    for (y in 1:L) {
      
      Inval<-as.period(interval(ymd(Temps[x]), ymd(Temps[y])),unit="days") #number of days between two dates
      matrixfauna.tmp[x,y]<-as.numeric(Inval, "days")
    }
  }
  
  colnames(matrixfauna.tmp) <- seq(1:L)
  rownames(matrixfauna.tmp) <- seq(1:L)
  matrixfauna.tmp<-abs(matrixfauna.tmp)
  matrixfauna.tmp[1:3,1:3]
  
  ##Working the time intervals matrix, we want values proportionnal to 15 days, so we round the real days interval to the closest 15 multiplier
  matrixfauna.tmp15<-matrix(0,ncol=L, nrow=L)
  
  for (A in 1:L) {
    
    for (B in 1:L) {
      INI<-matrixfauna.tmp[A,B]
      
      if ( is.integer(INI/15) == F ){
        Before<- (as.integer(INI/15))*15
        Diff1<-INI-Before
        After<- (as.integer(INI/15)+1)*15
        Diff2<-After-INI
        Tobb<-as.data.frame(rbind(c(Before,Diff1),c(After,Diff2))) #possible values before and after the real interval 
        colnames(Tobb)<-c("Val","Diff")
        matrixfauna.tmp15[A,B]<-Tobb$Val[Tobb$Diff <= 7] #Replacing in the matrix the real interval with the closest 15* interval 
      }
    }
  }
  matrixfauna.tmp15[1:3,1:3]
  
  ####Creating the final table and the plot ####
  Final<-data.frame(matrix(ncol = 2, nrow = 0))
  
  for (l in 0:(max(matrixfauna.tmp15)/15)) {
    
    N<-length(matrixfauna.sim2[which(matrixfauna.tmp15 == 15*l)])/2     # number of points used to calculate the variable below, /2 because we only consider using half the values of a full matrix (upper vs lower diagonals) 
    Sims <- mean(matrixfauna.sim2[which(matrixfauna.tmp15 == 15*l)]) #Mean/Median/MAx etc... of all the values of the same interval index (proportional to 15 days)
    Q1 <- quantile(matrixfauna.sim2[which(matrixfauna.tmp15 == 15*l)])[2] #Q25%
    Q2 <- quantile(matrixfauna.sim2[which(matrixfauna.tmp15 == 15*l)])[3] #Median
    Q3 <- quantile(matrixfauna.sim2[which(matrixfauna.tmp15 == 15*l)])[4] #Q75%
    D=15*l
    Y<-D/365
    
    Final<-rbind(Final,c(Sims,D,Y,Q1,Q2,Q3,N))
    
  }
  colnames(Final)<-c("MeanSims","Days","Year","Q1","Q2","Q3","N")
  Final=round(Final,digits=3)
  Final=na.omit(Final)
  
  
  #Not plotted this part simulate the stat_smooth function for the obtention of the equation values and the R²
  Object<-loess(Final$MeanSims~Final$Days) 
  PP<-summary(Object)$fitted  #Isolation of loess() resulting values 
  
  AmplBegin<- round(summary(loess(Final$Q3~Final$Days) )$fitted[1] - summary(loess(Final$Q1~Final$Days) )$fitted[1],digits=3)
  AmplEnd<- round(summary(loess(Final$Q3~Final$Days) )$fitted[length(Final$Days)] - summary(loess(Final$Q1~Final$Days) )$fitted[length(Final$Days)],digits=3)
  AmplLoss<- AmplBegin - AmplEnd
  
  Mk<-mkTrend(PP)
  
  if (o==1) {
    Ant=mkTrend(PP)
  } else if (o==2) {
    Cab=mkTrend(PP)
  } else if (o==3) {
    Gef=mkTrend(PP)
  } else if (o==4) {
    Don=mkTrend(PP)
  } else if (o==5) {
    At=mkTrend(PP)
  } else if (o==6) {
    St=mkTrend(PP)
  }
  
  
  #MValues from the mktrend
  trend<-Mk$`Corrected p.value`
  if ( trend < 0.01) {
    trend <- "< 0.01*"
  } else if (trend < 0.05){
    trend <- "< 0.05*"
  }  else {
    trend <- "> 0.05"
  }
  
  
  senslope<-round(Mk$`Sen's Slope`,digits=8)
  
  xfin<-(max(Final$Days)-round(max(Final$Days)*0.05))/365 #I hide the end of the chart, removing 5% to not plot the noise due to the lack of data between samples with 10 years apart
  
  Three<-Final$Q3
  One<-Final$Q1
  All<-Final$N
  Avera<-Final$MeanSims
  
  
  p<- ggplot(data=Final) + geom_line(aes(x=Year,y=MeanSims)) + geom_point(aes(x=Year,y=MeanSims))+ 
    #ggtitle(paste("Resilience in the water composition in", Sites$Code[o], "(2008-2018)"))+ 
    stat_smooth(aes(x=Year,y=MeanSims),show.legend=T,se = FALSE, fullrange = T) + 
    stat_smooth(aes(x=Year,y=Q1),lty=2,col="red",se = FALSE) +
    stat_smooth(aes(x=Year,y=Q3),lty=2,col="red",se = FALSE) +
    scale_color_manual(limits=c(4000,1),breaks=c("Q1,Q3","Mean"),values = c("red","blue"))+ 
    labs(y="Similarities between samples", x="Time Lapse between samples (Years)") +
    theme_classic()+
    theme_bw(base_line_size = 0)+
    ylim(c(0.2,0.8))+
    xlim(c(0,xfin))+ #I hide the end of the chart, removing 5% to not plot the noise due to the lack of data between samples with 10 years apart
    geom_text(x=2.5, y=0.8,col="Black", size = 5, label=paste(Sites$Code[o],"(2008-2019)"))+ #we take the mean value to make a point not the median
    geom_text(x=2.5, y=0.8,col="Black", size = 5,label=paste("\n\n\n Slope =",
                                                             format(senslope,scientific = TRUE,digits=3)," (p-value",trend,"), \n Amplitude loss ~ -",
                                                             AmplLoss,", N ~ ",round(mean(All))))+      #equation of lm(loess(data))
    geom_text(x=xfin+0.25, y=mean(Three)-0.03,col="Red",size = 6, label=paste("Q3"))+ 
    geom_text(x=xfin+0.25, y=mean(One)-0.03,col="Red",size = 6, label=paste("Q1"))+
    geom_text(x=xfin+0.25, y=mean(Avera)-0.03,col="Blue",size = 6, label=paste("µ"))
  
  Plotlist[[o]]<- p
  print(paste(Lieu,"fini"))
  
}

p1 <- Plotlist[[1]]
p2 <- Plotlist[[2]]
p3 <- Plotlist[[3]]
p4 <- Plotlist[[4]]
p5 <- Plotlist[[5]]
p6 <- Plotlist[[6]]

figbray <- cowplot::plot_grid(p1,p2,p3,p4,p5,p6, ncol=2, nrow=3)

png(file= paste("Figures/", format(Sys.Date(), "%d%m%Y"),Sites$Nom[6],"-BrayCurtis.png"),width = 210 , height = 150, units = "mm" , res = 300 , family = "Arial" , pointsize = 10 )
#pdf(paste("Figures/", format(Sys.Date(), "%d%m%Y"),"-BrayCurtis_legend.pdf"), width = 16 , height = 16)
#gridExtra::grid.arrange(figbray) #, bottom="Figure 2: Seasonal anomalies of environmental variables over the period of 2008-2019 \n(black line: median value between all stations; grey ribbon: Q1 and Q3 limits). ")
p6
dev.off()

# pdf(file= paste("Figures/", format(Sys.Date(), "%d%m%Y"),"-BrayCurtis2.pdf"),width = 18 , height = 18)
# print(fig)
# dev.off()


#### ****************
###### III) THE CTA- Community Trajectories Analysis ####
#### ****************

# library(vegclust)
# library(RColorBrewer)
# library(smacof)
# library(vegan)
# library(gridExtra)

##Importation
o<-c(1:6)
Lieu<-Sites$Num[o]
Tabsite<-Tot[Tot$Lieu.de.surveillance...Mnémonique %in% as.character(Lieu),]
Tabsite<-Tabsite[order(as.Date(Tabsite$Passage...Date, format="%d/%m/%Y")),]
Tabsite<-Tabsite[Tabsite$Passage...Année >=2008, ]
Tabsite<-Tabsite[Tabsite$Passage...Année <=2019, ]#period targetted

Tabsite<-Sp99(Tabsite,Sites)
Tabsite<-Tabsite[-which(Tabsite$Résultat...Nom.du.taxon %in% c("Pennées","Centriques","Peridiniales")) ,]

FirstY<- min(Tabsite$Passage...Année)
LastY<-max(Tabsite$Passage...Année)
Tabsite$Résultat...Valeur.de.la.mesure[Tabsite$Résultat...Valeur.de.la.mesure ==0]<-100
Tabsite<-Tabsite[Tabsite$Résultat...Valeur.de.la.mesure >100, ]
Tabsite$Résultat...Valeur.de.la.mesure<-round(log10(Tabsite$Résultat...Valeur.de.la.mesure),digits=2)

#***:Changing the winter's reference year for example : Winter 2007 = Dec 2007 + Jan et Feb 2008 and not Dec 2007 + Jan et Feb 2007 ***
Tabsite$Passage...Année[Tabsite$Passage...Mois %in% c(1,2)]<-Tabsite$Passage...Année[Tabsite$Passage...Mois %in% c(1,2)]-1
#*********

Save<-Tabsite
#
#Axes :
a<-1
b<-2

pdf(file= paste("Figures/", format(Sys.Date(), "%d%m%Y"),"CTA_legend with seasonal facets-Axes",a,b,".pdf"),width = 8 , height = 8)

par(mfrow = c(2,2),mar=c(7, 4, 1, 2) + 0.1) #default :  c(5, 4, 4, 2) + 0.1 (bottom, left, top, right)
#Legend on Spring only
{s<-"Spring"
  Tabsite<-Save[Save$Saison == s, ] #précision sur les données de départ
  
  # Analyses methods starts now :
  require(plyr)
  Tabsite$Résultat...Valeur.de.la.mesure<-round(log10(Tabsite$Résultat...Valeur.de.la.mesure),digits=2)
  
  colnames(Tabsite)[colnames(Tabsite)=="Résultat...Nom.du.taxon"]<-"Code"
  
  
  ##Thanks love :##
  Total <- ddply (Tabsite, .(Passage...Année,Code,Lieu.de.surveillance...Libellé), 
                  function(x)data.frame(Total = mean(x$Résultat...Valeur.de.la.mesure)))
  
  ##Binding year and season
  #
  colnames(Total)<-c("Passage...Date","Code","Lieu.de.surveillance...Libellé","Résultat...Valeur.de.la.mesure")
  
  #Binding this to the station name
  Total$group<-paste(Total$Passage...Date,Total$Lieu.de.surveillance...Libellé,sep="_")
  Total<-Total[,-c(1,3)]
  Total[1:5,]
  
  library(tidyr)
  library(ade4)
  matrixflora<-spread(Total,Code,Résultat...Valeur.de.la.mesure,fill=0)
  matrixflora[1:5,1:5]
  
  #Saving the first column apart for later
  Keep<-matrixflora[,1]
  matrixflora<-matrixflora[,-1]
  
  factor<-rep(unique(Tabsite$Lieu.de.surveillance...Libellé),each=15)
  substr(Keep,1,7)
  
  #####Bray-Curtis####
  matrixflora.diss <- vegdist(matrixflora,"bray") #renvoie une liste des distances
  matrixflora.diss2<- as.matrix(vegdist(matrixflora,"bray")) #créer la matrice de dissimilarité
  matrixflora.sim2<-1-matrixflora.diss2
  #hist(matrixflora.sim2)
  #write.table(matrixfauna.sim2,file = "Derived data sets/Matrice_sim_Antifer.csv",row.names=F,sep= ";", dec=".")
  
  
  #For seasonal Display
  trajectoryPCoA(matrixflora.diss, substr(Keep,6,36), substr(Keep,1,4), traj.colors = brewer.pal(6,"Dark2"), lwd = 2,length=0.1,axes=c(a,b))
  legend("bottomright", col=brewer.pal(8,"Dark2"),
         legend= paste(Sites$Code[order(Sites$Code[1:6])]), bty="n", lty=1, lwd = 2)
  legend("topright", legend= paste(s), bty="n")
  #title(paste("Community Trajectory Analysis",s),cex.main = 1.1)
  
  
  # ## To make a table of trajectories lenght :
  # P<-  trajectoryPCoA(matrixflora.diss, substr(Keep,6,36), substr(Keep,1,4), 
  #                     traj.colors = brewer.pal(6,"Dark2"), lwd = 2,length=0.1,axes=c(a,b))
  # 
  # Coord<- data.frame(Axe_1=P$points[,1],Axe_2=P$points[,2],site=substr(Keep,6,36),year=as.integer(substr(Keep,1,4)))
  # 
  # write.table(Coord,file = paste("Derived data sets/Matrice_CTA_",s,".csv"),row.names=F,sep= ";", dec=".")
  # 
} #Spring with station legend


for (s in (c("Summer","Autumn","Winter"))) {
  Tabsite<-Save[Save$Saison == s, ] #précision sur les données de départ
  
  # Analyses methods starts now :
  require(plyr)
  Tabsite$Résultat...Valeur.de.la.mesure<-round(log10(Tabsite$Résultat...Valeur.de.la.mesure),digits=2)
  
  colnames(Tabsite)[colnames(Tabsite)=="Résultat...Nom.du.taxon"]<-"Code"
  
  
  ##Preparing the dataset :##
  Total <- ddply (Tabsite, .(Passage...Année,Code,Lieu.de.surveillance...Libellé), 
                  function(x)data.frame(Total = mean(x$Résultat...Valeur.de.la.mesure)))
  
  ##Binding year and season
  #
  colnames(Total)<-c("Passage...Date","Code","Lieu.de.surveillance...Libellé","Résultat...Valeur.de.la.mesure")
  
  #Binding this to the station name
  Total$group<-paste(Total$Passage...Date,Total$Lieu.de.surveillance...Libellé,sep="_")
  Total<-Total[,-c(1,3)]
  Total[1:5,]
  
  library(tidyr)
  library(ade4)
  matrixflora<-spread(Total,Code,Résultat...Valeur.de.la.mesure,fill=0)
  matrixflora[1:5,1:5]
  
  #Saving the first column apart for later
  Keep<-matrixflora[,1]
  matrixflora<-matrixflora[,-1]
  
  factor<-rep(unique(Tabsite$Lieu.de.surveillance...Libellé),each=15)
  substr(Keep,1,7)
  
  #####Bray-Curtis####
  matrixflora.diss <- vegdist(matrixflora,"bray") #renvoie une liste des distances
  matrixflora.diss2<- as.matrix(vegdist(matrixflora,"bray")) #créer la matrice de dissimilarité
  matrixflora.sim2<-1-matrixflora.diss2
  #hist(matrixflora.sim2)
  #write.table(matrixfauna.sim2,file = "Derived data sets/Matrice_sim_Antifer.csv",row.names=F,sep= ";", dec=".")
  
  
  #For seasonal Display
  trajectoryPCoA(matrixflora.diss, substr(Keep,6,36), substr(Keep,1,4), traj.colors = brewer.pal(6,"Dark2"), lwd = 2,length=0.1,axes=c(a,b))
  #legend("bottomright", col=brewer.pal(8,"Dark2"),
  #legend= paste(Sites$Code[order(Sites$Code[1:6])]), bty="n", lty=1, lwd = 2)
  legend("topright", legend= paste(s), bty="n")
  #title(paste("Community Trajectory Analysis",s),cex.main = 1.1)
  
  ## To plot the coordinates on the axis 1 & 2:
  
  #Coord<- data.frame(Axe_1=P$points[,1],Axe_2=P$points[,2],site=substr(Keep,6,36),year=as.integer(substr(Keep,1,4)))
  
  # CTA_1 <- ggplot(Coord) + geom_line(aes(x=year,y=Axe_1,col=site),size=1) + 
  #   ggtitle(paste("Scores Axis 1",s)) + scale_color_manual(values=brewer.pal(6,"Dark2"))+
  #   ylab(paste("Axis 1")) + xlab(paste("Year"))
  # plot(CTA_1)
  # 
  # CTA_2<- ggplot(Coord) + geom_line(aes(x=year,y=Axe_2,col=site),size=1) + 
  #   ggtitle(paste("Scores Axis 2 ",s)) + scale_color_manual(values=brewer.pal(6,"Dark2"))+
  #   ylab(paste("Axis 2")) + xlab(paste("Year"))
  # plot(CTA_2)
  
  
  # ## To make a table of trajectories lenght :
  # P<-  trajectoryPCoA(matrixflora.diss, substr(Keep,6,36), substr(Keep,1,4), 
  #                     traj.colors = brewer.pal(6,"Dark2"), lwd = 2,length=0.1,axes=c(a,b))
  # 
  # Coord<- data.frame(Axe_1=P$points[,1],Axe_2=P$points[,2],site=substr(Keep,6,36),year=as.integer(substr(Keep,1,4)))
  # 
  # write.table(Coord,file = paste("Derived_datasets/Matrice_CTA_",s,".csv"),row.names=F,sep= ";", dec=".")
  # 
  
  
} # facets for seasons without station legend

dev.off()


#### ****************
###### IV) Costatis - Double Partial Triadic Analysis - Annual k-table = season + year - Slimani's plot ####
#### ****************

{
  ChosenSeason <- unique(Tot$Saison) #[3] ## order need to be the one of the table for later use
  #s<-ChosenSeason
  VisualName<- c("ANT","ATSO","CAB","DONV","GEF","STCA")
  
  #### ****************
  #### A - PTA 1 : Flore <<<<<< ####
  #### ****************
  {
    ##Préparation
    o<-c(1:6) 
    Lieu<-Sites$Num[order(Sites$Code)][o]
    Tabsite<-Tot[Tot$Lieu.de.surveillance...Mnémonique %in% as.character(Lieu),]
    Tabsite<-Tabsite[order(as.Date(Tabsite$Passage...Date, format="%d/%m/%Y")),]
    
    #***:Changing the winter's reference year for example : Winter 2007 = Dec 2007 + Jan et Feb 2008 and not Dec 2007 + Jan et Feb 2007 ***
    Tabsite$Passage...Année[Tabsite$Passage...Mois %in% c(1,2)]<-Tabsite$Passage...Année[Tabsite$Passage...Mois %in% c(1,2)]-1
    #*********
    
    FirstY<- min(Tabsite$Passage...Année)
    LastY<-max(Tabsite$Passage...Année)
    Tabsite$Résultat...Valeur.de.la.mesure[Tabsite$Résultat...Valeur.de.la.mesure ==0]<-100
    Tabsite<-Tabsite[Tabsite$Résultat...Valeur.de.la.mesure >100, ]
    
    Tabsite<-Tabsite[Tabsite$Passage...Année >=2008, ]
    
    Tabsite<-Sp99(Tabsite,Sites)
    Tabsite<-Tabsite[-which(Tabsite$Résultat...Nom.du.taxon %in% c("Pennées","Centriques","Peridiniales")) ,]
    
    # >>>>>>>> methods starts now :<<<<<<<<
    
    require(plyr)
    Tabsite$Résultat...Valeur.de.la.mesure<-round(log10(Tabsite$Résultat...Valeur.de.la.mesure),digits=2)
    
    ##Aggregation per season, year and site##
    colnames(Tabsite)[colnames(Tabsite)=="Résultat...Nom.du.taxon"]<-"Code"
    Total <- ddply (Tabsite, .(Passage...Année,Saison,Code,Station), 
                    function(x)data.frame(Total = mean(x$Résultat...Valeur.de.la.mesure)))
    
    ##For season selection
    Total<-Total[Total$Saison %in% ChosenSeason,]
    ###
    
    ##Binding year and season on one column
    Total$Passage...Année<-paste(Total$Passage...Année,Total$Saison,sep="-")
    #Total<-Total[-which(Total$Passage...Année %in% "2019-Winter"),]
    
    Total<-Total[,colnames(Total) %in% c( "Passage...Année","Station","Code","Total")]
    #****************
    #
    
    matrixflora<-spread(Total,Code,Total,fill=0)
    
    #Binding dates and station name
    matrixflora<-matrixflora[order(matrixflora$Passage...Année),]
    matrixflora$group<-paste(matrixflora$Passage...Année,matrixflora$Station,sep="_")
    matrixflora[1:5,]
    
    matrixflora[1:5,1:5]
    
    #Saving the first column apart for later
    Keep<-matrixflora[,c(1,2)]
    matrixflora<-matrixflora[,-c(1,2)]
    rownames(matrixflora) <- matrixflora$group
    matrixflora<-matrixflora[,-dim(matrixflora)[2]]
    #
    
    pca1 <- dudi.pca (matrixflora, scale = F, scan = F, nf = 3)
    wit1 <- wca (pca1, as.factor(as.character(Keep$Passage...Année)) , scan = F) # = wit1 is considering all sites independently
    ktaflora <- ktab.within (wit1,colnames = Keep$Station) # creating k-tables, one for every sites, so you have to put your variable of time in order (unique) then repeat for the number of sites. 
    
    
    kta4 <- t(ktaflora)
    pta0 <- pta (kta4, scan = F, 3)
    plot(pta0,option=c(1:4))
  } #Waiting for run
  
  
  pdf(file= paste("Figures/PTA flore inversed behind the scene_codes4",format(Sys.Date(), "%d%m%Y"),".pdf"))
  plot(pta0)
  dev.off()
  
  #### ****************
  #### B - PTA 2 : Environnement <<<<<< #####
  #### ****************
  
  DataEnv<-read.csv("Derived_datasets/DataEnv_Costatis-Nut-09-2020.csv",header = T,sep= ";", dec=".") #DataEnv_Costatis-Nut-04-2020
  DataEnv<-DataEnv[DataEnv$Lieu.de.surveillance...Libellé %in% Sites$Nom,]
  
  DataEnv<-DataEnv[DataEnv$Passage...Année >=2008, ]
  
  #***:Changing the winter's reference year for example : Winter 2007 = Dec 2007 + Jan et Feb 2008 and not Dec 2007 + Jan et Feb 2007 ***
  DataEnv$Passage...Année[DataEnv$Passage...Mois %in% c(1,2)]<-DataEnv$Passage...Année[DataEnv$Passage...Mois %in% c(1,2)]-1
  #*********
  
  for (i in (1:dim(DataEnv)[1])){
    DataEnv$Station[i]<-as.character(Sites$Code[Sites$Nom == as.character(DataEnv$Lieu.de.surveillance...Libellé[i])])
  } #Creating a column with sites' code
  
  {
    # all parameters available : c("NH4","NO3+NO2","PO4","SALI","SIOH","TEMP","Précipitations..mm.","Durée.d.insolation...min.Vitesse..ms.","Débit..m3s.")
    DataAdj<-DataEnv[DataEnv$Résultat...Code.paramètre %in% c("NH4","NO3+NO2","PO4"),]
    DataNsP<-ddply(DataAdj,.(Saison, Station, Passage...Année, Lieu.de.surveillance...Libellé),function(x)data.frame(
      Mean = sum(x$Mean[x$Résultat...Code.paramètre %in% c("NH4","NO3+NO2")]) / x$Mean[x$Résultat...Code.paramètre %in% c("PO4")],
      Med = sum(x$Med[x$Résultat...Code.paramètre %in% c("NH4","NO3+NO2")]) / x$Med[x$Résultat...Code.paramètre %in% c("PO4")],
      Qone = sum(x$Qone[x$Résultat...Code.paramètre %in% c("NH4","NO3+NO2")]) / x$Qone[x$Résultat...Code.paramètre %in% c("PO4")],
      Qthree = sum(x$Qthree[x$Résultat...Code.paramètre %in% c("NH4","NO3+NO2")]) / x$Qthree[x$Résultat...Code.paramètre %in% c("PO4")] ))
    
    
    DataNsP<-tibble::add_column(DataNsP, Résultat...Code.paramètre = rep("N/P",dim(DataNsP)[1]), .after = "Passage...Année")  
    
    DataEnv<-rbind(DataEnv,DataNsP)
  }# Creation of the N/P parameter and selection of the abiotic variables  
  
  { 
    Opti_SW_Wind <- 270
    DataEnv$Mean[DataEnv$Résultat...Code.paramètre == "Direction..degres."]<- 
      abs(100 - abs(DataEnv$Mean[DataEnv$Résultat...Code.paramètre == "Direction..degres."] - Opti_SW_Wind) / 180 *100)
    
    DataEnv$Med[DataEnv$Résultat...Code.paramètre == "Direction..degres."]<- 
      abs(100 - abs(DataEnv$Med[DataEnv$Résultat...Code.paramètre == "Direction..degres."] - Opti_SW_Wind) / 180 *100)
    
    DataEnv$Qone[DataEnv$Résultat...Code.paramètre == "Direction..degres."]<- 
      abs(100 - abs(DataEnv$Qone[DataEnv$Résultat...Code.paramètre == "Direction..degres."] - Opti_SW_Wind) / 180 *100)
    
    DataEnv$Qthree[DataEnv$Résultat...Code.paramètre == "Direction..degres."]<- 
      abs(100 - abs(DataEnv$Qthree[DataEnv$Résultat...Code.paramètre == "Direction..degres."] - Opti_SW_Wind) / 180 *100)
    
  }#transformation of the wind direction values so that it is a similarity to an optimal wind direction
  
  
  ChosenVar<- c("N/P","SALI","TEMP","Précipitations..mm.","Durée.d.insolation...min.", "Débit..m3s.","Vitesse..ms.","Direction..degres.","OXYGENE","TURB")
  
  DataEnv<-DataEnv[DataEnv$Résultat...Code.paramètre %in% ChosenVar,]
  
  Missing_value_replacement<-mean(DataEnv$Mean[DataEnv$Saison=="Winter" & DataEnv$Station == "STCA" & DataEnv$Résultat...Code.paramètre == "N/P" ])
  #
  {  
    DataEnv$Mean<-as.numeric(as.character(DataEnv$Mean)) 
    
    TabsiteEnv<-DataEnv
    FirstY<- min(TabsiteEnv$Passage...Année)
    LastY<-max(TabsiteEnv$Passage...Année)
    
    require(plyr)
    
    Total <- TabsiteEnv
    
    
    #If season selection  ********************
    Total<-Total[Total$Saison %in% ChosenSeason,]
    ###***************************************
    
    ##Binding year and season
    Total$Passage...Année<-paste(Total$Passage...Année,Total$Saison,sep="-")
    
    #Total<-Total[-which(Total$Passage...Année %in% "2019-Winter"),]
    
    Total <- Total[,colnames(Total) %in% c("Passage...Année","Résultat...Code.paramètre","Mean","Station")]
    #
    
    Total$Résultat...Code.paramètre<-as.character(Total$Résultat...Code.paramètre)
    
    Total$Résultat...Code.paramètre[Total$Résultat...Code.paramètre %in% "Débit..m3s."] <- "Flow"
    Total$Résultat...Code.paramètre[Total$Résultat...Code.paramètre %in% "Durée.d.insolation...min."] <- "Daylight"
    Total$Résultat...Code.paramètre[Total$Résultat...Code.paramètre %in% "Précipitations..mm."] <- "Rainfall"
    Total$Résultat...Code.paramètre[Total$Résultat...Code.paramètre %in% "Vitesse..ms."] <- "Wind speed"
    Total$Résultat...Code.paramètre[Total$Résultat...Code.paramètre %in% "Direction..degres."] <- "West Wind"
    Total$Résultat...Code.paramètre[Total$Résultat...Code.paramètre %in% "OXYGENE"] <- "Oxygen"
    Total$Résultat...Code.paramètre[Total$Résultat...Code.paramètre %in% "SALI"] <- "Salinity"
    Total$Résultat...Code.paramètre[Total$Résultat...Code.paramètre %in% "TEMP"] <- "Water temperature"
    Total$Résultat...Code.paramètre[Total$Résultat...Code.paramètre %in% "TURB"] <- "Turbidity"
    
    library(tidyr)
    library(ade4)
    matrixenv<-spread(Total,Résultat...Code.paramètre,Mean,fill=NA) #then you realise the world is bad... only turbidity is complete
    
    
    #Binding dates and station names
    matrixenv<-matrixenv[order(matrixenv$Passage...Année),]
    matrixenv$group<-paste(matrixenv$Passage...Année,matrixenv$Station,sep="_")
    
    matrixenv[1:5,1:5] 
    #>>>>>>>>>>>>>< Pour les données spécial costatis <<<<<<<<<<<<<<<<<
    
    
    library(tidyr)
    library(ade4)
    
    Keepenv<- matrixenv[,c(1,2)]
    #matrixenv<-scale(matrixenv[,c(2:length(colnames(matrixenv)))]) #already taken care of within PTA
    
    matrixenv <- matrixenv[,-c(1,2)]
    rownames(matrixenv)<-matrixenv$group
    matrixenv<-matrixenv[,-dim(matrixenv)[2]]
    
    matrixenv$`N/P`[which(is.na(matrixenv$`N/P`))]<- Missing_value_replacement
    
    
    
    
    ###>>>>> PTA application <<<<<
    pcaenv <- withinpca (matrixenv, as.factor(as.character(Keepenv$Passage...Année)), scaling = "partial", scannf = F, nf = 3)
    ktapenv <- ktab.within (pcaenv,colnames = Keepenv$Station) # creating k-tables, one for every sites, so you have to put your variable of time in order (unique) then repeat for the number of sites. 
    
    kta41 <- t(ktapenv)
    pta01 <- pta (kta41, scan = F, 3)
    plot(pta01,option=c(1:4))
    
  } #Waiting for run
  # 
  
  pdf(file= paste("Figures/PTA inversed env behind the scene_codes4",format(Sys.Date(), "%d%m%Y"),".pdf"))
  plot(pta01)
  dev.off()
  
  
  #### ****************
  #######  Costatis #####
  #### ****************
  cost1 <- costatis(ktapenv, ktaflora, scannf = FALSE)
  plot(cost1)
  round(cost1$RV,3)
  #
}

#### ****************
### Slimani's plot production ####
#### ****************
cost1 <- costatis(ktapenv, ktaflora, scannf = FALSE)
row.names(cost1$c1)[which(row.names(cost1$c1) == "N.P")]<-"N/P"


#Plot limits
{xlt <- c(min(cost1$supIX[,1], cost1$c1[,1]*5, cost1$supIY[,1],
              cost1$l1[,1]*7), max(cost1$supIX[,1], cost1$c1[,1]*5,
                                   cost1$supIY[,1], cost1$l1[,1]*7))
  ylt <- c(min(cost1$supIX[,2], cost1$c1[,2]*5, cost1$supIY[,2],
               cost1$l1[,2]*7), max(cost1$supIX[,2], cost1$c1[,2]*5,
                                    cost1$supIY[,2], cost1$l1[,2]*7))
  lim1 <- c(min(xlt, ylt), max(xlt, ylt))
}

#Plot Station vs Environnments axis
{sl1 <- s.label(cost1$c1*5, xlim=lim1, ylim=lim1, label =
                  row.names(cost1$c1), plabels = list(cex = 1.25, col = "brown",
                                                      optim = T), ppoints.cex=0, plabels.boxes = list(draw = F), plot = FALSE)
  sa1 <- s.arrow(cost1$c1*4, xlim=lim1, ylim=lim1, label =
                   row.names(cost1$c1), plabels.cex = 0, plabels.boxes.draw = FALSE, pline.col = "brown",
                 psub.cex = 0, plines.lwd = 1, plot = FALSE)
  
  sc1 <- s.class(cost1$supIX, xlim=lim1, ylim=lim1,
                 fac = ktaflora$TC[,2], ellipseSize = 0, starSize = 0.7,
                 plabels = list(cex=1, col = "darkgreen"), ppoints.cex = .5, ppoints.alpha =0.3, 
                 plines.lwd = 0.5, plines.col= "grey", ppoints.col="grey", plot = FALSE)
  
  ss1 <- superpose(superpose(sa1,sc1, plot = FALSE), sl1,
                   plot = FALSE)
}

#Calculation of greatness with coordinates
Good<-sqrt((cost1$l1$RS1)*(cost1$l1$RS1) + (cost1$l1$RS2)*(cost1$l1$RS2))
Selec<-cbind(cost1$l1,rownames(cost1$l1))
Selec<-cbind(Selec,Good)
colnames(Selec)[c(3,4)]<-c("Numb","Greatness")
Signi<-cost1$eig/sum(cost1$eig)

#Selection of the species to plot
Selec<-Selec[order(Selec$Greatness,decreasing = T),]#[c(1:50),]

ChgNom<- read.table(file = "Derived_datasets/Species Info-05-2020.csv",header=T,sep= ";", dec=".")
Selec<-merge(Selec,ChgNom[,c(1,2)],by.x = c("Numb"), by.y = c("Esp"),all.y=F,all.x= T)
rownames(cost1$l1)<-Selec$Nomination

Classy<-Tot[,colnames(Tot) %in% c("Résultat...Nom.du.taxon", "Class.worms2019")]
Classy<-unique(Classy)
Selec <- merge (Selec, Classy,  by.x = c("Numb"), by.y = c("Résultat...Nom.du.taxon"),all.y=F,all.x= T)

Selec$Nomination<-as.character(Selec$Nomination)

# Plot Station vs Species axis
{sl2dia <- s.label(cost1$l1[which(rownames(cost1$l1) %in% Selec$Nomination[which(Selec$Class.worms2019 == "Bacillariophyceae")]),]*15, 
                   xlim=lim1+c(-1,1), ylim=lim1+c(-1,1), label =
                     row.names(cost1$l1[which(rownames(cost1$l1) %in% Selec$Nomination[which(Selec$Class.worms2019 == "Bacillariophyceae")]),]), 
                   plabels = list(cex = 0.8, col = "darkblue", optim = T), plabels.boxes = list(draw = F), ppoints.cex=0, plot = FALSE)
  
  sa2dia <- s.arrow(cost1$l1[which(rownames(cost1$l1) %in% Selec$Nomination[which(Selec$Class.worms2019 == "Bacillariophyceae")]),]*13, 
                    xlim=lim1+c(-1,1), ylim=lim1+c(-1,1), label =
                      row.names(cost1$l1[which(rownames(cost1$l1) %in% Selec$Nomination[which(Selec$Class.worms2019 == "Bacillariophyceae")]),]),
                    plabels.cex = 0, psub.cex = 0, plines.col = "darkblue", plines.lwd = 0.2, plines.lty=2, plot = FALSE)
  
  sc2dia <- s.class(cost1$supIY, xlim=lim1, ylim=lim1,
                    fac = ktaflora$TC[,2], ellipseSize = 0, starSize = 0.7,
                    plabels = list(cex=1, col = "darkgreen"), ppoints.cex = .5,ppoints.alpha =0.3,
                    plines.lwd = 0.2, plines.col= "grey", ppoints.col="grey",plot = FALSE)
  
  ss2dia <- superpose(superpose(sl2dia, sa2dia, plot = FALSE), sc2dia,
                      plot = FALSE)
  
} #Diatoms
{sl2oth <- s.label(cost1$l1[which(rownames(cost1$l1) %in% Selec$Nomination[which(Selec$Class.worms2019 != "Bacillariophyceae")]),]*15, 
                   xlim=lim1+c(-1,1), ylim=lim1+c(-1,1), label =
                     row.names(cost1$l1[which(rownames(cost1$l1) %in% Selec$Nomination[which(Selec$Class.worms2019 != "Bacillariophyceae")]),]), 
                   plabels = list(cex = 0.8, col = "darkblue", optim = T), plabels.boxes = list(draw = F), ppoints.cex=0, plot = FALSE)
  
  sa2oth <- s.arrow(cost1$l1[which(rownames(cost1$l1) %in% Selec$Nomination[which(Selec$Class.worms2019 != "Bacillariophyceae")]),]*13, 
                    xlim=lim1+c(-1,1), ylim=lim1+c(-1,1), label =
                      row.names(cost1$l1[which(rownames(cost1$l1) %in% Selec$Nomination[which(Selec$Class.worms2019 != "Bacillariophyceae")]),]),
                    plabels.cex = 0, psub.cex = 0, plines.col = "darkblue", plines.lwd = 0.2,plines.lty=2,  plot = FALSE)
  
  sc2oth <- s.class(cost1$supIY, xlim=lim1, ylim=lim1,
                    fac = ktaflora$TC[,2], ellipseSize = 0, starSize = 0.7,
                    plabels = list(cex=1, col = "darkgreen"), ppoints.cex = .5,ppoints.alpha =0.3,
                    plines.lwd = 0.2, plines.col= "grey", ppoints.col="grey",plot = FALSE)
  
  ss2oth <- superpose(superpose(sl2oth, sa2oth, plot = FALSE), sc2oth,
                      plot = FALSE)
} #Others

mlay <- matrix(c(0,1,1,0, 0,1,1,0, 2,2,3,3, 2,2,3,3), byrow = T, nrow = 4)
st1 <- ADEgS(list(ss1,ss2dia,ss2oth), layout = mlay)


png(file= paste("Figures/Costatis-Slimani_Figures",format(Sys.Date(), "%d%m%Y"),".png"),width = 210 , height = 150, units = "mm" , res = 300 , family = "Arial" , pointsize = 10 )
pdf(file= paste("Figures/Costatis-Slimani_Figures",format(Sys.Date(), "%d%m%Y"),".pdf"),width = 10 , height = 10)
st1
dev.off()

pdf(file= paste("Figures/Costatis-Mainplot-Slimani_style2270",format(Sys.Date(), "%d%m%Y"),".pdf"))
plot(cost1)
dev.off()


