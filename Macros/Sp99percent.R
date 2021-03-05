# Nom               : Sp99percent
# Type              : 
# Objet             : Refonte du jeu de donnée avec selection uniquement des espèces composant plus de 1% des abondances des Sites selectionnés/toutes profondeurs considérées
# Input             : 
# Output            : 
# Auteur            : AL
# R version         : 3.5.1, (ggplot2 : V3.1.0)
# Date de création  : 23 JUILLET 2019
#_______________________________________________________________________________
#DataIn<-Dataset

Sp99<- function(Data,Sites) {
  
  DataIn<-Data
  Wanted <-NULL
  Unwanted <-NULL
  DataIn <- DataIn[DataIn$Lieu.de.surveillance...Mnémonique %in% Sites$Num,]
  DataIn<-DataIn[DataIn$Résultat...Valeur.de.la.mesure > 100,]
  
  if("Code" %in% colnames(DataIn)) {

    for (i in (1:length(unique(DataIn$Code)))) {
      
      MoyenTaxon <- mean(DataIn$Résultat...Valeur.de.la.mesure[DataIn$Code == unique(DataIn$Code)[i]])
      MoyenTotal <- mean(DataIn$Résultat...Valeur.de.la.mesure)
      
      if(MoyenTaxon/MoyenTotal > 0.01 ) { Wanted <- c(Wanted,as.character(unique(DataIn$Code)[i])) 
      } else { 
        Unwanted <- c(Unwanted,as.character(unique(DataIn$Code)[i]))
        }
    }
    DataOut<-Data[Data$Code %in% Wanted,]
  
  } else {
    
    for (i in (1:length(unique(DataIn$Résultat...Nom.du.taxon)))) {
      
      MoyenTaxon <- mean(DataIn$Résultat...Valeur.de.la.mesure[DataIn$Résultat...Nom.du.taxon == unique(DataIn$Résultat...Nom.du.taxon)[i]])
      MoyenTotal <- mean(DataIn$Résultat...Valeur.de.la.mesure)
      
      if(MoyenTaxon/MoyenTotal > 0.01 ) { Wanted <- c(Wanted,as.character(unique(DataIn$Résultat...Nom.du.taxon)[i])) 
      } else { 
        Unwanted <- c(Unwanted,as.character(unique(DataIn$Résultat...Nom.du.taxon)[i]))
      }
    }
    DataOut<-Data[Data$Résultat...Nom.du.taxon %in% Wanted,]
  }
  
  
  
  warning(paste(length(Unwanted),"were deleted from the data set as they represented less than 1% of the total abundance among them :"))
  warning(paste(Unwanted, " "))
  
  return(DataOut)
  
}


