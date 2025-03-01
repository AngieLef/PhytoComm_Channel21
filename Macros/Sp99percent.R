# Nom               : Sp99percent
# Type              : 
# Objet             : Refonte du jeu de donn�e avec selection uniquement des esp�ces composant plus de 1% des abondances des Sites selectionn�s/toutes profondeurs consid�r�es
# Input             : 
# Output            : 
# Auteur            : AL
# R version         : 3.5.1, (ggplot2 : V3.1.0)
# Date de cr�ation  : 23 JUILLET 2019
#_______________________________________________________________________________
#DataIn<-Dataset

Sp99<- function(Data,Sites) {
  
  DataIn<-Data
  Wanted <-NULL
  Unwanted <-NULL
  DataIn <- DataIn[DataIn$Lieu.de.surveillance...Mn�monique %in% Sites$Num,]
  DataIn<-DataIn[DataIn$R�sultat...Valeur.de.la.mesure > 100,]
  
  if("Code" %in% colnames(DataIn)) {

    for (i in (1:length(unique(DataIn$Code)))) {
      
      MoyenTaxon <- mean(DataIn$R�sultat...Valeur.de.la.mesure[DataIn$Code == unique(DataIn$Code)[i]])
      MoyenTotal <- mean(DataIn$R�sultat...Valeur.de.la.mesure)
      
      if(MoyenTaxon/MoyenTotal > 0.01 ) { Wanted <- c(Wanted,as.character(unique(DataIn$Code)[i])) 
      } else { 
        Unwanted <- c(Unwanted,as.character(unique(DataIn$Code)[i]))
        }
    }
    DataOut<-Data[Data$Code %in% Wanted,]
  
  } else {
    
    for (i in (1:length(unique(DataIn$R�sultat...Nom.du.taxon)))) {
      
      MoyenTaxon <- mean(DataIn$R�sultat...Valeur.de.la.mesure[DataIn$R�sultat...Nom.du.taxon == unique(DataIn$R�sultat...Nom.du.taxon)[i]])
      MoyenTotal <- mean(DataIn$R�sultat...Valeur.de.la.mesure)
      
      if(MoyenTaxon/MoyenTotal > 0.01 ) { Wanted <- c(Wanted,as.character(unique(DataIn$R�sultat...Nom.du.taxon)[i])) 
      } else { 
        Unwanted <- c(Unwanted,as.character(unique(DataIn$R�sultat...Nom.du.taxon)[i]))
      }
    }
    DataOut<-Data[Data$R�sultat...Nom.du.taxon %in% Wanted,]
  }
  
  
  
  warning(paste(length(Unwanted),"were deleted from the data set as they represented less than 1% of the total abundance among them :"))
  warning(paste(Unwanted, " "))
  
  return(DataOut)
  
}


