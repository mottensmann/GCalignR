RelativeAbundance = function(Chromatograms){
  # Computes the relative concentration of a substance
  # based on the chromatogram with merged rows and substracted
  # controls, if it was applied
  SumArea <- sum(Chromatograms$Area)
  SumHeight <- sum(Chromatograms$Height)
  for (i in 1:dim(Chromatograms)[1]){
    Chromatograms$AreaRel[i] <- (Chromatograms$Area[i]/SumArea)*100
    Chromatograms$HeightRel[i] <- (Chromatograms$Height[i]/SumHeight)*100
  }
  Chromatograms
}