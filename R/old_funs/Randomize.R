RandomizeOrder = function(Samples){
  # Shuffle the order of Chromatograms randomly
  N <- length(Samples) # How many samples are there
  RandOrder <- sample(1:N)
  ObjectNames <- character(0)
  Names <- names(Samples)
  NewOrder <- list()
  for (i in 1:N){
    NewOrder[i] <- (Samples[RandOrder[i]])
    ObjectNames[i] <- Names[RandOrder[i]]
  }
  names(NewOrder) <- ObjectNames
  NewOrder
}
