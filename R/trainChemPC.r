#train function for EI or GP
trainChemPC = function(trainData, targetVector){

	x = trainData
	y = targetVector
	if (is.null(trainData) ||  is.null(targetVector)){
	 	  print("Error : Input parameteres can not be null.")
	 	  return
	 }
		 
		  covfunc="covSum,covSEiso,covNoise"  	
		  loghyper= array(c(-1,-1,-1), dim=c(3,1))
		  loghyper = minimize(loghyper, 'gpr', -100, covfunc, x, y)
		  loghyper = loghyper[[1]]				     				
  
	return (loghyper)
}

