#predict function
predictChemPC = function( trainData, testData, targetData, method, loghyper = NULL){


	x = trainData
	xstar =testData
	y = targetData
	if (is.null(trainData) ||  is.null(testData) || is.null(method) || is.null(targetData) ){
	 	  print("Error : Input parameteres can not be null.")
	 	  return (NULL)
	 }
	simulation=0
	if (method== "random") simulation=1
	if (method== "EI") simulation=2
	if (method== "1NN") simulation=3
	if (method== "GP") simulation=4
	if(simulation ==0 ){

	 	  print("Error : function is called with wrong method Type parameter.")
	 	  return (NULL)
	}	
	if ((simulation==2 || simulation==4) && is.null(loghyper) ){

		  print("Error : function is called with null loghyper parameter.")
	 	  return (NULL)
	}

	if(is.na(dim(x)[2])){
		  print("Error : trainData must have more than one row.")
	 	  return (NULL)
	} 
 			
		        d1=dim(xstar)[1]
		            if (d1==0){
				print(paste("Error: Number of rows in the test data is ",d1,".",sep=""))
			        return (NULL)
			    }else if(d1== 1 ){		    		
			    	xstar= as.vector(xstar)
			    }
		        	#print(paste("Number of rows in the test data is ",d1,".",sep=""))
		        if (simulation==1){  # RA
		            index = round(runif(1)*d1) 
		            if(index== 0) index =1  #special case, so that mu[0] would not be meaning less, when d1 is 1 	
		        }else{ #  EI & NN & GP
				if (simulation==3) {# --------  1NN Selection (based on TanimotoCoef) ----------------
					 max_  = max(y ,na.rm=T)
			 		t_ =   y==	max_
			  		B= table(t_, 1:length(t_))   
			  		# B[names(B)==TRUE] 
			  	          row = as.integer(names(B[2,B[2,]==1])) 
					if (length(row)==0)  row=grep(TRUE,t_, fixed = TRUE)    # or row=match(0,B[1,])
		                		index=OneNN_Tc(x[row,], xstar)   #x is train data, xstar is test data
				}else { #EI & GP
				                #covfunc =c( "covSum","covSEiso","covNoise")
					      covfunc="covSum,covSEiso,covNoise"  						   
				               f.out= gpr(loghyper, covfunc, x, y, xstar , TRUE)
				               mu = f.out[1][[1]]
				               S2= f.out[2 ][[ 1 ]]
		              		     if (simulation==4 ){
		                   	  	  index=grep(mu [mu==max(mu)][1] ,mu, TRUE)[1]
		             		    }else{ # EI
				                    d1= dim(mu)[1]
				                    d2= dim(mu)[2]
				                    Qmax = max(y)
				                    ei = rep(1, d1)
				                    for (j in 1:d1){
					                        m = mu[j]
					                        var = S2[j]
					                        # max{sigma^2(x), squared difference of log-potency of nearest neighbor training compound to predicted log-potency} / Tc
					                        index_=OneNN_Tc(xstar[j,],x)
					                        s=sqrt(max( var, (y[ index_] - m ) ^ 2 ))
					                        #s=sqrt(max([var (y(index_)-m)^2]));
					                        u = m - Qmax
					                        ei[j] = (u * pnorm(u/s,0,1)) + (s * dnorm(u/s,0,1))
				                     }
		   	                		  #------------------------
			                    	  index = grep(ei[ei==max(ei)] [1], ei, TRUE)[1]#  %EI Selection
			     		  }#end if simulation ==4
				}#end if simulation ==3
    			  }#end if simulation ==1

	
return (index)
}


