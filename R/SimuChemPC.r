Normalize = function(X){
	 if(is.na(   dim(X)[2] )){
	 	mu= mean(X)	 
	 	norm= sweep(X,length(dim(X)) , mu, FUN="-")  #margin is dimention of X , we wanted to 
	 	sigma = sd(X, na.rm = FALSE)	
	 }else{
	 	 mu= colMeans(X) 	
	 	norm= sweep(X,length(dim(X)) , mu, FUN="-")  #margin is dimention of X , we wanted to 
	 	sigma = sapply(X,sd) #after R 2.14.0 only this works 
	 }	 
	 	 # use it for both trainX and trainY with different dimentions
	 sigma[sigma==0]=1	  
	 norm = sweep(norm, length(dim(X)), sigma, FUN="/" )
	 return ( norm) 	 
}

TanimotoCoef= function(X, Y){
	X= data.matrix(X)#R automatically converts matrixes to data frame, we have to convert them back
	Y=data.matrix	(Y)
	res = (X %*%t(Y ))/ ((X %*%t(X)) + (Y %*% t(Y) )- (X %*% t(Y)))
	#if(length(dim(x)) >1) res=res[1]#i added myself. some times more than  one feature is selected, and only one of them will be used in ONeNN
	# different fetures may lead us to different indexes
	return(res)
}

# 1NN based on Tanimoto Coefficient 
# X1 => vector / X2 => matrix
# returns index of 1NN row of X2 to X1
OneNN_Tc = function(X1,X2){
	ind=0
	tc=-1
    for (k in 1:dim(X2)[1]){
        K=TanimotoCoef(X1,X2[k,])
        #CheckTc(K, X1, X2[k,])
        	if ( K>tc){
            		tc=K
            	Ind=k
       	 }
    }
    index=Ind
    return( index )
}

SimuChemPC = function( dataFile, seedFile, method, repeatExperiment = 25)
{
	 if (is.null(dataFile) || is.null(seedFile) || is.null(method) ){
	 	  print("Error : Input parameteres can not be null.")
	 	  return
	 }
	simulation=0
	if (method== "random") simulation=1
	if (method== "EI") simulation=2
	if (method== "1NN") simulation=3
	if (method== "GP") simulation=4
	if(simulation ==0 ){
		 	  print("Error : function is called with wrong method parameter.")
	 	  return
	}	
	a<-read.table(seedFile)   # read data from file with tab delimiter
	Seed=a$V1 #convert the frame into list of integer values	
	DATA= read.table(dataFile)
	len=dim(DATA)[1]
	AttrCounter=rep(0,dim(DATA)[1]-2)
	TOTAL= repeatExperiment # number that the experiment repeats
	trainCompound=c()
		
	# 3D matrices
	D_3_train= c()        #array(1,TOTAL)# in matlab its a row of matrixes that their number is TOTAL
	D_3_test=  c()   #array(1,TOTAL)
	D_3_train_Normalized=  c()       # array(1,TOTAL)
	D_3_test_Normalized= c()    #array(1,TOTAL)

	# cell arrays(i.e. dynamic arrays with different dimensions)
	D_3_train_Normalized_afterFeatureSelection=array(1,TOTAL)
	D_3_test_Normalized_afterFeatureSelection=array(1,TOTAL)
		
 	pVals=c()
	hVals=array(1,TOTAL)
	adj_pVals=array(1,TOTAL)
	crit_pVals=array(1,TOTAL)
	pVals=array(1,TOTAL)
	no_feature_selected = 0		
  for(loop in 1:TOTAL) {
	    M=DATA
	    N=DATA [1:(len/2),1: dim(M)[2]]	            #array((len/2),dim(M)[2])
	    # do random partition (train & test data)
	    START=1
	    END=len
	    seed_i=1	
		for (i in 1:(len/2))       #take half of data  randomly and put in N as our train data
		    {
		        ran = round(START + (END-START) * Seed[seed_i])
		        seed_i= seed_i+1
		        N[i,]=M[ran,]
		        M = M[-ran,]
		        END=END-1
		       # print(ran)#to do :delete it
		  }
    	#============= Train and Test Data =====================================
    	end=dim(M)[2]
    	testX=M[,3:end-1]   # ignore first col (ID) as well as last one (potency)
    	test_ID = as.array(M[,1])         # extract first col (ID)
   	testY=as.array(log10(M[,end]) ) # extract last col (potency)

    	N_= N[,3:end-1] 	# ignore first col (ID) as well as last one (potency)
    	N__ = as.array(N[,1])         	 # extract first col (ID)
    	N___=as.array(log10(N[,end])) # extract last col (potency)
    	#===================== Get fractions of Train Data ====================
    	START=1
    	len_= dim(N_)[1]
    	END=len_
    	fraction=2/4
    	trainX=N_[1:(len_* fraction),1: dim(N_)[2]]	
    	trainID=N__[1: (dim(N__)*fraction)]
        	trainY=N___[1: (dim(N___)*fraction)]	
    	for (i in 1:(len_ * fraction))
    	{
    		ran = round(START + (END-START) * Seed[seed_i])
        		seed_i = seed_i+1
        		trainX[i,]=N_[ran,]
        		trainID[i]=N__[ran]
        		trainY[i]=N___[ran]
        		N_=N_[-ran,]
        		N__=N__[-ran]
        		N___=N___[-ran]
        		END=END-1
    	}
   #===================== CHECK =======================================
    D_3_train=c(D_3_train,list(trainX,trainY))#index is(loop*2)-1   like:  trainX=D_3_train[[(loop*2)-1]]
    D_3_test=c(D_3_test,list(testX,testY))
   #============= Normalization =====================================
    
    # normalization train data
    x =Normalize(trainX)
    sigmaX =  sapply(trainX,sd)
    sigmaX[sigmaX==0]=1	 
    muX = colMeans(trainX)    	
    y =Normalize(trainY)
    sigmaY = sd(trainY, na.rm = FALSE)
    sigmaY[sigmaY==0]=1	 
    muY = mean(trainY)   
    # normalization test data / use the same (mu & std) as train data
      testX_norm = sweep( testX,2, muX, FUN = "-")
      xstar = sweep( testX_norm,2 , sigmaX, FUN = "/") 
      	    
    testY_norm = sweep( testY, 1, muY, FUN = "-")
    MU2 = sweep( testY_norm, 1, sigmaY, FUN = "/")
    
    #===================== CHECK =======================================
    D_3_train_Normalized = c(D_3_train_Normalized,list(x,y))#loop*2 is the index for x    
    D_3_test_Normalized = c(D_3_test_Normalized,list(xstar,MU2))

    #==================== Feature Selection ==============================
    #apply multiple testing correction for feature selection
    op <- options(warn = (-1))
    rho =  cor(x, y,   method = "spearman") 	
    suppressWarnings(warning("Cannot compute exact p-value with ties"))
    warnings= NULL 	
    pVals_loop=x[1,]
    for(j in 1:dim(x)[2]){
    	pVals_loop[j]=  cor.test(x[,j], y,   method = "spearman")[3]$p.value
    }
    pVals=c(pVals,pVals_loop)
   t_adj_pVals= p.adjust(pVals_loop ,  method = "fdr") 	#hVals AND crit_pVals ARE NOT produced  [hVals{loop} crit_pVals{loop} adj_pVals{loop}]=fdr_bh(pVals{loop});      there is  MISTAKE IN MATLAB FUNCTION fdr_bh; IT DOES NOT INCLUDE NEGATIVE P-VALUES; 
   adj_pVals=c(adj_pVals, t_adj_pVals)# Adjusted values are not same with matlab 
   t_hVals = abs(t_adj_pVals)< .05 
  if (length( grep(TRUE,t_hVals))>0)
  {
  	  B= table(t_hVals, 1:length(t_hVals))
  	  Attr= as.integer(names(B[2,B[2,]==1])) # find them 
  	  if (length(Attr)==0)  Attr=grep(TRUE,t_hVals, fixed = TRUE)   #when there is only one answer
  	  hVals =	c(hVals, t_hVals)
  }else{	  
	# In case no feature is found to be significant, we select those features for 
	#which adjusted p-value is minimized
	 min_  = min(abs(t_adj_pVals),na.rm=T)
	 t_hVals = abs(t_adj_pVals)==	min_
	  B= table(t_hVals, 1:length(t_hVals))
  	  Attr= as.integer(names(B[2,B[2,]==1])) # find them 	
           if (length(Attr)==0)  Attr=grep(TRUE,t_hVals, fixed = TRUE)   
            no_feature_selected = no_feature_selected+1
  }
    
    AttrCounter[Attr]=AttrCounter[Attr]+1
    x=x[,Attr]  
    xstar=xstar[ , Attr]  
    
    D_3_train_Normalized_afterFeatureSelection = c(D_3_train_Normalized_afterFeatureSelection,x)
    D_3_test_Normalized_afterFeatureSelection = c( D_3_test_Normalized_afterFeatureSelection,xstar)	
    #==================== Remove duplication for test data =============
     print("Removing duplication for test data.")
    L2=duplicated(xstar)
    xstar = unique(xstar)  # and duplicated() returns index of those who are duplicated 
    p=1
    mu2=rep(1,length(test_ID)-length(L2[L2==TRUE]))
    testID=rep(1,length(test_ID)-length(L2[L2==TRUE]))
    for( i in 1:length(test_ID)){
	        if (L2[i]==FALSE){
	            mu2[p]=MU2[i]
	            testID[p]=test_ID[i]
	            p=p+1
	        }
     }
    #================= Initialization =================================
    print("Initialization.")
    trainCompound = c( trainCompound ,trainID)
    if( loop==1){
	        c.init = dim(test_ID)
	        TestCompound=rep(-1000,c.init *TOTAL)  #TestCompound / fill MATRIX with -1000
	        TestCompound = array(TestCompound , c(c.init ,TOTAL))
	        PotencyReal=rep(-1000,c.init * TOTAL) # Potency (Real) / fill MATRIX with -1000
	        PotencyReal =array(PotencyReal, c(c.init,TOTAL))
	}
    #================= Run simulation =================================
    print(paste("Running the method: ", method,".",sep=""))
     counter = 1 
    while (TRUE){    			
		        d1=dim(xstar)[1]
		        		    if (d1==0){
			            break
			    }else if(d1== 1 ){		    		
			    	xstar= as.vector(xstar)
			    }
		        	print(paste("Left number of rows in the test data: ",d1,".",sep=""))
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
		                		index=OneNN_Tc(x[row,], xstar)   #x is trian data, xstar is test data
				}else { #EI & GP
				                #covfunc =c( "covSum","covSEiso","covNoise")
					      covfunc="covSum,covSEiso,covNoise"  	
					       loghyper= array(c(-1,-1,-1), dim=c(3,1))
				                loghyper = minimize(loghyper, 'gpr', -100, covfunc, x, y)
				               loghyper = loghyper[[1]]
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
		           #========================================================================
			  PotencyReal [counter,loop] = mu2[index]
			  TestCompound[counter,loop] = testID[index]
			   print(paste("COUNTER -->",counter , sep = " "))
			   print(paste("Experiment number -->", loop, sep = " "))
			   counter = counter + 1
			   x= rbind(x , xstar[index,])  
			    yy= as.list(y) # this not nnot work y=rbind (y, mu2[index])
			    yy[length(yy)+1]=mu2[index]
			    y= (array(unlist(yy)))			   
			   xstar=xstar[-index,]
			   mu2=  mu2[-index]
			   testID=   testID[-index]			      		
		}#end while	    
	 }#end for loop
	 print(paste("storing the method result in file: " ,method, "_Feature Selection_Random.RData", sep=""))
	 dump(ls(),paste(method, "_Feature Selection_Random.RData", sep=""))	 	
	
}
