#' Poisson SVD with Offset estimation
#' @param Y:	data matrix
#' @param B:  the first K columns of the U from SVD(log(Y/T))
#' @param F:  the first K columns of the SV from SVD(log(Y/T))
#' @param K:  # of SVD components
#' @param  niter:  	# of maximum iterations for poissonSVD
#' @param  err:   	error bound on claiming convergence
#' @param  ln:    	link functions, "log", "sqrt", "identity"
#' @param  verbose: 	Log information
#' @return Returns PSVD results
#' 	a <- PSVDOS(Y,K=2)
#'
#' @references
#' [1] Lee, S., Chugh, P. E., Shen, H., Eberle, R., & Dittmer, D. P. (2013) Poisson factor models with applications to non-normalized microRNA profiling. Bioinformatics, 29(9), 1105-1111

PSVDOS <- function(Y,K=NULL, B = NULL, F = NULL, E = NULL,
		niter = 100, err = 0.0001,ln="log",verbose = 0,const=1,zerocount=exp(-3),
		rowoffset = NULL,colcenter=0,...)
{#, coloffset=NULL,...){

	n = dim(Y)[1]
	m = dim(Y)[2]
	rowoffsetupdate = 0 #No update for offset
	#	coloffsetupdate = 0
	if(is.null(K)){K=min(n,m)}
	if(is.null(B) | is.null(F)){
		Ytmp = Y;Ytmp[Y==0] = zerocount
		tmp = log(Ytmp)
		temp = svd(tmp-matrix(rep(apply(tmp,1,mean),m),n,m))#t(matrix(rep(apply(tmp,2,mean),n),m,n)))#+matrix(mean(tmp),n,m))

		if(is.null(rowoffset)){## the initial values of row offset value needs to be estimated.

			if(is.null(B)){
				B = matrix(temp$u[,1:K], nrow=n)
			}
			if(is.null(F)){
				tmpv = scale(temp$v[,1:K],scale=FALSE)
				F =tmpv%*%diag(temp$d[1:K], nrow=K)
			}
#			if(is.null(rowoffset)){
			rowoffset = rep(0,n)
			rowoffsetupdate = 1 #Estimate Offset vector iteratively.
			for(i in 1:n){
				temp=NULL
				temp <- glm(Y[i,]~F, family=poisson(link=ln))
				rowoffset[i] = exp(temp$coefficients[1])
				B[i,] = temp$coefficients[-1]
			}

#			}
			mu = mean(log(rowoffset))
#			print(mu)
			for(j in 1:m){
				temp=NULL
				temp <- glm(Y[,j]~-1+B, family=poisson(link=ln), offset=log(rowoffset))
				F[j,] = temp$coefficients#[-1]
			}
		}
		offtemp = matrix(rep(log(rowoffset),m),n,m)# + t(matrix(rep(log(coloffset),n),m,n))
		tmp = log(Ytmp)-offtemp
		temp = svd(tmp)
		B = matrix(temp$u[, 1:K], nrow=n)
		tmpv = matrix(scale(temp$v[,1:K],scale=FALSE, center=TRUE),m,K)
		F =tmpv%*%diag(temp$d[1:K], nrow=K)
		if(is.null(E)){
			E =tmpv
		}
	}

  	Bp = B;Fp = F;Ep = E;BF.d = 1;iter = 0;
	if (verbose==1) {print(paste("K = ",K, sep=""))}

	while(BF.d > err){
    		iter <- iter+1
		if (verbose==1) {print(c(iter, BF.d))}

		for(i in 1:n){
			#if (verbose==1) {print(paste("First Step, i=", i,sep=""))}
			temp=NULL
			temp<-glm(Y[i,]~Fp, family=poisson(link=ln))

			if(rowoffsetupdate==1){
				rowoffset[i] = exp(temp$coefficients[1])
			}
			Bp[i,] <- temp$coefficients[-1]
			#print(Bp[i,])
		}
		mu = mean(log(rowoffset))

		for(j in 1:m){
			#if (verbose==1) {print(paste("Second Step, j=", j,sep=""))}
			temp=NULL
			temp<-glm(Y[,j]~-1+Bp, family=poisson(link=ln), offset =log(rowoffset))
			Fp[j,] <- temp$coefficients
		}
		tp = Bp%*%t(Fp)
    		temp <- svd(tp-matrix(rep(apply(tp,1,mean),m),n,m))
		Bp <- matrix(temp$u[,1:K],nrow=n)
		Dp = temp$d
		Ep =  matrix(temp$v[,1:K], nrow=m)
		tmpv = matrix(scale(temp$v[,1:K],scale=FALSE),m,K)
		Fp =tmpv%*%diag(temp$d[1:K], nrow=K)
    		BF.d <- max(sqrt(sum((B-Bp)^2)), sqrt(sum((E-Ep)^2)))# sqrt(sum((F-Fp)^2)))
    		B <- Bp
    		F <- Fp
		E = Ep

    		if(iter > niter) {
      			print("Fail to converge! Increase the niter!")
      			break
		}
  	}

  	rowoffset = exp(log(rowoffset)-mu)
  	temp <- B%*%t(F)
  	if(ln=="log"){
		dev <- 2*sum(exp(temp)-Y*temp)
  	}

  	cat(iter, "\n")
  	return(list(B=B, F=F, E=E, D = Dp,iter=iter, BF.d=BF.d, deviance=dev, rowoffset = rowoffset, mu=mu))
}

