###cross validation by fixing eta
###x: features
###K: maximum number of latent variables allowed
cv.eta=function( x, t, d, foldid, K, eta, method ){

	n=length(t)
	nfold=max(foldid)
	cvraw=matrix(NA,nrow=nfold,ncol=K)

	for(i in 1:nfold){
		##training
   		o=which(foldid==i)
   		cox=coxph(Surv(t[-o],d[-o])~1)
     	res=residuals(cox,type="deviance")
     	mod=spls.cox(x=x[-o,],y=res,K=K,eta=eta,
       		kappa=0.5,scale.x=T,scale.y=F)

		for(k in 1:K){
			beta=NULL
			w=mod$wlist[[k]]
			A=mod$Alist[[k]]

			x_minus_i=scale(x[-o,A],mod$meanx[A],mod$normx[A])
			s_minus_i=x_minus_i%*%w
			fit=coxph(Surv(t[-o],d[-o])~s_minus_i)
		    beta=w%*%summary(fit)$coef[,1]

		    x_i=scale(x[o,A],mod$meanx[A],mod$normx[A])
		    s_i=x_i%*%w

		    xfull=scale(x[,A],mod$meanx[A],mod$normx[A])
			sfull=xfull%*%w

            if(method=="auc"){
            	xlp <- rep(NA, n)
       			xlp[-o] <- x_minus_i %*% beta
            	xlp[o] <- x_i %*% beta

            	AUCs <- getIndicCViAUCSurvROCTest(xlp[-o],xlp[o],
                  Surv.rsp=Surv(t[-o],d[-o]),
                  Surv.rsp.new=Surv(t[o],d[o]),
                  times.auc=seq(0,max(t),length.out=1000),
                  times.prederr=seq(0,max(t),length.out=1000)[-(990:1000)],
                  fit, plot.it=F)

               cvraw[i,k]=AUCs$AUC_survivalROC.test$iauc
            }

            if(method=="plik"){
            	pl_minus_i=-2*logplik(x=x_minus_i,time=t[-o],
            		status=d[-o],b=beta, return.all=F )

   			    plfull=-2*logplik(x=xfull,time=t,status=d,
   			    	b=beta, return.all=F )

     			cvraw[i,k]=plfull-pl_minus_i
            }
        }
	}

    if(method=="auc"){
    	cvm=apply(cvraw,2,mean,na.rm=TRUE)
		cvsd=sqrt(apply(cvraw,2,var,na.rm=TRUE))/nfold
	}

	if(method=="plik"){
		weights=rep(1,n)
		weights=as.vector( tapply(weights*d,foldid,sum) )
		cvraw=cvraw/weights
    	cvm=apply(cvraw,2,weighted.mean,w=weights)
    	cvsd=sqrt( apply( scale(cvraw,cvm,F)^2, 2, weighted.mean,
			 w=weights )/(nfold-1) )
	}

	result=cbind( cvm, cvsd, 1:K, rep(eta,K) )
	result
 }

#################################################
##-----------------------------------------------
##-----------------------------------------------
cv_splscox=function(x, t, d, foldid, K, eta.vec, method, parallel){

	if(parallel==T){
		cvmat=foreach(i=1:length(eta.vec),.combine='rbind')%dopar%{
			cv.eta(x=x,t=t,d=d,foldid=foldid,K=K,
			eta=eta.vec[i],method=method)}
	}

	if(parallel==F){
		cvmat=foreach(i=1:length(eta.vec),.combine='rbind')%do%{
			cv.eta(x=x,t=t,d=d,foldid=foldid,K=K,
			eta=eta.vec[i],method=method)}
	}

	opt.k=opt.eta=c(NA,NA)

	if(method=="auc"){
		idx=which.max(cvmat[,1])

		opt.k[1]=cvmat[idx,3]
 	 	opt.eta[1]=cvmat[idx,4]

 	 	se1=cvmat[idx,1]-cvmat[idx,2]
 	 	mat=cvmat[cvmat[,1]>se1,,drop=F]
 	 	##choose the most parsimonious model in terms of genes
   		##always choose the smallest k among those with largest eta
 	 	q=which.max(mat[,4])
		opt.k[2]=mat[q,3]
 	 	opt.eta[2]=mat[q,4]
 	 }

 	if(method=="plik"){
		idx=which.min(cvmat[,1])

		opt.k[1]=cvmat[idx,3]
 	 	opt.eta[1]=cvmat[idx,4]

 	 	se1=cvmat[idx,1]+cvmat[idx,2]
 	 	mat=cvmat[cvmat[,1]<se1,,drop=F]
 	 	##choose the most parsimonious model in terms of genes
   		##always choose the smallest k among those with largest eta
 	 	q=which.max(mat[,4])
		opt.k[2]=mat[q,3]
 	 	opt.eta[2]=mat[q,4]
 	 }
	list(cvmat=cvmat,opt.k=opt.k,opt.eta=opt.eta)
}
