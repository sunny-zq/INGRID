spls.cox <- function(x, y, K, eta, kappa=0.5, select="pls2", fit="regression", 
    scale.x=TRUE, scale.y=FALSE, eps=1e-04, maxstep=100) 
{   
    force(K)
    x <- as.matrix(x)
    n <- nrow(x)
    p <- ncol(x)
    ip <- c(1:p)
    y <- as.matrix(y)
    q <- ncol(y)
    one <- matrix(1, 1, n)
    
    mu <- one %*% y/n
    y <- scale(y, drop(mu), FALSE)
    meanx <- drop(one %*% x)/n
    x <- scale(x, meanx, FALSE)
    
    if (scale.x) {
        normx <- sqrt(drop(one %*% (x^2))/(n - 1))
        if (any(normx < .Machine$double.eps)) {
            stop("Some columns of the predictor matrix have 0 variance")
        }
        x <- scale(x, FALSE, normx)
    }else {normx <- rep(1, p)}
    
    if (scale.y) {
        normy <- sqrt(drop(one %*% (y^2))/(n - 1))
        if (any(normy < .Machine$double.eps)) {
            stop("Some columns of the response matrix have 0 variance")
        }
        y <- scale(y, FALSE, normy)
    }else {normy <- rep(1, q)}
    
    betahat <- matrix(0, p, q)
    
    ##objects for CV function
    blist <- wlist <- Alist <- list()
    
    x1 <- x
    y1 <- y
    
    type <- correctp.cox(x, y, eta, K, kappa, select, fit)
    eta <- type$eta
    K <- type$K
    kappa <- type$kappa
    select <- type$select
    fit <- type$fit
    
    if (is.null(colnames(x))) {
        xnames <- c(1:p)
    }else {xnames <- colnames(x)}
    
    new2As <- list()
           
    for (k in 1:K) {
        Z <- t(x1) %*% y1
        
        what <- spls.dv(Z, eta, kappa, eps, maxstep)
        A <- unique(ip[what != 0 | betahat[, 1] != 0])
        ##newly added variables
        new2A <- ip[what != 0 & betahat[, 1] == 0]
        xA <- x[, A, drop = FALSE]
        
        plsfit <- pls.cox(X=xA, Y=y, ncomp=min(k, length(A)), 
            mode=fit, scale.X=FALSE, scale.Y=FALSE)
		
		predplsfit <- predict.pls.cox(plsfit, newdata=xA, scale.X=F, 
            scale.Y=F)
            
		##update:
        betahat <- matrix(0, p, q)
        betahat[A, ] <- matrix(predplsfit$B.hat[, , plsfit$ncomp], 
            length(A), q)
        blist[[k]] <- betahat
        new2As[[k]] <- new2A
        Alist[[k]] <- A
        wlist[[k]] <- predplsfit$w
        
        if (select == "pls2") {
            y1 <- y - predplsfit$predict[, , plsfit$ncomp]
        }
    }
    
    if (!is.null(colnames(x))) {
        rownames(betahat) <- colnames(x)
    }
    if (q > 1 & !is.null(colnames(y))) {
        colnames(betahat) <- colnames(y)
    }
    
    object <- list(x=x, y=y, A=A, w=predplsfit$w, Alist=Alist,
    		wlist=wlist, blist=blist, betahat=betahat,
    		mu=mu, meanx=meanx, normx = normx, 
        	normy=normy, eta=eta, K=K, 
        	plsmod=plsfit)
    class(object) <- "spls"
    object
}

