pred_SOAP_step = function(betalist,ylist, tlist,spline_basis,nminus=1,sigma=0.1,weights = rep(1000,length(betalist)))
{
cfits = lapply(1:length(ylist), function(x){
	timei = tlist[[x]]
	xmat = lapply(1:length(betalist), function(i){
		pc_fit  = fd(betalist[[i]], spline_basis)
		eval.fd(timei, pc_fit)%>%as.numeric
	})%>%do.call(cbind,.)
	index_pc = length(betalist)
	yi = ylist[[x]]
	cfit = lm(yi~0+xmat)%>%coef
	if(length(cfit)<length(betalist)) cfit = c(cfit, rep(0,length(betalist)-length(cfit)))
	cfit
})%>%do.call(rbind,.)
yfits = lapply(1:nrow(cfits), function(i){
	cfit = cfits[i,]
	rowSums(mapply("*",cfit, betalist,SIMPLIFY=TRUE))
}
)%>%do.call(cbind,.)
yfitsfd = fd(yfits, spline_basis)
residuals = sapply(1:length(ylist), function(x){
	timei = tlist[[x]]
	resid = ylist[[x]] - as.numeric(eval.fd(timei,yfitsfd[x]))
	if (all(resid==0)) {
		return(NA)
	} else {
	return(mean(resid^2))
	}
})
sigmahat = mean(residuals[!is.na(residuals)])
residuals2 = lapply(1:length(ylist), function(x){
	timei = tlist[[x]]
	resid = ylist[[x]] - as.numeric(eval.fd(timei,yfitsfd[x]))
	resid
})
list(sigmahat = sigmahat, yfd_fit=  yfitsfd,sfit = cfits,residuals=residuals2)
}



select_gamma = function(previous_beta, gamma_pool, folds, observed, timepoints,spline_basis)
{
cv_res = lapply(gamma_pool, function(gam){
fit_error = c()
for (j in 1:length(folds)){
	test = folds[[j]]
	train = setdiff(1:length(observed),test)
pc_fit = third_FPC_conditional(runif(spline_basis$nbasis,-2,2), observed = observed[train], timepoints=timepoints[train],pc_index=length(previous_beta)+2, gamma = gam, betalist =previous_beta,threshold=1e-3)
betalist=  list()
betalist = previous_beta
betalist[[length(previous_beta)+1]]= pc_fit$beta
res = pred_SOAP(betalist,observed[test], timepoints[test], spline_basis)
res$sigmahat
fit_error  =c(fit_error,res$sigmahat)
}
fit_error 
})%>%do.call(rbind,.)
print(cv_res)
(gamma_select = gamma_pool[which.min(cv_res%>%rowMeans)])

pc_fit = third_FPC_conditional(runif(spline_basis$nbasis,-2,2), observed = observed, timepoints=timepoints,pc_index=length(previous_beta)+2, gamma = gamma_select, betalist =previous_beta,threshold=1e-3)
return(list(pc_fit=pc_fit, cv_res = cv_res))

}


first_FPC_soap  = function(beta1,observed,timepoints,threshold=1e-3,minit=3,gamma=0){
thresh = 1
it = 1
beta1_before = rep(0, length(beta1))
value = -1e14
	pc_fit = fd(beta1, spline_basis)
	##print(pc_fit)
	pc_fit = (1/sqrt(inprod(pc_fit, pc_fit))*pc_fit)
	#print(pc_fit)
	beta1 = coef(pc_fit)%>%as.numeric
	R = inprod(spline_basis, spline_basis,2,2)
	while (thresh >  threshold|it<minit){
	
	beta1_before = beta1
	value_before = value
alpha_fit = function(subj_index){
	(timei =timepoints[[subj_index]])
	xmati = eval.fd(timei, pc_fit)
	lm(observed[[subj_index]]~0+xmati)%>%coef%>%as.numeric
}

sfit= sapply(1:length(observed),alpha_fit)

residual_fit = function(subj_index){
	(timei =timepoints[[subj_index]])
	xmati = eval.fd(timei, pc_fit)
	lm(observed[[subj_index]]~xmati+0)%>%residuals
}

rfit0= lapply(1:length(observed),residual_fit)
rfit = lapply(rfit0, function(x) {x^2%>%mean})%>%do.call(c,.)
value = mean(rfit^2)+gamma*inprod(pc_fit,pc_fit,2,2)
# print(value)
if(abs(value_before - value/value_before) < threshold) break

yem = do.call(c,observed)


xalpha = lapply(1:length(observed),function(subj_index){

	(timei =timepoints[[subj_index]])
	xmati = eval.basis(timei, spline_basis)
	xmati*sfit[subj_index]

})%>%do.call(rbind,.)
yem[is.na(yem)] <- 0
beta1 = solve(t(xalpha%>%as.matrix)%*%(as.matrix(xalpha)) + gamma*R, t(xalpha%>%as.matrix)%*%as.matrix(yem))
pc_fit = fd(beta1, spline_basis)
pc_fit = (1/sqrt(inprod(pc_fit, pc_fit))*pc_fit)
# plot(pc_fit)
beta1 = coef(pc_fit)%>%as.numeric
thresh =  max(abs(beta1_before-beta1))
it = it+1
if(it%%30==0) {print(it);print(as.numeric(thresh))}
}
pc_fit = fd(beta1, spline_basis)
pc_fit = (1/sqrt(inprod(pc_fit, pc_fit))*pc_fit)
beta1 = coef(pc_fit)%>%as.numeric
rfit= lapply(1:length(observed),residual_fit)%>%do.call(c,.)
print("Done!")
cat(sprintf('The threshold is %s . \n',thresh))
return(list(beta=beta1,pc_fit = pc_fit,sfit = sfit,thresh = thresh,it=it,value=value,residual =rfit0 , gamma=gamma))
}




third_FPC_conditional_soap  = function(beta3, pc_index, observed, timepoints, betalist =previous_beta , threshold=1e-4,minit=1,gamma=0){
if(missing(pc_index)) pc_index= length(betalist)+3
thresh = 1
it = 1
R = inprod(spline_basis, spline_basis,2,2)
E = inprod(spline_basis, spline_basis,0,0)

pc_fit = fd(beta3, spline_basis)
pc_fit = (1/sqrt(inprod(pc_fit, pc_fit))*pc_fit)
beta3 = coef(pc_fit)%>%as.numeric

pc_fits_previous = lapply(1:length(betalist), function(x){
pc_fit1 = fd(betalist[[x]], spline_basis)
pc_fit1 = (1/sqrt(inprod(pc_fit1, pc_fit1))*pc_fit1)
pc_fit1
})
value = -1e14

observed2 =observed
observed2[which(sapply(observed,length)<=pc_index)]=NULL
timepoints2 = timepoints
timepoints2[which(sapply(observed,length)<=pc_index)]=NULL

beta3_before = rep(0, length(betalist[[1]]))
while (thresh >  threshold|it<minit){
	
	beta3_before = beta3
	value_before = value

alpha_fit = function(subj_index){
	(timei =timepoints2[[subj_index]])
	
	xmat_previous = lapply(pc_fits_previous, function(x){
	eval.fd(timei,x)
	})%>%do.call(cbind,.)
	xmati = eval.fd(timei, pc_fit)
	lm(observed2[[subj_index]]~xmat_previous+xmati+0)%>%coef%>%as.numeric
}

sfit= lapply(1:length(observed2),alpha_fit)%>%do.call(rbind,.)

residual_fit = function(subj_index){
	(timei =timepoints2[[subj_index]])
	xmat_previous = lapply(pc_fits_previous, function(x){
	eval.fd(timei,x)
	})%>%do.call(cbind,.)
	xmati = eval.fd(timei, pc_fit)
	lm(observed2[[subj_index]]~xmat_previous+xmati+0)%>%residuals
}

rfit= lapply(1:length(observed2),residual_fit)%>%do.call(c,.)

value = as.numeric(mean(rfit^2)+gamma*inprod(pc_fit,pc_fit,2,2))
if(abs(value_before - value/value_before) < threshold) break
N = sapply(observed2,length)%>%sum

yem = lapply(1:length(observed2), function(i){
(timei =timepoints2[[i]])
yfits_previous = lapply(1:length(pc_fits_previous), function(x){
	sfit[i,x]*eval.fd(timei,pc_fits_previous[[x]])%>%as.numeric
	})%>%do.call(rbind,.)%>%colSums

(observed2[[i]] - yfits_previous)/sqrt(length(timei))/sqrt(N)
})%>%do.call(c,.)

yem[is.na(yem)] <- 0
xalpha = lapply(1:length(observed2),function(subj_index){

	(timei =timepoints2[[subj_index]])
	xmati = eval.basis(timei, spline_basis)
	(xmati*sfit[subj_index,ncol(sfit)])/sqrt(length(timei))/sqrt(N)

})%>%do.call(rbind,.)


A = xalpha
qmat = 2*(t(A)%*%A+gamma*R)
pmat = as.numeric(-2*t(yem)%*%A)
betamat=  do.call(rbind,betalist)
cmat = rbind(betamat%*%E,-betamat%*%E)
# beta3  = lsei::qp(qmat, pmat, c=cmat, d=rep(0,length(betalist)))
beta3  = lsei::qp(qmat, pmat, e=cmat, f=c(rep(-1e-3,length(betalist)),rep(-1e-3,length(betalist))))
pc_fit = fd(beta3, spline_basis)
pc_fit = (1/sqrt(inprod(pc_fit, pc_fit))*pc_fit)
# plot(pc_fit)
beta3 = coef(pc_fit)%>%as.numeric

thresh =  max(abs(beta3_before-beta3))

loss_value= mean((yem - A%*%beta3)^2)+as.numeric(gamma*t(beta3)%*%R%*%beta3)
#print(value)
it = it+1
if(it%%10==0) {print(it);print(thresh)}
}
pc_fit = fd(beta3, spline_basis)
pc_fit = (1/sqrt(inprod(pc_fit, pc_fit))*pc_fit)
beta3 = coef(pc_fit)%>%as.numeric
print("Done!")
cat(sprintf('The threshold is %s . \n',thresh))
return(list(beta=beta3,pc_fit = pc_fit, sfit = sfit, previous_beta =betalist, thresh = thresh,it=it,value=value,gamma=gamma))
}



third_FPC_conditional2_soap  = function(beta3, observed, timepoints, betalist =previous_beta , threshold=1e-4,minit=1,gamma=0,range_score=1e5){
thresh = 1
it = 1
R = inprod(spline_basis, spline_basis,2,2)
E = inprod(spline_basis, spline_basis,0,0)

pc_fit = fd(beta3, spline_basis)
pc_fit = (1/sqrt(inprod(pc_fit, pc_fit))*pc_fit)
beta3 = coef(pc_fit)%>%as.numeric

pc_fits_previous = lapply(1:length(betalist), function(x){
pc_fit1 = fd(betalist[[x]], spline_basis)
pc_fit1 = (1/sqrt(inprod(pc_fit1, pc_fit1))*pc_fit1)
pc_fit1
})
value = -1e14

observed2 =observed
# observed2[which(sapply(observed,length)<=pc_index)]=NULL
timepoints2 = timepoints
# timepoints2[which(sapply(observed,length)<=pc_index)]=NULL

beta3_before = rep(0, length(betalist[[1]]))
while (thresh >  threshold|it<minit){
	
	beta3_before = beta3
	value_before = value

alpha_fit = function(subj_index){
	(timei =timepoints2[[subj_index]])
	xmati = eval.fd(timei, pc_fit)
	sf = lm(observed2[[subj_index]]~xmati+0)%>%coef%>%as.numeric
	if(sf>range_score) sf =  range_score
	if(sf < -range_score) sf= - range_score
	return(sf)
}

sfit= lapply(1:length(observed2),alpha_fit)%>%do.call(rbind,.)

residual_fit = function(subj_index){
	(timei =timepoints2[[subj_index]])
	xmati = eval.fd(timei, pc_fit)
	lm(observed2[[subj_index]]~xmati+0)%>%residuals

}


rfit0= lapply(1:length(observed2),residual_fit)
rfit  = lapply(rfit0, function(x) {x^2%>%mean})%>%do.call(c,.)

value = as.numeric(mean(rfit)+gamma*inprod(pc_fit,pc_fit,2,2))
# if(abs((value_before - value)/value_before) < threshold) break
N = sapply(observed2,length)%>%sum

yem = lapply(1:length(observed2), function(i){
(timei =timepoints2[[i]])
(observed2[[i]])/sqrt(length(timei))/sqrt(N)
})%>%do.call(c,.)


xalpha = lapply(1:length(observed2),function(subj_index){

	(timei =timepoints2[[subj_index]])
	xmati = eval.basis(timei, spline_basis)
	(xmati*sfit[subj_index,ncol(sfit)])/sqrt(length(timei))/sqrt(N)

})%>%do.call(rbind,.)


A = xalpha
qmat = 2*(t(A)%*%A+gamma*R)
pmat = as.numeric(-2*t(yem)%*%A)
betamat=  do.call(rbind,betalist)
cmat = rbind(betamat%*%E,-betamat%*%E)
# beta3  = lsei::qp(qmat, pmat, c=cmat, d=rep(0,length(betalist)))
beta3  = lsei::qp(qmat, pmat, e=cmat, f=c(rep(-1e-3,length(betalist)),rep(-1e-3,length(betalist))))
pc_fit = fd(beta3, spline_basis)
pc_fit = (1/sqrt(inprod(pc_fit, pc_fit))*pc_fit)
# plot(pc_fit)

if(inprod(pc_fit)<0) {
	pc_fit = - pc_fit
	
}
beta3 = coef(pc_fit)%>%as.numeric
thresh =  max(abs(beta3_before-beta3))

# loss_value= mean((yem - A%*%beta3)^2)+as.numeric(gamma*t(beta3)%*%R%*%beta3)
#print(value)
it = it+1
if(it%%10==0) {print(it);print(thresh)}
}
pc_fit = fd(beta3, spline_basis)
pc_fit = (1/sqrt(inprod(pc_fit, pc_fit))*pc_fit)
beta3 = coef(pc_fit)%>%as.numeric
print("Done!")
cat(sprintf('The threshold is %s . \n',thresh))
return(list(beta=beta3,pc_fit = pc_fit, sfit = sfit, previous_beta =betalist, residual= rfit0, thresh = thresh,it=it,value=value,gamma=gamma))
}


pred_SOAP = function(betalist,ylist, tlist,spline_basis,nminus=1){

cfits = lapply(1:length(ylist), function(x){
	print(x)
	timei = tlist[[x]]
	xmat = lapply(1:length(betalist), function(i){
		pc_fit  = fd(betalist[[i]], spline_basis)
		eval.fd(timei, pc_fit)%>%as.numeric
	})%>%do.call(cbind,.)
	index_pc = min(length(ylist[[x]])-nminus, length(betalist))
	if(index_pc<=1) {
	    index_pc=1
		cfit = rep(0,length(betalist))
		# cfit[which.min(abs(ylist[[x]]/xmat))]=(ylist[[x]]/xmat)[which.min(abs(ylist[[x]]/xmat))]

		cfit_temp  = sapply(1:length(betalist),function(j) {
			lm(ylist[[x]]~0+xmat[,j])%>%coef%>%as.numeric
		})
		cfit[which.min(abs(cfit_temp))] = cfit_temp[which.min(abs(cfit_temp))]
		
	} else {
			
			cfit = lm(ylist[[x]]~0+xmat[,1:index_pc])%>%coef%>%as.numeric
	}
	if(length(cfit)<length(betalist)) cfit = c(cfit, rep(0,length(betalist)-length(cfit)))
	cfit
})%>%do.call(rbind,.)

yfits = lapply(1:nrow(cfits), function(i){
	cfit = cfits[i,]
	rowSums(mapply("*",cfit, betalist,SIMPLIFY=TRUE))
}
)%>%do.call(cbind,.)

yfitsfd = fd(yfits, spline_basis)


residuals = sapply(1:length(ylist), function(x){
	timei = tlist[[x]]
	resid = ylist[[x]] - as.numeric(eval.fd(timei,yfitsfd[x]))
	if (all(resid==0)) {
		return(NA)
	} else {
	return(mean(resid^2))
	}
	
})
sigmahat = mean(residuals[!is.na(residuals)])
list(sigmahat = sigmahat, yfd_fit=  yfitsfd,sfit = cfits)

}