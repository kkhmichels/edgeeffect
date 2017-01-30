

library(lme4)
library(Hmisc)


# A function to get the rates of change that are
# not instant 

# cp is the changepoint of interest

# rate is the distance over which the value of the covariate
#      changes from 0 to 1

# dat is a data frame that contains the distance variable
#     indicating where data were observed

get.rate<-function(cp,rate,dat){
  if(rate=='instant'){
    a<-25
  }  
  if(rate=='fast'){
    a<-50
  }
  if(rate=='moderate'){
    a<-75
  }
  if(rate=='slow'){
    a<-100
  }
  if(rate=='veryslow'){
    a<-250
  }

  ratevar<- 0*(dat$dist2<(cp-a))+
            (1/(2*a)*(dat$dist2-cp)+1/2)*(dat$dist2>=(cp-a)&dat$dist2<(cp+a))+
            1*(dat$dist2>=(cp+a))
  return(ratevar)
}

# A function to plot the covariate

# x is a sequence of distances where the 
#   covariate values are evaluated

# cp is the changepoint of interest

# rate is half of the distance over which the value of the covariate
#      changes from 0 to 1

ratefunc<-function(x,cp,rate){
  if(rate=='instant'){a<-25}
  if(rate=='fast'){a<-50}
  if(rate=='moderate'){a<-75}
  if(rate=='slow'){a<-100}
  if(rate=='veryslow'){a<-250}
  val<- 0*(x<(cp-a))+
        (1/(2*a)*(x-cp)+1/2)*(x>=(cp-a)&x<(cp+a))+
        1*(x>=(cp+a))
  return(val)
}




########################
#                      #
# Changepoint analysis #
#                      #
########################

# "changepoint_analysis" performs the changepoint analysis for continuous and count variables.

# dat is a data frame containing the variable of interest ("varname") and
#     the variable "transect" to be used as a random effect

# cp is a vector of changepoints of interest

# varname is a character string containing the variable of interest

# vartype is a character string indicating the type of variable ('continuous' or 'count')

# ylab_varname is a character string for the y-axis label in the plot

# B is the number of bootstrap samples to generate confidence intervals
#   for the plot (default = 1000)



changepoint_analysis<-function(dat,cp,varname,vartype,ylab_varname,B=1000){

  # Getting all of the changepoint variables

  cp.vars<-matrix(0,length(dat[,1]),length(cp)*5)
  for(j in 1:length(cp)){
    cp.vars[,(5*j-4):(5*j)]<-cbind(get.rate(cp[j],'instant',master),
                              get.rate(cp[j],'fast',master),
                              get.rate(cp[j],'moderate',master),
                              get.rate(cp[j],'slow',master),
                              get.rate(cp[j],'veryslow',master))
  }

  cp.varnames<-paste(rep('changepoint',length(cp.vars[1,])),sort(rep(1:length(cp),5)),rep(c('instant','fast','moderate','slow','veryslow'),length(cp)),sep='')
  dimnames(cp.vars)[[2]]<-cp.varnames

  cp_dat<-data.frame(cp_name=paste(rep('changepoint',length(cp.vars[1,])),sort(rep(1:length(cp),5)),rep(c('instant','fast','moderate','slow','veryslow'),length(cp)),sep=''),
                     cp_val=sort(rep(1:length(cp),5)),
                     rate_val=rep(c('instant','fast','moderate','slow','veryslow'),length(cp)))


  dat2<-cbind(dat,cp.vars)

  # Fitting all of the models

  AICvals<-rep(0,length(cp)*5+1)

  for(i in 1:length(cp.varnames)){
    formula.fit<-paste(varname,'~',cp.varnames[i],'+(1|transect)',sep='')

    if(vartype=='continuous'){
      AICvals[i]<-summary(lmer(formula.fit,data=dat2))$AIC[1]
    }
    if(vartype=='count'){
      AICvals[i]<-summary(glmer(formula.fit,data=dat2,family='poisson'))$AIC[1]
    }

    if((i%%10)==0 | i==length(cp.varnames)){
      cat(paste(round(100*i/length(cp.varnames)),'% of changepoint models fit.\n',sep=''))
    }
  }

  # Fitting the null model with no changepoint

  formula.fit<-paste(varname,'~(1|transect)',sep='')

  if(vartype=='continuous'){
    AICvals[length(cp)*5+1]<-summary(lmer(formula.fit,data=dat2))$AIC[1]
  }
  if(vartype=='count'){
    AICvals[length(cp)*5+1]<-summary(glmer(formula.fit,data=dat2,family='poisson'))$AIC[1]
  }

  names(AICvals)<-c(cp.varnames,'NULL')


  # Determining the best combination of changepoint and rate simultaneously

  cp_select<-cp.list[cp_dat[cp_dat$cp_name==names(which.min(AICvals)),]$cp_val]
  rate_select<-as.character(cp_dat[cp_dat$cp_name==names(which.min(AICvals)),]$rate_val)
  if(rate_select=='fast'){a_select<-50}
  if(rate_select=='moderate'){a_select<-75}
  if(rate_select=='slow'){a_select<-100}
  if(rate_select=='veryslow'){a_select<-250}
  if(rate_select=='instant'){a_select<-25}


  # Comparing the selected model to the null model

  formula.NULL<-paste(varname,'~(1|transect)',sep='')
  formula.changepoint<-paste(varname,'~',names(which.min(AICvals)),'+(1|transect)',sep='')

  if(vartype=='continuous'){
    fit.NULL<-lmer(formula.NULL,data=dat2)
    fit.changepoint<-lmer(formula.changepoint,data=dat2)
  }
  if(vartype=='count'){
    fit.NULL<-glmer(formula.NULL,data=dat2,family='poisson')
    fit.changepoint<-glmer(formula.changepoint,data=dat2,family='poisson')
  }


  # Print summaries of parameter estimates and likelihood ratio test

  print(summary(fit.NULL))
  print(summary(fit.changepoint))
  print(anova(fit.NULL,fit.changepoint))

  # Pulling selected changepoint variable from data frame "dat2"

  cp_select_var<-as.matrix(dat2[dimnames(dat2)[[2]]==names(which.min(AICvals))])

  # Parametric bootstrap confidence intervals

  vals_line_plt<-seq(min(dat2$dist2),max(dat2$dist2),length=1000)
  cov_vals<-ratefunc(vals_line_plt,cp_select,rate_select)  

  if(vartype=='continuous'){
    regline.bs<-matrix(0,length(vals_line_plt),B)
    sig.error<-summary(fit.changepoint)$sigma
    sig.transect<-sqrt(summary(fit.changepoint)$varcor$transect[1,1])
    beta.fit<-summary(fit.changepoint)$coef[,1]
    formula.bs<-paste('response.bs~',names(which.min(AICvals)),'+(1|transect)',sep='')
    for(b in 1:B){
      re.bs<-rnorm(length(unique(dat2$transect)),0,sig.transect)
      re.bs<-re.bs[as.numeric(as.factor(dat2$transect))]
      dat2$response.bs<-rnorm(length(dat2[,1]),cbind(1,cp_select_var)%*%beta.fit,sig.error)+re.bs
      temp<-lmer(formula.bs,data=dat2)

      if((b%%100)==0 | b==B){
        cat(paste(round(100*b/B),'% of bootstrap samples processed.\n',sep=''))
      }

      regline.bs[,b]<-fixef(temp)[1]+fixef(temp)[2]*cov_vals
    }
  }
  if(vartype=='count'){
    regline.bs<-matrix(0,length(vals_line_plt),B)
    sig.transect<-sqrt(summary(fit.changepoint)$varcor$transect[1,1])
    beta.fit<-summary(fit.changepoint)$coef[,1]
    formula.bs<-paste('response.bs~',names(which.min(AICvals)),'+(1|transect)',sep='')
    for(b in 1:B){
      re.bs<-rnorm(length(unique(dat2$transect)),0,sig.transect)
      re.bs<-re.bs[as.numeric(as.factor(dat2$transect))]
      dat2$response.bs<-rpois(length(dat2[,1]),exp(cbind(1,cp_select_var)%*%beta.fit+re.bs))
      temp<-glmer(formula.bs,family='poisson',data=dat2)

      if((b%%100)==0 | b==B){
        cat(paste(round(100*b/B),'% of bootstrap samples processed.\n',sep=''))
      }

      regline.bs[,b]<-exp(fixef(temp)[1]+fixef(temp)[2]*cov_vals)
    }
  }

  #setting margins
  par(mar=c(6,5,4,5))
  #ord<-order(dat2$dist2)

  plot_dat<-data.frame(x=vals_line_plt,mean=rowMeans(regline.bs),L=unlist(apply(regline.bs,1,function(x){quantile(x,.025)})),U=unlist(apply(regline.bs,1,function(x){quantile(x,.975)})))

  plot(x=dat2$dist2,y=dat2[,dimnames(dat2)[[2]]==varname],pch=1,cex=0.5,lwd=1,xlab='',ylab=ylab_varname,axes=F,cex.lab=1.5)
  box()
  lines(plot_dat$mean~plot_dat$x,col='red',lwd=3)
  lines(plot_dat$L~plot_dat$x,col='blue',lwd=3,lty=2)
  lines(plot_dat$U~plot_dat$x,col='blue',lwd=3,lty=2)

  abline(v=cp_select-a_select) #end of change point rate
  abline(v=cp_select+a_select) #start of change point rate
  minor.tick(nx=10, ny=1, tick.ratio=0.5) #needs library Hmisc to add tick marks
  axis(1,at=seq(-500,1000,by=250), labels=T,cex.axis=1.5) #adds tick mark labels
  axis(2,labels=T,cex.axis=1.5,cex.lab=1.5)
  mtext("Managed",side=1,adj=0.075,padj=3.5,cex=1.5) #subtitle left
  mtext("Wilderness",side=1,adj=0.575,padj=3.5,cex=1.5) #subtitle right
  mtext("Distance (m)",side=1,adj=0.33,padj=5,cex=1.5)

  return(list(changepoint=cp_select,rate=rate_select,changepoint_region=c(cp_select-a_select,cp_select+a_select),variable=varname))

}




