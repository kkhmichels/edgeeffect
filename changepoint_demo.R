# set working directory
# for KM's computer
setwd("~/Desktop/2_Project/_Drafts/Chapters/1_Modern/App Veg Sci/GitHub")

source('changepoint_functions.R')

master<-read.csv('modern.csv')
dimnames(master)[[2]]<-c('no','transect','ogplot','newplot','dist','slope',
                         'aspect','height','windfalls')

#setting up negative values for managed plots and positive values for wilderness plots
master$dist2<-master$dist*I(master$newplot<=12)+(-1)*master$dist*I(master$newplot>12)

head(master)

cp.list<-c(seq(-475,-75,50),0,seq(75,500,50),650,900)

#pdf('plots.and.changepoints.pdf')
plot(master$transect~master$dist2,pch=20,main='Possible changepoints in red',ylab='Transect',xlab='Distance (m)')
abline(v=cp.list,col='red',lty=2)
#dev.off()




##################################

# Analysis for CONTINUOUS variables

# Example: canopy height

fit_height<-changepoint_analysis(master,cp.list,'height','continuous','Canopy Height (m)')
print(fit_height)


################################## 

# Analysis for COUNT variables

# Example: windfalls

fit_windfalls<-changepoint_analysis(master,cp.list,'windfalls','count','Number of Windfalls')
print(fit_windfalls)




