library(fitdistrplus)

args <- commandArgs(trailingOnly = TRUE)

if (length(args)<2) {
       stop("Please input the information of arguments:
        args[1]         input file
        args[2]		CI
	args[3]         output file
  Usage:
        Rscript GetCutoff.R input 0.995 output"
        )
}

input <- args[1]
CI<- args[2]
output <- args[3]

ci<-(1-as.numeric(CI))/2+as.numeric(CI);
dat = read.csv(input,header = F,sep="\t",stringsAsFactors = F)
dat$cutoff<-NA
dat$average<- NA
dat$sd<- NA
#colnames(dat) = c('Marker',"DupRatio","SampleNumber","KLD")
for (i in 1:nrow(dat)){
	data<-unlist(strsplit(dat[i,4],split=","))
	fit<-fitdist(as.numeric(data), distr = "gamma", method = "mle")
	shape<-as.numeric(summary(fit)[[1]][1])
	rate<-as.numeric(summary(fit)[[1]][2])
	dat$average[i]<-shape/rate
	dat$sd[i]<-shape/(rate*rate)
	dat$cutoff[i]<-qgamma(as.numeric(ci),shape, rate=rate, lower.tail=TRUE, log.p = FALSE)
	}
write.table(dat,file=output,row.names=FALSE,quote=F,sep='\t')

