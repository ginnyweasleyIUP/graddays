### functions needed to generate Figures for GMST recon paper in Nature Geoscience. RN 2019/03/29


### read a table containing a time series, first column = years
read.ts<-function(filename,sep=";",header=T,skip=0){
  ind<-as.matrix(read.table(filename,sep=sep,header=header,skip=skip))
  data<-ts(ind[,-1],start=ind[1,1])
  return(data)
}

### apply a function to a timeseries-matrix. the outcome is a time series with the same tsp
tsapply<-function(x,dim,fun){
  out<-ts(apply(x,dim,fun),start=start(x)[1])
}

### anomalien machen für vektor/matrix
anomalies.period<-function(data,start,end){
  if(is.null(dim(data))==T){
    data<-data-mean(window(data,start=start,end=end),na.rm=T)
  }else{ 
    sy<-which(time(data)==start)
    ey<-which(time(data)==end)
    data<-ts(apply(data,2,function(x) x-mean(x[sy:ey],na.rm=T)),start=start(data)[1])
  }
  return(data)
}

###filter time series with custom filter:
tsfilt<-function(x,width=31,method="loess",cut.end=T){
  #possibilities
  #"loess" --> loessfilt function
  #"spline2" --> splinesmoother2 from Dave
  # "gauss" --> gaussfilter
  #"rm" --> running mean
  #"hamming" --> hamming
  #"bw" --> butterworth
  if(method=="loess"){
    filtered<-loessfilt(x,width,cut.end=cut.end)
  }
  if(method=="spline2"){
    if(length(which(is.na(x)))>0){
      if(is.null(dim(x))==T){
        filtered<-splinesmoother2.nas(x,width,cut.end=cut.end)
      }else{ 
        filtered<-ts(apply(x,2,function(y) splinesmoother2.nas(y,width,cut.end)),start=start(x)[1])
      }
    }else{
      filtered<-splinesmoother2(x,width)
      if(cut.end==T){
        sx<-fy(x)+floor(width/2)-1
        ex<-ly(x)-floor((width-0.5)/2)+1
        if(is.null(dim(x))==T){
          window(filtered,end=sx)<-NA
          window(filtered,start=ex)<-NA 
        }else{
          for (i in 1:dim(x)[2]){
            window(filtered[,i],end=sx[i])<-NA
            window(filtered[,i],start=ex[i])<-NA
          }
        }
      }
    }
  }
  if(method=="gauss"){
    if(is.null(dim(x))==T){
      filtered<-gauss.na(x,width)
    }else{ 
      filtered<-ts(apply(x,2,function(y) gauss.na(y,width)),start=start(x)[1])
    }
  }
  if(method=="rm"){
    if(is.null(dim(x))==T){
      filtered<-rollmean.na(x,width)
    }else{ 
      if(length(which(is.na(x)))>0){
        filtered<-ts(apply(x,2,function(y) rollmean.na(y,width)),start=start(x)[1])
      }else{
        filtered<-rollmean(x,width)
      }
    } 
  }
  if(method=="hamming"){
    require(oce)
    filtered<-stats::filter(x,makeFilter("hamming", 50, asKernel=FALSE))
  }
  if(method=="bw"){
    if(is.null(dim(x))==T){
      filtered<-butterfilt.na(x,width)
    }else{ 
      filtered<-tsapply(x,2,function(y) butterfilt.na(y,width))
    }
    if(cut.end==T){
      sx<-start(x)[1]+floor(width/2)-1
      ex<-end(x)[1]-floor((width-0.5)/2)+1
      window(filtered,end=sx)<-NA
      window(filtered,start=ex)<-NA
    }
  }
  
  
  filtered
}

#now low pass butterfilt only
butterfilt.na<-function(y,tsc,type="low",order=2){
  #require(signal)
  # dt is expected to be equal to 1 sampling unit!!
  #first get rid of NAs at begining and end
  sy<-min(which(!is.na(y)))
  ey<-max(which(!is.na(y)))
  x<-y[sy:ey]
  mx<-mean(x)
  #anomalies to full peroid mean may help to reduce end effects in some cases
  x<-x-mx
  nx<-tsc*order
  x2<-c(rep(mean(x[1:tsc]),nx),x,rep(mean(x[(length(x)-tsc+1):length(x)]),nx))
  bf <- signal:::butter(order, 1/tsc, type=type)
  b1 <- signal:::filtfilt(bf, x2)
  b1 <- b1[-c(1:nx,(length(b1)-nx+1):length(b1))]
  #pl.mts(cbind(x2,b1))
  b1<-b1+mx
  z<-y
  z[sy:ey]<-b1
  return(z)
  #detach(package:signal,unload=TRUE,force=TRUE)
  #unloadNamespace("signal")
  #  freqz(b1)
  #  zplane(bf)
}


## same but allowing NAs at the beginning and end of series
## There are end effect with this filter, 2 options to deal with this
## are implemented:
## end.m="allmean":
##  subtract the mean because of end effects of bw filter, where it is biased towards 0
##  note that the end effects remain if the ts is far away from 0 at one end.
## end.m="pad":
##  add the mean over the last tsc.low years tsc.up*2 times at the end, analogous at the beginning.
##  Tthis may also lead to artifacts depending on time scales chosen.
## use end.m="n" or so to not do any correction
bandpass.tsc.na<-function(y,tsc.low,tsc.up,cut.end=T,end.m="allmean"){
  #require(signal)
  # dt is expected to be equal to 1 sampling unit!!
  sy<-min(which(!is.na(y)))
  ey<-max(which(!is.na(y)))
  x<-y[sy:ey]
  
  if(end.m=="allmean"){
    mx<-mean(x)
    x<-x-mx 
  }
  if(tsc.low>tsc.up) warning("Timescales messed up")
  if(end.m=="pad"){
    nx<-tsc.up*2
    nx2<-tsc.low
    x<-c(rep(mean(x[1:nx2]),nx),x,rep(mean(x[(length(x)-nx2+1):length(x)]),nx))
  }
  bf <- signal:::butter(2, c(1/tsc.up,1/tsc.low), type="pass")
  b1 <- signal:::filtfilt(bf, x)
  if(end.m=="pad") b1 <- b1[-c(1:nx,(length(b1)-nx+1):length(b1))]
  #b1<-b1+mx
  if(cut.end==T){
    b1[1:floor(tsc.up/2)]<-NA
    b1[(length(b1)-floor((tsc.up-0.5)/2)+1):length(b1)]<-NA
  }
  z<-y
  z[sy:ey]<-b1
  return(z)
  #detach(package:signal,unload=TRUE,force=TRUE)
  #unloadNamespace("signal")
  #  freqz(b1)
  #  zplane(bf)
}


###plot time series with shading based on quantiles
# dolab= should the legend with the shading be plotted on the right?
# addmedian: should the line of the ensenbmle median be plotted? 
# inq is the increase between the quantiles to plot.
# quants= pre-defined quantiles to plot shaded.
# mcol=color to plot the median,
# scol=color(rgb form, only the first 3 values, not the alpha value) for the quantile shading.
# sepfac col and maxcol define some increase factor for darkness of the colours and the darkest value to use (third value in the rgb call)
# xl=xlim,yl=ylim
plot.quantile.densities<-function(x,dolab=T,addmedian=T,inq=0.01,quants=seq(inq,1-inq,by=inq),mcol="darkblue",scol=c(0,0,1),sepfac.col=21,maxcol=0.05,xl=NA,yl=NA){
  quantsx<-(0:100)/(length(quants)+1)
  mq<-which(quants==0.5)
  if(is.na(xl)) xl<-tsp(x)[1:2]
  if(is.na(yl)) yl<-c(min(x,na.rm=T),max(x,na.rm=T))
  xq<-qnorm(quants)
  if(sepfac.col<(1/maxcol)) sepfac.col<-1/maxcol+1
  colsy2a<-(xq-min(xq)+1)/(max(xq-min(xq)+1)*sepfac.col)
  colsy2<-colsy2a
  colsy2<-colsy2/max(colsy2)*maxcol
  labs<-c(50,"40 / 60", "30 / 70","20 / 80", "10 / 90", "0 / 100")
  qts<-ts(t(apply(x,1,quantile,probs=quants,na.rm=T)),start=start(x))
  plot(qts[,mq],yaxs="i",xaxs="i",col=mcol,lwd=2,ylim=yl,xlim=xl,xlab="",ylab="",tick=F,labels=F,bty="n",t="n")
  for(i in 1:((length(quants)-1)/2)){
    tspolygon(qts[,i],qts[,length(quants)-i+1],col=rgb(scol[1],scol[2],scol[3],colsy2[i]))
  }
  if(addmedian){
    par(new=T)
    plot(qts[,mq],yaxs="i",xaxs="i",lwd=2,col=mcol,ylim=yl,xlim=xl,xlab="",ylab="",tick=F,labels=F,bty="n")
  }
  xats.l<-c(xl[2]+(xl[2]-xl[1])*0.01,xl[2]+(xl[2]-xl[1])*0.01+(xl[2]-xl[1])*0.01)
  if(dolab){
    for(i in 1:((length(quantsx)-1)/2)){
      rect(xats.l[1],yl[2]-(quantsx[length(quantsx)-i+1]-quantsx[i])*(yl[2]-yl[1]),xats.l[2],yl[2],col=rgb(scol[1],scol[2],scol[3],colsy2[i]),border=NA,xpd=T)
    }
    axis(4,line=1,at=seq(yl[2],yl[1],by=-(yl[2]-yl[1])/5),labels=labs,cex.axis=1,las=1,hadj=0.3,cex.axis=1,mgp=c(3,1.5,0))
  }
}

#add polygon to a timeseries plot
#lower and upper are the timeseries for the lower and upper bounds resp.
tspolygon<-function(lower,upper,col=rgb(0,0,0,0.4),border=NA){
  if(is.na(min(lower))){
    lower<-na.omit(lower)
  }
  if(is.na(min(upper))){
    upper<-na.omit(upper)
  }
  if(is.na(border)) border<-col
  xx<-c(tsp(lower)[1]:tsp(lower)[2],tsp(lower)[2]:tsp(lower)[1])
  polygon(xx,c(lower,rev(upper)),border=border,col=col)
}


#make boxplots using quantile-based whiskers instead of fractions of the interquartile range. end of whiskers will be at quantiles defined by min and max.
#call bxp.q to allow for multiple datasets; plot via bxp()
boxplot.quantiles<-function(x,min=0.1,max=0.9){
  bx<-boxplot(x,plot=F)
  if(is.list(x)==T){
    minx<-unlist(lapply(x,function(x) quantile(x,probs=min,na.rm=T)))
    maxx<-unlist(lapply(x,function(x) quantile(x,probs=max,na.rm=T)))
    bx$stats[1,]<-minx
    bx$stats[5,]<-maxx
  }else{
    minx<-quantile(x,probs=min,na.rm=T)
    maxx<-quantile(x,probs=max,na.rm=T) 
    bx$stats[1]<-minx
    bx$stats[5]<-maxx
  }
  
  bx$stats[1]<-minx
  bx$stats[5]<-maxx
  if(is.list(x)==T){
    bx$out<-unlist(lapply(x,function(x) x[which(x<minx | x>maxx)]))
  }else{
    bx$out<-x[which(x<minx | x>maxx)]
  }
  bx$group<-rep(1,length(bx$out))
  return(bx)
}


bxp.q<-function(x,min=0.1,max=0.9){
  if(is.list(x)==T){
    
    for(n in 1:length(x)){
      bx<-boxplot.quantiles(x[[n]],min,max)
      if(n==1){
        out<-bx
      }else{
        out$stats<-cbind(out$stats,bx$stats)
        out$n<-c(out$n,bx$n)
        out$conf<-cbind(out$conf,bx$conf)
        out$out<-c(out$out,bx$out)
        out$group<-c(out$group,rep(n,length(bx$out)))
        out$names<-c(out$names,n)
      }
    }
  }else{
    if(length(dim(x)[2])>0){
      
      for(n in 1:dim(x)[2]){
        bx<-boxplot.quantiles(x[,n],min,max)
        if(n==1){
          out<-bx
        }else{
          out$stats<-cbind(out$stats,bx$stats)
          out$n<-c(out$n,bx$n)
          out$conf<-cbind(out$conf,bx$conf)
          out$out<-c(out$out,bx$out)
          out$group<-c(out$group,rep(n,length(bx$out)))
          out$names<-c(out$names,n)
        }
      }
    }else{
      out<-boxplot.quantiles(x,min,max)
    }
  }
  return(out)
}

##fast linear trend and detrend functions
fasttrend<-function(y){
  l<-length(y)
  x<-1:l
  tx<-x-(l/2+0.5)
  sum((y-mean(y))*(tx))/(sum((tx)^2))
}

#"manual" linear detrending faster than with lm
fastdetrend<-function(y,nas=F){
  l<-length(y)
  x<-1:l
  tx<-x-(l/2+0.5)
  if(nas){
    tr<-sum((y-mean(y,na.rm=T))*(tx),na.rm=T)/(sum((tx)^2))
  }else{
    tr<-sum((y-mean(y))*(tx))/(sum((tx)^2))
  }
  z<-y-(x-1)*tr
}


