
###############################################################Script to analze the FLY CLEAR-CLIP data
#------------------------------motif enriched locus distribution
library(RColorBrewer)
cols  = brewer.pal(8, "Dark2")
num = read.table("~/lab309/Fly_clash/1_seedDistribution/site.txt")
numfseq = num[num$V1=="fseq" & num$V3>0,]
numrseq = num[num$V1=="rseq" & num$V3>0,]

fuecdf = function(x,n){
  ee = ecdf(x)
  xx = unique(sort(c(seq(0,200, length=201), knots(ee))))
  result = data.frame(xvalue = xx,yvlaue = ee(xx)*length(x)/n)
  result
}

fmer6 = numfseq[numfseq$V2=="6mer",]
rmer6 = numrseq[numrseq$V2=="6mer",]
FFresult1 = fuecdf(fmer6$V3,nrow(numfseq))
RRresult1 = fuecdf(rmer6$V3,nrow(numrseq))

fmer8 = numfseq[numfseq$V2=="8mer",]
rmer8 = numrseq[numrseq$V2=="8mer",]
FFresult2 = fuecdf(fmer8$V3,nrow(numfseq))
RRresult2 = fuecdf(rmer8$V3,nrow(numrseq))

fmer7A1 = numfseq[numfseq$V2=="7merA1",]
rmer7A1 = numrseq[numrseq$V2=="7merA1",]
FFresult3 = fuecdf(fmer7A1$V3,nrow(numfseq))
RRresult3 = fuecdf(rmer7A1$V3,nrow(numrseq))

fmer7m8 = numfseq[numfseq$V2=="7merM8",]
rmer7m8 = numrseq[numrseq$V2=="7merM8",]
FFresult4 = fuecdf(fmer7m8$V3,nrow(numfseq))
RRresult4 = fuecdf(rmer7m8$V3,nrow(numrseq))

fmer6off = numfseq[numfseq$V2=="6merOffset",]
rmer6off = numrseq[numrseq$V2=="6merOffset",]
FFresult5 = fuecdf(fmer6off$V3,nrow(numfseq))
RRresult5 = fuecdf(rmer6off$V3,nrow(numrseq))

write.table(cbind(FFresult1,FFresult2,FFresult3,FFresult4,FFresult5),"forward.txt")
write.table(cbind(RRresult1,RRresult2,RRresult3,RRresult4,RRresult5),"reverse.txt")


linewidth = 3
tiff("~/lab309/CLASH/CLASHALL/CLASH-1111/seedDistribution/seedDistribution_R.tiff",width = 6, height =6, units = "in",res=600,pointsize =12,
     type="cairo",compression = "lzw")

par(mar=c(5.1, 4.1, 4.1, 2.1), mgp=c(2, 0.8, 0), las=0)
cols = colorRampPalette(brewer.pal(6,"Reds"))(6)
fmer6 = numfseq[numfseq$V2=="6mer",]
rmer6 = numrseq[numrseq$V2=="6mer",]
FFresult = fuecdf(fmer6$V3,nrow(numfseq))
RRresult = fuecdf(rmer6$V3,nrow(numrseq))
plot(FFresult,col = cols[3],lwd=linewidth,xlim=c(-200,200),type="l",ylab="CDF: seed type presence",
     xlab = "Distance from ligation point (nts)",ylim=c(0,0.7),las=1,xaxt="n")
lines(-RRresult$xvalue,RRresult$yvlaue,col = cols[3],lwd=linewidth)
axis(1, at=c(-200,-150,-100,-50,0,50,100,150,200), labels=c(-200,-150,-100,-50,0,50,100,150,200)) 
legend("topleft",c("8mer", "7mer-A1", "7mer-m8","6mer","6mer-offset","Total"), 
       col =c(cols[6:2],1),lwd=linewidth,bty="n")


fmer8 = numfseq[numfseq$V2=="8mer",]
rmer8 = numrseq[numrseq$V2=="8mer",]
FFresult = fuecdf(fmer8$V3,nrow(numfseq))
RRresult = fuecdf(rmer8$V3,nrow(numrseq))
lines(FFresult,col = cols[6],lwd=linewidth)
lines(-RRresult$xvalue,RRresult$yvlaue,col = cols[6],lwd=linewidth)

fmer7A1 = numfseq[numfseq$V2=="7merA1",]
rmer7A1 = numrseq[numrseq$V2=="7merA1",]
FFresult = fuecdf(fmer7A1$V3,nrow(numfseq))
RRresult = fuecdf(rmer7A1$V3,nrow(numrseq))
lines(FFresult,col = cols[5],lwd=linewidth)
lines(-RRresult$xvalue,RRresult$yvlaue,col = cols[5],lwd=linewidth)

fmer7m8 = numfseq[numfseq$V2=="7merM8",]
rmer7m8 = numrseq[numrseq$V2=="7merM8",]
FFresult = fuecdf(fmer7m8$V3,nrow(numfseq))
RRresult = fuecdf(rmer7m8$V3,nrow(numrseq))
lines(FFresult,col = cols[4],lwd=linewidth)
lines(-RRresult$xvalue,RRresult$yvlaue,col = cols[4],lwd=linewidth)

fmer6off = numfseq[numfseq$V2=="6merOffset",]
rmer6off = numrseq[numrseq$V2=="6merOffset",]
FFresult = fuecdf(fmer6off$V3,nrow(numfseq))
RRresult = fuecdf(rmer6off$V3,nrow(numrseq))
lines(FFresult,col = cols[2],lwd=linewidth)
lines(-RRresult$xvalue,RRresult$yvlaue,col = cols[2],lwd=linewidth)


fmerTotal = c(fmer6$V3,fmer7A1$V3,fmer7m8$V3,fmer8$V3,fmer6off$V3)
rmerTotal = c(rmer6$V3,rmer7A1$V3,rmer7m8$V3,rmer8$V3,rmer6off$V3)
FFresult = fuecdf(fmerTotal,nrow(numfseq))
RRresult = fuecdf(rmerTotal,nrow(numrseq))
lines(FFresult,col = 1,lwd=linewidth)
lines(-RRresult$xvalue,RRresult$yvlaue,col = 1,lwd=linewidth)

box()
dev.off()



#------------------------------statistic analysis of binding free energy
setwd("~/lab309/Fly_clash/2_freeenergy/")
energy_orig = read.table("freeenergy_clash.txt",sep="\t")
energy_shuffinter = read.table("freeenergy_shuffleI.txt",sep="\t")


library(RColorBrewer)
cols  = brewer.pal(8, "Blues")

plot(density(energy_shuffinter$V2),col=cols[6],lwd=2.5,lty=1,las=1,ylab="Probability density",xlab="Free energy (kcal/mol)",main="",xlim=c(-35,0))
lines(density(energy_orig$V2),col=2,lwd=2.5,lty=1)
legend("topleft",c("CLASH chimeras","Shuffled interaction","Shuffled mRNAs","Random mRNAs")
       ,lwd=2.5,col=c(2,cols[4],cols[6],cols[8]),bty="n", lty=c(1,1,1,1),cex=.8)




#-----------------------------statistic analysis of miRNA abundance and chimeras abundance 
setwd("~/lab309/Fly_clash/3_freq/")
test = read.table("freq_flyclash.txt",sep="\t")

colg2 = colorRampPalette(brewer.pal(6,"Blues"))(6)

testFW0 = test[test$V1=='F' & test$V2=='W0',c(3,4,5)]
testFWH = test[test$V1=='F' & test$V2=='WH',c(3,4,5)]
testFC0 = test[test$V1=='F' & test$V2=='C0',c(3,4,5)]
testFCH = test[test$V1=='F' & test$V2=='CH',c(3,4,5)]

plot(log2(testFFF[testFFF$V2=='P01' | testFFF$V2=='P30',]$V5),log2(testFFF[testFFF$V2=='P01' | testFFF$V2=='P30',]$V4),cex=.8,pch=19,
     col=c(rep(cols[1],74),rep(cols[2],67)),
     xlim=c(0,20),ylim=c(0,10), xlab = "miRNA abundance",ylab="Chimera abundance",las=1)
legend("topleft",c("3h PE     r=0.78", "30h PE   r=0.76"), 
       col =c(cols[1],cols[2]),pch=19,bty="n")



#---------------------------interaction site conservation analysis
setwd("~/lab309/Fly_clash/5_conservation/")

score = read.table("conscore.txt")

plot(198:398,apply(score,2,mean)[198:398],type="l",ylim=c(1.5,2.1),lwd=2,col="gray20")

lines(300:306,apply(score,2,mean)[300:306],col="red",lwd=3)


#-------------------------Homer miRNA target motif analysis
library(pheatmap)
cols = brewer.pal(9, "Blues")
setwd("~/lab309/Fly_clash/4_motif/")

homer = read.table("uHomer_motif.txt",row.names = 1,sep="\t")
breaks = seq(0, 1, by = 0.01)

pheatmap(homer,cluster_cols =F,clustering_distance_rows = "correlation",border_color = "gray90",
         color = colorRampPalette(c("white", "brown4","brown4"))(length(breaks)),breaks = breaks)

pheatmap(homer,cluster_cols =F,clustering_distance_rows = "correlation",border_color = "white",
         color = colorRampPalette(c("gray90", "orange3","orange3"))(length(breaks)),breaks = breaks)

pheatmap(homer,cluster_cols =F,clustering_distance_rows = "correlation",border_color = "NA",
         colorRampPalette(brewer.pal(n = 7, name ="Reds"))(100),breaks = breaks)




#---------------------------expression miR-34 null
miR34d3 = read.table("~/lab309/Fly_clash/6_expAnalysis/miR34/miR34null20d.txt")

assayT = miR34d3[miR34d3$V4>0,]
assayC = miR34d3[miR34d3$V4==0,]
assayO = miR34d3[miR34d3$V4==2,]

fuecdf = function(x,n){
  ee = ecdf(x)
  xx = unique(sort(c(seq(0,n, length=n+1), knots(ee))))
  result = data.frame(xvalue = xx,yvalue = ee(xx)*length(x)/n)
  result
}
library(RColorBrewer)
cols  = brewer.pal(8, "Set1")

background = fuecdf(log2(assayC$V2/assayC$V3),nrow(assayC))
target = fuecdf(log2(assayT$V2/assayT$V3),nrow(assayT))
clipSupport =  fuecdf(log2(assayO$V2/assayO$V3),nrow(assayO))

plot(c(0,1),xlim=c(-1,1),col=1,type="n",las=1,ylab="CDF",xlab="Fold change (log2)")
lines(background$xvalue,background$yvalue,col=1,lwd=2)
lines(target$xvalue,target$yvalue,col=cols[1],lwd=2)
lines(clipSupport$xvalue,clipSupport$yvalue,col=cols[2],lwd=2)

text(-0.6,1,"Non-miR-34-5p chim. (3,898)",col=1,cex=1)
text(-0.6,0.95,"miR-34-5p chrim. (1,314) p = 9.996e-08",col=cols[1],cex=1)
text(-0.6,0.9,"miR-34-5p chrim.+peak (871) p = 1.676e-09",col=cols[2],cex=1)




#--------------------------expression miR-305 null 
miR34d3 = read.table("~/lab309/Fly_clash/6_expAnalysis/miR305/miR305over.txt")

assayT = miR34d3[miR34d3$V4>0,]
assayC = miR34d3[miR34d3$V4==0,]
assayO = miR34d3[miR34d3$V4==2,]

fuecdf = function(x,n){
  ee = ecdf(x)
  xx = unique(sort(c(seq(0,n, length=n+1), knots(ee))))
  result = data.frame(xvalue = xx,yvalue = ee(xx)*length(x)/n)
  result
}
library(RColorBrewer)
cols  = brewer.pal(8, "Set1")

background = fuecdf(log2(assayC$V2/assayC$V3),nrow(assayC))
target = fuecdf(log2(assayT$V2/assayT$V3),nrow(assayT))
clipSupport =  fuecdf(log2(assayO$V2/assayO$V3),nrow(assayO))

plot(c(0,1),xlim=c(-1,1),col=1,type="n",las=1,ylab="CDF",xlab="Fold change (log2)")
lines(background$xvalue,background$yvalue,col=1,lwd=2)
lines(target$xvalue,target$yvalue,col=cols[1],lwd=2)
lines(clipSupport$xvalue,clipSupport$yvalue,col=cols[2],lwd=2)

text(-0.6,1,"Non-miR-305-5p chim. (3,696)",col=1,cex=1)
text(-0.6,0.95,"miR-305-5p chrim. (447) p = 0.005618",col=cols[1],cex=1)
text(-0.6,0.9,"miR-305-5p chrim.+peak (265) p = 0.0191",col=cols[2],cex=1)





#----------------------chimera frequency comparison

cctable = read.table('~/lab508/flysequencing/combineCLASHclip-m.txt',header=F)

#filterwithPeaks
cchim = cctable[(cctable$V10+cctable$V11+cctable$V12+cctable$V13)>0.4,]

plot(log2(cchim$V10),log2(cchim$V12),xlim=c(-3,10),ylim=c(-3,10),type="n")
for (i in 1:nrow(cchim)){
  a = cchim[i,10]
  b = cchim[i,12]
  cccolor = 'gray50'
  if(a==0.1 & b==0.1){
    cccolor = 'gray50'
  }
  if(a==0.1 & b>0.1){
    cccolor = 'darkorange'
  }
  if(a>0 & b==0.1){
    cccolor ='deepskyblue'
  }
  if(a>0.1 & b>0.1){
    if(a/b>=2){
      cccolor='deepskyblue'
    }
    if(b/a>=2){
      cccolor='darkorange'
    }
  }
  points(log2(a),log2(b),col=cccolor,cex=.5,pch=19)
}

for (i in 1:nrow(cchim)){
  a = cchim[i,10]
  b = cchim[i,12]
  if(cchim[i,14]==1){
    cccolor='red'
    points(log2(a),log2(b),col=cccolor,cex=.5,pch=17)
  }
}



#--------------------------------------miRNA abundance
miRNA = read.table('~/Dropbox/DME_CLASH/CLASH_data_analysis/flyAgobound.txt',header=T,row.names = 1)

plot(miRNA[,9:10],pch=19,cex=.7,col="gray60",las=1)
for(i in 1:nrow(miRNA)){
  a = miRNA[i,9]
  b = miRNA[i,10]
  if(abs(a)>=1 | abs(b)>=1){
    points(a,b,col="gray20",cex=0.7,pch=19)
  }
}

for(i in 1:nrow(miRNA)){
  a = miRNA[i,9]
  b = miRNA[i,10]
  if(miRNA[i,12]==1){
    cccolor='red'
    points(a,b,col=cccolor,cex=0.7,pch=19)
  } 
}

abline(h=0,lty=2,col='gray50')
abline(v=0,lty=2,col="gray50")




#------------------------------------diff interaction analysis

clkclip = read.table('/media/xiaonan/backupVT/Fly_clash/clipcluster/clipmatrix_clk.txt',header=F)
#clkclip=clkclip[clkclip$V15=='R',]
c0vsw0 = log2((clkclip$V8+1)/(clkclip$V6+1))
chvswh = log2((clkclip$V9+1)/(clkclip$V7+1))

library("RColorBrewer")
setcolor = brewer.pal(n = 12, name = "Paired")


thresh = 0
plot(c0vsw0,chvswh,pch=19,cex=.4,type="n",ylim=c(-10,10))
for(i in 1:length(c0vsw0)){
  colp = 'gray50'
  cpch = 19
  clkclip[i,11]=0
  if(c0vsw0[i]<=(-thresh) & chvswh[i]>(-thresh) & chvswh[i]<(thresh)){
    colp = setcolor[1]
    clkclip[i,11]=1
  }
  if(c0vsw0[i]>=(thresh) & chvswh[i]>(-thresh) & chvswh[i]<(thresh)){
    colp = setcolor[2]
    clkclip[i,11]=2
  }
  
  if(chvswh[i]<=(-thresh) & c0vsw0[i]>(-thresh) & c0vsw0[i]<(thresh)){
    colp = setcolor[3]
    clkclip[i,11]=3
  }
  if(chvswh[i]>=(thresh) & c0vsw0[i]>(-thresh) & c0vsw0[i]<(thresh)){
    colp = setcolor[4]
    clkclip[i,11]=4
  }
  #if(clkclip[i,15]=='R'){
  #  colp=setcolor[6]
  #}
  points(c0vsw0[i],chvswh[i],pch=cpch,cex=.2,col=colp)
}

#---------------------------------------------plot the clip peak along the gene


clkclip = read.table('/media/xiaonan/backupVT/Fly_clash/clipcluster/clipmatrix_clk.txt',header=F)

#Pdp1 FBgn0016694
#Vri FBgn0016076
#tim FBgn0014396
#cwo FBgn0259938
#Per FBgn0003068

plotclip = clkclip[clkclip$V1=='FBgn0016694',]

xmin = min(plotclip$V3)
xmax = 7825000#max(plotclip$V4)

par(mfrow=c(4,1),mar=c(4,4,1,1))

for(j in 6:9){
  ymin= 0
  ymax = 400#max(plotclip[,j])
  plot(c(xmin,xmax),c(ymin,ymax),type="n")
  for(i in 1:nrow(plotclip)){
    peak = plotclip[i,j]
    if(peak>0.1){
      pos =  (plotclip[i,3]+plotclip[i,4])/2
      pd = peak*cos(seq(-pi/2,pi/2,by=0.01))
      ld = seq(pos-49,pos+50,by=99.1/315)
      #ld = seq(pos-24,pos+25,by=49/157)
      polygon(ld,pd,col=2,border=NA)
    }
  }
}


