suppressMessages(library(corrplot))

to.list <- function( x, y ) { lapply( (unique( y )), function( i ) { which( x == i ) } ) }
newto.list <- function( x, y, z ) { lapply(z , function( i ) { x[which( y == i )] } ) }

args <- commandArgs(TRUE)

a = as.matrix(read.table(paste(args[1],"/sitesData.txt",sep=""))) #sitesData

info = read.table(paste(args[1],"/info.txt",sep="")) #infoBest

b = info[,3]
newb <- sort(unique(b),decreasing=T)

wlines = 12
h = 1000
w = 600

if(length(b) > 20000) { wlines = 40}
newa = matrix(nrow=dim(a)[1] + wlines * (length(newb)-1) , ncol = dim(a)[2])
newk = 1
k = 1
for (i in newb) {
    newa[newk:(newk+length(which(b==i)) -1),] = a[b==i,]
    newk = newk +length(which(b==i))

    if(i > min(newb)) {
        for(l in newk:(newk+wlines-1)) {
            newa[l,] =  a[1,] -  a[1,] + 5
        }
        newk = newk + wlines
    }
    
}

png(paste(args[1],"/fullPartition.png",sep=""),h=h/100,width = w/100,units="in",res=600)
image(1:dim(newa)[2],1:(newk-1),t(newa),col=c("green","blue","orange","red","black","white"),xaxt='n',yaxt='n',ylab=paste(length(b),"sequences clustered into",length(newb), "modules",sep=" "),xlab="motif occurrence",useRaster = TRUE)

starts = c(0,which(colSums(newa) == 5 * dim(newa)[1]),dim(newa)[2])

for(i in 1:(length(starts)-1)) {
    text((starts[i]+starts[i+1])/2,dim(newa)[1]+2,paste("sites",i,sep =" "),srt=90,col="darkgreen",xpd=NA,adj=0,cex=0.8)
}
legend(grconvertX(0.5, "device"), grconvertY(1, "device"), c("A","C","G","T"), fill = c("green","blue","orange","red"),xpd = NA)
garbage = dev.off()




sites = info[,4:dim(info)[2]]
sites[!is.na(sites)] = 1
sites[is.na(sites)] = 0

k = matrix(nrow=length(unique(b)),ncol=dim(sites)[2])

for(i in 1:length(unique(b))) {
    k[i,] = apply(sites[b==i,],2,mean)
}

png(paste(args[1],"/circlePlot.png",sep=""),h=dim(k)[1]/2 , w = dim(k)[2]/2 ,units="in",res=300)
par(xpd=TRUE)
newk = k
newk[newk == 0] = NA
rownames(newk) = paste("module",c(1:dim(newk)[1]))
colnames(newk) = paste("motif",c(1:dim(newk)[2]))
corrplot(newk,is.corr=F,na.label.col="white",cl.pos="b",cl.ratio=0.2,cl.length=3,cl.cex=0.5,outline=T,tl.cex=0.7,tl.col="darkgreen")#,col=col2(20))
points(-10, 7, pch = 8,col="red")
garbage = dev.off()

