


BID <- tpose(fread(lf("BID.dat")))
colnames(BID) <- seq(0,18,by=1) # FH
colnames(BID) <- c(0,4,5,6,seq(8,19,by=1))

#pdf("CorrVisual.pdf", width=18, height=9)

par(mfrow=c(1,2))
BID <- tpose(fread(lf("r0--BID.dat")))[-20,]
colnames(BID) <- seq(0,18,by=1) # FH
corrplot::corrplot(cor(BID), method = "color", type="lower", 
		    tl.srt = 0, tl.offset = 1)
trm <- c(2,3,4,8)
BID <- tpose(fread(lf("r1--BID.dat")))[-trm,]
colnames(BID) <- c(0,4,5,6,seq(8,19,by=1))
fwrite("IH--BID.dat", x=BID, sep=" ")
corrplot::corrplot(cor(BID), method = "color", type="lower", 
		    tl.srt = 0, tl.offset = 1)
dev.off()
