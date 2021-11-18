#First must have dada2 and Bioconductor installed: https://benjjneb.github.io/dada2/dada-installation.html 
#Must also have phyloseq installed: http://joey711.github.io/phyloseq/install.html
dyn.load(paste("RPluMA", .Platform$dynlib.ext, sep=""))
source("RPluMA.R")

library(ShortRead)
library(dada2); 
packageVersion("dada2")


input <- function(inputfile) {
  pfix <<- prefix()
  if (length(pfix) != 0) {
     pfix <<- paste(pfix, "/", sep="")
  }

  parameters <- read.table(inputfile, as.is=T);
  rownames(parameters) <- parameters[,1]
  fastq <<- toString(parameters["FASTQ", 2])
  filtered <<- paste(pfix, toString(parameters["filtered", 2]), sep="")
  dadapfx <<- paste(pfix, toString(parameters["DADA", 2]), sep="")
  chimpfx <<- paste(pfix, toString(parameters["merged", 2]), sep="")
}

run <- function() {
	dadaFs <<- readRDS(paste(dadapfx, "forward", "rds", sep="."))
	dadaRs <<- readRDS(paste(dadapfx, "reverse", "rds", sep="."))
        out <<- read.csv(filtered, header=TRUE)
        seqtab.nochim <<- readRDS(paste(chimpfx, "seqtab", "rds", sep="."))
        merger1.nochim <<- readRDS(paste(chimpfx, "merger1", "rds", sep="."))
}

output <- function(outputfile) {
fnFs <- sort(list.files(fastq, pattern="_R1.fastq", full.names = TRUE))
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
print(fnFs)
print(sample.names)
print(out)
out = out[,-1]
avgin<-mean(out[,1]) #average number of input reads per sanmple

mean(out[,2]) #average number of reads per sample after filterAndTrim

b = mean(out[,2]) - mean(out[,1]) #filtered reads - #input reads

perclost<-b/mean(out[,1]) #percentage change between reads in and reads filtered

print("Track Reads...")


getN <- function(x) sum(getUniques(x))

track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(merger1.nochim, getN), rowSums(seqtab.nochim))

print("Tracked")


# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)

colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
print(sample.names)
rownames(track) <- sample.names

head(track)



x = (mean(track[,3]) + mean(track[,4]))/2 #mean denoised reads per sample

x

percdenoise<-(x - mean(track[,1]))/mean(track[,1]) #percent change between input reads and denoised reads

percdloss<-(percdenoise-perclost) #reads lost at denoising


readsprior<-mean(track[,5]) #mean merged reads per sample

totalpercloss<-(mean(track[,5])-mean(track[,1]))/mean(track[,1]) #percentage change between input reads and merged reads

mergepercloss<-(totalpercloss-percdenoise) #reads lost at merging

#create single column table with the data of interest. 
#In order from top to bottom: average number of input reads, percent reads lost at initial filtering, percent reads lost at denoising, % lost at merging, total % lost, average number of reads prior to taxa assignment


etable<-cbind(abs(avgin),abs(perclost)*100,abs(percdloss)*100,abs(mergepercloss)*100,abs(totalpercloss)*100,readsprior) #make a table for the values i want

etablemat<-matrix(etable, nrow=6, ncol=1)

print(track)
print(etablemat)
write.csv(track, paste(outputfile, "track", "csv", sep="."))
write.csv(etablemat, paste(outputfile, "etable", "csv", sep="."))

}


