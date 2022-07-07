####
library(data.table)
library(GenomicRanges)
library(ggplot2)
library(regioneR)
library(plotly)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(hrbrthemes)
####
setwd('/my_dir/R_eprint_project/files/')

# Read eclip
temp = list.files(pattern="*.bed.gz")
myfiles = lapply(temp, function(x) fread(x,
                                         sep="\t",
                                         header=F,
                                         colClasses = c('character',rep('integer',2),
                                                        rep('character',3),rep('numeric',4)),
                                         col.names = c("chr", "start", "end", "name", "score",
                                                       "strand","signalValue","pValue","qValue","peak_point_source")))
names(myfiles) <- make.names(gsub("*.bed.gz$", "", temp))
# get only hg19 IDR files
sbset = lapply(myfiles, function(x) grepl('IDR',unique(x$name)))
myfiles = myfiles[unlist(sbset)]
# filter chromosome names
validchrname = c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10",  
                 "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", 
                 "chr20", "chr21", "chr22", "chrX", "chrY")

lapply(myfiles,function(x) setkey(x,chr))
myfiles <- lapply(myfiles, function(x) x[chr%in%validchrname])
# view
lapply(myfiles[1:2], head)
# Make Grange object
myfiles <- lapply(myfiles, function(x) makeGRangesFromDataFrame(x,
                                                                keep.extra.columns=TRUE,
                                                                ignore.strand=FALSE,
                                                                seqinfo=NULL,
                                                                seqnames.field="chr",
                                                                start.field="start",
                                                                end.field="end",
                                                                strand.field="strand",
                                                                starts.in.df.are.0based=TRUE))
# Make a Grange list object
myfiles <- GRangesList(myfiles)


#################################################################
# Read eprint peaks
pk.df = fread("fus.merge.peak.counts.inputfiltered.fusvsneg.deseq.sfglobal.batch1.wald.txt", sep="\t", header=T, quote="")
# drop last column, it's a duplicate of name
pk.df <- pk.df[,-32]
# paste chr in front of geneChr
pk.df[,geneChr:=paste0('chr',geneChr)]
# change geneChr 23 to X
pk.df[geneChr=='chr23',geneChr:='chrX']
# change strand encode
pk.df[,geneStrand := fcase(geneStrand==1,'+',
                           geneStrand==2,'-')]

# tt = unique(pk.df,by='geneId')
# tt[,count_median:=.(rowMedians(as.matrix(.SD))),.SDcols=c("sample5","sample6","sample7","sample8")]

# Median
pk.df[,count_median:=.(rowMedians(as.matrix(.SD))),.SDcols=c("sample5","sample6","sample7","sample8")]
qt = pk.df[,quantile(count_median)]
qt

pk.df[,qt:=
           fcase(count_median>=qt[1] & count_median<=qt[2],'A',
                 count_median>qt[2] & count_median<=qt[3],'B',
                 count_median>qt[3] & count_median<=qt[4],'C',
                 count_median>qt[4],'D'
           )
]

# frequency of quantiles
fq = pk.df[,.(.N),by=qt]
frq = fq[,N]
frq

# create peak grange
pk <- makeGRangesFromDataFrame(pk.df,
                               keep.extra.columns=TRUE,
                               ignore.strand=FALSE,
                               seqinfo=NULL,
                               seqnames.field="seqnames",
                               start.field="start",
                               end.field="end",
                               strand.field="strand",
                               starts.in.df.are.0based=TRUE)

extend_grange <- function(pk,bp) {
  tmp_start = ifelse(strand(pk)=='+',start(pk)-(bp+1),end(pk)-(bp))
  tmp_end   = ifelse(strand(pk)=='-',end(pk)+(bp),start(pk)+(bp-1))
  start(pk) = tmp_start
  end(pk)   = tmp_end
  return(pk)
}

pk <- extend_grange(pk,50)

##################################################################################################

# Read gene set
gene.set = fread("fus.input.gene.counts.txt",sep = '\t',header = T,quote = "")
setnames(gene.set,"GeneID","geneId")

# median for input
gene.set[,count_median:=.(rowMedians(as.matrix(.SD))),.SDcols=c("sample5","sample6","sample7","sample8")]

gene.set[,qt:=
        fcase(count_median>qt[1] & count_median<=qt[2],'A',
              count_median>qt[2] & count_median<=qt[3],'B',
              count_median>qt[3] & count_median<=qt[4],'C',
              count_median>qt[4],'D'
        )
]

dt = gene.set[!is.na(qt)]

# Is there any gene coordinate with less than 100bp ?
i = which(dt[,End-Start<100])
# increase by 101
dt[i,End:=Start+100]

# Start randomization
d = data.table()

start.time = Sys.time()
for (i in 1:1000) {
caseA = sample(dt[qt=='A',geneId],size = frq[1],replace = T)
caseB = sample(dt[qt=='B',geneId],size = frq[2],replace = T)
caseC = sample(dt[qt=='C',geneId],size = frq[3],replace = T)
caseD = sample(dt[qt=='D',geneId],size = frq[4],replace = T)
allcases = c(caseA,caseB,caseC,caseD)
# match
rdt = dt[match(allcases,dt[,geneId])]
rdt[,newstart:=round(runif(.N,min = Start,max = End-100))]
rdt[,newend:=newstart+100]
rdt = makeGRangesFromDataFrame(rdt,
                         keep.extra.columns=TRUE,
                         ignore.strand=FALSE,
                         seqinfo=NULL,
                         seqnames.field="Chr",
                         start.field="newstart",
                         end.field="newend",
                         strand.field="Strand",
                         starts.in.df.are.0based=TRUE)

v <- countOverlaps(rdt,myfiles)
d[,paste0('random_gen',i):=v]
}
end.time = Sys.time()
end.time - start.time


saveRDS(d,file = 'd.rds')
# load d
d = readRDS(file = 'd.rds')
## Read d
d[1:5,1:5]
# save number of peaks
n_peaks = dim(d)[1]
# melt columns so we can count by peaks/overlaps groups
d <- melt(d,measure.vars = colnames(d),variable.name = 'peaks',value.name = 'overlaps')
# groupby and count
d <- d[,.N,.(peaks,overlaps)]
# get pct of mean
d = d[,.(pct=(mean(N)/n_peaks)*100,pct_sd=(sd(N)/n_peaks)*100),keyby=overlaps]
# get full pct
full_pct = d[-1,sum(pct)]
# subset
d = d[2:13]
# cumsum
d[,cm_pct:=full_pct-cumsum(c(0,pct[-length(pct)]))]


##########################################################################

# Count Overlaps
result <- as.data.table(countOverlaps(pk,myfiles))
# change column name
setnames(result,'V1','overlaps')
setkey(result,overlaps)
n_peaks = dim(result)[1]
# count how many by key
result <- result[,.(.N),by=overlaps]
#order
result <- result[order(overlaps)]
# add pct column 
result[,pct:=(N/n_peaks)*100]
# pct remaining after 0 overlaps
remainder = result[-1,sum(pct)]
# only 12 peaks from first
result = result[2:13]
# sum of vector pct
result[,cm_pct:=remainder-cumsum(c(0,pct[-length(pct)]))]

###########################################

# combine
dat <- rbind(result,d,fill=T)
dat[,peaks:=rep(factor(c('eprint_peaks','mean_random_peaks')),times = c(nrow(result),nrow(d)))]

# \u2265 for greater than or equal

ggplot(data = dat,aes(x=factor(overlaps),y=cm_pct)) + geom_col(aes(fill=peaks),position = 'dodge2',colour="NA") +
  geom_errorbar(aes(ymin=cm_pct-pct_sd,ymax=cm_pct+pct_sd),
                position=position_dodge2(padding = 0.5)) +
  scale_x_discrete(name="overlaps",breaks=seq_along(1:12),labels=paste0('\u2265',seq_along(1:12))) +
  scale_fill_ft() +
  theme_classic()
