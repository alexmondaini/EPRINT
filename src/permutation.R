####
library(data.table)
library(GenomicRanges)
library(ggplot2)
library(regioneR)
library(plotly)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(hrbrthemes)
library(ggplot2)
####
setwd('/my_dir/EPRINT/')

# Read eclip
temp = list.files(path='data/eclip-accessions',pattern="*.bed.gz",full.names = T)
myfiles = lapply(temp, function(x) fread(x,
                                         sep="\t",
                                         header=F,
                                         colClasses = c('character',rep('integer',2),
                                                        rep('character',3),rep('numeric',4)),
                                         col.names = c("chr", "start", "end", "name", "score",
                                                       "strand","signalValue","pValue","qValue","peak_point_source")))
names(myfiles) <- gsub("*.bed.gz$","",basename(temp))
# get only hg19 IDR files
sbset = lapply(myfiles, function(x) grepl('IDR',unique(x$name)))
myfiles = myfiles[unlist(sbset)]
# filter chromosome names
validchrname = c(paste0("chr",seq(22)),"chrX","chrY")

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
pk.df = fread("data/fus-eprint/fus.merge.peak.counts.inputfiltered.fusvsneg.deseq.sfglobal.batch1.wald.txt", sep="\t", header=T, quote="")
# drop last column, it's a duplicate of name
pk.df <- pk.df[,-32]
# paste chr in front of geneChr
pk.df[,geneChr:=paste0('chr',geneChr)]
# change geneChr 23 to X
pk.df[geneChr=='chr23',geneChr:='chrX']
# change strand encode
pk.df[,geneStrand := fcase(geneStrand==1,'+',
                           geneStrand==2,'-')]

# Over the unique geneIDs

# tt = unique(pk.df,by='geneId')
# tt[,count_median:=.(rowMedians(as.matrix(.SD))),.SDcols=c("sample5","sample6","sample7","sample8")]

# Median
pk.df[,count_median:=.(rowMedians(as.matrix(.SD))),.SDcols=c("sample5","sample6","sample7","sample8")]
qt.pk = pk.df[,quantile(count_median)]
qt.pk

pk.df[,qt.pk:=
           fcase(count_median>=qt.pk[1] & count_median<=qt.pk[2],'A',
                 count_median>qt.pk[2] & count_median<=qt.pk[3],'B',
                 count_median>qt.pk[3] & count_median<=qt.pk[4],'C',
                 count_median>qt.pk[4],'D'
           )
]

# frequency of peak quantiles
fq.pk = pk.df[,.(.N),by=qt.pk]
setkey(fq.pk,qt.pk)
frq = fq.pk[,N]
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
gene.set = fread("data/fus-eprint/fus.input.gene.counts.txt",sep = '\t',header = T,quote = "")
setnames(gene.set,"GeneID","geneId")

# median for input
gene.set[,count_median:=.(rowMedians(as.matrix(.SD))),.SDcols=c("sample5","sample6","sample7","sample8")]

# qt.gene = gene.set[,quantile(count_median)]
# 
# # frequency of gene quantiles
# fq.gene = gene.set[,.(.N),by=qt.gene]
# setkey(fq.gene,qt.pk)
# frq.gene = fq.gene[,N]
# frq.gene

gene.set[,qt.pk:=
        fcase(count_median>qt.pk[1] & count_median<=qt.pk[2],'A',
              count_median>qt.pk[2] & count_median<=qt.pk[3],'B',
              count_median>qt.pk[3] & count_median<=qt.pk[4],'C',
              count_median>qt.pk[4],'D'
        )
]

# get only genes that are present in peaks
dt = gene.set[gene.set$geneId %in% pk.df$geneId,]
dt = dt[!is.na(qt.pk)]

# Is there any gene coordinate with less than 100bp ?
i = which(dt[,End-Start<100])
# increase by 101
dt[i,End:=Start+100]

# Start randomization
d = data.table()

start.time = Sys.time()
for (i in 1:1000) {
caseA = sample(dt[qt.pk=='A',geneId],size = frq[1],replace = T)
caseB = sample(dt[qt.pk=='B',geneId],size = frq[2],replace = T)
caseC = sample(dt[qt.pk=='C',geneId],size = frq[3],replace = T)
caseD = sample(dt[qt.pk=='D',geneId],size = frq[4],replace = T)
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


#########################################################################
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


# produce dotplot
ggplot(mapping = aes(x=factor(overlaps),y=cm_pct),data = dat) +
  geom_dotplot(aes(fill=peaks),colour="NA",binaxis = "y", stackdir = "centerwhole",binwidth = 1.2) +
  scale_x_discrete(name="overlaps",breaks=seq_along(1:12),labels=paste0('\u2265',seq_along(1:12))) +
  theme_classic(base_size = 13) + 
  scale_fill_grey(labels=c("eprint_peaks"="eprint\npeaks","mean_random_peaks"="random\npeaks")) +
  labs(y='cumulative percentage',fill=NULL) +
  theme(legend.position = "top",
        legend.key.size = unit(1,'cm'),
        )


# the best one with poinrange
ggplot(mapping = aes(x=factor(overlaps),y=cm_pct,colour=peaks),data = dat) +
  geom_pointrange(aes(ymin=cm_pct-pct_sd-.5,ymax=cm_pct+pct_sd+.5),fatten = 5) +
  scale_x_discrete(name="overlaps",breaks=seq_along(1:12),labels=paste0('\u2265',seq_along(1:12))) +
  theme_classic(base_size = 13) + 
  scale_colour_grey(labels=c("eprint_peaks"="eprint\npeaks","mean_random_peaks"="random\npeaks")) +
  labs(y='cumulative percentage',fill=NULL) +
  theme(legend.position = "top",
        legend.key.size = unit(1,'cm'),
  )+
  guides(color = guide_legend(override.aes=list(linetype = c("blank", "solid")))) 
