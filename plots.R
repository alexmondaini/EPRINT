####
library(data.table)
library(GenomicRanges)
library(ggplot2)
library(regioneR)
library(plotly)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(hrbrthemes)
####
#Read all Eclip files
setwd('/my_dir/R_eprint_project/files/')
temp = list.files(pattern="*.bed.gz")
myfiles = lapply(temp, function(x) fread(x,
                                         sep="\t",
                                         header=F,
                                         colClasses = c('character',rep('integer',2),
                                                        rep('character',3),rep('numeric',4)),
                                         col.names = c("chr", "start", "end", "name", "score", "strand","signalValue","pValue","qValue","peak_point_source")))
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

################################
# Density
require(scales)
number_of_peaks = lapply(myfiles, function(x) nrow(x))
df = do.call(rbind.data.frame, number_of_peaks)
names(df) <- 'number_of_peaks'
# PLot Density of # of peaks
ggplot(df,aes(number_of_peaks)) + 
  geom_area(stat = "bin",bins=100) +
  scale_x_continuous(breaks = seq(0,30000,2000),labels = comma) + 
  scale_y_continuous(breaks = seq(0,30,2)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle("Distribution of number of peaks over all accessions")

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

# Make Grange list
myfiles <- GRangesList(myfiles)

# TXDB
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
transcriptome <-  transcripts(txdb)
transcriptome <-  transcriptome[seqnames(transcriptome)%in%validchrname]

# assign lenghts to chromosomes
len = seqlengths(transcriptome)[1:24]
lapply(myfiles, function(x) seqlengths(x) <- len)

# widths
ls_width = lapply(myfiles, function(x) width(x))
n.obs <- sapply(ls_width, length)
seq.max <- seq_len(max(n.obs))
mat <- sapply(ls_width, "[", i = seq.max)
df = as.data.table(mat)
# Melt for plotting
df = melt(df,variable.name = 'accession',value.name = 'width')
# Plot width of each eclip dataset
ggplot(df,aes(width)) + geom_density() + facet_wrap(vars(accession))

# Coverage
# Get the percentage of bases in each element of the list of ranges which intersect with transcriptome
per = lapply(myfiles,function(x){
  (sum(width(GenomeInfoDb::intersect(reduce(x), reduce(transcriptome), ignore.strand = F))) / sum(width(reduce(x))))*100
})

# Buil df
setDT(per)
per = melt(per,variable.name = 'accession',value.name = 'percent')
# plot barplot
p = ggplot(per,aes(x=accession,y=percent,colour=accession)) + 
  geom_col(width=0.4, position = position_dodge(width=0.5))  +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        plot.title = element_text(hjust = 0.5)) +
  ggtitle("Plot of % coverage of each accession \n by Txdb transcriptome")

ggsave(file='../images/coverage.png',p, width = 1920/72, height = 1080/72, dpi = 72)


# Read eprint

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

# make a Grange object
pk <- makeGRangesFromDataFrame(pk.df,
                               keep.extra.columns=TRUE,
                               ignore.strand=FALSE,
                               seqinfo=NULL,
                               seqnames.field="seqnames",
                               start.field="start",
                               end.field="end",
                               strand.field="strand",
                               starts.in.df.are.0based=TRUE)

# defnining lenghts for each chromosmoe
seqlengths(pk) = head(len,-1)

# Distance to the closest peak
dst_near = lapply(myfiles, function(x) distanceToNearest(pk,x))
dst_near = lapply(dst_near, function(x) mcols(x))
dst_near = lapply(dst_near, function(x) as.data.table(x))
dst_near = rbindlist(dst_near,idcol = "accession")

# PLot
ggplot(data=dst_near, aes(x=distance, group=accession, fill=accession)) +
  geom_density(adjust=1.5,position = 'fill') +
  theme_ipsum() +
  # facet_wrap(~accession) +
  theme(
    legend.position="none",
    panel.spacing = unit(0.1, "lines"),
    axis.ticks.x=element_blank(),
    plot.title = element_text(hjust = 0.5)
  ) +
  scale_x_continuous(breaks = seq(0,250000000,25000000),labels = comma) +
  ggtitle("Density of Distance from EPRINT to Nearest Region in \n in each ENCODE/ECLIP accession")

# Transform to plot correlation
cor_dt = pk.df[,.(.N,mean_read_count=mean(c(sample5,sample6,sample7,sample8))),by=geneId]

ggplot(data = cor_dt,aes(x=log2(mean_read_count),y=N)) + geom_point() + geom_smooth(method = 'lm')

# or with correlation coefficient in the graph
library(ggpubr)
cor_dt[,log2_mean_read_count:=log2(mean_read_count)]
setnames(cor_dt,"N","peaks_per_gene")

sp <- ggscatter(cor_dt, x = "log2_mean_read_count", y = "peaks_per_gene",
                add = "reg.line",  # Add regressin line
                add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                conf.int = TRUE # Add confidence interval
)
# Add correlation coefficient
sp + stat_cor(method = "pearson", label.x = 9, label.y = 100) +
  ggtitle("Eprint peaks mean read count per gene ")

