#############
library(regioneR)
library(DESeq2)
#############
setwd('C:/Users/alexa/Documents/HKU/r_projects/ak_scripts/')

pk.df = read.table("fus.merge.peak.eprint.counts.gene.counts.txt", sep="\t", header=T, quote="")

deseq.obj = readRDS('fus.merge.peak.counts.fusvsneg.deseq.lrt.dds')
res.df = read.table(file="fus.merge.peak.counts.fusvsneg.deseq.lrt.txt",sep="\t", header=T)


# Some insights 
ggplot(pk.df, aes(width)) +
  geom_histogram(bins = 50,binwidth = 10) +
  ggtitle('DE peaks | Our peaks') +
  annotate('text',x=200,y=40000,label=paste('Total number of peaks',length(pk.df$width),sep = ': ')) +
  scale_x_continuous(limits=c(0, 300),breaks = seq(0,300,by=25),"WIDTH") +
  ylab('Count')


ggplot(pk.df, aes(width)) +
  geom_step(stat = 'ecdf',col='darkgreen',size=1.2) +
  ggtitle('DE peaks | Our peaks') +
  #annotate('text',x=200,y=40000,label=paste('Total number of peaks',length(pk.df$width),sep = ': ')) +
  scale_x_continuous(limits=c(0, 300),breaks = seq(0,300,by=25),"WIDTH") +
  ylab('Cumulative Density')
##############


# Subset pk.df by the peaks that already passed these filters:
# 1. include only intragenic peaks i.e. column intragenic should be 1
# 2. sum(eprint counts) > 10
# 3. sum(input counts) > 10

index_pk.df <- row.names(deseq.obj)
pk.df <- subset(pk.df,name%in%index_pk.df)
pk.df <- cbind(pk.df,log2FoldChange = res.df$log2FoldChange,padj = res.df$padj)

# Extend peak start by 100 bp if strand +, else extended peak end by -100bp.

extension <- function(peak_df) {
  idx <- peak_df$strand == '+'
  peak_df$end[idx] = peak_df$start[idx] + 100
  peak_df$start[!idx] = peak_df$end[!idx] - 100
  return(peak_df)
}

pk.df <- extension(pk.df)

# padj < 0.05 and log2foldchange > 0 (upregulated peaks)
# padj < 0.05 and log2foldchange < 0 (downregulated peaks)
# padj < 0.05 (all.peaks)
# all peaks no matter which padj (universe == pk.df)

up_peaks <- subset(pk.df,log2FoldChange > 0 & padj < .05)
down_peaks <- subset(pk.df,log2FoldChange < 0 & padj < .05)
all.de.peaks <- subset(pk.df,padj < .05)

list_up_down_all <- list(up_peaks = up_peaks,down_peaks = down_peaks,all.de.peaks = all.de.peaks,
                         universe.peaks = pk.df)

# Transform all sorts of peak sets into granges objects

keys <- names(list_up_down_all)
de_grange_list = list()
j = 1
for (df in list_up_down_all) {
  de_grange_list[[keys[j]]] <- makeGRangesFromDataFrame(df,
                                                        keep.extra.columns=TRUE,
                                                        ignore.strand=FALSE,
                                                        seqinfo=NULL,
                                                        seqnames.field="seqnames",
                                                        start.field="start",
                                                        end.field="end",
                                                        strand.field="strand",
                                                        starts.in.df.are.0based=TRUE)
  j <- j+1
  
}

# Remove universe from list because the other loops will not need it explicitly
pk.df <- de_grange_list[['universe.peaks']]
de_grange_list[['universe.peaks']] <- NULL

# Start Permutation tests


# assign as many labels as we have tests
pt_labels <- vector()
for (i in names(de_grange_list)) {
  for (j in names(grange_list)) {
    pt_labels <- c(pt_labels,paste(i,j,sep = '_'))
  }
}
pt_labels

# create empty list append to list the results of the tests
pt <- list()
k <- 1
for (i in de_grange_list) {
  for (j in grange_list) {
    pt[[pt_labels[k]]] <- permTest(A=i,
                                   ntimes=500,
                                   randomize.function=resampleRegions,
                                   universe=pk.df,
                                   evaluate.function=numOverlaps,
                                   B=j,
                                   verbose=FALSE)
    k <- k+1
  }
}

# Display all plots

x <- c(11.71,58.00,57.00,7.11,31.38,26.90,16.90,83.00,66.00)
y <- c(0.21,0.04,0.04,0.29,0.07,0.06,0.18,0.03,0.04)

counter = 1
for (i in pt_labels) {
  plot(pt[[i]])
  text(x=x[counter],y=y[counter],i)
  counter <- counter + 1
}


