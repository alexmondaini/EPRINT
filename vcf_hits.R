#
library(GenomicRanges)
library(data.table)
library(ggplot2)
#

# Load
df = fread(file = '/my_dir/cromwell-executions/SelectVariants/output/combined_vcfs.txt')

df[,experiment:=fifelse(SAMPLE%in%c('10249_sample1','10249_sample1_variant_filtered',
                                  '10249_sample2','10249_sample2_variant_filtered',
                                  '10249_sample3','10249_sample3_variant_filtered',
                                  '10249_sample4','10249_sample4_variant_filtered'),
                        "eprint",
                        "input")]

# Make Granges
gr = makeGRangesFromDataFrame(df,
                         keep.extra.columns=TRUE,
                         ignore.strand=TRUE,
                         seqinfo=NULL,
                         seqnames.field="CHR",
                         start.field="START",
                         end.field="END",
                         starts.in.df.are.0based=FALSE)
genome(gr) <- "hg19"
gr

# Load peaks
pk.df = fread("R_eprint_project/files/fus.merge.peak.counts.inputfiltered.fusvsneg.deseq.sfglobal.batch1.wald.txt", sep="\t", header=T, quote="")
# drop last column, it's a duplicate of name
pk.df <- pk.df[,-32]
# paste chr in front of geneChr
pk.df[,geneChr:=paste0('chr',geneChr)]
# change geneChr 23 to X
pk.df[geneChr=='chr23',geneChr:='chrX']
# change strand encode
pk.df[,geneStrand := fcase(geneStrand==1,'+',
                           geneStrand==2,'-')]
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
genome(pk) <- "hg19"
pk

# Find overlaps
hits <- queryHits(findOverlaps(gr,pk))
hit_df <- df[hits,]

# Plot
ggplot(hit_df[SIZE_DELS<10,],aes(SIZE_DELS)) + geom_histogram(binwidth = 1,fill="white",color=4) + facet_wrap(vars(experiment))+
  scale_x_continuous(breaks = seq(1,10,1))

ggplot(hit_df[SIZE_DELS>10,],aes(SIZE_DELS)) + geom_histogram(binwidth = 10,fill="white",color=4) + facet_wrap(vars(experiment)) +
  scale_x_continuous(breaks = seq(10,200,by = 10)) 

# Number of overlaps for eprint and input
hit_df[experiment=="input",.N]
hit_df[experiment=="eprint",.N]

# Around 2 times more deletion overlaps in eprint 

