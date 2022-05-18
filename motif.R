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

df[,ID_:=fifelse(ID=='.',"unknown","known")]

# remove doubtful SNPs
df = df[!(VARIANT=='SNP'&SIZE_DELS!=0),]

# Make Granges
variants = makeGRangesFromDataFrame(df,
                                    keep.extra.columns=TRUE,
                                    ignore.strand=TRUE,
                                    seqinfo=NULL,
                                    seqnames.field="CHR",
                                    start.field="START",
                                    end.field="END",
                                    starts.in.df.are.0based=FALSE)
genome(variants) <- "hg19"
variants

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