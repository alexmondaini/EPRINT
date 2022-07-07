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

# look at combinations of SNP
df[VARIANT=='SNP',SNP_subs:=paste(REF,ALT,sep = '_')]

# condense variants
df[,
   SNP_no_strand:=
     fcase(SNP_subs=='A_G'|SNP_subs=='T_C','A_G|T_C',
           SNP_subs=='G_T'|SNP_subs=='C_A','G_T|C_A',
           SNP_subs=='G_A'|SNP_subs=='C_T','G_A|C_T',
           SNP_subs=='A_C'|SNP_subs=='T_G','A_C|T_C',
           SNP_subs=='G_C'|SNP_subs=='C_G','G_C|C_G',
           SNP_subs=='A_T'|SNP_subs=='T_A','A_T|T_A')
]

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
hits <- queryHits(findOverlaps(variants,pk))
hit_df <- df[hits,]

# Find subjecthits
i <- subjectHits(findOverlaps(variants,pk))
hit_strand <- as.vector(strand(pk[i,]))
# create new column in hit_df
hit_df[,hit_strand:=hit_strand]

# Counts of intersectioned SNPs and peaks
p <- hit_df[,.N,by=.(ID_,experiment,VARIANT)]
p
# Counts of non intersectioned SNPs and peaks
df[,.N,by=.(ID_,experiment,VARIANT)]
# Percentage of intersected SNPs 
p[,N:=hit_df[,.N,by=.(ID_,experiment,VARIANT)][,N]/df[,.N,by=.(ID_,experiment,VARIANT)][,N]*100]
p

#color
cbPalette <- c("#999999", "#E69F00","#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
#Plot
ggplot(p,aes(x=experiment,y=N)) +
  geom_col(aes(fill=ID_)) +
  scale_fill_manual(values=cbPalette) +
  facet_wrap(vars(VARIANT)) +
  scale_y_continuous(breaks = seq(0,12,2)) +
  labs(y='pct of variants under peaks',fill='Variant type')

# Plot by SNP substitution type and strand

ggplot(hit_df[VARIANT=='SNP',],aes(x=experiment)) +
  geom_bar(aes(fill=ID_)) +
  scale_fill_manual(values=cbPalette) +
  facet_wrap(vars(SNP_subs,hit_strand),scale='free') +
  labs(y='Count',fill='Variant type')

# Grouping SNP subs
p <- hit_df[VARIANT=='SNP',.N,keyby=.(ID_,experiment,SNP_no_strand)]
p
p <- p[df[VARIANT=='SNP',.N,keyby=.(ID_,experiment)]]
# this is the total
df[VARIANT=='SNP',.N,keyby=.(ID_,experiment)]
p[,pct:=N/i.N*100]

# plot
ggplot(p,aes(x=experiment,y=pct)) +
  geom_col(aes(fill=ID_)) +
  scale_fill_manual(values=cbPalette) +
  facet_wrap(vars(SNP_no_strand)) +
  labs(y='Count',fill='Variant type')


