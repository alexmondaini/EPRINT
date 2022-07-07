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

# Find overlaps
hits <- queryHits(findOverlaps(variants,pk))
hit_df <- df[hits,]

# Find subjecthits
i <- subjectHits(findOverlaps(variants,pk))
hit_strand <- as.vector(strand(pk[i,]))
# create new column in hit_df
hit_df[,hit_strand:=hit_strand]
# condense variants
hit_df[,fcase(SNP_subs=='')]

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
#plot
ggplot(p,aes(x=experiment,y=N)) +
  geom_bar(aes(fill=ID_),stat='identity') +
  scale_fill_manual(values=cbPalette) +
  facet_wrap(vars(VARIANT)) +
  scale_y_continuous(breaks = seq(0,12,2)) +
  labs(y='pct of variants under peaks',fill='Variant type')

# Plot SNPs known and unknown - old
ggplot(hit_df[SIZE_DELS==0,],aes(y=SIZE_DELS/df[,.N,by=.(ID_,experiment)][,N]*100)) +
  geom_bar(aes(fill=ID_),position = position_stack(reverse = TRUE)) +
  facet_wrap(vars(experiment)) +
  theme(legend.position = "top",axis.text.y = element_blank()) +
  labs(y='SNP',fill='SNP_type')

# Plot by SNP substitution type
p <- hit_df[VARIANT=='SNP',.N,keyby=.(ID_,experiment,SNP_subs)]
p
df[SIZE_DELS==0,.N,keyby=.(ID_,experiment,SNP_subs)]
p[,N:=hit_df[SIZE_DELS==0,.N,keyby=.(ID_,experiment,SNP_subs)][,N]/df[SIZE_DELS==0,.N,keyby=.(ID_,experiment,SNP_subs)][,N]*100]
p

# Plot SNPs by experiment and substitution type
ggplot(p,aes(x=experiment,y=N)) +
  geom_bar(aes(fill=ID_),stat = 'identity') +
  facet_wrap(vars(SNP_subs)) +
  scale_fill_manual(values=cbPalette) +
  labs(y='pct of SNP type',fill='SNP type')

# Plot SNPs by experiment and substitution type
ggplot(hit_df[SIZE_DELS==0,],aes(y=SIZE_DELS)) +
  geom_bar(aes(fill=ID_),position = position_stack(reverse = TRUE)) +
  facet_wrap(vars(experiment,SNP_subs)) +
  theme(legend.position = "top",axis.text.y = element_blank()) +
  labs(y='SNP',fill='SNP_type')

# Plot DELs < 10
p <- hit_df[SIZE_DELS>=1&SIZE_DELS<=10,.N,keyby=.(experiment,ID_)]
p[,N:=hit_df[SIZE_DELS>=1&SIZE_DELS<=10,.N,keyby=.(experiment,ID_)][,N]/df[SIZE_DELS>=1&SIZE_DELS<=10,.N,keyby=.(experiment,ID_)][,N]*100]
p

# Plot DELs by type
ggplot(p,aes(x=experiment,y=N)) +
  geom_bar(aes(fill=ID_),stat='identity') +
  labs(y='pct of DELs under peaks',fill='SNP type') +
  scale_fill_manual(values = cbPalette)

# Plot DELs by type
ggplot(hit_df[SIZE_DELS>=1&SIZE_DELS<=10,],aes(x=SIZE_DELS)) +
  geom_histogram(aes(fill=ID_),binwidth = 1) +
  facet_wrap(vars(experiment)) +
  scale_x_continuous(breaks = seq(1,10,1)) +
  labs(fill='DEL_type')

# Plot DELs > 10
p <- hit_df[SIZE_DELS>=10,.N,keyby=.(experiment,ID_)]
p[,N:=hit_df[SIZE_DELS>=10,.N,keyby=.(experiment,ID_)][,N]/df[SIZE_DELS>=10,.N,keyby=.(experiment,ID_)][,N]*100]
p

# Plot Dels > 10
ggplot(p,aes(x=experiment,y=N)) +
  geom_bar(aes(fill=ID_),stat='identity') +
  labs(y='pct of DELs under peaks',fill='SNP type') +
  scale_fill_manual(values = cbPalette)

# Plot DELs >10 
ggplot(hit_df[SIZE_DELS>=10,],aes(x=SIZE_DELS)) +
  geom_histogram(aes(fill=ID_),binwidth = 10) +
  facet_wrap(vars(experiment)) +
  labs(fill='DEL_type')

# Overall representation
p <- hit_df[,.N,keyby=.(experiment,ID_,VARIANT)]
p[,N:=hit_df[,.N,keyby=.(experiment,ID_,VARIANT)][,N]/df[,.N,keyby=.(experiment,ID_,VARIANT)][,N]*100]

# Overall representation
ggplot(p,aes(x=experiment,y=N)) +
  geom_bar(aes(fill=ID_),stat='identity') +
  facet_wrap(vars(VARIANT)) +
  scale_fill_manual(values = cbPalette) +
  labs(y='pct of variants peaks under',fill='Variant type')

# Overall representation
ggplot(hit_df,aes(SIZE_DELS)) +
  geom_histogram(binwidth = 1,fill="white",color=4) + 
  facet_wrap(vars(experiment,VARIANT,ID_),scales = 'free',nrow = 2)

# Number of overlaps for eprint and input
hit_df[experiment=="input",.N]
hit_df[experiment=="eprint",.N]

# Subsetbyoverlaps both ways
variants_sbset <- subsetByOverlaps(variants,pk)
pk_sbset <- subsetByOverlaps(pk,variants)
# Get start point of eprint peak subset
strt_pk = GRanges(seqnames = seqnames(pk_sbset),ranges = start(pk_sbset),strand = strand(pk_sbset))


# Distance to nearest
dist <- distanceToNearest(variants_sbset,strt_pk)
mcols(dist) <-  c(mcols(dist),mcols(variants_sbset))


# Plot by experiment vartiant_type and ID
dis = setDT(data.frame(dist))
plt = ggplot(dis)

# dist to start of peak
plt +
  geom_bar(aes(distance,fill=ID_)) +
  scale_fill_manual(values = cbPalette) +
  scale_x_binned(breaks = seq(0,412,by = 10)) + 
  facet_wrap(vars(VARIANT,experiment),scales='free_y') + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(fill='Variant type')
#old
plt + 
  geom_histogram(binwidth = 10,aes(x=distance,fill=..x..)) + 
  facet_wrap(vars(experiment,ID_,VARIANT),ncol = 2) +
  scale_color_continuous() +
  scale_x_continuous(breaks = seq(0,200,by = 10)) +
  labs(title='Distance of SNPs in base pair to eprint peak midpoint') +
  guides(fill="none")

# PLot SNP subtypes
plt +
  geom_bar(data=dis[VARIANT=='SNP'],aes(distance,fill=ID_)) +
  scale_x_binned(breaks = seq(0,412,by = 30)) + 
  scale_fill_manual(values = cbPalette) +
  facet_wrap(vars(experiment,SNP_subs),scales='free_y') + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(fill='SNP type')

# old
plt +
  geom_histogram(binwidth = 10,aes(x=distance,fill=..x..)) + 
  facet_wrap(vars(experiment,SNP_subs)) +
  scale_color_continuous() +
  scale_x_continuous(breaks = seq(0,200,by = 10)) +
  labs(title='Distance of SNPs in base pair to eprint peak start') +
  guides(fill="none") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# The whole plot
dist_tot <- distanceToNearest(variants[seqnames(variants)!='chrY'],pk)
mcols(dist_tot) <- c(mcols(dist_tot),mcols(variants[seqnames(variants)!='chrY']))

# Plot by experiment vartiant_type and ID
dis = setDT(data.frame(dist_tot))
plt = ggplot(dis)

# dist to start of peak
plt +
  geom_boxplot(data=dis[VARIANT=='SNP'],aes(x=experiment,y=distance,fill=ID_)) +
  scale_y_continuous(limits = c(0,10000),breaks = seq(0,10000,1000)) +
  scale_fill_manual(values = cbPalette) +
  labs(fill='SNP type')
