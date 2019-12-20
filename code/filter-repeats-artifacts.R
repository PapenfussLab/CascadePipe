library(reshape2)
library(GenomicRanges)
library(stringr)
library(data.table)

args=commandArgs(trailingOnly=T)

mat<-read.delim(args[1],sep="\t",header=1)

mat.cat1.variants<-unique(as.vector(mat[which(mat$category==1),"variant"]))

print(dim(mat))
#keep only variants regarded as category 1 in at least one sample. 
mat<-mat[which(mat$variant %in% mat.cat1.variants),]

print(tail(mat))
rownames(mat)<-1:nrow(mat)
print(tail(mat))
print(dim(mat))

gr_mat<-makeGRangesFromDataFrame(mat,start.field="coordinate",end.field="coordinate",ignore.strand=TRUE)

print("Reading repeatmasker file...")
repeatmasker_df<-fread(file="../external_db/hg19.fa.tab")
#repeatmasker_df<-read.delim(file="../external_db/hg19.fa.tab",sep="\t",header=F)
colnames(repeatmasker_df)<-c("score1","score2","score3","score4","chr","start","end","reminder","strand","repeat","family","score5","score6","score7","score8")

repeatmasker_df$chr<-str_replace_all(repeatmasker_df$chr,"^chr([[:alnum:]]+)","\\1")
gr_repeatmasker<-makeGRangesFromDataFrame(repeatmasker_df,ignore.strand=TRUE)


print("Reading low complexity regions file...")
lc_df<-read.delim(file="../external_db/btu356_LCR-hs37d5.bed",sep="\t",header=F)
colnames(lc_df)<-c("chr","start","end")

lc_df$chr<-str_replace_all(lc_df$chr,"^chr([[:alnum:]]+)","\\1")
lc_df$start <- lc_df$start + 1
gr_lc<-makeGRangesFromDataFrame(lc_df,ignore.strand=TRUE)


print("Reading artifactual regions file...")
art_df<-read.delim(file="../external_db/all.gene.symbols.with.coordinates.mainChrom.txt",sep="\t",header=1)
gr_art<-makeGRangesFromDataFrame(art_df,ignore.strand=TRUE)


print("Reading centromere and telomere file...")
cen_tel_df<-read.delim(file="../external_db/hg19_centromere_telomere",sep="\t",header=F)

colnames(cen_tel_df)<-c("bin","chr","start","end","ix","n","size","type","bridge")
cen_tel_df$chr<-str_replace_all(cen_tel_df$chr,"^chr([[:alnum:]]+)","\\1")
cen_tel_df$start <- cen_tel_df$start + 1
gr_cen_tel<-makeGRangesFromDataFrame(cen_tel_df,ignore.strand=TRUE)


print(gr_repeatmasker)
print(length(gr_repeatmasker))
print(gr_lc)
print(length(gr_lc))
print(gr_art)
print(length(gr_art))
print(gr_cen_tel)
print(length(gr_cen_tel))

gr_all<-c(gr_repeatmasker,gr_lc,gr_art,gr_cen_tel)

print(gr_all)
print(length(gr_all))

hits_reps<-findOverlaps(gr_all,gr_mat)
print(subjectHits(hits_reps))


mat.filtered<-mat[setdiff(rownames(mat),rownames(mat)[subjectHits(hits_reps)]),]
print(dim(mat.filtered))

write.table(x=mat.filtered,file=args[2],sep="\t",col.names=T,row.names=F,quote=F)

