library(GenomicRanges)
library(stringr)


args=commandArgs(trailingOnly=T)

mat<-read.delim(args[1],sep="\t",header=1)

gr_mat<-makeGRangesFromDataFrame(mat,start.field="coordinate",end.field="coordinate",ignore.strand=TRUE)

bed_file<-read.delim(args[2],header=0,sep="\t")
colnames(bed_file)<-c("chr","start","end")

print(head(bed_file))

bed_file$start <- bed_file$start + 1

gr_BED<-makeGRangesFromDataFrame(bed_file,ignore.strand=TRUE)

print(gr_BED)
print(gr_mat)

hits<-findOverlaps(gr_BED,gr_mat)

print(hits)

write.table(mat[rownames(mat)[subjectHits(hits)],],file=args[3],sep="\t",quote=FALSE,row.names = F,col.names=T)

