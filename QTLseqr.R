#!R 3.50

library("ggplot2")
library("QTLseqr")
library("argparse")
################################
parser = ArgumentParser()
parser$add_argument("--data", help="the SNP.table", type="character",required=TRUE)
parser$add_argument("--highbulk", help='mut bulk name',required=TRUE)
parser$add_argument("--lowbulk", help='wt bulk name', required=TRUE)
parser$add_argument("--q_value", help='plot using q value',type="double", default=0.01)
parser$add_argument("--windowsize", help='size of window',type="double", default=1e6)
parser$add_argument("--bulksize", help='size of bulk',type="integer", default=30 )
parser$add_argument("--popstruc", help='size of bulk',default="F2")
args <- parser$parse_args()
data = args$data
highbulk = args$highbulk
lowbulk = args$lowbulk
q_value = args$q_value
windowsize = args$windowsize
bulksize = args$bulksize
popstruc = args$popstruc
#import data
chromfile <- read.table('chromlist.txt',stringsAsFactors=F)
Chroms <- as.vector(chromfile[,1]) 

df <-
     importFromGATK(
     file = data,
     highBulk = highbulk,
     lowBulk = lowbulk,
     chromList = Chroms
     )
#Filter SNPs based on some criteria
ggplot(data = df) + geom_histogram(aes(x = DP.HIGH + DP.LOW)) + xlim(0,500)
ggsave( file = "Depth.png", width = 5, height = 6, type = "cairo", dpi = 600)

ggplot(data = df) + geom_histogram(aes(x = REF_FRQ))
ggsave( file = "ref_frequency.png", width = 5, height = 6, type = "cairo", dpi = 600)

ggplot(data = df) +
    geom_histogram(aes(x = SNPindex.HIGH))
ggsave( file = "High_index.png", width = 5, height = 6, type = "cairo", dpi = 600)


ggplot(data = df) +
    geom_histogram(aes(x = SNPindex.LOW))
ggsave( file = "Low_index.png", width = 5, height = 6, type = "cairo", dpi = 600)

df_filt <-
          filterSNPs(
          SNPset = df,
          refAlleleFreq = 0.20,
          minTotalDepth = 100,
          maxTotalDepth = 400,
          minSampleDepth = 40,
          minGQ = 99
          )

#Run G' analysis
df_filt <- runGprimeAnalysis(SNPset = df_filt,
          windowSize = windowsize,
          outlierFilter = "deltaSNP")

#Run QTLseq analysis
df_filt <- runQTLseqAnalysis(
          SNPset = df_filt,
          windowSize = windowsize,
          popStruc = popstruc,
          bulkSize = c(bulksize,bulksize),
          replications = 10000,
          intervals = c(95, 99)
          )
#G’ distribution plots

png("hampel.png",width=2800,height=600)
plotGprimeDist(SNPset = df_filt, outlierFilter = "Hampel")
dev.off()

png("delta_SNP.png",width=2800,height=600)
plotGprimeDist(SNPset =df_filt, outlierFilter = "deltaSNP", filterThreshold = 0.1)
dev.off()

#QTL analysis plots
png("SNPs.png",width=2800,height=600)
plotQTLStats(SNPset = df_filt, var = "nSNPs")
dev.off()

png("DeltaSNP.png",width=2800,height=600)
plotQTLStats(SNPset = df_filt, var = "deltaSNP", plotIntervals = TRUE)
dev.off()

png("Gprime.png",width=2800,height=600)
plotQTLStats(SNPset = df_filt, var = "Gprime", plotThreshold = TRUE, q = q_value)
dev.off()

png("negLog10Pval.png",width=2800,height=600)
plotQTLStats(SNPset = df_filt, var = "negLog10Pval", plotThreshold = TRUE, q = q_value)
dev.off()

#Extracting data
#QTL <- getSigRegions(SNPset = df_fit, alpha = 0.01)
#write.table(QTL,"QTL.table")

Gprime_results <- getQTLTable(SNPset = df_filt, method = "Gprime",alpha = 0.01, export = TRUE)
QTLseq_results <- getQTLTable(SNPset = df_filt, method = "QTLseq",alpha = 0.01, export = TRUE)
