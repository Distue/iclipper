# R file to generate visual output from iclipper outputs
# Author: Thomas Schwarzl <schwarzl@embl.de>
# Licence: MIT

suppressMessages(require(ggplot2))
suppressMessages(require(pheatmap))
suppressMessages(require(seqLogo))

# --------------------------------------------------------------
# Functions to count the GC Content
# --------------------------------------------------------------
getGCcontent <-  function(seq){
   return(length(gregexpr('[GCgc]', seq)[[1]])) #/nchar(seq))
}

# --------------------------------------------------------------
# Run
# --------------------------------------------------------------



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Adapters
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

adapters <- read.delim("data_adapters.txt", sep="", header=T)

frequency <- colSums(adapters)
sum(frequency)

# remove the unknown columns and rows
adapters <- adapters[rownames(adapters) != "unknown",colnames(adapters) != "unknown" ]

# calculates sums and gc contents
adapters.df <- data.frame(row.names = rownames(adapters))
adapters.df[,'rowSum'] <- rowSums(adapters)
adapters.df[,'gcContent'] <- unlist(lapply(rownames(adapters), getGCcontent))
adapters.table <- as.data.frame(table(adapters.df['rowSum']))
adapters.df[,'freq'] <- unlist(lapply(adapters.df[,'rowSum'], function(x) adapters.table[adapters.table[,1] == x,2]))
adapters.df <- adapters.df[order(adapters.df['rowSum'], decreasing = T),]

# calculates frequencies
adapters.fd <- data.frame(row.names = substr(colnames(adapters), 2, 1e9))
adapters.fd[,'length'] <- rownames(adapters.fd)
adapters.fd[,'freq'] <- colSums(adapters)
adapters.fd[,'gc'] <- unlist(lapply(colnames(adapters), function(ylength) {
   gclist <- unlist(lapply(rownames(adapters), function(xadapters) {
      val <- adapters[xadapters, ylength]
      if (val > 0) {
         return(rep(adapters.df[xadapters, 'gcContent'], val))
      } else {
         return(NA)
      }
   }))

   return(median(gclist, na.rm = T))
}))

# Plot the histogram for adapter distribution
ggplot(adapters.df, aes(x=rowSum)) +
   # Histogram with density instead of count on y-axis
   geom_histogram(aes(y=..density..),
                  binwidth=.5,
                  colour="black", fill="white") +
   # Overlay with transparent density plot
   geom_density(alpha=.2, fill="#FF6666")  +
   ggtitle("Distribution of adapter usage with density curve")


# Plot the histogram for adapter distribution
ggplot(adapters.df, aes(x=rowSum)) +
   # Histogram
   geom_histogram(binwidth=1, colour="black", fill="white") +
   ggtitle("Histogram of adapter usage")

# Make each dot partially transparent, with 1/4 opacity
# For heavy overplotting, try using smaller values
ggplot(adapters.df, aes(x=gcContent, y=rowSum, size=freq)) +
   geom_point(shape=19) +      # Use solid circles
   #              alpha=1/4)    +  # 1/4 opacity
   scale_size_continuous(range = c(2,15)) +
   ggtitle("Frequency of Adapter Occurances=") +
   xlab("GC nucleotides") +
   ylab("# reads")

# Jitter the points
ggplot(adapters.df, aes(x=gcContent, y=rowSum))+
   geom_point(shape=1,      # Use hollow circles
              position=position_jitter(width=0.1,height=.1))


# plot the GC content versus the occurance, aggregated
adapters.agg <- aggregate(adapters.df$gcContent, by=list(adapters.df$rowSum), FUN=median)
colnames(adapters.agg) <- c("occurance", "medianGC")

ggplot(adapters.agg, aes(x=occurance, y=medianGC)) +
   geom_point(shape=1)  + geom_smooth() +
   ggtitle("Adapter occurance")


# Plot the histogram for adapter distribution
ggplot(adapters.fd, aes(x=length, y=freq, size=gc)) +
   geom_point(shape=19) +
   scale_size_continuous(range = c(2,15)) +
   xlab("Read length") +
   ylab("Occurance")

ggplot(adapters.fd, aes(x=length, y=gc, size=freq)) +
   stat_boxplot(geom ='errorbar') +
   geom_boxplot() + # shorthand for  stat_boxplot(geom='boxplot')
   scale_size_continuous(range = c(0.01,4))

tiff("adapters-gc-rowSum.tiff")
dev.off()


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Read distribution
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
read.dist <- as.data.frame(colSums(adapters))
read.dist[,"length"] <- as.numeric(sub('X', '', rownames(read.dist)))
colnames(read.dist)[1] <- c("reads")
read.dist[,'logreads'] <- log(read.dist[,'reads'], base=10)
read.dist[,'perc']  <- read.dist[,'reads']/sum(read.dist[,'reads']) * 100


ggplot(read.dist, aes(x=length , y=reads, group=1)) +
  #   geom_line() +
  # geom_point(size=3, fill="white") +
   xlab("Read length") + ylab("log10 Read Count") +
   ggtitle("Read distribution over length") +
   geom_area() +
   coord_trans(ytrans = 'log10', limx = c(floor(min(read.dist[,'length'])),ceiling(max(read.dist[,'length']))), limy =c(floor(min(read.dist[,'logreads'])),ceiling(max(read.dist[,'logreads'])))) +
   annotation_logticks(scaled=FALSE) +
   scale_x_continuous(breaks = trans_breaks("log10", function(x) 10^x),
                      labels = trans_format("log10", math_format(10^.x))) +
   scale_y_continuous(breaks = trans_breaks("log10", function(x) 10^x),
                      labels = trans_format("log10", math_format(10^.x)))


ggplot(read.dist, aes(x=length , y=logreads, group=1)) +
   geom_line() +
   geom_point(size=3, fill="white")


ggplot(read.dist, aes(x=length , y=perc, group=1)) +
   geom_line() +
   geom_point(size=3, fill="white")



ggplot(read.dist, aes(x=length , y=logreads, group=1)) +
   geom_area() +
   geom_point(size=3, fill="white")


write.table(read.dist, file="read.dist.txt", sep="\t")



  # ylim(c(floor(min(read.dist[,'logreads'])),ceiling(max(read.dist[,'logreads'])))) +
   # geom_area()

   # geom_line() +
   # geom_point(size=3, fill="white") +


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Frequency of deletions
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

deletions <- t(read.delim("data_deletions.txt", sep="", header=T))
deletions <- deletions[rownames(deletions) != "unknown" ,colnames(deletions) != "unknown"]

# Replace
#is.nan.data.frame <- function(x)
#   do.call(cbind, lapply(x, is.nan))
#deletions[is.nan(deletions)] <- 0

pheatmap(log(deletions + 1),
         cluster_rows = F,
         cluster_cols = F,
         main = "Deletion Frequency",
         fontsize_rows = 20,
         fontsize_cols = 20,
         xlab = "log deletion frequency",
         ylab = "read length")

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Deletion sites
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

deletion.sites <- t(read.delim("data_deletion-sites.txt", sep="", header=T))
deletion.sites <- deletion.sites[rownames(deletion.sites) != "unknown", colnames(deletion.sites) != "unknown"]

pheatmap(log(deletion.sites + 1),
         cluster_rows = F,
         cluster_cols = F,
         main = "Deletion Sites",
         fontsize_rows = 20,
         fontsize_cols = 20,
       #  scale="row",
         xlab = "deletion sites",
         ylab = "read length")




deletion.sites.sum <- as.data.frame(rowSums(t(deletion.sites)) / sum(deletion.sites) )
deletion.sites.sum[,"position"] <- as.numeric(rownames(deletion.sites.sum))
colnames(deletion.sites.sum)[1] <- c("freq")


ggplot(deletion.sites.sum, aes(x=position, y=freq) ) +
   geom_area(aes(colour = "red", fill= "red")) +
  # geom_ribbon(data=deletion.sites.sum, ) +
   ggtitle("Frequency of deletions") +

   geom_area(aes(colour = "red", fill= "red"))
   #  geom_ribbon(data=deletion.sites.sum, )

# geom_point(size=1)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Deletions per Base
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
deletion.bases <- read.delim("data_deletion-sites-perbase.txt", sep="", header=T)
deletion.bases <- deletion.bases[rownames(deletion.bases) != "unknown", colnames(deletion.bases) != "unknown"]
deletion.bases.ratio <- sweep(deletion.bases,2,colSums(deletion.bases),`/`)

pheatmap(deletion.bases,
         cluster_rows = F,
         cluster_cols = F,
         main = "Deletion Sites per Base",
         fontsize_rows = 20,
         fontsize_cols = 20,
         xlab = "bases",
         ylab = "read length",
         scale = "column",
         legend = F)

seqLogo(deletion.bases.ratio, ic.scale=FALSE)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Frequency of insertions
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

insertions <- t(read.delim("data_insertions.txt", sep="", header=T))
insertions <- insertions[rownames(insertions) != "unknown" ,colnames(insertions) != "unknown"]

# Replace
#is.nan.data.frame <- function(x)
#   do.call(cbind, lapply(x, is.nan))
#deletions[is.nan(deletions)] <- 0

pheatmap(log(insertions + 1),
         cluster_rows = F,
         cluster_cols = F,
         main = "Insertions Frequency",
         fontsize_rows = 20,
         fontsize_cols = 20,
         xlab = "log insertions frequency",
         ylab = "read length")


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Insertion sites
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

insertion.sites <- t(read.delim("data_insertion-sites.txt", sep="", header=T))
insertion.sites <- insertion.sites[rownames(insertion.sites) != "unknown", colnames(insertion.sites) != "unknown"]

pheatmap(log(insertion.sites + 1),
         cluster_rows = F,
         cluster_cols = F,
         main = "Insertion Sites",
         fontsize_rows = 20,
         fontsize_cols = 20,
         xlab = "insertion position",
         ylab = "read length")


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Insertion per Base
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
insertion.bases <- read.delim("data_insertion-sitesperbase.txt", sep="", header=T)
insertion.bases <- insertion.bases[rownames(insertion.bases) != "unknown", colnames(insertion.bases) != "unknown"]

insertion.bases.ratio <- sweep(insertion.bases,2,colSums(insertion.bases),`/`)

pheatmap(insertion.bases,
         cluster_rows = F,
         cluster_cols = F,
         main = "Insertion Sites per Base",
         fontsize_rows = 20,
         fontsize_cols = 20,
         xlab = "bases",
         ylab = "read length",
         scale = "column",
         legend = F)

seqLogo(insertion.bases.ratio, ic.scale=FALSE)



intron.exon.to.deletion <- read.delim("data_intron-exon-to-deletion.txt", sep="", header=T)
intron.exon.to.deletion <- intron.exon.to.deletion[rownames(intron.exon.to.deletion) != "unknown", colnames(intron.exon.to.deletion) != "unknown"]
rownames(intron.exon.to.deletion) <- as.numeric(rownames(intron.exon.to.deletion))
intron.exon.to.deletion <- intron.exon.to.deletion[order(as.numeric(rownames(intron.exon.to.deletion))),]
intron.exon.to.deletion <- t(intron.exon.to.deletion)
mat2 <- intron.exon.to.deletion / rowSums(intron.exon.to.deletion) * 100
colnames(mat2) <- paste("X", colnames(mat2), sep="")
png("data_intron-exon-to-deletion.png")
pheatmap(mat2,
         cluster_rows = F,
         cluster_cols = F,
         main = "Closest Deletion to Intron-Exon Junction",
         xlab = "distance to closest deletion",
         ylab = "read length",
         scale = "none")
dev.off()


intron.exon.to.deletion.all <- colSums(intron.exon.to.deletion)
intron.exon.to.deletion.all.df <- as.data.frame(cbind(as.numeric(names(intron.exon.to.deletion.all)), as.numeric(intron.exon.to.deletion.all)))
colnames(intron.exon.to.deletion.all.df) <- c("dist", "del")
ggplot(intron.exon.to.deletion.all.df, aes(x=dist, y=del) ) +
   geom_area(aes(colour = "red", fill= "red")) +
 #  geom_ribbon(data=deletion.sites.sum, ) +
   ggtitle("Deletions to Intron-Exon Junctions")
  # geom_point(size=1)





exon.intron.to.deletion <- read.delim("data_exon-intron-to-deletion.txt", sep="", header=T)
exon.intron.to.deletion <- exon.intron.to.deletion[rownames(exon.intron.to.deletion) != "unknown", colnames(exon.intron.to.deletion) != "unknown"]
rownames(exon.intron.to.deletion) <- as.numeric(rownames(exon.intron.to.deletion))
exon.intron.to.deletion <- exon.intron.to.deletion[order(as.numeric(rownames(exon.intron.to.deletion))),]
exon.intron.to.deletion <- t(exon.intron.to.deletion)

exon.intron.to.deletion.all <- colSums(exon.intron.to.deletion)
exon.intron.to.deletion.all.df <- as.data.frame(cbind(as.numeric(names(exon.intron.to.deletion.all)), as.numeric(exon.intron.to.deletion.all)))
colnames(exon.intron.to.deletion.all.df) <- c("dist", "del")


ggplot(exon.intron.to.deletion.all.df, aes(x=dist, y=del) ) +
   geom_area(aes(colour = "red", fill= "red")) +
   xlim(c(-100,100)) +
   #  geom_ribbon(data=deletion.sites.sum, ) +
   ggtitle("Deletions to Exon-Intron Junctions")
  #geom_point(size=1)



exon.intron.to.deletion.all.df

exon.intron.to.deletion.all.df



exon.intron.to.deletion.all.df <- exon.intron.to.deletion.all.df[exon.intron.to.deletion.all.df[,1] >= -100 & exon.intron.to.deletion.all.df[,1] < 0,]
intron.exon.to.deletion.all.df <- intron.exon.to.deletion.all.df[intron.exon.to.deletion.all.df[,1] <=  100 & intron.exon.to.deletion.all.df[,1] > 0,]

mat <- rbind(exon.intron.to.deletion.all.df, intron.exon.to.deletion.all.df)
ggplot(mat, aes(x=dist, y=del) ) +
   geom_area(aes(colour = "red", fill= "red")) +
   #  geom_ribbon(data=deletion.sites.sum, ) +
   ggtitle("Deletions to Exon-Exon Junctions")
#geom_point(size=1)




a <- exon.intron.to.deletion[,as.numeric(colnames(exon.intron.to.deletion)) < 0 &
                                as.numeric(colnames(exon.intron.to.deletion)) > -100]

pheatmap(a[!rownames(a) %in% "X42",],
         cluster_rows = F,
         cluster_cols = F,
         main = "Exon left",
         xlab = "distance to closest deletion",
         ylab = "read length",
         scale = "none")






# --------

intron.exon.to.start <- read.delim("data_intron-exon-to-start.txt", sep="", header=T)
intron.exon.to.start <- intron.exon.to.start[rownames(intron.exon.to.start) != "unknown", colnames(intron.exon.to.start) != "unknown"]
rownames(intron.exon.to.start) <- as.numeric(rownames(intron.exon.to.start))
intron.exon.to.start <- intron.exon.to.start[order(as.numeric(rownames(intron.exon.to.start))),]
intron.exon.to.start <- t(intron.exon.to.start)
rownames(intron.exon.to.start) <- as.numeric(sub("X", "", rownames(intron.exon.to.start)))
a <- intron.exon.to.start / rowSums(intron.exon.to.start) * 100

pheatmap(a[rownames(a) > 18 & rownames(a) < 38, as.numeric(colnames(a)) > -50 & as.numeric(colnames(a)) < 50  ] ,
         cluster_rows = F,
         cluster_cols = F,
         main = "Intron-Exon to Start",
         xlab = "distance to closest start",
         ylab = "read length",
         scale = "none")

#!rownames(a) %in% c("X15", "X16", "X40", "X42")


#-------

readcount <- read.delim("data_readcount.txt", sep="", header=T)
readcount.cum <- cumsum(readcount)
readcount.cum <- readcount.cum/max(readcount.cum)
readcount.cum[,"length"] <- as.numeric(rownames(readcount.cum))

g <- ggplot(readcount.cum, aes(x=length, y=counts) ) +
   #geom_area(aes(colour = "red", fill= "red", guide = FALSE)) +
   geom_line(aes(group=1)) +
   #geom_point(size=2) +
   ylim(c(0,1)) +
   ylab("") +
   xlab("read length") +
   ggtitle("Read Counts") +
   theme_bw() + # black and white theme
   theme(axis.text  = element_text(size = 15),
         axis.title = element_text(size = 16),
         plot.title = element_text(size = 20)) +
   guides(colour = FALSE)

png(filename = "plot_readcounts-cummulative.png", width = 300, height = 300, units = "px")
g
dev.off()

pdf("plot_readcounts-cummulative.pdf")
g
dev.off()

tiff(filename = "plot_readcounts-cummulative.tiff", width = 300, height = 300, units = "px")
g
dev.off()



# CITS / CIMS
cims.dist.cits.to.del <- read.table("bed/cims-dist-CITS-to-deletion-n10.p0.05.bed", header = F)
colnames(cims.dist.cits.to.del) <- c( "dist" )

plot(density(cims.dist.cits.to.del[cims.dist.cits.to.del[,1] < 50,1]))
hist(cims.dist.cits.to.del[cims.dist.cits.to.del[,1] < 50,1], breaks=20)

sum(cims.dist.cits.to.del[,1] > 10)
sum(cims.dist.cits.to.del[,1] <= 10)



cims.dist.cits.to.del.cum <- cumsum(cims.dist.cits.to.del)
cims.dist.cits.to.del.cum <- cims.dist.cits.to.del.cum/max(cims.dist.cits.to.del.cum)
cims.dist.cits.to.del.cum[,"count"] <- as.numeric(rownames(cims.dist.cits.to.del.cum))


g <- ggplot(cims.dist.cits.to.del.cum, aes(x=count, y=dist) ) +
   #geom_area(aes(colour = "red", fill= "red", guide = FALSE)) +
   geom_line(aes(group=1)) +
   #geom_point(size=2) +
   #ylim(c(0,1)) +
   ylab("Cummulative distance to truncation site") +
   xlab("Deletions") +
   ggtitle("Cumm Dist") +
   theme_bw() + # black and white theme
   theme(axis.text  = element_text(size = 15),
         axis.title = element_text(size = 16),
         plot.title = element_text(size = 20)) +
   guides(colour = FALSE)
g
