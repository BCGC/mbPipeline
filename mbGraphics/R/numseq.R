#' Plots number of sequences
#' 
#' This function generates files, plotting each independent variable,
#' whether it be categorical or continuous, versus the calculated
#' number of sequences. 
#'
#' @param meta.frame A data frame containing all data. Two columns must be titled 'adiv' and
#' 'nseqs' and must contain the alpha diversity values and number of sequences respectively. All
#' other information be metadata or categorical/continuous variables that describe the samples.
#' @param data.types A list containing text indicators of each meta.frame column data type
#' in order. Input text for each list element are 'cat', 'cont', and 'meta'.
#'
#' @return No output is returned. Graph file generated: "NumSequences.pdf"
#'
#' @keywords keywords
#'
#' @export
#' 
#' @examples
#' meta.frame <- read.table("IECG57L01.graphics_data.txt", skip=1, header=TRUE, stringsAsFactors=FALSE)
#' data.types <- strsplit(readLines("IECG57L01.graphics_data.txt")[1], '\t')[[1]]
#' numseq(meta.frame, data.types)
numseq <- function(meta.frame, data.types)
{
  library(gdata)
  library(MASS)
    pdf(file="NumSequences.pdf")

    xbound <- 0

    for(i in 1:length(data.types)){
        if (data.types[i] == 'cat')
        {
            xbound <- xbound + (length(unique(meta.frame[,i]))) + 2
        }
    }

    xcounter <- -3

    for(i in 1:length(data.types))
    {
        column <- meta.frame[,i]

        if (data.types[i] == "cat") #grabs categorical info rows
        {
            xcounter <- xcounter + 3
            groups <- unique(column) #grabs different types of the category

            xcounter.beforeplot <- xcounter+1

            for(j in 1:length(groups))
            {

                index <- which(meta.frame[,i] == groups[j])

                attach(meta.frame[index,])
                #nseqs is now extrated for all of groups[j]
                max <- max(nseqs)
                min <- min(nseqs)

                if(j == 1)
                {
                    if(xcounter > 0)
                    {
                        par(new='T')#determines whether to start a new graph
                    }
                    xcounter <- xcounter + 1
                    plot(rep(xcounter, length(nseqs)), nseqs, pch='-', cex=2, xlim=c(0, xbound), ylim=c(min(meta.frame$nseqs), max(meta.frame$nseqs)), xaxt='n', xlab = '', ylab = 'n Sequences', main = 'Number of Sequences per Sample', bty='n')
                    lines(c(xcounter, xcounter), c(min, max), lwd = 5, col = rgb(0,0,0,.7))
                    mtext(groups[j], side=1, at=xcounter)
                    par(new='F')
                }
                else
                {
                    xcounter <- xcounter + 1
                    par(new='T')
                    plot(rep(xcounter, length(nseqs)), nseqs, pch='-', cex=2, xlim=c(0, xbound), ylim=c(min(meta.frame$nseqs), max(meta.frame$nseqs)), xaxt='n', xlab = '', ylab = '', main = '', bty='n')
                    lines(c(xcounter, xcounter), c(min, max), lwd = 5, col = rgb(0,0,0,.7))
                    mtext(groups[j], side=1, at=xcounter)
                    par(new='F')
                }

                detach()
            }
            #following plots the group label in the center of the category section
            lines(c(xcounter.beforeplot, xcounter), rep(par('usr')[3], 2), lwd = 2, xpd = TRUE)
            mtext(colnames(meta.frame)[i], side = 1, line = 2, at=(xcounter.beforeplot+xcounter)/2)
        }
    }

    for(i in 1:length(data.types))
    {
        if(data.types[i] == "cont") #to plot data (continuous format)
        {
            column <- meta.frame[,i]
            plot(column, meta.frame$nseqs, xlab = colnames(meta.frame)[i], ylab = 'n Sequences', main = 'Number of Sequences per Sample')
        }
    }

    dev.off()
}

#meta.frame <- read.table("IECG57L01.graphics_data.txt", skip=1, header=TRUE, stringsAsFactors=FALSE)
#data.types <- strsplit(readLines("IECG57L01.graphics_data.txt")[1], '\t')[[1]]

#numseq(meta.frame, data.types)
