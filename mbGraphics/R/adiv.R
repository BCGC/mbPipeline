#' Plots alpha diversity calculations
#' 
#' This function generates files, plotting each independent variable,
#' whether it be categorical or continuous, versus the calculated
#' inverse simpson diversity index. 
#'
#' @param meta.frame A data frame containing all data. Two columns must be titled 'adiv' and
#' 'nseqs' and must contain the alpha diversity values and number of sequences respectively. All
#' other information be metadata or categorical/continuous variables that describe the samples.
#' @param data.types A list containing text indicators of each meta.frame column data type
#' in order. Input text for each list element are 'cat', 'cont', and 'meta'.
#'
#' @return No output is returned. Graph file generated: "AphaDiersity.pdf"
#'
#' @keywords keywords
#'
#' @export
#' 
#' @examples
#' meta.frame <- read.table("IECG57L01.graphics_data.txt", skip=1, header=TRUE, stringsAsFactors=FALSE)
#' data.types <- strsplit(readLines("IECG57L01.graphics_data.txt")[1], '\t')[[1]]
#' adiv(meta.frame, data.types)
adiv <- function(meta.frame, data.types)
{
  library(gdata)
  library(MASS)
    pdf(file="AlphaDiversity.pdf")

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

                sub <- meta.frame[index,]
                #adiv is now extrated for all of groups[j]
                max <- max(sub$adiv)
                min <- min(sub$adiv)

                if(j == 1)
                {
                    if(xcounter > 0)
                    {
                        par(new='T')#determines whether to start a new graph
                    }
                    xcounter <- xcounter + 1
                    plot(rep(xcounter, length(sub$adiv)), sub$adiv, pch='-', cex=2, xlim=c(0, xbound), ylim=c(min(meta.frame$adiv), max(meta.frame$adiv)), xaxt='n', xlab = '', ylab = 'Inverse Simpson Diversity Index', main = 'Alpha Diversity per Sample', bty='n')
                    lines(c(xcounter, xcounter), c(min, max), lwd = 5, col = rgb(0,0,0,.7))
                    mtext(groups[j], side=1, at=xcounter)
                    par(new='F')
                }
                else
                {
                    xcounter <- xcounter + 1
                    par(new='T')
                    plot(rep(xcounter, length(sub$adiv)), sub$adiv, pch='-', cex=2, xlim=c(0, xbound), ylim=c(min(meta.frame$adiv), max(meta.frame$adiv)), xaxt='n', xlab = '', ylab = '', main = '', bty='n')
                    lines(c(xcounter, xcounter), c(min, max), lwd = 5, col = rgb(0,0,0,.7))
                    mtext(groups[j], side=1, at=xcounter)
                    par(new='F')
                }
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
            plot(column, meta.frame$adiv, xlab = colnames(meta.frame)[i], ylab = 'Inverse Simpson Diversity Index', main = 'Alpha Diversity per Sample')
        }
    }

    dev.off()
}

#meta.frame <- read.table("IECG57L01.graphics_data.txt", skip=1, header=TRUE, stringsAsFactors=FALSE)
#data.types <- strsplit(readLines("IECG57L01.graphics_data.txt")[1], '\t')[[1]]

#adiv(meta.frame, data.types)
