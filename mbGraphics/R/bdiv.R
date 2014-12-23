#' Plots beta diversity calculations
#' 
#' This function generates files, plotting between-sample diversity,
#' i.e. beta diversity, that plots the thetayc values of all catergorical
#' variable combinations. 
#'
#' @param meta.frame A data frame containing all data. Two columns must be titled 'adiv' and
#' 'nseqs' and must contain the alpha diversity values and number of sequences respectively. All
#' other information be metadata or categorical/continuous variables that describe the samples.
#' @param data.types A list containing text indicators of each meta.frame column data type
#' in order. Input text for each list element are 'cat', 'cont', and 'meta'.
#' @param beta.frame A data frame containing all beta diversity calculations for each sample
#' combination. Data frame should be formatted as "sample1", "sample2", "bdiv", "cmin", "cmax".
#' The cmin and cmax values represent the confidence range for the beta diversity value.
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
#' beta.frame <- read.table("IECG57L01.beta_data.out", header=TRUE, stringsAsFactors=FALSE)
#' bdiv(meta.frame, data.types, beta.frame)

bdiv <- function(meta.frame, data.types, beta.frame)
{
  library(gdata)
  library(MASS)
    pdf(file="BetaDiversity.pdf")

    clist <- c(rgb(27, 158, 119, max=255), rgb(217, 95, 2, max=255), rgb(117, 112, 179, max=255), rgb(231, 41, 138, max=255), rgb(102, 166, 30, max=255), rgb(230, 171, 2, max=255))
    cindex <- 1

    #This for loop appends all category information (2 columns for 2 samples)


    labelcounter <- 6
   # sample1types <- vector()
   # sample2types <- vector()
    for(i in 1:length(data.types))
    {
        if(data.types[i] == 'cat')
        {
            sample1types <- vector()
            sample2types <- vector()
            for(j in 1:nrow(beta.frame))
            {
                sample1 <- beta.frame$sample1[j]
                sample1types <- append(sample1types, meta.frame[,i][which(meta.frame[,1]==sample1)[1]])
               # print(which(meta.frame[,1]==sample1))
               # print(sample1)
               # print("")

                sample2 <- beta.frame$sample2[j]
                sample2types <- append(sample2types, meta.frame[,i][which(meta.frame[,1]==sample2)[1]])
            }

            label1 <- paste('cat', as.character(i), 'sam1', sep='')
            label2 <- paste('cat', as.character(i), 'sam2', sep='')
            beta.frame = cbind(beta.frame, sample1types, stringsAsFactors=FALSE)
            beta.frame = cbind(beta.frame, sample2types, stringsAsFactors=FALSE)
            colnames(beta.frame)[labelcounter] <- label1
            colnames(beta.frame)[labelcounter+1] <- label2
            labelcounter <- labelcounter + 2
        }
    }

    numgraphs <- 0
    for(i in 1:length(data.types))
    {
        if(data.types[i] == 'cat')
        {
            numgraphs <- numgraphs + 1
        }
    }

    colcounter <- 6
    for(i in 1:numgraphs)
    {
        s1 <- colcounter
        colcounter <- colcounter + 1
        s2 <- colcounter
        colcounter <- colcounter + 1

        #This next part gets all the possible combinations of types
        types <- unique(append(beta.frame[,s1], beta.frame[,s2]))
        combinations <- data.frame(type1=rep(NA, factorial(length(types))), type2=rep(NA, factorial(length(types))), stringsAsFactors=FALSE)
        counter <- 0
        for(j in 1:length(types))
        {
            counter <- counter + 1
            combinations[counter,] <- c(types[j], types[j])
        }
        for(j in 1:(length(types)-1))
        {
            for(k in (j+1):length(types))
            {
                counter <- counter +1
                combinations[counter,] <- c(types[j], types[k])
            }
        }
        plotting <- FALSE
        xcounter <- 0
        #This next part gets each combination group of Within and plots
        for(j in 1:length(types))
        {
            valuestoplot <- vector()
            for(k in 1:nrow(beta.frame))
            {
                if((beta.frame[k,s1] == combinations[j,1]) & (beta.frame[k,s2] == combinations[j,2]))
                {
                    valuestoplot <- append(valuestoplot, k)
                }
            }
            #Plots the withins valuestoplot - for each within combination
            c <- clist[cindex]
            cindex <- cindex + 1
            if(cindex == 6){
                cindex <- 1
            }
            if(length(valuestoplot) != 0)
            {
                tmp <- xcounter + 1
                for(k in 1:length(valuestoplot))
                {
                    if(!plotting)
                    {
                        par(new ='F')
                        xcounter <- xcounter + 1
                        plot(xcounter, beta.frame$bdiv[valuestoplot[k]], pch='-', cex=2, xlim=c(0, (nrow(beta.frame)+5*(nrow(combinations))-5)), ylim=c(min(beta.frame$cmin), max(beta.frame$cmax)), xaxt = 'n', xlab = '', ylab = 'thetayc', main = 'Between - Sample (Beta) Diversity', bty='n', col = c)
                        lines(c(xcounter, xcounter), c(beta.frame$cmin[valuestoplot[k]], beta.frame$cmax[valuestoplot[k]]), lwd = 5, col = c)
                        plotting <- TRUE
                    }
                    else
                    {
                        xcounter <- xcounter + 1
                        par(new='T')
                        plot(xcounter, beta.frame$bdiv[valuestoplot[k]], pch='-', cex=2, xlim=c(0, (nrow(beta.frame)+5*(nrow(combinations))-5)), ylim=c(min(beta.frame$cmin), max(beta.frame$cmax)), xaxt = 'n', xlab = '', ylab = '', main = '', bty='n', col = c)
                        lines(c(xcounter, xcounter), c(beta.frame$cmin[valuestoplot[k]], beta.frame$cmax[valuestoplot[k]]), lwd = 5, col = c)
                    }
                }
                mtext(paste(types[j],types[j]), side=1, at=(tmp+xcounter)/2)
                lines(c(tmp-1, xcounter+1), rep(par('usr')[3], 2), lwd = 5, xpd = TRUE, col=c)
                xcounter <- xcounter + 5
            }
            #Changes the color type
            #TODO
        }
        mtext('Within', side=1, line=2, at=(xcounter-4)/2)
        lines(c(1-1, xcounter-5+1), rep(par('usr')[3], 2), lwd = 3, xpd = TRUE, col = rgb(0,0,0,.6))

        #This next part gets each combination group of Between and plots
        btwnstart <- xcounter + 1
        for(j in (length(types)+1):nrow(combinations))
        {
            valuestoplot <- vector()
            for(k in 1:nrow(beta.frame))
            {
                if((beta.frame[k,s1] == combinations[j,1]) & (beta.frame[k,s2] == combinations[j,2]))
                {
                    valuestoplot <- append(valuestoplot, k)
                }
            }
            #Plots the between valuestoplot - for each between combination
            c <- clist[cindex]
            cindex <- cindex + 1
            if(cindex == 6){
                cindex <- 1
            }
            if(length(valuestoplot) != 0)
            {
                tmp <- xcounter + 1
                for(k in 1:length(valuestoplot))
                {
                    if(!plotting)
                    {
                        par(new ='F')
                        xcounter <- xcounter + 1
                        plot(xcounter, beta.frame$bdiv[valuestoplot[k]], pch='-', cex=2, xlim=c(0, (nrow(beta.frame)+5*(nrow(combinations))-5)), ylim=c(min(beta.frame$cmin), max(beta.frame$cmax)), xaxt = 'n', xlab = '', ylab = 'thetayc', main = 'Between - Sample (Beta) Diversity', bty='n', col = c)
                        lines(c(xcounter, xcounter), c(beta.frame$cmin[valuestoplot[k]], beta.frame$cmax[valuestoplot[k]]), lwd = 5, col = c)
                        plotting <- TRUE
                    }
                    else
                    {
                        xcounter <- xcounter + 1
                        par(new='T')
                        plot(xcounter, beta.frame$bdiv[valuestoplot[k]], pch='-', cex=2, xlim=c(0, (nrow(beta.frame)+5*(nrow(combinations))-5)), ylim=c(min(beta.frame$cmin), max(beta.frame$cmax)), xaxt = 'n', xlab = '', ylab = '', main = '', bty='n', col = c)
                        lines(c(xcounter, xcounter), c(beta.frame$cmin[valuestoplot[k]], beta.frame$cmax[valuestoplot[k]]), lwd = 5, col = c)
                    }
                }
                mtext(paste(combinations[j,1], combinations[j,2]), side=1, at=(tmp+xcounter)/2)
                lines(c(tmp-1, xcounter+1), rep(par('usr')[3], 2), lwd = 5, xpd = TRUE, col=c)
            }
            xcounter <- xcounter + 5

        }
        mtext('Between', side=1, line=2, at=(btwnstart+xcounter-5)/2)
        lines(c(btwnstart-1, xcounter-5+1), rep(par('usr')[3], 2), lwd = 3, xpd = TRUE, col = rgb(0,0,0,.6))

        #subtitle
        subtitlename <- colnames(meta.frame)[as.numeric(substring(colnames(beta.frame)[s1], 4, 4))]
        mtext(paste('by', subtitlename))
    }
    dev.off()
}

#meta.frame <- read.table("IECG57L01.graphics_data.txt", skip=1, header=TRUE, stringsAsFactors=FALSE)
#beta.frame <- read.table("IECG57L01.beta_data.out", header=TRUE, stringsAsFactors=FALSE)
#data.types <- strsplit(readLines("IECG57L01.graphics_data.txt")[1], '\t')[[1]]

#bdiv(meta.frame, data.types, beta.frame)
