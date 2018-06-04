#' Create scatter plots in three views.
#'
#' @description Visualize data in X-Y plane, X-Z plane, and Y-Z
#' plane. The input object (matrix or data.frame) must contain at least 3
#' columns.
#'
#' @param X A matrix or a data.frame that contains at least 3 columns of numeric
#' data. If there are more than 3 columns in X, only the first 3 columns will be
#' used.
#' @param labels A vector containing row labels of X for display. All vector
#' elements should be of type "character" (as.character). The length of
#' vector equals the number of rows in X.
#' @param col.pat.table A data.frame that associates colors and patterns to row
#' labels (see labels before). It needs to contain 3 columns which represent
#' unique labels (column 1) and associated patterns (column 2) and colors
#' (column 3). Default = NA.
#' @param only.row A vector that contains subset of row numbers that are
#' selected to be plotted. Default = NA.
#'
#' @details col.pat.table needs to be set properly. For example, given the
#' vector labels consists of 5 elements as in c('pop1','pop2','pop2','pop2',
#' 'pop1'), the data.frame col.pat.table must contain 2 rows corresponding to
#' 'pop1' and 'pop2', and 3 columns to describe their characteristics (i.e.
#' patterns and colors). See the section "Examples" for more details.
#'
#' @export
#'
#' @importFrom grDevices rgb
#' @importFrom graphics axis legend mtext par plot points
#'
#' @examples
#'
#' #Load simulated dataset
#' data(example_SNP)
#'
#' PCs <- cal.pc.linear(simsnp$snp, no.pc = 3)
#' plot3views( PCs$PC, sample_labels)
#'
#' #To change colors and patterns using symbols
#' all.labels <- unique(sample_labels)
#' my.colors <- c('pink', 'yellow', 'cyan', 'green')
#' my.patterns<- c(0,1,2,3)
#' my.table <- data.frame(all.labels, my.patterns, my.colors)
#' plot3views(PCs$PC, labels = sample_labels, col.pat.table = my.table)
#'
#' #To change colors and patterns using characters
#' my.patterns<- c('o', 'x', '&', '#')
#' #Use cbind to keep my.patterns as 'character'
#' my.table <- cbind(all.labels, my.patterns, my.colors)
#' plot3views(PCs$PC, labels = sample_labels, col.pat.table = my.table)


plot3views <- function(X, labels, col.pat.table = NA, only.row = NA){
  map.color = c("red", rgb(0, 68, 27, maxColorValue = 255),
                "blue", rgb(231, 41, 138, maxColorValue = 255),
                "darkorange", "black")
  map.color = c(map.color, rgb(102, 37, 6, maxColorValue = 255),
                rgb(63, 0, 125 , maxColorValue = 255), "green")
  map.color = c(map.color, "cyan", rgb(250, 159, 181, maxColorValue = 255),
                "yellow", "darkgrey")
  map.color = c(map.color, rgb(116, 196, 118, maxColorValue = 255))

  map.pch = c(1,0,2:18,35:38,60:64,94,126)
  map.pch = c(map.pch,33:34,42,45,47,40,91,123,41,92,93,125)
  map.pch = c(map.pch,49:57,97:107,109:110,112:119,121:122)
  map.pch = c(map.pch,65:78,81:82,84:85,89)

  map.pattern = c()
  for (i in 1:length(map.pch))
    for (j in 1:length(map.color)){
      tmp = c(i,j)
      map.pattern = rbind(map.pattern,tmp)
    }

  if (class(labels) == "data.frame"){
    labels = labels[,1]
    u.label = sort(unique(labels))
  }else{
    u.label = sort(unique(labels))
  }

  par(mfrow=c(2,2))

  #Top-Left
  par(mar=c(4, 1, 1, 2))
  plot(c(min(X[,1]),max(X[,1])),c(min(X[,2]),max(X[,2])),type="n",xlab="",ylab="",main="",axes=FALSE)
  axis(side=1, labels=TRUE, line=2)
  axis(side=4, labels=TRUE, line=2)
  mtext("PC1", side=1, line=0.5)
  mtext("PC2", side=4, line=0.5)
  set_legend = NULL
  set_pch = NULL
  set_col = NULL
  for (k in 1:length(u.label)){
    if (anyNA(col.pat.table)){
      spch = map.pch[map.pattern[k,1]]
      scolor = map.color[map.pattern[k,2]]
    }else{
      idx = which(col.pat.table[,1] == u.label[k])
      spch = col.pat.table[idx,2]
      scolor = col.pat.table[idx,3]
    }

    if (anyNA(only.row)){
      points(X[labels %in% u.label[k],1],X[labels %in% u.label[k],2],col=scolor,pch=spch)
    }else{
      tmp.idx1=which(labels %in% u.label[k])
      tmp.idx2=intersect(tmp.idx1,only.row)
      if (length(tmp.idx2)>0){
        points(X[tmp.idx2,1],X[tmp.idx2,2],col=scolor,pch=spch)
      }
    }
    set_pch = c(set_pch,spch)
    set_col = c(set_col, scolor)
  }

  #Top-Right
  par(mar=c(4, 3.5, 1, 1))
  plot(c(min(X[,3]),max(X[,3])),c(min(X[,2]),max(X[,2])),type="n",xlab="",ylab="",main="",axes=FALSE)
  axis(side=1, labels=TRUE, line=2)
  mtext("PC3", side=1, line=0.5)
  for (k in 1:length(u.label)){
    if (anyNA(col.pat.table)){
      spch = map.pch[map.pattern[k,1]]
      scolor = map.color[map.pattern[k,2]]
    }else{
      idx = which(col.pat.table[,1] == u.label[k])
      spch = col.pat.table[idx,2]
      scolor = col.pat.table[idx,3]
    }
    if (anyNA(only.row)){
      points(X[labels %in% u.label[k],3],X[labels %in% u.label[k],2],col=scolor,pch=spch)
    }else{
      tmp.idx1=which(labels %in% u.label[k])
      tmp.idx2=intersect(tmp.idx1,only.row)
      if (length(tmp.idx2)>0){
        points(X[tmp.idx2,3],X[tmp.idx2,2],col=scolor,pch=spch)
      }
    }
  }

  #Bottom-Left
  par(mar=c(1, 1, 0.5, 2))
  plot(c(min(X[,1]),max(X[,1])),c(min(X[,3]),max(X[,3])),type="n",xlab="",ylab="",main="",axes=FALSE)
  axis(side=4, labels=TRUE, line=2)
  mtext("PC3", side=4, line=0.5)
  for (k in 1:length(u.label)){
    if (anyNA(col.pat.table)){
      spch = map.pch[map.pattern[k,1]]
      scolor = map.color[map.pattern[k,2]]
    }else{
      idx = which(col.pat.table[,1] == u.label[k])
      spch = col.pat.table[idx,2]
      scolor = col.pat.table[idx,3]
    }

    if (anyNA(only.row)){
      points(X[labels %in% u.label[k],1],X[labels %in% u.label[k],3],col=scolor,pch=spch)
    }else{
      tmp.idx1=which(labels %in% u.label[k])
      tmp.idx2=intersect(tmp.idx1,only.row)
      if (length(tmp.idx2)>0){
        points(X[tmp.idx2,1],X[tmp.idx2,3],col=scolor,pch=spch)
      }
    }
  }

  #Bottom-Right
  plot(1, type = "n", axes=FALSE, xlab="", ylab="")
  if (length(u.label)>54){
    legend('center', inset=0, legend=u.label, pch=set_pch, col=set_col, ncol=4)
  }else if (length(u.label)>36){
    legend('center', inset=0, legend=u.label, pch=set_pch, col=set_col, ncol=3)
  }else if (length(u.label)>18){
    legend('center', inset=0, legend=u.label, pch=set_pch, col=set_col, ncol=2)
  }else{
    legend('center', inset=0, legend=u.label, pch=set_pch, col=set_col, ncol=1)
  }

  invisible(NULL)
}


