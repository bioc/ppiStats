##need Biobase, RColorBrewer, lattice, grid

createSummaryTables <- function(dataGraphs){

  dNames <- names(dataGraphs)
  deg <- lapply(dataGraphs, degree)
  indeg <- lapply(deg, function(x) x$inDegree)
  outdeg <- lapply(deg, function(x) x$outDegree)
  vb <- lapply(outdeg, function(x) names(x[x>0]))
  vp <- lapply(indeg, function(x) names(x[x>0]))
  vbp <- mapply(intersect, vb, vp)

  numVB <- listLen(vb)
  numVP <- listLen(vp)
  numVBP <- listLen(vbp)

  numInteractions <- sapply(indeg, sum)

  EDA <- data.frame(
                    
                    VB=numVB,
                    VP=numVP,
                    VBP=numVBP,
                    "VBP/VB"=round(numVBP/numVB, digits=2),
                    "VP/VB"=round(numVP/numVB,digits=2),
                    "TI"=numInteractions,
                    "TI/VB"=round(numInteractions/numVB,digits=2), 
                    check.names=FALSE,  row.names=dNames)

  return(EDA)

}

viabilityCharts <- function(dataGraphs, total=6466){
tot <- total

EDA <- createSummaryTables(dataGraphs)

EDAsub <- 
    with(EDA, 
         data.frame(Expt = rownames(EDA),
                    VB   = VB - VBP,
                    VBP  = VBP,
                    VP   = VP - VBP,
                    unt  = tot - (VB + VP - VBP)))

bcol <- c(brewer.pal(9, "Pastel1")[9], brewer.pal(12, "Paired")[1:3])[c(2, 3, 4, 1)]

bc <- barchart(reorder(Expt, -unt) ~ VB + VBP + VP + unt, 
             data = EDAsub, stack = TRUE,
             auto.key = 
             list(text = c("Viable Bait only", "Both Viable Prey and Bait",
                    "Viable Prey only", "Absent"),columns = 1, adj = 1), 
             xlab = "Number of proteins",
             par.settings = list(superpose.polygon = list(border = "transparent", col = bcol )))

plot(bc)
  

}


inOutScatterCharts <- function(dataGraphs, pThresh=0.01, pLevels=1e-4){

  bpRed = new.env(parent=globalenv(), hash=FALSE)
  ##bpMat <- lapply(dataGraphs, function(x) {as(x, "matrix")})

  if(!exists("bpMat")) {
    bpMat = new.env(parent=globalenv(), hash=FALSE)
    expNames <- names(dataGraphs)
    for(p in 1:length(dataGraphs)) {
      
      m = as(dataGraphs[[p]], "matrix")
      ## delete self-edges
      diag(m) = 0  
      
      stopifnot(identical(rownames(m), colnames(m)))
      vbp = rownames(m)[ (rowSums(m)>0) & (colSums(m)>0) ]
      
      m = m[vbp, vbp]
      
      if(nrow(m)>1) {
        assign(expNames[[p]], m, envir=bpMat)
      } else {
        cat(sprintf("Omitting %s, there is nothing much to do.\n", expNames[[p]]))
      }
    }
  }
  
  out = file("unrecipInOutDistribIncludes.tex", open="wt")
  
  for(name in ls(bpMat)) {
    f = assessSymmetry(bpMat[[name]])
    sel = (f$p>=pThresh)
    assign(name, bpMat[[name]][sel, sel], envir=bpRed)
    
    myPDF = function(x, ch, ...) {
      fn = sprintf("scp-%s-%s.pdf", name, ch)
      pdf(file=fn, ...)
      x
      dev.off()
      return(fn)
    }
    myEPS = function(x, ch, ...) {
      fn = sprintf("scp-%s-%s.eps", name, ch)
      postscript(file=fn, horizontal = FALSE, onefile = FALSE, paper = "special", ...)
      x
      dev.off()
      return(fn)
    }
    
    f1=myPDF(scpFun(f, "identity", pLevels = pLevels), ch="ident", width=4, height=4)
    f2=myPDF(scpFun(f, "sqrt", pLevels=pLevels), ch="sqrt", width=4, height=4)
    f2e=myEPS(scpFun(f, "sqrt", pLevels=pLevels), ch="sqrt", width=4, height=4)
    f3=myPDF(hist(f$p, main=name, xlab='p', col="skyblue", breaks=seq(0, 1, by=0.01)), 
      ch="hist", width=6, height=2.1)
    
    cat("\\begin{figure}[tp]\\begin{center}\n", file=out)
    cat(sprintf("\\includegraphics[width=0.5\\textwidth]{%s}\n", f1), file=out)
    cat(sprintf("\\includegraphics[width=0.5\\textwidth]{%s}\\\n", f2), file=out)
    cat(sprintf("\\includegraphics[width=0.8\\textwidth]{%s}\n", f3), file=out)
    cat(sprintf("\\caption{Scatterplots of in- and out-degree and symmetry $p$-values for %s}\n",
                name), file=out)
    cat("\\end{center}\\end{figure}\n", file=out)
  }
  
  close(out)
}

##don't export
pvalColors <- function(x, pLevels) brewer.pal(3, "Paired")[1+(x<pLevels)]

scpFun <- function(f, what, pLevels) {
  switch(what,
         identity = {
           trsf = function(x) x
           xlab = expression(n['out'])
           ylab = expression(n['in'])
         },
         sqrt = {
           trsf = function(x) sign(x)*sqrt(abs(x))
           xlab = expression(sqrt(n['out']))
           ylab = expression(sqrt(n['in']))
         })

  plx   = trsf(jitter(f$deg[, c('no', 'ni')]))
  axlim = c(0, max(plx))
  par(mai=c(0.9, 0.9, 0.01, 0.01))
  plot(plx, xlim=axlim, ylim=axlim,
       xlab=xlab, ylab=ylab, pch=20, main="",
       col=pvalColors(f$p, pLevels))

  for(k in 1:ncol(f$contours)) {
    px = f$contours[,k]
    py = (0:(length(px)-1)) - px
    lines(trsf(px), trsf(py), col="#808080")
    lines(trsf(py), trsf(px), col="#808080")
  }
}
