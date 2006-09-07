buildGOhtml <- function(hgTestOutPut, direc = "Over", Onto="CC", pvalueCut=0.001,
                        pop="Universe", uni = "Universe", title = "scratch", table.center = TRUE){

    filename <- paste(title, Onto, direc,pop,uni, ".html", sep="")
    outfile <- file(filename, "w")
    type <- "text/css"
    cat("<html>", "<head>", "<TITLE>GO Hypergeometric Test</TITLE>",
        "</head>", "<body bgcolor=#FFFFFF >", "<H1 ALIGN=CENTER >GO Hypergeometic Test</H1>",
        ## CSS to allow reasonable spacing for multiple links per cell
        paste("<style type=", type,">",sep = ""), "p{ margin-top: 1px; margin-bottom: 1px; padding-left: 10px; text-indent: -10px }",
        "</style>", file = outfile, sep = "\n")
    if (!missing(title)){
        title <- paste(title, Onto, ": ", pop, "/", uni, "-", direc, sep="")
        cat("<CENTER><H1 ALIGN=\"CENTER\">", title, " </H1></CENTER>\n",
            file = outfile, sep = "\n")}
    if (table.center)
      cat("<CENTER> \n", file = outfile)
    gt <- as.list(GOTERM)
    yG2P <- as.list(YEASTGO2ALLPROBES)
    
    goTerm2html <- function(GOID,othernames, table.head,
                            compSize = NULL){

        
        
 	GOID <- sort(GOID)
    	goClass <- gt[GOID]
	numGenes <- sapply(GOID, function(x) length(yG2P[[x]]))
        linkName <- sapply(goClass, function(x) x@Term)
        goDef <- sapply(goClass, function(x) x@Definition)
        print(length(goDef))
        print(length(linkName))
        linkName <- mapply(function(x,y,q) paste(x,":",y, "-", "Number of Genes annotated:", q,  sep=" "), GOID, linkName, numGenes)
	
        urlVect <- sapply(GOID, function(x) {paste("http://www.ebi.ac.uk/ego/DisplayGoTerm?id=", x, sep="")})
        
    	
        
        
        cat("<TABLE BORDER=4>", file = outfile, sep = "\n")
        if (!missing(table.head)) {
            headout <- paste("<TH>", table.head, "</TH>")
            cat("<TR>", headout, "</TR>", file = outfile, sep = "\n")
        }
        
        ##Here is where the url should go:
        
        for ( i in 1:length(urlVect)){
            cat("<TR> ", "<TD> ", "<a href =", urlVect[i], " onmouseover=", "\"return escape('",goDef[i],"')\"", "> ", linkName[i], " </a>", " </TD>", " </TR>"," \n",  file=outfile, sep = "")
        }
        
        cat("</TABLE>", file = outfile)
        cat("<br>", file = outfile)
        
        
        
    }
    
    
    
    somefunc <- function(z){

        
        
	mapply(function(x,y) {gT=names(which(pvalues(y)>pvalueCut));
                              goTerm2html(GOID=gT, table.head=x)},names(z),z)
    }
    
    for(i in 1:length(hgTestOutPut)){
        somefunc(hgTestOutPut[i])
    }
    if (table.center)
          cat("</CENTER> \n", file = outfile)
    cat(paste("<script language=", "JavaScript ", "type=", "text/javascript ", "src=",
              "wz_tooltip.js", "></script>", sep=""), file=outfile)
    cat("</body>", "</html>", sep = "\n", file = outfile)
    close(outfile)
}
