buildGOhtml <- function(hgTestOutPut, direc = "Over", Onto="CC", pvalueCut=0.001){

goTerm2html <- function(GOID, title, othernames, table.head,
       	                table.center = TRUE, compSize = NULL, direction=direc,
               	        Ont = Onto){
  	gt <- as.list(GOTERM)
    	goClass <- gt[GOID]	
	    linkName <- sapply(goClass, function(x) x@Term)
	    linkName <- mapply(function(x,y) paste(x,"-",y, sep=" "), GOID, linkName)
	
	    urlVect <- sapply(GOID, function(x) {paste("http://www.ebi.ac.uk/ego/DisplayGoTerm?id=", x, sep="")})
	    filename <- paste(title, Ont, direction, ".html", sep="")
    	outfile <- file(filename, "w")
    	type <- "text/css"
	    cat("<html>", "<head>", "<TITLE>GO Hypergeometric Test</TITLE>",
        "</head>", "<body bgcolor=#FFFFFF >", "<H1 ALIGN=CENTER >GO Hypergeometic Test</H1>",
        ## CSS to allow reasonable spacing for multiple links per cell
        paste("<style type=", type,">",sep = ""), "p{ margin-top: 1px; margin-bottom: 1px; padding-left: 10px; text-indent: -10px }",
        "</style>", file = outfile, sep = "\n")
    if (!missing(title)){
        title <- paste(title, "-", Ont, ": ", direction, sep="")
        cat("<CENTER><H1 ALIGN=\"CENTER\">", title, " </H1></CENTER>\n",
            file = outfile, sep = "\n")}
    if (table.center)
      cat("<CENTER> \n", file = outfile)
    cat("<TABLE BORDER=4>", file = outfile, sep = "\n")
    if (!missing(table.head)) {
        headout <- paste("<TH>", table.head, "</TH>")
        cat("<TR>", headout, "</TR>", file = outfile, sep = "\n")
    }
    
    ##Here is where the url should go:
    
    for ( i in 1:length(urlVect)){
        cat("<TR>", "<TD>", "<a href =", urlVect[i], ">", linkName[i], "</a>", "</TD>", "</TR>","\n",  file=outfile, sep = " ")
    }
    
    cat("</TABLE>", file = outfile)
    if (table.center)
      cat("</CENTER> \n", file = outfile)
    cat("</body>", "</html>", sep = "\n", file = outfile)
    close(outfile)
}


mapply(function(x,y) {gT=names(which(pvalues(y)>pvalueCut));
                  goTerm2html(GOID=gT, title=x, direction=direc,
                              Ont=Onto)},names(hgTestOutPut),hgTestOutPut)

}