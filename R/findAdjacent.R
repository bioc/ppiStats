#a function that will return a list of adjacent nodes for a specified bait
# and NULL if that bait is not a part of the graph

findAdjacent = function(x,bait){
	if(bait %in% nodes(x)){
		ans = unname(unlist(adj(x,bait)))
		if(length(ans)>0) ans = unique(c(bait,ans))
	} else ans=NULL
	ans
}

