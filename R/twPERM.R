twowayPERM = function( v1, v2, NPERM, stat, seed = 123) {
    e1 = new.env()
    set.seed(seed)
    e1$seed = seed
    environment(stat) = e1
    nv2 = !v2
    ans = rep(NA, NPERM)
    for(i in 1:NPERM) {
       vs = sample(v1)
       e1$n11 = sum(vs & v2)
       e1$n12 = sum(vs & nv2)
       e1$n21 = sum(!vs & v2)
       ans[i] = stat()
    }
    ans
  }

 ##x is a vector of length 4 with n11, n12, n21 and n22
 ##and we want to generate two binary vectors that have 
 ##length equal to the sum, and with stats equal to those given
 makeBinVect = function(n11, n12, n21, n22) {
    num = n11+n12+n21+n22
    n1o = n11+n12
    n1z = n21+n22
    n20 = n11+n21
    n2z = n12+n22
    v1 = v2 = rep(FALSE, num)
    v1[1:n11] = v2[1:n11] = TRUE 
    v1[(n11+1):(n11+n12)] = TRUE
    v2[(n11+n12+1):(n11+n12+n21)] = TRUE
    return(list(v1, v2))
 }


