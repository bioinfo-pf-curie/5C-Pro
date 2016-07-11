# Nicolas Servant
rm(list=ls())
require(HiTC)

args <- commandArgs(TRUE)
la <- length(args)
if (la > 0){
    for (i in 1:la)
          eval(parse(text=args[[i]]))
  }

message("##input", input)
message("##fbed", fbed)
message("##rbed", rbed)



## load data
d <- importC(input, ygi.bed=fbed, xgi.bed=rbed, allPairwise=TRUE)


## raw maps
HiTC:::saveContactMaps(HTClist(d$chrXchrX), con=gsub(".matrix", ".png", file.path(outdir,basename(input))), trim.range=.95, col.pos=c("white", "orange", "red", "black"))

