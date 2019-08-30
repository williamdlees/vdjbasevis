genloc <- geno_db[s,]$ALLELES
names(genloc) <- geno_db[s,]$GENE_LOC
tmp<-tapply(genloc,names(genloc),function(x){
  x<-sort(x,decreasing=FALSE)
  L<-length(x)
  rep(x,each=12/L)
})
matrix(unlist(tmp),genes_n,12,byrow=TRUE)