rm(list=ls())

getwd()

suppressMessages(library(vegan))


otu=read.table(file.choose(), sep = "\t", header = TRUE, row.names=1) 
env.dat=read.table(file.choose(), sep = "\t", header = TRUE, row.names=1) 

env.st = decostand(env.dat, method="standardize", MARGIN=2)


samp.fg = colnames(otu)
samp.env= rownames(env.st)
my.env = match(samp.fg, samp.env)
env.st2 = na.omit(env.st[my.env, ]) 
samp.env= rownames(env.st2)
my.fg = match(samp.env, samp.fg)
otu = otu[, my.fg]

env.st3=env.st2

otu = t(otu)
otu <- decostand(otu, method = 'standardize', MARGIN =2)

C.whole=rda(otu,env.st3);
# for env selection by CCA inflation factors


inf_factor = vif.cca(C.whole)


na_env = which(is.na(inf_factor)) 
if(isTRUE(length(na_env) > "0") ){ 
  inf_factor = inf_factor[-na_env]
  env.st3=env.st3[,-na_env] 
}
inf_factor

max_env = which(inf_factor == max(inf_factor))
env.st4 = env.st3
while ( inf_factor[max_env] > 20)
{
  env.st4 = env.st4[,-max_env]
  C.reduced = cca(otu, env.st4)
  inf_factor = vif.cca(C.reduced)
  max_env = which(inf_factor == max(inf_factor))
}
output2 = inf_factor ;output2
env.st4

# for F and p values
ind.p = array(0,dim=c(1,ncol(env.st4))) 
ind.F = array(0,dim=c(1,ncol(env.st4))) 
for(j in 1:ncol(env.st4)){
  ind.cca = rda(otu, env.st4[,j]) #ind.cca = cca(otu, env.st[,j], env.st[,-j])  #
  ind.sig = anova(ind.cca,step=1000)
  ind.p[1,j] = ind.sig$Pr[1]
  ind.F[1,j] = ind.sig$F[1]
}
ind.p
ind.F

colnames(ind.p) = colnames(env.st4) 

inf_Fp=rbind(output2,ind.F,ind.p)
row.names(inf_Fp)=c("inf_factor","F","p") 
inf_Fp 

#
C.whole =   rda(otu, env.st4)##cca(otu, env.st4)


x.sig = anova(C.whole)


str(x.sig)

x.p = x.sig$Pr[1] ;
x.F = x.sig$F[1]  ;


output1 = summary(C.whole)
output1
str(output1)
a=output1$sites;a  #
b=output1$cont$importance;b ##eigenvals(C.whole)
c=output1$biplot;c  #

##plot
ca1=round(b[2,1],2);ca1
ca2=round(b[2,2],2);ca2


plot(C.whole,dis=c('wa','cn'),xlab=paste("RDA1=",ca1),ylab=paste("RDA2=",ca2), main = paste("(F= ",round(x.F,2)," , p=",x.p,")",sep=''))
points(C.whole, pch=21, col="gray", bg="gray", cex=0.7)

write.table(a,file="cca_site.txt",sep="\t",col.names=NA)
write.table(b,file="cca_evale.txt",sep="\t",col.names=NA)
write.table(c,file="cca_env.txt",sep="\t",col.names=NA)
write.table(inf_Fp,file="cca_inf_Fp.txt",sep="\t",col.names=NA)
write.table(A,file="rdasy.txt",sep="\t",col.names = NA)



