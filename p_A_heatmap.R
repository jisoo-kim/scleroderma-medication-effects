load("Data.RData") #load d
mv = c("Pred","MTX","MMF","CTX","IVIG","AZA","Rituximab","Tocilizumab","HCQ","TNF","LEF")
pd = d[!is.na(med) & med == 1, mv,with=F]

mvs = c("Pred","MTX","MMF","CTX","IVIG","AZA","RTX","TCZ","HCQ","TNF","LEF")

library(ggplot2)
#library(hrbrthemes)
require(gridExtra)

if(diag_cmb == "all"){
  lcm0 = lcm1 = expand.grid(X = mvs, Y = mvs); lcm0$count = lcm1$count = NA
  k = 1
  for( i in 1:length(mvs)){
    for(j in 1:length(mvs)){
      if(i==j){
        lcm1[k,3] = sum(pd[,i,with=F]==1)
        tmp = apply(pd[,-i,with=F],1,sum)
        lcm0[k,3] = sum(pd[,i,with=F]==1 & tmp==0)
        
      } else {
        tmp = apply(pd[,c(i,j),with=F],1,sum)
        lcm0[k,3] = lcm1[k,3] = sum(tmp == 2)
      }
      k = k+1
    }
  }
  lcm1$ctlab = as.character(lcm1$count); lcm1$ctlab[lcm1$count == 0] = ""
  lcm0$ctlab = as.character(lcm0$count); lcm0$ctlab[lcm0$count == 0] = ""
  p1 = ggplot(lcm1, aes(X, Y, fill= count)) + 
    geom_tile(aes(fill = count)) + geom_text(aes(label=ctlab)) +
    theme(axis.title.x=element_blank(),axis.title.y=element_blank())+
    scale_fill_gradient(low="white", high="blue") + ggtitle("Diag include Combination")
  p0 = ggplot(lcm0, aes(X, Y, fill= count)) + 
    geom_tile(aes(fill = count)) + geom_text(aes(label=ctlab)) +
    theme(axis.title.x=element_blank(),axis.title.y=element_blank())+
    scale_fill_gradient(low="white", high="blue")  + ggtitle("Diag one Med")
  grid.arrange(p0, p1, ncol=2)
  
  remove(lcm0,lcm1)
} else {
  lcm = expand.grid(X = mvs, Y = mvs); lcm$count = NA
  k = 1
  for( i in 1:length(mvs)){
    for(j in 1:length(mvs)){
      if(i==j){
        if(diag_cmb == TRUE){
          lcm[k,3] = sum(pd[,i,with=F]==1)  
        } else {
          tmp = apply(pd[,-i,with=F],1,sum)
          lcm[k,3] = sum(pd[,i,with=F]==1 & tmp==0)
        }
      } else {
        tmp = apply(pd[,c(i,j),with=F],1,sum)
        lcm[k,3] = sum(tmp == 2)
      }
      k = k+1
    }
  }
  lcm$ctlab = as.character(lcm$count); lcm$ctlab[lcm$count == 0] = ""
  p = ggplot(lcm, aes(X, Y, fill= count)) + 
    geom_tile(aes(fill = count)) + geom_text(aes(label=ctlab)) +
    theme(axis.title.x=element_blank(),axis.title.y=element_blank())+
    scale_fill_gradient(low="white", high="blue") 
  plot(p)
  
  remove(lcm)
}

remove(pd, d)
rm(list = ls.str(mode = 'character'))
rm(list = ls.str(mode = 'numeric'))