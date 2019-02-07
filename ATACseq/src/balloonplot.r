
library(ggpubr)

isnull <- function(x)
  {
    if(is.null(x)) return(TRUE) else{ if(nrow(x)==0) return(TRUE) else return(FALSE)}
  }


for(i in 1:length(all.path.res))
  {    
pathwaymat=NULL
for(j in 1:length(all.path.res[[i]][[2]]))
  {
### indicate whether genes/peaks are up/down-regulated
    x=NULL
    if(!isnull(all.path.res[[i]][[2]][[j]]) & !isnull(all.path.res[[i]][[3]][[j]]))
      {
        x=rbind(cbind(all.path.res[[i]][[2]][[j]],"up"),cbind(all.path.res[[i]][[3]][[j]],"down"))
      }else{
        if(!isnull(all.path.res[[i]][[2]][[j]]) )
          x=cbind(all.path.res[[i]][[2]][[j]],"up")
        
        if(!isnull(all.path.res[[i]][[3]][[j]]))
          x=cbind(all.path.res[[i]][[3]][[j]],"down")
      }
### choose between up/down p-value that are more significant 
    if(!is.null(x))
      {
        tablex=table(x[,1])
        double=names(tablex)[tablex==2]

        if(length(double)>0)
          {
            y= x[!(x[,1] %in% double),]

            for(l in 1:length(double))
              {
                tempx= x[x[,1] %in% double[l],]
                y=rbind(y,tempx[as.numeric(tempx[,2])==min(as.numeric(tempx[,2])),])
              }
          }else y=x
        temp.name=strsplit(tissue.gender.type[j]," ")[[1]]
        y=cbind(y,temp.name[1],paste(temp.name[-1],collapse=" "))
        pathwaymat=rbind(pathwaymat,y)
      }
  }


if(!is.null(pathwaymat)) balloonplot[[i]]=rbind(balloonplot[[i]],pathwaymat)
}




###




