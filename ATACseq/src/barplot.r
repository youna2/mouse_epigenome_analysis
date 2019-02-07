
for(i in 1:2)
  {    
    pathwaymat=NULL
    for(j in 1:length(all.path.res2[[i]][[2]]))
      {
### indicate whether genes/peaks are up/down-regulated
        x=NULL
        if(!isnull(all.path.res2[[i]][[2]][[j]]) & !isnull(all.path.res2[[i]][[3]][[j]]))
          {
            x=rbind(cbind(all.path.res2[[i]][[2]][[j]],"+"),cbind(all.path.res2[[i]][[3]][[j]],"-"))
          }else{
            if(!isnull(all.path.res2[[i]][[2]][[j]]) )
              x=cbind(all.path.res2[[i]][[2]][[j]],"+")
            
            if(!isnull(all.path.res2[[i]][[3]][[j]]))
              x=cbind(all.path.res2[[i]][[3]][[j]],"-")
          }
### choose between up/down p-value that are more significant 

        temp.name=strsplit(tissue.gender.type[j]," ")[[1]]
        x=cbind(x,temp.name[1],paste(temp.name[-1],collapse=" "))
        pathwaymat=rbind(pathwaymat,x)
      }
    if(!is.null(pathwaymat)) barplot[[i]]=rbind(barplot[[i]],pathwaymat)
  }


