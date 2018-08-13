


multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                    ncol = cols, nrow = ceiling(numPlots/cols))
  }

 if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}



                                        #load("/data/youna/mouse_aging/enrichpathway.Rdata")
q1=q2=q3=q4=NA
pathwaytype=c("(cell type)","(immune module)")

for(n in 1:length(all.path.res))
  {
    enrichpath=all.path.res[[n]]

    m=ceiling(length(enrichpath[[2]])/4)
    for(l in 1:m-1)
      {
        for(k in (1+l*4):(min(4+l*4,length(enrichpath[[1]]))))    
          {
            p0=p1=NULL

            if(!is.null(enrichpath[[1]][[k]]))
              if( nrow(enrichpath[[1]][[k]])>0)
                {
                  logp=c(-log(as.numeric(rbind(enrichpath[[1]][[k]])[,2])))
                  pthresh=min(abs(logp))
                  df=data.frame(x=enrichpath[[1]][[k]][,1],y= -log(as.numeric(enrichpath[[1]][[k]][,2])))

                  p0<-ggplot(data=df, aes(x=x, y=y))+  geom_bar(stat = "identity", position="identity",width=0.7,size=0.5,color="white") +
                    geom_hline(color="black",size=0.25,yintercept = 0) +
                      geom_hline(color="firebrick4",size=0.25,linetype=2,yintercept = pthresh)  + coord_flip()+labs(x=paste("Enriched gene set",pathwaytype[n]), y="Enrichment, -log10 P") +  theme_minimal(base_size = 8) +  theme(axis.text.y = element_text(size=8),aspect.ratio = 1,panel.grid.major.y = element_line(color="honeydew2"))+ggtitle(paste(tissue.gender.type[k],"both direction"))

                }
            


            module=c(rbind(enrichpath[[2]][[k]])[,1],rbind(enrichpath[[3]][[k]])[,1])
            if(length(module)>0)
              {
                logp=c(-log(as.numeric(rbind(enrichpath[[2]][[k]])[,2])),log(as.numeric(rbind(enrichpath[[3]][[k]])[,2])))
                z=rep("-",length(logp));z[logp>0]="+"
                df=data.frame(x=module,y=logp,z=z)
                pthresh=min(abs(df$y))

                p1=ggplot(df, aes(x,y,fill=z)) +
                  geom_bar(stat = "identity", position="identity",width=0.7,size=0.5,color="white") +
                    geom_hline(color="black",size=0.25,yintercept = 0) +
                      geom_hline(color="firebrick4",size=0.25,linetype=2,yintercept = c(-pthresh,pthresh)) +
                        scale_fill_manual(values = c("-"="dodgerblue","+"="firebrick1"),guide=guide_legend(title = NULL)) +
                          coord_flip() +
                            scale_y_continuous(labels = c(15,10,5,0,5,10,15),breaks = c(-15,-10,-5,0,5,10,15)) +
                                        # facet_wrap(~Contrast,nrow=1) +
                              labs(x=paste("Enriched gene set",pathwaytype[n]), y="Enrichment, -log10 P") +
                                theme_minimal(base_size = 8) +
                                  theme(axis.text.y = element_text(size=8),aspect.ratio = 1,panel.grid.major.y = element_line(color="honeydew2"))+ggtitle(tissue.gender.type[k])

              }
           if(k %% 4==1) {q1=p1}
           if(k %% 4==2) {q2=p1}
           if(k %% 4==3) {q3=p1}
           if(k %% 4==0) {q4=p1}
            
          }
       multiplot(q1,q2,q3,q4,cols=2)
        q1=q2=q3=q4=NA
        
      }
  }

