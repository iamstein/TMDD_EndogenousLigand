library(ggplot2)
library(dplyr)
library(grid)
library(gridExtra)

#graphics functions that are used

#saves the plots, appending the filename to the bottom and adding draft if draft.flag=TRUE
saveplot = function(width,height,dirs,fpref,draft.flag,filex.flag=TRUE,p=last_plot(),draft.x=Inf,draft.y=Inf,filetype="png") {
  filedir        = file.path(dirs$results)
  dirs$filename  = paste0(dirs$outpref,fpref,".",filetype)         #get the full filename
  p = draft(draft.flag,x=draft.x,y=draft.y)
  if (filex.flag) p = filex(dirs)
  ggsave(plot=p,width=width,height=height,file.path(filedir,dirs$filename))
  return(p)
}

#add draft annotation
#p = plot object
#draft.flag = TRUE or FALSE - indicates whether to add draft
draft = function(draft.flag,p=last_plot(),x=Inf,y=Inf,label="DRAFT",color="plum",hjust=1.2,vjust=1.2) {
  if (draft.flag)
    p = p + annotate("text",x=x,y=y,label=label,color=color,hjust=hjust,vjust=vjust)
  return(p)
}

#append filenames to bottom of the plot by greating a Grob
filex = function(dirs,p=last_plot(),fontsize=8) {
  bottom.txt = paste0(dirs$topz,"\n",dirs$pgm.local,dirs$Rname,"\n",dirs$results.local,dirs$filename)
  p = arrangeGrob(p, bottom = textGrob(bottom.txt,gp=gpar(fontsize=fontsize)))
  return(p)
}

#save table as a csv file
savetable = function(tab,dirs,fpref) {
  filedir   = file.path(dirs$results)
  filename  = paste0(dirs$outpref,fpref,".csv") 
  write.csv(tab,file=file.path(filedir,filename),row.names=FALSE,quote=FALSE)
  print(kable(tab))
}

#change x scale to days, weeks, months or years, from days (by default)
#the first letter of scale_str can be d, w, m, y 
#the remaining characters of scale_str are numbers which indicate the end of the x-axis
#increment is how often to space the xticks
#t.start is where to start from
xscale = function(scale_str,p=last_plot(),increment=1,t.start=0) {
  s1 = substr(scale_str,1,1)
  n  = as.numeric(substr(scale_str,2,nchar(scale_str)))
  
  if (s1=="d") {
    xlabel = "Day"
    sc     = 1
  } else if (s1=="w") {
    xlabel = "Week"
    sc     = 7
  } else if (s1=="m") {
    xlabel = "Month"
    sc     = 30.4375
  } else if (s1=="y") {
    xlabel = "Year"
    sc     = 365.25
  } else {
    error("invalid initial string")
  }
  
  breaks = seq(t.start*sc,n*sc,sc*increment)
  limits = c(t.start*sc,n*sc)
  labels = breaks/sc
  
  #p = p + coord_cartesian(xlim=limits)
  p = p + scale_x_continuous(breaks=breaks,labels=labels,lim=limits)
  p = p + xlab(xlabel)
}

#the functions below are for having nicer logscale labels 
#and for easily spacing by 10s (dx=1) or 3s(dx=.5)
log.breaks = function(dx=1,x1=-10,x2=10) {
  breaks = signif(10^seq(x1,x2,dx),1) }

log.labels = function(dx=1,x1=-10,x2=10) {
  labels = sprintf("%1g",log.breaks(dx,x1,x2)) }

log.breaks    = function(dx=1,x1=-10,x2=10)     {signif(10^seq(x1,x2,dx),1) }
log.labels    = function(dx=1,x1=-10,x2=10)     {sprintf("%1g",log.breaks(dx,x1,x2)) }

scale.x.log10 = function(decade.spacing=1,x1=-10,x2=10,...) {
  if (length(intersect(c("breaks","labels"),names(list(...))))==0) {
    scale_command = scale_x_log10(breaks=log.breaks(decade.spacing,x1,x2),
                                  labels=log.labels(decade.spacing,x1,x2),
                                  ...)
  } else {    
    scale_command = scale_x_log10(...)
  }
  return(list(
    scale_command,
    annotation_logticks(base = 10, sides = "b", color = "grey50"))
  )
}

scale.y.log10 = function(decade.spacing=1,x1=-10,x2=10,...) {
  if (length(intersect(c("breaks","labels"),names(list(...))))==0) {
    scale_command = scale_y_log10(breaks=log.breaks(decade.spacing,x1,x2),
                                  labels=log.labels(decade.spacing,x1,x2),
                                  ...)
  } else {    
    scale_command = scale_y_log10(...)
  }
  return(list(
    scale_command,
    annotation_logticks(base = 10, sides = "l", color = "grey50"))
  )}
#remove legend from plot
no.legend = function() {
  x = theme(legend.position="none") }

plot.dose0 = function() {
  x = geom_vline(aes(xintercept=0),color="grey50") }


