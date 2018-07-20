# inputs for multi.solve
# model    = model file
# p0       = baseline parameters
# event    = eventTable
# explore  = explore data.frame
multi.solve = function(model,p0,event,explore,explore2=NA){
  mod = model
  ev  = event
  ex  = explore
  lseq= function(from,to,length.out){exp(seq(log(from), log(to), length.out = length.out))}
  
  #error check - make sure p is a vector or data.frame
    if (class(p0) == "numeric") {}
    else if (class(p0)=="data.frame" && nrow(p0)==1) {
      name = names(p0)
      p0    = as.numeric(p0)
      names(p0) = name
    } else {stop("unrecognized data type for p0")}

  #error check - make sure p contains a dose call
    if (!("dose" %in% names(p0)))
        stop("dose in p0 is required")
    
  if (!is.na(explore2[1])) {
    ex2 = explore2
    ex2$flag = TRUE #flag for searching through second parameter simultaneously
  } else {
    ex2 = data.frame(flag=FALSE,fold.min=1,fold.max=1,fold.n=1,param=NA,units=NA) #dummy ex2 for no searching
  }
  
  outnum = 0
  OUT = list()
  for (par.name in explore$param) {
    ex = explore[par.name,]
    
    if (!("fold.range" %in% names(explore))) {
      fold.range = lseq( ex$fold.min, ex$fold.max, ex$fold.n)
    } else {
      fold.range = explore[par.name,"fold.range"][[1]]
    }
    
    for (foldpar2 in lseq(ex2$fold.min,ex2$fold.max,ex2$fold.n)) {
    for (foldpar  in fold.range) {
      outnum = outnum+1
      #set up parameters
      p      = p0
      if (ex2$flag) {
        par2.name    = ex2$param
        par2.val.ex2 = p0[par2.name]*foldpar2 #keep this so as not to overwrite, for label.group2 below
        p[par2.name] = par2.val.ex2
      } else {
        par2.name    = "none" #this will only appear in the group.ind below
        par2.val.ex2 = NA
      }
      par.val     = p0[par.name]*foldpar
      p[par.name] = par.val
      
      #set up dosing (where dose and dosing interval are parameters that can vary)
      event  = eventTable()
      event$add.sampling(ev$t.sample)
      
      if (is.element("tau",names(p)))
        event$add.dosing(dose=p["dose"], nbr.doses=ev$n.dose, dosing.interval=p["tau"], dosing.to=ev$cmt)
      else
        event$add.dosing(dose=p["dose"], nbr.doses=1, dosing.to=ev$cmt)
      
      
      #solve ode
      out      = mod$rxode$solve(mod$repar(p), event, mod$init(p))
      out      = mod$rxout(out,p)
      
      #compute Cnonlin for plotting
      Cnonlin = mod$fun.Cnonlin(mod$repar(p))
      out$Cnonlin = Cnonlin
      #stopifnot(!is.na(Cnonlin))
      #stopifnot(!any(is.na(out$D)),!any(is.na(out$time)))
      if (!is.na(Cnonlin)) {
        tnonlin = approx(x=out$D,y=out$time,xout=Cnonlin)
        out$tnonlin = tnonlin$y
      } else {
        tnonlin = NA
      }
      
      #dose normalized
      out$Dnorm        = out$D/p["dose"]
      out$Cnonlin.norm = Cnonlin/p["dose"]
      
      #add useful key quantities for plotting
      out$par  = par.name
      out$pval = signif(par.val,1)
      out$pfold= signif(foldpar,1)
      #out$order= ex$order
      
      out$label.group= paste0(par.name," (",ex$units,")\n",
                              signif(ex$title.scale*p0[par.name]*ex$fold.min,1),"-",signif(ex$title.scale*p0[par.name]*ex$fold.max,1))

      #parameter2
        out$par2       = par2.name
        out$pval2.ex2  = signif(par2.val.ex2,2)
        out$pfold2     = signif(foldpar2,2)
        out$label.group2=paste0(par2.name,"=",signif(par2.val.ex2,2)," ",ex2$units)
      
      out$label.ind  = outnum

      #create giant output list
      OUT[[outnum]] = out
    } }
  }
  OUT = bind_rows(OUT)
  OUT = mutate(OUT, label.group = factor(label.group ,levels=unique(label.group)),
                    label.group2= factor(label.group2,levels=unique(label.group2)))
  OUT$model = as.character(model$name)
     
  
  return(OUT)
}
