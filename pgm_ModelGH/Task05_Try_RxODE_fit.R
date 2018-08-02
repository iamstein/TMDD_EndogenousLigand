library(dplyr)
library(RxODE)
library(nlmixr)
library(ggplot2)

ode <- "
C2 = centr/V2;
C3 = peri/V3;
d/dt(centr) = - CL*C2 - Q*C2 + Q*C3;
d/dt(peri)  =   Q*C2  - Q*C3;
"
sys1    = RxODE(model = ode)
param   = c(CL=.0793, V2=.64, Q=.292, V3=9.63)
init    = c(0,0)
time    = sort(c(1:23,seq(24,240,by=24)))	#note that if your pred=0, prop error will blow in your face


ev <- eventTable() %>%
  add.dosing(1000, 1, 24, rate=5000) %>%
  add.sampling(time)

dat = sys1$solve(param,ev,init) %>%
  as.data.frame()

dat.fit = dat %>% 
  select(time,cp=C2) %>%
  mutate(cp = cp*exp(.1*rnorm(nrow(dat))))

ggplot(dat.fit,aes(x=time,y=cp)) + geom_point()


foo = function(th)
{
  pars   = th[1:4]
  names(pars) = names(param)
  pred = sys1$solve(pars,ev,init)[,"C2"]
  sig = th[5]
  res = log(dat.fit[,"cp"]) - log(pred)
  llik = res^2/sig^2 + log(sig^2)
  sum(llik)
}

param.init = param*exp(rnorm(length(param)))
fit = nlminb(c(param.init, .1), foo, control = list(trace=T))

ev.fit = ev %>%
  add.sampling(seq(0,240,.1))

model.prediction.final = sys1$solve(fit$par[1:4],ev.fit,init) %>%
  as.data.frame() %>%
  mutate(type = "final fit")
model.prediction.init = sys1$solve(param.init,ev.fit,init) %>%
  as.data.frame() %>%
  mutate(type = "initial guess")

model.prediction = bind_rows(model.prediction.final,model.prediction.init)

g = ggplot(dat.fit,aes(x=time,y=cp))
g = g + geom_point()
g = g + geom_line(data=model.prediction,aes(y=C2,color=type,group=type))
g = g + scale_y_log10()
print(g)
           
