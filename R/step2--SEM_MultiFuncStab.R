# Please run step1--stability_calcuation.R

library(lavaan)
library(semPlot)
library(vegan)
library(cowplot)
source("R_sources/other_source.R")
source("R_sources/permSEM_source.R")

data_sem=list(cbind(invariability,univariate_stability$space.stab),
              cbind(resist_to.D,univariate_stability$resist.D),
              cbind(resist_to.I,univariate_stability$resist.I))
data=NULL
for (i in 1:3){
  data[[i]]=as.data.frame(data_sem[[i]])
  data[[i]]=as.data.frame(apply(data[[i]],2,stand))
}

###################  buid initial sem model
model.1 <-
'
SF + LF  ~ D + I + PR
Pl~D + I + PR
Multifunc ~ SF + LF + Pl + D + I + PR
SF + LF ~ Pl
#SF~~LF
'

model.1=add.variance(model.1)

Group_Effects=NULL
# reps=10000
sample.n.dots=min(reps,40)
########################stability
df=df1=data[[1]]

#df$Pl[!is.na(df$Pl)]=resid(lm(df$Pl~df$PR))

names(df)
apply(is.na(df),2,sum)
apply(df,2,sd,na.rm=T)

saturated_fit=lavaan(model = model.1 ,data = df, std.ov=TRUE, auto.var=TRUE, auto.fix.first=TRUE)
#model.1=delete.variables.in.model(saturated_fit,delete.pars = c("NPP","","Pll" , "SF.mass" , "LF.div" , "LF.mass","PR.dff"))
saturated_fit=lavaan(model = model.1 ,data = df, std.ov=TRUE, auto.var=TRUE, auto.fix.first=TRUE)
fit=fit1=lavaan(model = model.1, data = df, ordered = NULL,
                sampling.weights   = NULL,
                sample.cov = NULL, sample.mean = NULL, sample.th = NULL,
                sample.nobs = NULL,
                group = NULL, cluster = NULL)
summary(fit1)

fitMeasures(fit,c("chisq","df","pvalue", "cfi", "rmsea","rmsea.pvalue"))
modificationIndices(fit,sort.=TRUE,minimum.value= 1)

semPaths(fit1,layout = "circle2",what = "est",whatLabels = "std.all",style="ram",residuals=T,nCharNodes=6,
         sizeMan=10,sizeMan2 = 5,sizeLat=20,sizeInt=10)
semPaths(fit1,layout = "circle2",whatLabels = "std.all")
delete.paths=list(list(lhs="Multifunc",rhs="SF",op="~"),list(lhs="Multifunc",rhs="LF",op="~"),list(lhs="Multifunc",rhs="Pl",op="~"),
                  list(lhs="SF",rhs="D",op="~"),list(lhs="SF",rhs="I",op="~"),list(lhs="SF",rhs="PR",op="~"),
                  list(lhs="LF",rhs="D",op="~"),list(lhs="LF",rhs="I",op="~"),list(lhs="LF",rhs="PR",op="~"),
                  list(lhs="Pl",rhs="D",op="~"),list(lhs="Pl",rhs="I",op="~"),list(lhs="Pl",rhs="PR",op="~"))

fitM=NULL;nam_fitM=c("original","Pl~PR","LF~PR","SF~PR","Pl~I","LF~I","SF~I","Pl~D","LF~D","SF~D","Multifunc~Pl","Multifunc~LF","Multifunc~SF")
fitM$invariability=fitmeasures(fit1,c("chisq","df","pvalue", "cfi", "rmsea","rmsea.pvalue"),output = "matrix");colnames(fitM$invariability)="original"
for(delete.i in 1:length(delete.paths)){
  md=delete.variables.in.model(saturated_fit,delete.paths=delete.paths[[delete.i]])
  fit.delete.path=lavaan(model = md, data = df)
  temp=fitmeasures(fit.delete.path,c("chisq","df","pvalue", "cfi", "rmsea","rmsea.pvalue"),output = "matrix");colnames(temp)=paste(delete.paths[[delete.i]][c(1,3,2)],collapse = "")
  fitM$invariability=cbind(fitM$invariability,temp)
}
round(fitM$invariability[,nam_fitM],4)
#MD=stepwise.model(fit,data = df,op="~~")
#MD1=MD$best.fit.Model
#reg=MD$Stepwise.regression
#bestfit1 = lavaan(model = MD$best.fit.Model, data = df, std.ov=TRUE, auto.var=TRUE, auto.fix.first=TRUE)
#summary(bestfit1)
#fitMeasures(bestfit1,c("chisq","df","pvalue", "cfi", "rmsea","rmsea.pvalue"))
#semPaths(bestfit1,layout = "circle2",what = "std.all",whatLabels = "std.all",style="ram",residuals=T,nCharNodes=6,sizeMan=10,sizeMan2 = 5,sizeLat=20,sizeInt=10)
#semPaths(bestfit1,layout = "circle2",whatLabels = "std.all")

effects.stat=sem.effects.stat(fit = fit1,data = df)
effects.stat$std.effects.summary
effects.stat$effects.summary
set.seed(1000)
Group_Effects$Invariability=
  lavaan.path.stat(fit,data=df,effect.type=c("total.effects", "direct.effects", "indirect.effects"),
                                   group.rhs=c(D="D",
                                               I="I",
                                               PR="PR",
                                               Pl="Pl",LF="LF",SF="SF"),
                                   group.lhs=c(Pl="Pl",LF="LF",SF="SF",
                                               Func.Stab="Multifunc"),
                                   via=c(comm="Pl+LF+SF"),
                                   rh=NULL,lh=NULL,standardized = TRUE,
                                   bootstrap=T,reps=reps,return.boot=T,boot.depth = NULL,permutation = T,
                                   ind.in.dist = invariability$IND,dist.list = data_SEM.list,boot.perm.simultaneous = F,
                                   exclude.paths=c("Pl1~PR"))

eff=Group_Effects$Invariability
eff$std.effects.summary
eff$std_effects

rown=rownames(eff$std_effects$total.effects)
coln=colnames(eff$std_effects$total.effects)
data.boot=databarplot=NULL


for(independ in c("D","I","PR")){
  for(depend in c("Func.Stab","Pl","SF","LF")){
    independ.i=which(coln==independ)
    depend.i=which(rown==depend)
    for(effect.i in c("total.effects", "indirect.effects","direct.effects","indirect.by_Pl","indirect.by_LF","indirect.by_SF","indirect.by_interac")){
      
      temp0=eff$std.effetc.boot[[effect.i]][depend.i,independ.i,]
      
      tempbar=data.frame(effect.type=effect.i,
                         effects=eff$std_effects[[effect.i]][depend,independ],
                         p.value=eff$p.std_effects[[effect.i]][depend,independ],
                         se=eff$se.std_effects[[effect.i]][depend,independ],
                         ci=eff$se.std_effects[[effect.i]][depend,independ]*1,
                         group=paste0(depend,"~",independ),
                         depend=depend,
                         independ=independ)
      
      databarplot=rbind(databarplot,tempbar)
      
      temp=data.frame(effect.type=effect.i,
                      group=paste0(depend,"~",independ),
                      effects=eff$std.effetc.boot[[effect.i]][depend.i,independ.i,1:sample.n.dots],
                      group=paste0(depend,"~",independ),
                      depend=depend,
                      independ=independ)
      data.boot=rbind(data.boot,temp)
    }
  }
}

#################################################################################################### barplot
plot_bar1=NULL
for(i in c("total.effects", "indirect.effects","direct.effects","indirect.by_Pl","indirect.by_LF","indirect.by_SF","indirect.by_interac")){
  temp.data=data.boot[data.boot$effect.type==i&!data.boot$depend%in%c("Func.Stab"),];temp.data$color=temp.data$depend
  dataforjiter=temp.data
  dataforjiter[dataforjiter$depend%in%c("Pl","SF","LF"),]$depend="Comm.Stab"
  
  #dataforjiter[dataforjiter$depend%in%c("Pl","SF","LF"),]$color="Comm.Stab"
  #dataforjiter[dataforjiter$color%in%c("MS1","MS2","MS3"),]$color=dataforjiter[dataforjiter$color%in%c("MF1","MF2","MF3"),]$color
  #dataforjiter[dataforjiter$depend%in%c("MF1","MF2","MF3"),]$depend="Multi.Func"
  
  
  temp.data.bar=databarplot[databarplot$effect.type==i,]
  log.func=temp.data.bar$depend%in%c("Func.Stab","Comm.Stab")
  data.bar.func=temp.data.bar[log.func,]
  data.bar.comm=temp.data.bar[!log.func,]
  
  #temp.slop.p=slope.p[slope.p$group==i&slope.p$effect.type!="direct.effects",]
  max.gap=max(c(0,temp.data.bar$effects+temp.data.bar$ci),na.rm=T)-min(c(temp.data.bar$effects-temp.data.bar$ci,0),na.rm=T)
  bk.y=max(c(round(max.gap/4,digits = 1),0.1),na.rm=T)
  
  plot_bar1[[i]]=ggplot(data.bar.func,aes(x = independ,y = effects))+
    facet_wrap(~depend)+
    #geom_smooth(data=temp.data,aes(x = PR,y = effects),method = "lm",se=F)+
    geom_jitter(data = dataforjiter,na.rm = T,
                aes(x = independ,y = effects,col=color),height=0,size=0.1,shape=16,alpha=0.8)+
    geom_bar(aes(x = independ,y = effects),stat='identity',data = data.bar.func,
             position='stack',fill=NA,color="black",size=0.25)+
    geom_errorbar(aes(ymin = effects-ci, ymax = effects + ci),width = 0.5,size=0.375)+
    geom_text( aes(x=independ, y=effects + sign(effects)*(ci+0.025*max.gap),
                   label=as.symbol.pvalue(p.value,sig.level=c(0.001,0.01,0.05))),size=4) +
    scale_color_manual(values = c("#FF0000","#00FFFF","#0000FF","#FF00FF","#00FF00", "#FFFF00")[c(1,2,4,6,5,3)])+
    labs(x=element_blank(), y = paste0(i," (SEM 1)"))+
    scale_y_continuous(breaks = sort(c(-seq(bk.y,5,bk.y),seq(0,5,bk.y))))+
    theme_bw()+
    theme(text=element_text(size = 8, family = "sans", face = "plain"),axis.title.x = element_blank(),
          panel.grid.major.x = element_blank(),    panel.grid.minor.x =  element_blank(),  panel.grid.minor.y =  element_blank(),
          panel.background = element_rect(color = 'black', fill = 'transparent'),legend.position = "none")
}
pp <- plot_grid(plot_bar1[[1]],plot_bar1[[2]],
                plot_bar1[[3]],plot_bar1[[4]],
                plot_bar1[[5]],plot_bar1[[6]],
                plot_bar1[[7]],
                hjust = 0,vjust = 0,
                label_fontfamily = "sans", label_size = 12,  label_x =c(0.15,0.15), label_y =1, ncol = 4, nrow=2)
#+ theme(plot.margin = unit(c(0.8,0,0,0), "cm"))
pp
#mkdir("SEM/H1_4/barplot")
#ggsave("SEM/H1_4/barplot/barplot_invariability.jpg",pp ,width = 6, height = 2)
#ggsave("SEM/H1_4/barplot/barplot_invariability.pdf",pp ,width = 6, height = 2)



####################################################resist_to.D
model.to.D=delete.variables.in.model(saturated_fit,delete.pars = "D")
#model.to.D=paste(model.to.D,"; LF~~SF")
df=df2=data[[2]]
names(df)

fit = fit2 = lavaan(model = model.to.D ,data = df,auto.var=TRUE, auto.fix.first=TRUE,auto.cov.lv.x=TRUE)#,test ="bootstrap",verbose=TRUE)
summary(fit2)

fitMeasures(fit,c("chisq","df","pvalue", "cfi", "rmsea","rmsea.pvalue"))
modificationIndices(fit,sort.=TRUE,minimum.value= 0)

semPaths(fit,layout = "circle2",what = "est",whatLabels = "std.all",style="ram",residuals=T,nCharNodes=6,
         sizeMan=10,sizeMan2 = 5,sizeLat=20,sizeInt=10)
semPaths(fit,layout = "circle2",whatLabels = "std.all")

fitM$resist.D=fitmeasures(fit2,c("chisq","df","pvalue", "cfi", "rmsea","rmsea.pvalue"),output = "matrix");colnames(fitM$resist.D)="original"
for(delete.i in 1:length(delete.paths)){
  md=delete.variables.in.model(fit2,delete.paths=delete.paths[[delete.i]])
  fit.delete.path=lavaan(model = md, data = df)
  temp=fitmeasures(fit.delete.path,c("chisq","df","pvalue", "cfi", "rmsea","rmsea.pvalue"),output = "matrix");colnames(temp)=paste(delete.paths[[delete.i]][c(1,3,2)],collapse = "")
  fitM$resist.D=cbind(fitM$resist.D,temp)
  temp=NA
  if(delete.paths[[delete.i]]$rhs=="D"){
    fitM$resist.D[,delete.i+1]=NA
  }
}
round(fitM$resist.D[,nam_fitM],3)


  
#MD=stepwise.model(fit,data = df,op="~~")
#MD2=MD$best.fit.Model
#reg=MD$Stepwise.regression
# bestfit2=lavaan(model = MD$best.fit.Model, data = df)
#summary(bestfit2)
#fitMeasures(bestfit2,c("chisq","df","pvalue", "cfi", "rmsea","rmsea.pvalue"))
#semPaths(bestfit2,layout = "circle2",whatLabels = "std.all")

effects.stat=sem.effects.stat(fit = fit2,data = df)
effects.stat$std.effects.summary
set.seed(1000)
Group_Effects$Resistance_to.D=
  lavaan.path.stat(fit,data=df2,effect.type=c("total.effects", "direct.effects", "indirect.effects"),
                                          group.rhs=c(D="D",
                                                      I="I",
                                                      PR="PR",
                                                      Pl="Pl",LF="LF",SF="SF"),
                                          group.lhs=c(Pl="Pl",LF="LF",SF="SF",
                                                      Func.Stab="Multifunc"),
                                          via=c(comm="Pl+LF+SF"),
                                          rh=NULL,lh=NULL,standardized = TRUE,
                                   bootstrap=T,reps=reps,return.boot=T,boot.depth = NULL,permutation = T,
                                   ind.in.dist = resist_to.D$IND,dist.list = data_SEM.list,boot.perm.simultaneous = F,
                                   exclude.paths=c("Pl1~PR"))


eff=Group_Effects$Resistance_to.D
eff$std.effects.summary
eff$std_effects

rown=rownames(eff$std_effects$total.effects)
coln=colnames(eff$std_effects$total.effects)
data.boot=databarplot=NULL


for(independ in c("D","I","PR")){
  for(depend in c("Func.Stab","Pl","SF","LF")){
    independ.i=which(coln==independ)
    depend.i=which(rown==depend)
    for(effect.i in c("total.effects", "indirect.effects","direct.effects","indirect.by_Pl","indirect.by_LF","indirect.by_SF","indirect.by_interac")){
      
      temp0=eff$std.effetc.boot[[effect.i]][depend.i,independ.i,]
      
      tempbar=data.frame(effect.type=effect.i,
                         effects=eff$std_effects[[effect.i]][depend,independ],
                         p.value=eff$p.std_effects[[effect.i]][depend,independ],
                         se=eff$se.std_effects[[effect.i]][depend,independ],
                         ci=eff$se.std_effects[[effect.i]][depend,independ]*1,
                         group=paste0(depend,"~",independ),
                         depend=depend,
                         independ=independ)
      
      databarplot=rbind(databarplot,tempbar)
      
      temp=data.frame(effect.type=effect.i,
                      group=paste0(depend,"~",independ),
                      effects=eff$std.effetc.boot[[effect.i]][depend.i,independ.i,1:sample.n.dots],
                      group=paste0(depend,"~",independ),
                      depend=depend,
                      independ=independ)
      data.boot=rbind(data.boot,temp)
    }
  }
}

#################################################################################################### barplot
plot_bar2=NULL
for(i in c("total.effects", "indirect.effects","direct.effects","indirect.by_Pl","indirect.by_LF","indirect.by_SF","indirect.by_interac")){
  temp.data=data.boot[data.boot$effect.type==i&!data.boot$depend%in%c("Func.Stab"),];temp.data$color=temp.data$depend
  dataforjiter=temp.data
  dataforjiter[dataforjiter$depend%in%c("Pl","SF","LF"),]$depend="Comm.Stab"
  
  #dataforjiter[dataforjiter$depend%in%c("Pl","SF","LF"),]$color="Comm.Stab"
  #dataforjiter[dataforjiter$color%in%c("MS1","MS2","MS3"),]$color=dataforjiter[dataforjiter$color%in%c("MF1","MF2","MF3"),]$color
  #dataforjiter[dataforjiter$depend%in%c("MF1","MF2","MF3"),]$depend="Multi.Func"
  
  
  temp.data.bar=databarplot[databarplot$effect.type==i,]
  log.func=temp.data.bar$depend%in%c("Func.Stab","Comm.Stab")
  data.bar.func=temp.data.bar[log.func,]
  data.bar.comm=temp.data.bar[!log.func,]
  
  #temp.slop.p=slope.p[slope.p$group==i&slope.p$effect.type!="direct.effects",]
  max.gap=max(c(0,temp.data.bar$effects+temp.data.bar$ci),na.rm=T)-min(c(temp.data.bar$effects-temp.data.bar$ci,0),na.rm=T)
  bk.y=max(c(round(max.gap/4,digits = 1),0.1),na.rm=T)
  
  plot_bar2[[i]]=ggplot(data.bar.func,aes(x = independ,y = effects))+
    facet_wrap(~depend)+
    #geom_smooth(data=temp.data,aes(x = PR,y = effects),method = "lm",se=F)+
    geom_jitter(data = dataforjiter,na.rm = T,
                aes(x = independ,y = effects,col=color),height=0,size=0.1,shape=16,alpha=0.8)+
    geom_bar(aes(x = independ,y = effects),stat='identity',data = data.bar.func,
             position='stack',fill=NA,color="black",size=0.25)+
    geom_errorbar(aes(ymin = effects-ci, ymax = effects + ci),width = 0.5,size=0.375)+
    geom_text( aes(x=independ, y=effects + sign(effects)*(ci+0.025*max.gap),
                   label=as.symbol.pvalue(p.value,sig.level=c(0.001,0.01,0.05))),size=4) +
    scale_color_manual(values = c("#FF0000","#00FFFF","#0000FF","#FF00FF","#00FF00", "#FFFF00")[c(1,2,4,6,5,3)])+
    labs(x=element_blank(), y = paste0(i," (SEM 2)"))+
    scale_y_continuous(breaks = sort(c(-seq(bk.y,5,bk.y),seq(0,5,bk.y))))+
    theme_bw()+
    theme(text=element_text(size = 8, family = "sans", face = "plain"),axis.title.x = element_blank(),
          panel.grid.major.x = element_blank(),    panel.grid.minor.x =  element_blank(),  panel.grid.minor.y =  element_blank(),
          panel.background = element_rect(color = 'black', fill = 'transparent'),legend.position = "none")
}
pp <- plot_grid(plot_bar2[[1]],plot_bar2[[2]],
                plot_bar2[[3]],plot_bar2[[4]],
                plot_bar2[[5]],plot_bar2[[6]],
                plot_bar2[[7]],
                hjust = 0,vjust = 0,
                label_fontfamily = "sans", label_size = 12,  label_x =c(0.15,0.15), label_y =1, ncol = 4, nrow=2)
#+ theme(plot.margin = unit(c(0.8,0,0,0), "cm"))
pp
#mkdir("SEM/H1_4/barplot")
#ggsave("SEM/H1_4/barplot/barplot_resist.toD.jpg",pp ,width = 6, height = 2)
#ggsave("SEM/H1_4/barplot/barplot_resist.toD.pdf",pp ,width = 6, height = 2)

#################################################################
####################################################resist_to.I
model.to.I=delete.variables.in.model(saturated_fit,delete.pars = c("Pl.de","I"))

df=df3=data[[3]]
names(df)
fit = fit3 = lavaan(model = model.to.I ,data = df,auto.var=TRUE, auto.fix.first=TRUE,auto.cov.lv.x=TRUE) #,test ="bootstrap",verbose=TRUE)

#fit = fit3 = lavaan(model = model.to.I ,data = df,auto.var=TRUE, auto.fix.first=TRUE,auto.cov.lv.x=TRUE,test ="bootstrap",verbose=TRUE) #,test ="bootstrap",verbose=TRUE)
summary(fit3)

fitMeasures(fit,c("chisq","df","pvalue", "cfi", "rmsea","rmsea.pvalue"))
modificationIndices(fit,sort.=TRUE,minimum.value= 0)

semPaths(fit,layout = "circle2",what = "est",whatLabels = "std.all",style="ram",residuals=T,nCharNodes=6,
         sizeMan=10,sizeMan2 = 5,sizeLat=20,sizeInt=10)
semPaths(fit,layout = "circle2",whatLabels = "std.all")
fitM$resist.I=fitmeasures(fit3,c("chisq","df","pvalue", "cfi", "rmsea","rmsea.pvalue"),output = "matrix");colnames(fitM$resist.I)="original"
for(delete.i in 1:length(delete.paths)){
  md=delete.variables.in.model(fit3,delete.paths=delete.paths[[delete.i]])
  fit.delete.path=lavaan(model = md, data = df)
  temp=fitmeasures(fit.delete.path,c("chisq","df","pvalue", "cfi", "rmsea","rmsea.pvalue"),output = "matrix");colnames(temp)=paste(delete.paths[[delete.i]][c(1,3,2)],collapse = "")
  fitM$resist.I=cbind(fitM$resist.I,temp)
  if(delete.paths[[delete.i]]$rhs=="I"){
    fitM$resist.I[,delete.i+1]=NA
  }
}
round(fitM$resist.I[,nam_fitM],3)

#MD=stepwise.model(fit,data = df,op="~~")
#MD3=MD$best.fit.Model
#reg=MD$Stepwise.regression
#fit=bestfit3=lavaan(model = MD$best.fit.Model, data = df)
#summary(fit)
#fitMeasures(fit,c("chisq","df","pvalue", "cfi", "rmsea","rmsea.pvalue"))

effects.stat=sem.effects.stat(fit = fit3,data = df)
effects.stat$std.effects.summary

Group_Effects$Resistance_to.I=
  lavaan.path.stat(fit = fit,data=df,effect.type=c("total.effects", "direct.effects", "indirect.effects"),
                                          group.rhs=c(D="D",
                                                      I="I",
                                                      PR="PR",
                                                      Pl="Pl",LF="LF",SF="SF"),
                                          group.lhs=c(Pl="Pl",LF="LF",SF="SF",
                                                      Func.Stab="Multifunc"),
                                          via=c(comm="Pl+LF+SF"),
                                          rh=NULL,lh=NULL,standardized = TRUE,
                                   bootstrap=T,reps=reps,return.boot=T,boot.depth = NULL,permutation = T,
                                   ind.in.dist = resist_to.I$IND,dist.list = data_SEM.list,boot.perm.simultaneous = F,
                                   exclude.paths=c("Pl1~PR"))


eff=Group_Effects$Resistance_to.I
eff$std.effects.summary
eff$std_effects

rown=rownames(eff$std_effects$total.effects)
coln=colnames(eff$std_effects$total.effects)
data.boot=databarplot=NULL


for(independ in c("D","I","PR")){
  for(depend in c("Func.Stab","Pl","SF","LF")){
    independ.i=which(coln==independ)
    depend.i=which(rown==depend)
    for(effect.i in c("total.effects", "indirect.effects","direct.effects","indirect.by_Pl","indirect.by_LF","indirect.by_SF","indirect.by_interac")){
      
      temp0=eff$std.effetc.boot[[effect.i]][depend.i,independ.i,]
      
      tempbar=data.frame(effect.type=effect.i,
                         effects=eff$std_effects[[effect.i]][depend,independ],
                         p.value=eff$p.std_effects[[effect.i]][depend,independ],
                         se=eff$se.std_effects[[effect.i]][depend,independ],
                         ci=eff$se.std_effects[[effect.i]][depend,independ]*1,
                         group=paste0(depend,"~",independ),
                         depend=depend,
                         independ=independ)
      
      databarplot=rbind(databarplot,tempbar)
      
      temp=data.frame(effect.type=effect.i,
                      group=paste0(depend,"~",independ),
                      effects=eff$std.effetc.boot[[effect.i]][depend.i,independ.i,1:sample.n.dots],
                      group=paste0(depend,"~",independ),
                      depend=depend,
                      independ=independ)
      data.boot=rbind(data.boot,temp)
    }
  }
}

#################################################################################################### barplot
plot_bar3=NULL
for(i in c("total.effects", "indirect.effects","direct.effects","indirect.by_Pl","indirect.by_LF","indirect.by_SF","indirect.by_interac")){
  temp.data=data.boot[data.boot$effect.type==i&!data.boot$depend%in%c("Func.Stab"),];temp.data$color=temp.data$depend
  dataforjiter=temp.data
  dataforjiter[dataforjiter$depend%in%c("Pl","SF","LF"),]$depend="Comm.Stab"
  
  #dataforjiter[dataforjiter$depend%in%c("Pl","SF","LF"),]$color="Comm.Stab"
  #dataforjiter[dataforjiter$color%in%c("MS1","MS2","MS3"),]$color=dataforjiter[dataforjiter$color%in%c("MF1","MF2","MF3"),]$color
  #dataforjiter[dataforjiter$depend%in%c("MF1","MF2","MF3"),]$depend="Multi.Func"
  
  
  temp.data.bar=databarplot[databarplot$effect.type==i,]
  log.func=temp.data.bar$depend%in%c("Func.Stab","Comm.Stab")
  data.bar.func=temp.data.bar[log.func,]
  data.bar.comm=temp.data.bar[!log.func,]
  
  #temp.slop.p=slope.p[slope.p$group==i&slope.p$effect.type!="direct.effects",]
  max.gap=max(c(0,temp.data.bar$effects+temp.data.bar$ci),na.rm=T)-min(c(temp.data.bar$effects-temp.data.bar$ci,0),na.rm=T)
  bk.y=max(c(round(max.gap/4,digits = 1),0.1),na.rm=T)
  
  plot_bar3[[i]]=ggplot(data.bar.func,aes(x = independ,y = effects))+
    facet_wrap(~depend)+
    #geom_smooth(data=temp.data,aes(x = PR,y = effects),method = "lm",se=F)+
    geom_jitter(data = dataforjiter,na.rm = T,
                aes(x = independ,y = effects,col=color),height=0,size=0.1,shape=16,alpha=0.8)+
    geom_bar(aes(x = independ,y = effects),stat='identity',data = data.bar.func,
             position='stack',fill=NA,color="black",size=0.25)+
    geom_errorbar(aes(ymin = effects-ci, ymax = effects + ci),width = 0.5,size=0.375)+
    geom_text( aes(x=independ, y=effects + sign(effects)*(ci+0.025*max.gap),
                   label=as.symbol.pvalue(p.value,sig.level=c(0.001,0.01,0.05))),size=4) +
    scale_color_manual(values = c("#FF0000","#00FFFF","#0000FF","#FF00FF","#00FF00", "#FFFF00")[c(1,2,4,6,5,3)])+
    labs(x=element_blank(), y = paste0(i," (SEM 3)"))+
    scale_y_continuous(breaks = sort(c(-seq(bk.y,5,bk.y),seq(0,5,bk.y))))+
    theme_bw()+
    theme(text=element_text(size = 8, family = "sans", face = "plain"),axis.title.x = element_blank(),
          panel.grid.major.x = element_blank(),    panel.grid.minor.x =  element_blank(),  panel.grid.minor.y =  element_blank(),
          panel.background = element_rect(color = 'black', fill = 'transparent'),legend.position = "none")
}
pp <- plot_grid(plot_bar3[[1]],plot_bar3[[2]],
                plot_bar3[[3]],plot_bar3[[4]],
                plot_bar3[[5]],plot_bar3[[6]],
                plot_bar3[[7]],
                hjust = 0,vjust = 0,
                label_fontfamily = "sans", label_size = 12,  label_x =c(0.15,0.15), label_y =1, ncol = 4, nrow=2)
#+ theme(plot.margin = unit(c(0.8,0,0,0), "cm"))
pp
#mkdir("SEM/H1_4/barplot")
#ggsave("SEM/H1_4/barplot/barplot_resist.toI.jpg",pp ,width = 6, height = 2)
#ggsave("SEM/H1_4/barplot/barplot_resist.toI.pdf",pp ,width = 6, height = 2)

##########################################stability

pp <- plot_grid(plot_bar1[[1]],plot_bar1[[2]],
                plot_bar2[[1]],plot_bar2[[2]],
                plot_bar3[[1]],plot_bar3[[2]],
                #labels = c( "PR residuals", "Invasion residuals"),
                hjust = 0,vjust = 0,
                label_fontfamily = "sans", label_size = 12,  label_x =c(0.15,0.15), label_y =1, ncol = 2, nrow=3)
#+  theme(plot.margin = unit(c(0.8,0,0,0), "cm"))
pp
#mkdir("SEM/H1_4/barplot/")
#ggsave("SEM/H1_4/barplot/stability_effect.jpg",pp ,width = 6, height = 4.5)
#ggsave("SEM/H1_4/barplot/stability_effect.pdf",pp ,width = 6, height = 4.5)

########################################################  Summary: Invariability + Resistance to Drought + Resistance to Invasion

Summary_effects=NULL
#for(i in names(Summary_effects)){Summary_effects[[i]]=NA}
for(eff.i in names(Group_Effects$Invariability$std.effetc.boot)){
  #for(eff.i in "indirect.effects"){
  temp.effects=NULL
  temp.effects$Invariability = Group_Effects$Invariability$std.effetc.boot[[eff.i]]
  rownames(temp.effects$Invariability)=rownames(Group_Effects$Invariability$std_effects[[eff.i]])
  colnames(temp.effects$Invariability)=colnames(Group_Effects$Invariability$std_effects[[eff.i]])
  temp.effects$Resistance_to.D = Group_Effects$Resistance_to.D$std.effetc.boot[[eff.i]]
  rownames(temp.effects$Resistance_to.D)=rownames(Group_Effects$Resistance_to.D$std_effects[[eff.i]])
  colnames(temp.effects$Resistance_to.D)=colnames(Group_Effects$Resistance_to.D$std_effects[[eff.i]])
  temp.effects$Resistance_to.I = Group_Effects$Resistance_to.I$std.effetc.boot[[eff.i]]
  rownames(temp.effects$Resistance_to.I)=rownames(Group_Effects$Resistance_to.I$std_effects[[eff.i]])
  colnames(temp.effects$Resistance_to.I)=colnames(Group_Effects$Resistance_to.I$std_effects[[eff.i]])
  
  temp_effect.array=array(NA,dim=c(dim(temp.effects$Invariability),3),
                          dimnames = list(rownames(temp.effects$Invariability),
                                          colnames(temp.effects$Invariability),
                                          1:reps,c("Invariability","Resistance_to.D","Resistance_to.I")))
  for(i in names(temp.effects)){
    for(rown in rownames(temp.effects[[i]])){
      for(coln in colnames(temp.effects[[i]])){
        if(!is.null(temp.effects[[i]][rown,coln,])){
          temp_effect.array[rown,coln,,i]=temp.effects[[i]][rown,coln,]
        }
      }
    }
  }
  Summary_effects$std.effetc.boot[[eff.i]] =temp_effect.array.meanboot=apply(temp_effect.array,MARGIN = c(1,2,3),mean,na.rm=T)
  Summary_effects$se.std_effects[[eff.i]] =temp_effect.array.se=apply(temp_effect.array.meanboot,MARGIN = c(1,2),sd,na.rm=T)
  Summary_effects$std_effects[[eff.i]] =temp_effect.array.mean=apply(temp_effect.array.meanboot,MARGIN = c(1,2),mean,na.rm=T)
  Summary_effects$p.std_effects[[eff.i]] =
    pt(-abs(temp_effect.array.mean)/(temp_effect.array.se+10^-26),df = reps-1)
  
  summary.temp=temp_effect.array.mean;summary.temp[]=""
  
  summary.temp[]=paste(char.align(temp_effect.array.mean,align = "right",digits = 2),
                       as.symbol.pvalue(Summary_effects$p.std_effects[[eff.i]],align = "left",digits = 2,sig.level = c(0.001,0.01,0.05)))
  Summary_effects$std.effects.summary[[eff.i]]=summary.temp
  
}

Group_Effects$Stability=Summary_effects
#save(list = "Group_Effects",file = "SEM/H1_4/Group_Effects.RData")
eff=Group_Effects$Stability
eff$std.effects.summary
eff$std_effects

rown=rownames(eff$std_effects$total.effects)
coln=colnames(eff$std_effects$total.effects)
data.boot=databarplot=NULL


for(independ in c("D","I","PR")){
  for(depend in c("Func.Stab","Pl","SF","LF")){
    independ.i=which(coln==independ)
    depend.i=which(rown==depend)
    for(effect.i in c("total.effects", "indirect.effects","direct.effects","indirect.by_Pl","indirect.by_LF","indirect.by_SF","indirect.by_interac")){
      
      temp0=eff$std.effetc.boot[[effect.i]][depend.i,independ.i,]
      
      tempbar=data.frame(effect.type=effect.i,
                         effects=eff$std_effects[[effect.i]][depend,independ],
                         p.value=eff$p.std_effects[[effect.i]][depend,independ],
                         se=eff$se.std_effects[[effect.i]][depend,independ],
                         ci=eff$se.std_effects[[effect.i]][depend,independ]*1,
                         group=paste0(depend,"~",independ),
                         depend=depend,
                         independ=independ)
      
      databarplot=rbind(databarplot,tempbar)
      
      temp=data.frame(effect.type=effect.i,
                      group=paste0(depend,"~",independ),
                      effects=eff$std.effetc.boot[[effect.i]][depend.i,independ.i,1:sample.n.dots],
                      group=paste0(depend,"~",independ),
                      depend=depend,
                      independ=independ)
      data.boot=rbind(data.boot,temp)
    }
  }
}

#################################################################################################### barplot
plot_bar=NULL
for(i in names(Group_Effects$Invariability$std.effetc.boot)){
  temp.data=data.boot[data.boot$effect.type==i&!data.boot$depend%in%c("Func.Stab","Comm.Stab"),];temp.data$color=temp.data$depend
  dataforjiter=temp.data
  dataforjiter[dataforjiter$depend%in%c("Pl","SF","LF"),]$depend="Comm.Stab"
  
  #dataforjiter[dataforjiter$depend%in%c("Pl","SF","LF"),]$color="Comm.Stab"
  #dataforjiter[dataforjiter$color%in%c("MS1","MS2","MS3"),]$color=dataforjiter[dataforjiter$color%in%c("MF1","MF2","MF3"),]$color
  #dataforjiter[dataforjiter$depend%in%c("MF1","MF2","MF3"),]$depend="Multi.Func"
  
  
  temp.data.bar=databarplot[databarplot$effect.type==i,]
  log.func=temp.data.bar$depend%in%c("Func.Stab","Comm.Stab")
  data.bar.func=temp.data.bar[log.func,]
  data.bar.comm=temp.data.bar[!log.func,]
  
  #temp.slop.p=slope.p[slope.p$group==i&slope.p$effect.type!="direct.effects",]
  max.gap=max(c(0,temp.data.bar$effects+temp.data.bar$ci),na.rm=T)-min(c(temp.data.bar$effects-temp.data.bar$ci,0),na.rm=T)
  bk.y=max(c(round(max.gap/4,digits = 1),0.1),na.rm=T)
  
  plot_bar[[i]]=ggplot(data.bar.func,aes(x = independ,y = effects))+
    facet_wrap(~depend)+
    #geom_smooth(data=temp.data,aes(x = PR,y = effects),method = "lm",se=F)+
    geom_jitter(data = dataforjiter,na.rm = T,
                aes(x = independ,y = effects,col=color),height=0,size=0.1,shape=16,alpha=0.8)+
    geom_bar(aes(x = independ,y = effects),stat='identity',data = data.bar.func,
             position='stack',fill=NA,color="black",size=0.25)+
    geom_errorbar(aes(ymin = effects-ci, ymax = effects + ci),width = 0.5,size=0.375)+
    geom_text( aes(x=independ, y=effects + sign(effects)*(ci+0.025*max.gap),
                   label=as.symbol.pvalue(p.value,sig.level=c(0.001,0.01,0.05))),size=4) +
    scale_color_manual(values = c("#FF0000","#00FFFF","#0000FF","#FF00FF","#00FF00", "#FFFF00")[c(1,2,4,6,5,3)])+
    labs(x=element_blank(), y = paste0(i," (SEM 1+2+3)"))+
    scale_y_continuous(breaks = sort(c(-seq(bk.y,5,bk.y),seq(0,5,bk.y))))+
    theme_bw()+
    theme(text=element_text(size = 8, family = "sans", face = "plain"),axis.title.x = element_blank(),
          panel.grid.major.x = element_blank(),    panel.grid.minor.x =  element_blank(),  panel.grid.minor.y =  element_blank(),
          panel.background = element_rect(color = 'black', fill = 'transparent'),legend.position = "none")
}
pp <- plot_grid(plot_bar[[1]],plot_bar[[2]],
                plot_bar[[3]],plot_bar[[4]],
                plot_bar[[5]],plot_bar[[6]],
                plot_bar[[7]],
                hjust = 0,vjust = 0,
                label_fontfamily = "sans", label_size = 12,  label_x =c(0.15,0.15), label_y =1, ncol = 4, nrow=2)
#+ theme(plot.margin = unit(c(0.8,0,0,0), "cm"))
pp
#mkdir("SEM/H1_4/barplot")
#ggsave("SEM/H1_4/barplot/barplot_Stability.jpg",pp ,width = 6, height = 1.5)
#ggsave("SEM/H1_4/barplot/barplot_Stability.pdf",pp ,width = 6, height = 1.5)

pp <- plot_grid(plot_bar1[[1]],plot_bar1[[2]],
                plot_bar2[[1]],plot_bar2[[2]],
                plot_bar3[[1]],plot_bar3[[2]],
                plot_bar[[1]],plot_bar[[2]],
                #labels = c( "PR residuals", "Invasion residuals"),
                hjust = 0,vjust = 0,
                label_fontfamily = "sans", label_size = 12,  label_x =c(0.15,0.15), label_y =1, ncol = 2, nrow=4)
#+  theme(plot.margin = unit(c(0.8,0,0,0), "cm"))
pp
#mkdir("SEM/H1_4/barplot/")
#ggsave("SEM/H1_4/barplot/stability_effect_all.jpg",pp ,width = 6, height = 6)
#ggsave("SEM/H1_4/barplot/stability_effect_all.pdf",pp ,width = 6, height = 6)



################################################################################################
data.boot=databarplot=NULL
for(iii in 1:3){
  eff=Group_Effects[[iii]]
  eff$std.effects.summary
  eff$std_effects
  
  rown=rownames(eff$std_effects$total.effects)
  coln=colnames(eff$std_effects$total.effects)
  
  
  
  for(independ in c("D","I","PR")){
    for(depend in c("Func.Stab","Pl","SF","LF")){
      independ.i=which(coln==independ)
      depend.i=which(rown==depend)
      for(effect.i in c("total.effects", "indirect.effects","direct.effects","indirect.by_Pl","indirect.by_LF","indirect.by_SF","indirect.by_interac")){
        
        temp0=eff$std.effetc.boot[[effect.i]][depend.i,independ.i,]
        
        tempbar=data.frame(effect.type=effect.i,
                           effects=eff$std_effects[[effect.i]][depend,independ],
                           p.value=eff$p.std_effects[[effect.i]][depend,independ],
                           se=eff$se.std_effects[[effect.i]][depend,independ],
                           ci=eff$se.std_effects[[effect.i]][depend,independ]*1,
                           group=paste0(depend,"~",independ),
                           depend=depend,
                           independ=independ,independ.SEMi=paste(independ,iii))
        
        databarplot=rbind(databarplot,tempbar)
        
        temp=data.frame(effect.type=effect.i,
                        group=paste0(depend,"~",independ),
                        effects=eff$std.effetc.boot[[effect.i]][depend.i,independ.i,1:sample.n.dots],
                        group=paste0(depend,"~",independ),
                        depend=depend,
                        independ=independ,independ.SEMi=paste(independ,iii))
        data.boot=rbind(data.boot,temp)
      }
    }
  }
}


#################################################################################################### barplot for effects on communities
plot_bar_community1=NULL
temp.data=data.boot[!data.boot$depend%in%c("Func.Stab","Comm.Stab"),];temp.data$color=temp.data$depend

dataforjiter=temp.data
dataforjiter[dataforjiter$depend%in%c("Pl","SF","LF"),]$depend="Comm.Stab"

#dataforjiter[dataforjiter$depend%in%c("Pl","SF","LF"),]$color="Comm.Stab"
#dataforjiter[dataforjiter$color%in%c("MS1","MS2","MS3"),]$color=dataforjiter[dataforjiter$color%in%c("MF1","MF2","MF3"),]$color
#dataforjiter[dataforjiter$depend%in%c("MF1","MF2","MF3"),]$depend="Multi.Func"


#temp.data.bar=databarplot[databarplot$effect.type==i,]
temp.data.bar=databarplot[!databarplot$independ.SEMi%in%c("D 2","I 3"),]
log.func=temp.data.bar$depend%in%c("Comm.Stab")
data.bar.func=temp.data.bar[log.func,]
data.bar.comm=temp.data.bar[!log.func,]
data.bar.comm$depend=factor(data.bar.comm$depend,levels = c("Pl","SF","LF","Func.Stab"))
#data.bar.comm$depend=paste(as.numeric(data.bar.comm$depend),data.bar.comm$depend)

data.bar.comm$effect.type=factor(data.bar.comm$effect.type,levels = c("total.effects","indirect.effects","direct.effects",
                                                                      "indirect.by_Pl","indirect.by_LF","indirect.by_SF","indirect.by_interac"))
#temp.slop.p=slope.p[slope.p$group==i&slope.p$effect.type!="direct.effects",]
for(j in names(Group_Effects$Invariability$std.effetc.boot)){
  judge_bar=data.bar.comm$depend=="Func.Stab"&data.bar.comm$effect.type==j
  max.gap=max(c(0,data.bar.comm[judge_bar,]$effects+data.bar.comm[judge_bar,]$ci),na.rm=T)-min(c(data.bar.comm[judge_bar,]$effects-data.bar.comm[judge_bar,]$ci,0),na.rm=T)
  bk.y=max(c(round(max.gap/4,digits = 2),0.001),na.rm=T)
  
  
  plot_bar_community1[[j]] = ggplot(data.bar.comm[data.bar.comm$depend=="Func.Stab",],aes(x = independ.SEMi,y = effects))+
    #facet_wrap(~paste(independ),nrow=1)+
    #geom_smooth(data=temp.data,aes(x = PR,y = effects),method = "lm",se=F)+
    #geom_jitter(data = dataforjiter[dataforjiter$independ=="Func.Stab",],na.rm = T,
    #aes(x = independ.SEMi,y = effects,col=color),height=0,width = 0.2,size=0.1,shape=16,alpha=0.8)+

    geom_bar(aes(x = independ.SEMi,y = effects),stat='identity',position='stack',fill=NA,color="black",size=0.25,
             data = data.bar.comm[judge_bar,])+
    geom_errorbar(aes(ymin = effects-ci, ymax = effects + ci),width = 0.5,size=0.375,
                  data = data.bar.comm[judge_bar,]) +
    geom_text( aes(x=independ.SEMi, y=effects + sign(effects)*(ci+0.1*max.gap),
                   label=as.symbol.pvalue(p.value,sig.level=c(0.001,0.01,0.05))),size=4,
               data = data.bar.comm[judge_bar,])+
    #scale_color_manual(values = c("#FF0000","#00FFFF","#0000FF","#FF00FF","#00FF00", "#FFFF00")[c(1,5,3,2,4,6)])+
    labs(x=element_blank(), y = paste(j,"on Func.Stab"))+
    scale_y_continuous(breaks = sort(c(-seq(bk.y,5,bk.y),seq(0,5,bk.y))))+
    theme_bw()+
    theme(text=element_text(size = 8, family = "sans", face = "plain"),axis.title.x = element_blank(),
          panel.grid.major.x = element_blank(),    panel.grid.minor.x =  element_blank(),  panel.grid.minor.y =  element_blank(),
          panel.background = element_rect(color = 'black', fill = 'transparent'),legend.position = "none")

  
}
pp <- plot_grid(plot_bar_community1[[1]],
                plot_bar_community1[[2]],
                plot_bar_community1[[3]],
                plot_bar_community1[[4]],
                plot_bar_community1[[5]],
                plot_bar_community1[[6]],
                plot_bar_community1[[7]],
                hjust = 0,vjust = 0,
                label_fontfamily = "sans", label_size = 12,  label_x =c(0.15,0.15), label_y =1, ncol = 4, nrow=2)
#+ theme(plot.margin = unit(c(0.8,0,0,0), "cm"))
pp
############################################################################################################
plot_bar_cummulation=NULL
temp.data=data.boot[!data.boot$depend%in%c("Func.Stab"),];temp.data$color=temp.data$depend

dataforjiter=temp.data
dataforjiter[dataforjiter$depend%in%c("Pl","SF","LF"),]$depend="Comm.Stab"

#dataforjiter[dataforjiter$depend%in%c("Pl","SF","LF"),]$color="Comm.Stab"
#dataforjiter[dataforjiter$color%in%c("MS1","MS2","MS3"),]$color=dataforjiter[dataforjiter$color%in%c("MF1","MF2","MF3"),]$color
#dataforjiter[dataforjiter$depend%in%c("MF1","MF2","MF3"),]$depend="Multi.Func"


#temp.data.bar=databarplot[databarplot$effect.type==i,]
temp.data.bar=databarplot[!databarplot$independ.SEMi%in%c("D 2","I 3"),]
log.func=temp.data.bar$depend%in%c("Comm.Stab")
data.bar.func=temp.data.bar[log.func,]
data.bar.comm=temp.data.bar[!log.func,]
data.bar.comm$depend=factor(data.bar.comm$depend,levels = c("Pl","SF","LF","Func.Stab"))
#data.bar.comm$depend=paste(as.numeric(data.bar.comm$depend),data.bar.comm$depend)
data.bar.comm$indirect.by=data.bar.comm$effect.type
data.bar.comm$effect.type=factor(data.bar.comm$effect.type,levels = c("total.effects","indirect.effects","direct.effects",
                                                                      "indirect.by_Pl","indirect.by_LF","indirect.by_SF","indirect.by_interac"))

data.bar.comm$indirect.by[data.bar.comm$indirect.by%in%c("indirect.by_Pl","indirect.by_LF","indirect.by_SF","indirect.by_interac")]="indirect.by"
data.bar.comm$indirect.by=factor(data.bar.comm$indirect.by,levels = c("total.effects","indirect.effects","direct.effects","indirect.by"))

#temp.slop.p=slope.p[slope.p$group==i&slope.p$effect.type!="direct.effects",]
for(j in levels(data.bar.comm$indirect.by)){
  judge_bar=data.bar.comm$depend=="Func.Stab"&data.bar.comm$indirect.by==j
  max.gap=max(c(0,data.bar.comm[judge_bar,]$effects+data.bar.comm[judge_bar,]$ci),na.rm=T)-min(c(data.bar.comm[judge_bar,]$effects-data.bar.comm[judge_bar,]$ci,0),na.rm=T)
  bk.y=max(c(round(max.gap/4,digits = 2),0.05),na.rm=T)
  
  
  plot_bar_cummulation[[j]] = ggplot(data.bar.comm[data.bar.comm$depend=="Func.Stab",],aes(x = independ.SEMi,y = effects))+
    #facet_wrap(~paste(independ),nrow=1)+
    #geom_smooth(data=temp.data,aes(x = PR,y = effects),method = "lm",se=F)+
    #geom_jitter(data = dataforjiter[dataforjiter$independ=="Func.Stab",],na.rm = T,
    #aes(x = independ.SEMi,y = effects,col=color),height=0,width = 0.2,size=0.1,shape=16,alpha=0.8)+
    
    geom_bar(aes(group=indirect.by, x = independ.SEMi,y = effects,fill=effect.type),stat='identity',position='stack',color="black",size=0.25,
             data = data.bar.comm[judge_bar,],alpha=0.5)+
    #geom_errorbar(aes(ymin = effects-ci, ymax = effects + ci,color=effect.type),width = 0.5,size=0.375,
                  #data = data.bar.comm[judge_bar,]) +
    geom_text( aes(x=independ.SEMi, y=effects + sign(effects)*(0*ci+0.0*max.gap),color=effect.type,
                   label=as.symbol.pvalue(p.value,sig.level=c(0.001,0.01,0.05))),size=4,
               data = data.bar.comm[judge_bar,])+
    #scale_color_manual(values = c("#FF0000","#00FFFF","#0000FF","#FF00FF","#00FF00", "#FFFF00")[c(1,5,3,2,4,6)])+
    labs(x=element_blank(), y = paste(j,"on Func.Stab"))+
    scale_y_continuous(breaks = sort(c(-seq(bk.y,5,bk.y),seq(0,5,bk.y))))+
    theme_bw()+
    theme(text=element_text(size = 8, family = "sans", face = "plain"),axis.title.x = element_blank(),
          panel.grid.major.x = element_blank(),    panel.grid.minor.x =  element_blank(),  panel.grid.minor.y =  element_blank(),
          panel.background = element_rect(color = 'black', fill = 'transparent'),legend.position = "none")
  
  
}
pp <- plot_grid(plot_bar_cummulation[[1]],
                plot_bar_cummulation[[3]],
                plot_bar_cummulation[[2]],
                plot_bar_cummulation[[4]],
                hjust = 0,vjust = 0,
                label_fontfamily = "sans", label_size = 12,  label_x =c(0.15,0.15), label_y =1, ncol = 2, nrow=2)
#+ theme(plot.margin = unit(c(0.8,0,0,0), "cm"))
pp

temp = data.bar.comm[data.bar.comm$depend=="Func.Stab",]
partition_total=partition_indirect=effect=NULL
for(i in unique(temp$independ.SEMi)){
  temp1=temp[temp$independ.SEMi==i,]
  to.toal=temp1$effects/temp1[temp1$effect.type=="total.effects",]$effects
  to.indirect=temp1$effects/temp1[temp1$effect.type=="indirect.effects",]$effects
  effect.i=as.numeric(temp1$effects)
  names(to.toal)=names(to.indirect)=names(effect.i)=temp1$effect.type
  partition_total[[i]]=to.toal
  partition_indirect[[i]]=to.indirect
  effect[[i]]=effect.i
}
partition_total=as.data.frame(partition_total)
partition_indirect=as.data.frame(partition_indirect)
effect=as.data.frame(effect)
partition_total
partition_indirect
effect
############################################################################################################×ÝÏò±È½Ï
plot_bar_community2=plot_bar_community3=NULL
temp.data=data.boot[!data.boot$depend%in%c("Func.Stab"),];temp.data$color=temp.data$depend

dataforjiter=temp.data
dataforjiter[dataforjiter$depend%in%c("Pl","SF","LF"),]$depend="Comm.Stab"

#dataforjiter[dataforjiter$depend%in%c("Pl","SF","LF"),]$color="Comm.Stab"
#dataforjiter[dataforjiter$color%in%c("MS1","MS2","MS3"),]$color=dataforjiter[dataforjiter$color%in%c("MF1","MF2","MF3"),]$color
#dataforjiter[dataforjiter$depend%in%c("MF1","MF2","MF3"),]$depend="Multi.Func"


#temp.data.bar=databarplot[databarplot$effect.type==i,]
temp.data.bar=databarplot[!databarplot$independ.SEMi%in%c("D 2","I 3"),]
log.func=temp.data.bar$depend%in%c("Comm.Stab")
data.bar.func=temp.data.bar[log.func,]
data.bar.comm=temp.data.bar[!log.func,]
data.bar.comm$depend=factor(data.bar.comm$depend,levels = c("Pl","SF","LF","Func.Stab"))
#data.bar.comm$depend=paste(as.numeric(data.bar.comm$depend),data.bar.comm$depend)

data.bar.comm$effect.type=factor(data.bar.comm$effect.type,levels = c("total.effects","indirect.effects","direct.effects"))
#temp.slop.p=slope.p[slope.p$group==i&slope.p$effect.type!="direct.effects",]
for(j in sort(unique(data.bar.comm$depend))){
  judge_total=data.bar.comm$depend==j&data.bar.comm$effect.type=="total.effects"
  max.gap=max(c(0,data.bar.comm[judge_total,]$effects+data.bar.comm[judge_total,]$ci),na.rm=T)-min(c(data.bar.comm[judge_total,]$effects-data.bar.comm[judge_total,]$ci,0),na.rm=T)
  bk.y=max(c(round(2*max.gap/4,digits = 1),0.02)/2,na.rm=T)
  
  
  plot_bar_community2[[j]] = ggplot(data.bar.comm[data.bar.comm$depend==j,],aes(x = independ.SEMi,y = effects))+
    #facet_wrap(~paste(independ),nrow=1)+
    #geom_smooth(data=temp.data,aes(x = PR,y = effects),method = "lm",se=F)+
    #geom_jitter(data = dataforjiter[dataforjiter$independ==j,],na.rm = T,
    #aes(x = independ.SEMi,y = effects,col=color),height=0,width = 0.2,size=0.1,shape=16,alpha=0.8)+
    geom_bar(aes(x = independ.SEMi,y = effects,fill=effect.type),stat='identity',position='stack',color=NA,size=0.25,alpha=0.25,
             data = data.bar.comm[data.bar.comm$depend==j&data.bar.comm$effect.type%in%c("indirect.effects","direct.effects"),])+
    geom_bar(aes(x = independ.SEMi,y = effects),stat='identity',position='stack',fill=NA,color="black",size=0.25,
             data = data.bar.comm[judge_total,])+
    #geom_errorbar(aes(ymin = effects-ci, ymax = effects + ci),width = 0.5,size=0.375,
                  #data = data.bar.comm[judge_total,])+
    geom_text( aes(x=independ.SEMi, y=effects,
                   label=as.symbol.pvalue(p.value,sig.level=c(0.001,0.01,0.05))),size=4,
               data = data.bar.comm[data.bar.comm$depend==j&
                                      (judge_total|data.bar.comm$effect.type%in%c("indirect.effects","direct.effects")),])+
    scale_color_manual(values = c("#FF0000","#00FFFF","#0000FF","#FF00FF","#00FF00", "#FFFF00")[c(1,5,3,2,4,6)])+
    labs(x=element_blank(), y = paste(" Effects on",j))+
    scale_y_continuous(breaks = sort(c(-seq(bk.y,5,bk.y),seq(0,5,bk.y))))+
    theme_bw()+
    theme(text=element_text(size = 8, family = "sans", face = "plain"),axis.title.x = element_blank(),
          panel.grid.major.x = element_blank(),    panel.grid.minor.x =  element_blank(),  panel.grid.minor.y =  element_blank(),
          panel.background = element_rect(color = 'black', fill = 'transparent'),legend.position = "none")
  
  plot_bar_community3[[j]] = ggplot(data.bar.comm[data.bar.comm$depend==j,],aes(y = independ.SEMi,x = effects))+
    #facet_wrap(~paste(independ),nrow=1)+
    #geom_smooth(data=temp.data,aes(y = PR,x = effects),method = "lm",se=F)+
    #geom_jitter(data = dataforjiter[dataforjiter$independ==j,],na.rm = T,
    #aes(y = independ.SEMi,x = effects,col=color),height=0,width = 0.2,size=0.1,shape=16,alpha=0.8)+
    geom_bar(aes(y = independ.SEMi,x = effects,fill=effect.type),stat='identity',position='stack',color=NA,size=0.25,alpha=0.25,
             data = data.bar.comm[data.bar.comm$depend==j&!judge_total,])+
    geom_bar(aes(y = independ.SEMi,x = effects),stat='identity',position='stack',fill=NA,color="black",size=0.25,
             data = data.bar.comm[judge_total,])+
    geom_errorbar(aes(xmin = effects-ci, xmax = effects + ci),width = 0.5,size=0.375,
                  data = data.bar.comm[judge_total,])+
    geom_text( aes(y=independ.SEMi, x=effects + sign(effects)*(ci+0.1*max.gap),
                   label=as.symbol.pvalue(p.value,sig.level=c(0.001,0.01,0.05))),size=4,
               data = data.bar.comm[judge_total,])+
    scale_color_manual(values = c("#FF0000","#00FFFF","#0000FF","#FF00FF","#00FF00", "#FFFF00")[c(1,5,3,2,4,6)])+
    labs(y=element_blank(), x = paste(" Effects on",j))+
    scale_x_continuous(breaks = sort(c(-seq(bk.y,5,bk.y),seq(0,5,bk.y))))+
    theme_bw()+
    theme(text=element_text(size = 8, family = "sans", face = "plain"),axis.title.y = element_blank(),
          panel.grid.major.y = element_blank(),    panel.grid.minor.y =  element_blank(),  panel.grid.minor.x =  element_blank(),
          panel.background = element_rect(color = 'black', fill = 'transparent'),legend.position = "none")
  
}

#}
pp <- plot_grid(plot_bar_community2[[1]],
                plot_bar_community2[[2]],
                plot_bar_community2[[3]],
                plot_bar_community2[[4]],


                hjust = 0,vjust = 0,
                label_fontfamily = "sans", label_size = 12,  label_x =c(0.15,0.15), label_y =1, ncol = 3, nrow=2)
#+ theme(plot.margin = unit(c(0.8,0,0,0), "cm"))
pp
pp <- plot_grid(plot_bar_community2[[4]],
                plot_bar_cummulation[[4]],
                
                hjust = 0,vjust = 0,
                label_fontfamily = "sans", label_size = 12,  label_x =c(0.15,0.15), label_y =1, ncol = 1, nrow=2)
#+ theme(plot.margin = unit(c(0.8,0,0,0), "cm"))
pp

#ggsave("SEM/H1_4/barplot/effect.final.jpg",pp ,width = 3, height = 6)
#ggsave("SEM/H1_4/barplot/effect.final.pdf",pp ,width = 3, height = 6)


pp <- plot_grid(plot_bar_community3[[1]],
                plot_bar_community3[[2]],
                plot_bar_community3[[3]],
                plot_bar_community3[[4]],

                hjust = 0,vjust = 0,
                label_fontfamily = "sans", label_size = 12,  label_x =c(0.15,0.15), label_y =1, ncol = 3, nrow=2)
#+ theme(plot.margin = unit(c(0.8,0,0,0), "cm"))
pp
#ggsave("SEM/H1_4/barplot/effect_summary_data2.jpg",pp ,width = 6, height = 6)
#ggsave("SEM/H1_4/barplot/effect_summary_data2.pdf",pp ,width = 6, height = 6)
#################################################################################################### barplot for effects of communities on function

temp=data.bar.comm[data.bar.comm$depend=="Func.Stab",]
re=NULL
for(i in unique(temp$independ.SEMi)){
  temp1=temp[temp$independ.SEMi==i,]
  indirect_percent=temp1[1:3,]$effects/sum(temp1[1:3,]$effects)*2
  names(indirect_percent)=c("total","indirect","direct")
  re[[i]]=indirect_percent
}
re=as.data.frame(re)
re

temp=data.bar.comm[data.bar.comm$depend=="Func.Stab",]
re=NULL
for(i in unique(temp$independ.SEMi)){
  temp1=temp[temp$independ.SEMi==i,]
  indirect_percent=temp1$effects/sum(temp1[1:3,]$effects)*2
  names(indirect_percent)=c("total.effects","indirect.effects","direct.effects","indirect.by_Pl","indirect.by_LF","indirect.by_SF","indirect.by_interac")
  re[[i]]=indirect_percent
}
re=as.data.frame(re)
re
