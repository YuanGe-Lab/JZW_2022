# Please run step1--stability_calcuation.R

library(lavaan)
library(semPlot)
library(vegan)
source("R_sources/other_source.R")
source("R_sources/permSEM_source.R")

reps=1000
RESULT_unique=NULL
for(iiii in c("Pl.mass","light","Pl.do","SF.mass","LF.mass","Soil.C","Soil.P","Soil.N","Glomalin","k","Glucosidase","Protease","Nitratase","Dehydrogenase")){
  data_sem=list(cbind(invariability,univariate_stability$space.stab),
                cbind(resist_to.D,univariate_stability$resist.D),
                cbind(resist_to.I,univariate_stability$resist.I))
  data=NULL
  for (datai in 1:3){
    data[[datai]]=as.data.frame(data_sem[[datai]])
    #data[[datai]]$PR=2^data[[datai]]$PR
    data[[datai]]=as.data.frame(apply(data[[datai]],2,stand))
    #data[[datai]]$NPP=data[[datai]]$Pl.mass
    data[[datai]]$Multifunc=data[[datai]][[iiii]]
  }
  
  ###################  buid initial sem model
  model.1 <-
    '
SF + LF + Multifunc ~ D + I + PR
Pl~D+I+PR
Multifunc ~ SF + LF + Pl
SF + LF ~ Pl
'

  model.1=add.variance(model.1)

  Group_Effects=NULL
  # reps = 10000
  sample.n.dots=min(reps,40)
  ########################stability
  df=df1=data[[1]]
  df$LF.div=apply(df[,c("LF.PD","LF.SR")],1,mean)
  df$SF.div=apply(df[,c("SF.PD","SF.SR")],1,mean)
  df$Pl.div=apply(df[,c("Pl.PD","Pl.SR")],1,mean)
  #df$Pl[!is.na(df$Pl)]=resid(lm(df$Pl~df$PR))
  
  
  names(df)
  apply(is.na(df),2,sum)
  apply(df,2,sd,na.rm=T)
  #df$SF=df$SF*10
  #df$LF=df$LF*10
  
  saturated_fit=lavaan(model = model.1 ,data = df)
  model.1=delete.variables.in.model(saturated_fit,delete.pars = c("NPP","","Pll" , "SF.mass" , "LF.div" , "LF.mass","PR.dff"))
  saturated_fit=lavaan(model = model.1 ,data = df)
  fit=fit1=lavaan(model = model.1, data = df, ordered = NULL,
                  sampling.weights   = NULL,
                  sample.cov = NULL, sample.mean = NULL, sample.th = NULL,
                  sample.nobs = NULL,
                  group = NULL, cluster = NULL)
  summary(fit)
  fitMeasures(fit)
  fitMeasures(fit,c("chisq","df","pvalue ", "cfi", "rmsea","rmsea.pvalue"))
  modificationIndices(fit,sort.=TRUE,minimum.value= 1)
  
  #semPaths(fit,layout = "circle2",what = "est",whatLabels = "std.all",style="ram",residuals=T,nCharNodes=6,
  #sizeMan=10,sizeMan2 = 5,sizeLat=20,sizeInt=10)
  semPaths(fit,layout = "circle2",whatLabels = "std.all")
  
  #MD=find.best.fit.model(fit,data = df,op="~~")
  #MD1=MD$best.fit.Model
  #reg=MD$Stepwise.regression
  #fit = bestfit1 = lavaan(model = MD$best.fit.Model, data = df)
  #summary(fit)
  #fitMeasures(fit,c("chisq","df","pvalue ", "cfi", "rmsea","rmsea.pvalue"))
  #semPaths(fit,layout = "circle2",whatLabels = "std.all")
  
  effects.stat=sem.effects.stat(fit = fit,data = df)
  effects.stat$std.effects.summary
  effects.stat$effects.summary
  
  Group_Effects$Invariability=
    grouped.sem.effects.stat.permute(fit,data=df,effect.type=c("total.effects", "direct.effects", "indirect.effects"),
                                     group.rhs=c(D="D",
                                                 I="I",
                                                 PR="PR",
                                                 #Comm.Stab="Pl + LF + SF",
                                                 Pl="Pl",LF="LF",SF="SF"),
                                     group.lhs=c(Pl="Pl",LF="LF",SF="SF",NPP="NPP",
                                                 Comm.Stab="Pl + LF + SF",
                                                 Func.Stab="Multifunc"),
                                     via=NULL,
                                     rh=NULL,lh=NULL,standardized = TRUE,
                                     bootstrap=T,reps=reps,return.boot=T,boot.depth = NULL,permutation = T,
                                     ind.in.dist = invariability$IND,dist.list = data_SEM.list,boot.perm.simultaneous = F)
  
  eff=Group_Effects$Invariability
  eff$std.effects.summary
  eff$std_effects
  
  rown=rownames(eff$std_effects$total.effects)
  coln=colnames(eff$std_effects$total.effects)
  data.boot=databarplot=NULL
  
  
  for(independ in c("D","I","PR")){
    for(depend in c("Func.Stab","Comm.Stab","Pl","SF","LF","NPP")){
      independ.i=which(coln==independ)
      depend.i=which(rown==depend)
      for(effect.i in c("total.effects", "indirect.effects","direct.effects")){
        
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
  for(i in c("total.effects", "indirect.effects")){
    temp.data=data.boot[data.boot$effect.type==i&!data.boot$depend%in%c("Func.Stab","Comm.Stab","NPP"),];temp.data$color=temp.data$depend
    dataforjiter=temp.data
    dataforjiter[dataforjiter$depend%in%c("Pl","SF","LF"),]$depend="Comm.Stab"
    
    #dataforjiter[dataforjiter$depend%in%c("Pl","SF","LF"),]$color="Comm.Stab"
    #dataforjiter[dataforjiter$color%in%c("MS1","MS2","MS3"),]$color=dataforjiter[dataforjiter$color%in%c("MF1","MF2","MF3"),]$color
    #dataforjiter[dataforjiter$depend%in%c("MF1","MF2","MF3"),]$depend="Multi.Func"
    
    
    temp.data.bar=databarplot[databarplot$effect.type==i,]
    log.func=temp.data.bar$depend%in%c("Func.Stab","Comm.Stab","NPP")
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
                  hjust = 0,vjust = 0,
                  label_fontfamily = "sans", label_size = 12,  label_x =c(0.15,0.15), label_y =1, ncol = 2, nrow=1)
  #+ theme(plot.margin = unit(c(0.8,0,0,0), "cm"))
  pp
  #mkdir("SEM/H1_4/barplot")
  #ggsave("SEM/H1_4/barplot/barplot_invariability.jpg",pp ,width = 6, height = 2)
  #ggsave("SEM/H1_4/barplot/barplot_invariability.pdf",pp ,width = 6, height = 2)
  
  
  
  ####################################################resist_to.D
  model.to.D=delete.variables.in.model(saturated_fit,delete.pars = "D")
  
  df=df2=data[[2]]
  names(df)
  
  fit = fit2 = lavaan(model = model.to.D ,data = df)#,test ="bootstrap",verbose=TRUE)
  summary(fit)
  fitMeasures(fit)
  fitMeasures(fit,c("chisq","df","pvalue ", "cfi", "rmsea","rmsea.pvalue"))
  modificationIndices(fit,sort.=TRUE,minimum.value= 0)
  
  semPaths(fit,layout = "circle2",what = "est",whatLabels = "std.all",style="ram",residuals=T,nCharNodes=6,
           sizeMan=10,sizeMan2 = 5,sizeLat=20,sizeInt=10)
  semPaths(fit,layout = "circle2",whatLabels = "std.all")
  
  #MD=find.best.fit.model(fit,data = df,op="~~")
  #MD2=MD$best.fit.Model
  #reg=MD$Stepwise.regression
  #fit=bestfit2=lavaan(model = MD$best.fit.Model, data = df)
  #summary(fit)
  #fitMeasures(fit,c("chisq","df","pvalue ", "cfi", "rmsea","rmsea.pvalue"))
  #semPaths(fit,layout = "circle2",whatLabels = "std.all")
  
  effects.stat=sem.effects.stat(fit = fit,data = df)
  effects.stat$std.effects.summary
  
  Group_Effects$Resistance_to.D=
    grouped.sem.effects.stat.permute(fit,data=df,effect.type=c("total.effects", "direct.effects", "indirect.effects"),
                                     group.rhs=c(D="D",
                                                 I="I",
                                                 PR="PR",
                                                 #Comm.Stab="Pl + LF + SF",
                                                 Pl="Pl",LF="LF",SF="SF"),
                                     group.lhs=c(Pl="Pl",LF="LF",SF="SF",NPP="NPP",
                                                 Comm.Stab="Pl + LF + SF",
                                                 Func.Stab="Multifunc" ),
                                     via=NULL,
                                     rh=NULL,lh=NULL,standardized = TRUE,
                                     bootstrap=T,reps=reps,return.boot=T,boot.depth = NULL,permutation = T,
                                     ind.in.dist = resist_to.D$IND,dist.list = data_SEM.list,boot.perm.simultaneous = F)
  
  
  eff=Group_Effects$Resistance_to.D
  eff$std.effects.summary
  eff$std_effects
  
  rown=rownames(eff$std_effects$total.effects)
  coln=colnames(eff$std_effects$total.effects)
  data.boot=databarplot=NULL
  
  
  for(independ in c("D","I","PR")){
    for(depend in c("Func.Stab","Comm.Stab","Pl","SF","LF","NPP")){
      independ.i=which(coln==independ)
      depend.i=which(rown==depend)
      for(effect.i in c("total.effects", "indirect.effects","direct.effects")){
        
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
  for(i in c("total.effects", "indirect.effects")){
    temp.data=data.boot[data.boot$effect.type==i&!data.boot$depend%in%c("Func.Stab","Comm.Stab","NPP"),];temp.data$color=temp.data$depend
    dataforjiter=temp.data
    dataforjiter[dataforjiter$depend%in%c("Pl","SF","LF"),]$depend="Comm.Stab"
    
    #dataforjiter[dataforjiter$depend%in%c("Pl","SF","LF"),]$color="Comm.Stab"
    #dataforjiter[dataforjiter$color%in%c("MS1","MS2","MS3"),]$color=dataforjiter[dataforjiter$color%in%c("MF1","MF2","MF3"),]$color
    #dataforjiter[dataforjiter$depend%in%c("MF1","MF2","MF3"),]$depend="Multi.Func"
    
    
    temp.data.bar=databarplot[databarplot$effect.type==i,]
    log.func=temp.data.bar$depend%in%c("Func.Stab","Comm.Stab","NPP")
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
      labs(x=element_blank(), y = paste0(i," (SEM 1)"))+
      scale_y_continuous(breaks = sort(c(-seq(bk.y,5,bk.y),seq(0,5,bk.y))))+
      theme_bw()+
      theme(text=element_text(size = 8, family = "sans", face = "plain"),axis.title.x = element_blank(),
            panel.grid.major.x = element_blank(),    panel.grid.minor.x =  element_blank(),  panel.grid.minor.y =  element_blank(),
            panel.background = element_rect(color = 'black', fill = 'transparent'),legend.position = "none")
  }
  pp <- plot_grid(plot_bar2[[1]],plot_bar2[[2]],
                  hjust = 0,vjust = 0,
                  label_fontfamily = "sans", label_size = 12,  label_x =c(0.15,0.15), label_y =1, ncol = 2, nrow=1)
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
  fit = fit3 = lavaan(model = model.to.I ,data = df) #,test ="bootstrap",verbose=TRUE)
  summary(fit)
  fitMeasures(fit)
  fitMeasures(fit,c("chisq","df","pvalue ", "cfi", "rmsea","rmsea.pvalue"))
  modificationIndices(fit,sort.=TRUE,minimum.value= 0)
  
  semPaths(fit,layout = "circle2",what = "est",whatLabels = "std.all",style="ram",residuals=T,nCharNodes=6,
           sizeMan=10,sizeMan2 = 5,sizeLat=20,sizeInt=10)
  semPaths(fit,layout = "circle2",whatLabels = "std.all")
  
  #MD=find.best.fit.model(fit,data = df,op="~~")
  #MD3=MD$best.fit.Model
  #reg=MD$Stepwise.regression
  #fit=bestfit3=lavaan(model = MD$best.fit.Model, data = df)
  #summary(fit)
  #fitMeasures(fit,c("chisq","df","pvalue ", "cfi", "rmsea","rmsea.pvalue"))
  
  effects.stat=sem.effects.stat(fit = fit,data = df)
  effects.stat$std.effects.summary
  
  Group_Effects$Resistance_to.I=
    grouped.sem.effects.stat.permute(fit = fit,data=df,effect.type=c("total.effects", "direct.effects", "indirect.effects"),
                                     group.rhs=c(D="D",
                                                 I="I",
                                                 PR="PR",
                                                 #Comm.Stab="Pl + LF + SF",
                                                 Pl="Pl",LF="LF",SF="SF"),
                                     group.lhs=c(Pl="Pl",LF="LF",SF="SF",NPP="NPP",
                                                 Comm.Stab="Pl + LF + SF",
                                                 Func.Stab="Multifunc"),
                                     via=NULL,
                                     rh=NULL,lh=NULL,standardized = TRUE,
                                     bootstrap=T,reps=reps,return.boot=T,boot.depth = NULL,permutation = T,
                                     ind.in.dist = resist_to.I$IND,dist.list = data_SEM.list,boot.perm.simultaneous = F)
  
  
  eff=Group_Effects$Resistance_to.I
  eff$std.effects.summary
  eff$std_effects
  
  rown=rownames(eff$std_effects$total.effects)
  coln=colnames(eff$std_effects$total.effects)
  data.boot=databarplot=NULL
  
  
  for(independ in c("D","I","PR")){
    for(depend in c("Func.Stab","Comm.Stab","Pl","SF","LF","NPP")){
      independ.i=which(coln==independ)
      depend.i=which(rown==depend)
      for(effect.i in c("total.effects", "indirect.effects","direct.effects")){
        
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
  for(i in c("total.effects", "indirect.effects")){
    temp.data=data.boot[data.boot$effect.type==i&!data.boot$depend%in%c("Func.Stab","Comm.Stab","NPP"),];temp.data$color=temp.data$depend
    dataforjiter=temp.data
    dataforjiter[dataforjiter$depend%in%c("Pl","SF","LF"),]$depend="Comm.Stab"
    
    #dataforjiter[dataforjiter$depend%in%c("Pl","SF","LF"),]$color="Comm.Stab"
    #dataforjiter[dataforjiter$color%in%c("MS1","MS2","MS3"),]$color=dataforjiter[dataforjiter$color%in%c("MF1","MF2","MF3"),]$color
    #dataforjiter[dataforjiter$depend%in%c("MF1","MF2","MF3"),]$depend="Multi.Func"
    
    
    temp.data.bar=databarplot[databarplot$effect.type==i,]
    log.func=temp.data.bar$depend%in%c("Func.Stab","Comm.Stab","NPP")
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
                  hjust = 0,vjust = 0,
                  label_fontfamily = "sans", label_size = 12,  label_x =c(0.15,0.15), label_y =1, ncol = 2, nrow=1)
  #+ theme(plot.margin = unit(c(0.8,0,0,0), "cm"))
  pp
  mkdir("SEM/H1_4/barplot")
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
  mkdir("SEM/H1_4/barplot/")
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
  1.195963e-02
  
  for(independ in c("D","I","PR")){
    for(depend in c("Func.Stab","Comm.Stab","Pl","SF","LF","NPP")){
      independ.i=which(coln==independ)
      depend.i=which(rown==depend)
      for(effect.i in c("total.effects", "indirect.effects","direct.effects")){
        
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
  for(i in c("total.effects", "indirect.effects")){
    temp.data=data.boot[data.boot$effect.type==i&!data.boot$depend%in%c("Func.Stab","Comm.Stab","NPP"),];temp.data$color=temp.data$depend
    dataforjiter=temp.data
    dataforjiter[dataforjiter$depend%in%c("Pl","SF","LF"),]$depend="Comm.Stab"
    
    #dataforjiter[dataforjiter$depend%in%c("Pl","SF","LF"),]$color="Comm.Stab"
    #dataforjiter[dataforjiter$color%in%c("MS1","MS2","MS3"),]$color=dataforjiter[dataforjiter$color%in%c("MF1","MF2","MF3"),]$color
    #dataforjiter[dataforjiter$depend%in%c("MF1","MF2","MF3"),]$depend="Multi.Func"
    
    
    temp.data.bar=databarplot[databarplot$effect.type==i,]
    log.func=temp.data.bar$depend%in%c("Func.Stab","Comm.Stab","NPP")
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
                  hjust = 0,vjust = 0,
                  label_fontfamily = "sans", label_size = 12,  label_x =c(0.15,0.15), label_y =1, ncol = 2, nrow=1)
  #+ theme(plot.margin = unit(c(0.8,0,0,0), "cm"))
  pp
  mkdir("SEM/H1_4/barplot")
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
  mkdir("SEM/H1_4/barplot/")
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
      for(depend in c("Func.Stab","Comm.Stab","Pl","SF","LF","NPP")){
        independ.i=which(coln==independ)
        depend.i=which(rown==depend)
        for(effect.i in c("total.effects", "indirect.effects","direct.effects")){
          
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
  #for(i in c("total.effects", "indirect.effects","direct.effects")){
  #temp.data=data.boot[data.boot$effect.type==i&!data.boot$depend%in%c("Func.Stab","Comm.Stab","NPP"),];temp.data$color=temp.data$depend
  temp.data=data.boot[!data.boot$depend%in%c("Func.Stab","Comm.Stab","NPP"),];temp.data$color=temp.data$depend
  
  dataforjiter=temp.data
  dataforjiter[dataforjiter$depend%in%c("Pl","SF","LF"),]$depend="Comm.Stab"
  
  #dataforjiter[dataforjiter$depend%in%c("Pl","SF","LF"),]$color="Comm.Stab"
  #dataforjiter[dataforjiter$color%in%c("MS1","MS2","MS3"),]$color=dataforjiter[dataforjiter$color%in%c("MF1","MF2","MF3"),]$color
  #dataforjiter[dataforjiter$depend%in%c("MF1","MF2","MF3"),]$depend="Multi.Func"
  
  
  #temp.data.bar=databarplot[databarplot$effect.type==i,]
  temp.data.bar=databarplot
  log.func=temp.data.bar$depend%in%c("Comm.Stab")
  data.bar.func=temp.data.bar[log.func,]
  data.bar.comm=temp.data.bar[!log.func,]
  data.bar.comm$depend=factor(data.bar.comm$depend,levels = c("Pl","SF","LF","Func.Stab","NPP" ))
  data.bar.comm$depend=paste(as.numeric(data.bar.comm$depend),data.bar.comm$depend)
  
  data.bar.comm$effect.type=factor(data.bar.comm$effect.type,levels = c("total.effects","indirect.effects","direct.effects"))
  #temp.slop.p=slope.p[slope.p$group==i&slope.p$effect.type!="direct.effects",]
  for(j in c("D","I","PR")){
    judge_total=data.bar.comm$independ==j&data.bar.comm$effect.type=="total.effects"
    max.gap=max(c(0,data.bar.comm[judge_total,]$effects+data.bar.comm[judge_total,]$ci),na.rm=T)-min(c(data.bar.comm[judge_total,]$effects-data.bar.comm[judge_total,]$ci,0),na.rm=T)
    bk.y=max(c(round(max.gap/4,digits = 1),0.02),na.rm=T)
    
    
    plot_bar_community1[[j]] = ggplot(data.bar.comm[data.bar.comm$independ==j,],aes(x = independ.SEMi,y = effects))+
      facet_wrap(~paste(depend),nrow=1)+
      #geom_smooth(data=temp.data,aes(x = PR,y = effects),method = "lm",se=F)+
      #geom_jitter(data = dataforjiter[dataforjiter$independ==j,],na.rm = T,
      #aes(x = independ.SEMi,y = effects,col=color),height=0,width = 0.2,size=0.1,shape=16,alpha=0.8)+
      geom_bar(aes(x = independ.SEMi,y = effects,fill=effect.type),stat='identity',position='stack',color=NA,size=0.25,alpha=0.25,
               data = data.bar.comm[data.bar.comm$independ==j&!judge_total,])+
      geom_bar(aes(x = independ.SEMi,y = effects),stat='identity',position='stack',fill=NA,color="black",size=0.25,
               data = data.bar.comm[judge_total,])+
      geom_errorbar(aes(ymin = effects-ci, ymax = effects + ci),width = 0.5,size=0.375,
                    data = data.bar.comm[judge_total,])+
      geom_text( aes(x=independ.SEMi, y=effects + sign(effects)*(ci+0.025*max.gap),
                     label=as.symbol.pvalue(p.value,sig.level=c(0.001,0.01,0.05))),size=4,
                 data = data.bar.comm[judge_total,])+
      scale_color_manual(values = c("#FF0000","#00FFFF","#0000FF","#FF00FF","#00FF00", "#FFFF00")[c(1,5,3,2,4,6)])+
      labs(x=element_blank(), y = paste0(j," effects"))+
      scale_y_continuous(breaks = sort(c(-seq(bk.y,5,bk.y),seq(0,5,bk.y))))+
      theme_bw()+
      theme(text=element_text(size = 8, family = "sans", face = "plain"),axis.title.x = element_blank(),
            panel.grid.major.x = element_blank(),    panel.grid.minor.x =  element_blank(),  panel.grid.minor.y =  element_blank(),
            panel.background = element_rect(color = 'black', fill = 'transparent'),legend.position = "none")
  }
  
  #}
  pp <- plot_grid(plot_bar_community1[[1]],
                  plot_bar_community1[[2]],
                  plot_bar_community1[[3]],
                  hjust = 0,vjust = 0,
                  label_fontfamily = "sans", label_size = 12,  label_x =c(0.15,0.15), label_y =1, ncol = 1, nrow=3)
  #+ theme(plot.margin = unit(c(0.8,0,0,0), "cm"))
  pp
  ############################################################################################################×ÝÏò±È½Ï
  plot_bar_community2=plot_bar_community3=NULL
  temp.data=data.boot[!data.boot$depend%in%c("Func.Stab","Comm.Stab","NPP"),];temp.data$color=temp.data$depend
  
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
  data.bar.comm$depend=factor(data.bar.comm$depend,levels = c("Pl","SF","LF","Func.Stab","NPP"))
  #data.bar.comm$depend=paste(as.numeric(data.bar.comm$depend),data.bar.comm$depend)
  
  data.bar.comm$effect.type=factor(data.bar.comm$effect.type,levels = c("total.effects","indirect.effects","direct.effects"))
  #temp.slop.p=slope.p[slope.p$group==i&slope.p$effect.type!="direct.effects",]
  for(j in sort(unique(data.bar.comm$depend))){
    judge_total=data.bar.comm$depend==j&data.bar.comm$effect.type=="total.effects"
    max.gap=max(c(0,data.bar.comm[judge_total,]$effects+data.bar.comm[judge_total,]$ci),na.rm=T)-min(c(data.bar.comm[judge_total,]$effects-data.bar.comm[judge_total,]$ci,0),na.rm=T)
    bk.y=max(c(round(max.gap/4,digits = 1),0.02),na.rm=T)
    
    
    plot_bar_community2[[j]] = ggplot(data.bar.comm[data.bar.comm$depend==j,],aes(x = independ.SEMi,y = effects))+
      #facet_wrap(~paste(independ),nrow=1)+
      #geom_smooth(data=temp.data,aes(x = PR,y = effects),method = "lm",se=F)+
      #geom_jitter(data = dataforjiter[dataforjiter$independ==j,],na.rm = T,
      #aes(x = independ.SEMi,y = effects,col=color),height=0,width = 0.2,size=0.1,shape=16,alpha=0.8)+
      geom_bar(aes(x = independ.SEMi,y = effects,fill=effect.type),stat='identity',position='stack',color=NA,size=0.25,alpha=0.25,
               data = data.bar.comm[data.bar.comm$depend==j&!judge_total,])+
      geom_bar(aes(x = independ.SEMi,y = effects),stat='identity',position='stack',fill=NA,color="black",size=0.25,
               data = data.bar.comm[judge_total,])+
      geom_errorbar(aes(ymin = effects-ci, ymax = effects + ci),width = 0.5,size=0.375,
                    data = data.bar.comm[judge_total,])+
      geom_text( aes(x=independ.SEMi, y=effects + sign(effects)*(ci+0.5*max.gap),
                     label=as.symbol.pvalue(p.value,sig.level=c(0.001,0.01,0.05))),size=4,
                 data = data.bar.comm[judge_total,])+
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
                  plot_bar_community2[[5]],
                  plot_bar_community2[[5]],
                  hjust = 0,vjust = 0,
                  label_fontfamily = "sans", label_size = 12,  label_x =c(0.15,0.15), label_y =1, ncol = 3, nrow=2)
  #+ theme(plot.margin = unit(c(0.8,0,0,0), "cm"))
  pp
  pp2 <- plot_grid(plot_bar_community3[[1]],
                  plot_bar_community3[[2]],
                  plot_bar_community3[[3]],
                  plot_bar_community3[[4]],
                  
                  hjust = 0,vjust = 0,
                  label_fontfamily = "sans", label_size = 12,  label_x =c(0.15,0.15), label_y =1, ncol = 3, nrow=2)
  #+ theme(plot.margin = unit(c(0.8,0,0,0), "cm"))
  pp2
  mkdir("SEM/H1_4/barplot/unique/")
  ggsave(paste("SEM/H1_4/barplot/unique/effect_summary",iiii,".jpg"),pp ,width = 6, height = 6)
  ggsave(paste("SEM/H1_4/barplot/unique/effect_summary",iiii,".pdf"),pp ,width = 6, height = 6)
  
  RESULT_unique[[iiii]]$models=list(df1=df1,fit1=fit1,df2=df2,fit2=fit2,df3=df3,fit3=fit3)                       
  RESULT_unique[[iiii]]$permre=Group_Effects
  RESULT_unique[[iiii]]$pp=pp
}
save(list = "RESULT_unique",file = "SEM/H1_4/barplot/unique/RESULT_unique.RData",compress = T)
direct=NULL
for(iiii in c("Pl.mass","Pl.do","light","SF.mass","LF.mass","Soil.C","Soil.P","Soil.N","Glomalin","k","Glucosidase","Protease","Nitratase","Dehydrogenase")){
  for(jjjj in names(RESULT_unique[[iiii]]$permre)){
    temp=RESULT_unique[[iiii]]$permre[[jjjj]]$std.effects.summary
    temp=t(as.matrix(temp$direct.effects[6,]));rownames(temp)=iiii
    direct[[jjjj]]= rbind(direct[[jjjj]],temp) 
  }
}
direct
write.csv("UniFuncStab_direct.effects",paste0("SEM/H1_4/barplot/unique/UniFuncStab_direct.effects.csv"))
#write.table(FM,"UniFuncStab_direct.effects.csv",sep = ",",append = T, col.names = NA)
for(i in names(direct)){
  write.table(i,"SEM/H1_4/barplot/unique/UniFuncStab_direct.effects.csv",sep = ",",append = T)
  write.table(direct[[i]],"SEM/H1_4/barplot/unique/UniFuncStab_direct.effects.csv",sep = ",",append = T, col.names = NA)
}
for(i in names(direct)){
  write.csv(direct[[i]],paste0("SEM/H1_4/barplot/unique/paths.",i,".csv"))
}
