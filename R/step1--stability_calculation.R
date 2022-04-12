
library(vegan)
library(GUniFrac)
library(picante)
library(phangorn)
library(ggpmisc)
library(cowplot)
source("R_sources/other_source.R")
source("R_sources/permSEM_source.R")
order_fauna=F
same_PS=F
SQRT_biomass=F
std.Method="no"

map=read.csv("data_files/map.csv",row.names = 1)
nrow(map)
tabs = list(Pl=NULL,SF=NULL,LF=NULL)
tabs$Pl=read.csv("data_files/plant_comm.csv", row.names = 1, header=T)
tabs$SF = read.csv("data_files/soil_fauna_comm.csv", row.names = 1, header=T)
tabs$LF = read.csv("data_files/litter_fauna_comm.csv", row.names = 1, header=T)

ncol(tabs$Pl)
ncol(tabs$SF)
ncol(tabs$LF)

trees=NULL
trees$Pl=midpoint(read.tree("data_files/cds_plant.tre"))
trees$SF=midpoint(read.tree("data_files/cds_fauna.tre"))
trees$LF=trees$SF

a=NULL
for(i in names(tabs)){
  tabs[[i]]=tabs[[i]][,rownames(map)]
  ma=match.phylo.comm(trees[[i]],t(tabs[[i]]))
  trees[[i]] = ma$phy; tabs[[i]]=t(ma$comm)
  a=c(a,all(apply(tabs[[i]],2,sum) > 0 ))
}
a

dim(tabs$Pl)
dim(tabs$SF)
dim(tabs$LF)


alpha=NULL
alpha_all=read.csv("data_files/alpha_all.csv",row.names = 1)
alpha_all=alpha_all[rownames(map),]
for(i in names(tabs)){
  alpha[[i]]=picante::pd(samp = t(tabs[[i]]), tree =  trees[[i]], include.root = T)
  colnames(alpha[[i]])=paste(i, colnames(alpha[[i]]), sep=".")
  for(j in colnames(alpha[[i]])){
    alpha_all[[j]] = 0
    alpha_all[rownames(alpha[[i]]),j] = alpha[[i]][,j]
  }
  alpha_all[[paste(i, "mass", sep=".")]] = NA
  alpha_all[rownames(alpha[[i]]),paste(i, "mass", sep=".")] = apply(tabs[[i]],2,sum)
}



names(alpha_all)
stand=function(x) {x=as.matrix(x);return((x-mean(x,na.rm=T))/sd(x,na.rm=T))}
stand.range=function(x) {x=as.matrix(x);(x-min(x,na.rm=T))/sd(x,na.rm=T)}
alpha_all=alpha_all[,c("Pl.mass","Pl.do","In.up","In.do","SF.mass","LF.mass","Pl.PD","Pl.SR","SF.PD","SF.SR","LF.PD","LF.SR")]

alpha_all$Pl.tt=alpha_all$Pl.mass+alpha_all$Pl.do

all(rownames(alpha_all) == rownames(map))

ENV=read.csv("data_files/ENV_variables.csv",row.names = 1)
ENV=ENV[rownames(map),]
all(rownames(alpha_all) ==rownames(map))
alpha_all=cbind(alpha_all,as.data.frame(ENV))

uniq_functions = c("Pl.mass","light","Pl.do","SF.mass","LF.mass","Soil.C","Soil.P","Soil.N","Glomalin","k","Glucosidase","Protease","Nitratase","Dehydrogenase")
  
all(uniq_functions %in% names(alpha_all))

if(SQRT_biomass){alpha_all[,c("Pl.tt","Pl.do","Pl.mass","SF.mass","LF.mass")]=
  sqrt(alpha_all[,c("Pl.tt","Pl.do","Pl.mass","SF.mass","LF.mass")])}
if(!std.Method == "no"){
  alpha_all=decostand(alpha_all,method =std.Method,MARGIN = 2)
}

apply(alpha_all[,uniq_functions],2,sd,na.rm=T)
apply(alpha_all[,uniq_functions],2,mean,na.rm=T)


similarities=NULL
mat_temp=matrix(0,nrow(map),nrow(map),dimnames = list(rownames(map),rownames(map)))
similarities=NULL
for(i in names(tabs)){
  unifracs=UniFrac(if(SQRT_biomass){sqrt(tabs[[i]])}else{tabs[[i]]},trees[[i]],total.normalize = F)
  similarities[[i]]=1-as.matrix(unifracs$w_unifrac[,,1])
  mat_temp[]=0;mat_temp[rownames(similarities[[i]]),colnames(similarities[[i]])]=similarities[[i]][]
  mat_temp[!rownames(mat_temp)%in%colnames(tabs[[i]]),!rownames(mat_temp)%in%colnames(tabs[[i]])]=0
  similarities[[i]]=mat_temp
}

multi.variate.mat=list(  Multifunc=alpha_all[,uniq_functions] )
apply(multi.variate.mat$Multifunc,2,mean,na.rm=T)

stability.fun=function(df,method=c("global","univariate.mean"),
                       na.action=c("na.rm","resample"),reps=1000,dist.method="bray",
                       stand.method=c("stand.range","mean.as.weight","z.score"),weight=NULL){
  if(is.null(weight)){weight=1}
  if(!is.null(stand.method)){
    if(stand.method[[1]]=="stand.range"){
      df=t(t(apply(df,2,stand.range))*weight)
    }else if(stand.method[[1]]=="mean.as.weight"){
      df=t(t(apply(df,2,function(x) {if(all(is.na(x))){NA}else{x/mean(x,na.rm=TRUE)}}))*weight)
    }else if(stand.method[[1]]=="z.score"){
      df=t(t(apply(df,2,stand))*weight)
      dist.method="euclidean"
    }
  }
  univariate.stab <- function(univariates){
    temp_sim_mat <- matrix(NA,nrow(univariates),nrow(univariates),dimnames = list(rownames(univariates),rownames(univariates)))
    univariate_similarity <- array(NA,dim = c(nrow(univariates),nrow(univariates),ncol(univariates)),
                                   dimnames = list(rownames(univariates),rownames(univariates),colnames(univariates)))
    for (i in colnames(univariates)){
      temp_sim_mat[]<-NA
      
      temp1 <- univariates[,i];names(temp1) <- rownames(univariates)
      
      temp1<-temp1[rownames(univariates)]
      log.na <- is.na(temp1)
      log.0 <- !log.na&temp1 == 0
      temp_sim_mat[log.0,!log.na]<-temp_sim_mat[!log.na,log.0]<-0
      temp_sim_mat[log.0,log.0] <- 0
      
      temp2 <- temp1[!log.na&!log.0]
      square_mat <- row.mat(temp2)
      temp_sim_mat[!log.na&!log.0,!log.na&!log.0] <- 1-abs(square_mat-t(square_mat))/(square_mat+t(square_mat))
      univariate_similarity[,,i] <- temp_sim_mat
    }
    return(univariate_similarity)
  }
  if(na.action[1]=="resample"){
    rand.df=array(NA,dim=c(nrow(df),ncol(df),reps),dimnames = list(rownames(df),colnames(df),1:reps))
    rand.stab=array(NA,dim=c(nrow(df),nrow(df),reps),dimnames = list(rownames(df),rownames(df),1:reps))
    for(i in 1:reps){
      rand.df[,,i]=apply(df,2,missing.data.process,method="resample.self")
    }
    
    if(method[1]=="global"){
      for(i in 1:reps){
        rand.stab[,,i]=1-as.matrix(vegdist(rand.df[,,i],method=dist.method))
      }
      stab=apply(rand.stab,c(1,2),mean)
    }else if(method[1]=="univariate.mean"){
      for(i in 1:reps){
        print(i)
        rand.stab[,,i]=apply(univariate.stab(rand.df[,,i]),c(1,2),mean,na.rm=T)
      }
      stab=apply(rand.stab,c(1,2),mean,na.rm=TRUE)
    }
  }else if(na.action[1]=="na.rm"){
    if(method[1]=="global"){
      stab=1-as.matrix(vegdist(df,method=dist.method,na.rm=TRUE))
    }else if(method[1]=="univariate.mean"){
      stab=apply(univariate.stab(df),c(1,2),mean,na.rm=T)
    }
  }
  return(stab)
}

for(i in names(multi.variate.mat)){
  similarities[[i]]=stability.fun(multi.variate.mat[[i]],na.action="resample",reps = 100,
                                        stand.method="stand.range",
                                        weight=NULL)
}
#####################################################calculating multifunction
data_SEM.list=c(similarities)
names(data_SEM.list)

invariability=resist_to.I=resist_to.D=stab=NULL
abs.PR = 0.5 # if (abs.PR < 1 ); then similarity with same PSR level will be calculated
same_PS = FALSE # if(same_PS); then similarity will be calculated in univariate (one-to-one) correspondent way; otherwise in multivariate way
for (i in names(data_SEM.list)){
  substability_mat=data_SEM.list[[i]]
  D=map$D
  I=map$I
  PR=map$PR
  
  judge_low=lower.tri(substability_mat)
  I.row=row.mat(I)
  D.row=row.mat(D)
  PR.row=row.mat(PR)
  sort.row=row.mat(1:nrow(map))
  judge.low.up=!sort.row==t(sort.row)
  PS_row=row.mat(map$PS)
  IND=matrix(1:length(data_SEM.list[[1]]),nrow = nrow(data_SEM.list[[1]]))
  
  substability=data.frame(dist=substability_mat[judge.low.up],
                          PR.row=PR.row[judge.low.up],PR.col=t(PR.row)[judge.low.up],
                          D.row=D.row[judge.low.up],D.col=t(D.row)[judge.low.up],
                          I.row=I.row[judge.low.up],I.col=t(I.row)[judge.low.up],
                          n.row=sort.row[judge.low.up],n.col=t(sort.row)[judge.low.up],
                          PS.row=PS_row[judge.low.up],PS.col=t(PS_row)[judge.low.up],
                          IND=IND[judge.low.up])

  judge_PS=if(same_PS){substability$PS.row==substability$PS.col}else{TRUE}
  
  within=substability[abs(substability$PR.row-substability$PR.col)<abs.PR&
                        substability$I.row==substability$I.col&
                        substability$D.row==substability$D.col&
                        substability$n.col>substability$n.row,]
  
  bet.I=substability[abs(substability$PR.row-substability$PR.col)<abs.PR&
                       (substability$I.row-substability$I.col)==-1&
                       substability$D.row==substability$D.col&
                       judge_PS,]
  
  bet.D=substability[abs(substability$PR.row-substability$PR.col)<abs.PR&
                       substability$I.row==substability$I.col&
                       (substability$D.row==0&substability$D.col>0)&
                       judge_PS,]
  
  invariability[[i]]=within$dist; invariability.IND=data.frame(D=within$D.col,I=within$I.col,PR=(within$PR.col+within$PR.row)/2,
                                                               IND=within$IND,PR.dff=abs(within$PR.col-within$PR.row))
  resist_to.I[[i]]=bet.I$dist; resist_to.I.IND=data.frame(D=bet.I$D.col,I=abs(bet.I$I.col-bet.I$I.row),PR=(bet.I$PR.col+bet.I$PR.row)/2,
                                                          IND=bet.I$IND,PR.dff=abs(bet.I$PR.col-bet.I$PR.row))
  resist_to.D[[i]]=bet.D$dist; resist_to.D.IND=data.frame(D=abs(bet.D$D.col-bet.D$D.row),I=bet.D$I.col,PR=(bet.D$PR.col+bet.D$PR.row)/2,
                                                          IND=bet.D$IND,PR.dff=abs(bet.D$PR.col-bet.D$PR.row))
  stab[[i]]=substability$dist; stab.IND=data.frame(D.dff=abs(substability$D.col-substability$D.row),I.dff=abs(substability$I.col-substability$I.row),PR.dff=abs(substability$PR.col-substability$PR.row),
                                                   D=abs(substability$D.col+substability$D.row)/2,I=abs(substability$I.col+substability$I.row)/2,PR=abs(substability$PR.col+substability$PR.row)/2)
}

invariability=cbind(as.data.frame(invariability),invariability.IND)
resist_to.I=cbind(as.data.frame(resist_to.I),resist_to.I.IND)
resist_to.D=cbind(as.data.frame(resist_to.D),resist_to.D.IND)
stab=cbind(as.data.frame(stab),stab.IND)

univariates=alpha_all
names(univariates)

sum((!is.na(univariates)&univariates<0))
sum((!is.na(univariates)&univariates==0))
sum(is.na(univariates))
apply(is.na(univariates),2,sum)

mean_univarriate_list=NULL
univariate_similarity_list=NULL
temp_sim_mat=matrix(NA,nrow(map),nrow(map),dimnames = list(rownames(map),rownames(map)))

for (i in names(univariates)){
  temp_sim_mat[]<-NA
  temp_mean.mat<-temp_sim_mat
  temp1=univariates[[i]];names(temp1)=rownames(univariates)
  
  temp1<-temp1[rownames(map)]
  log.na=is.na(temp1)
  log.0=!log.na&temp1==0
  temp_sim_mat[log.0,!log.na]<-temp_sim_mat[!log.na,log.0]<-0
  temp_sim_mat[log.0,log.0]=NA
  
  temp2<-temp1[!log.na&!log.0]
  square_mat=row.mat(temp2)
  if(!i%in%c("In.do.fail","In.up.fail")){
    temp_sim_mat[!log.na&!log.0,!log.na&!log.0]=1-abs(square_mat-t(square_mat))/(square_mat+t(square_mat))
    univariate_similarity_list[[i]]=temp_sim_mat
  }

  ###################
  if(!i%in%c("In.do","In.up")){
    temp_mean.mat[]<-(row.mat(temp1)+t(row.mat(temp1)))/2
    mean_univarriate_list[[i]]=temp_mean.mat
  }

}

#######################invariability
univariate_stability=NULL
univariate_mean=NULL
for (i in names(univariate_similarity_list)){
  univariate_stability[["space.stab"]][[i]]=univariate_similarity_list[[i]][invariability.IND$IND]
  univariate_stability[["resist.D"]][[i]]=univariate_similarity_list[[i]][resist_to.D.IND$IND]
  univariate_stability[["resist.I"]][[i]]=univariate_similarity_list[[i]][resist_to.I.IND$IND]
}
for(i in  names(univariate_similarity_list)){
  univariate_mean[["space.stab"]][[i]]=mean_univarriate_list[[i]][invariability.IND$IND]
  univariate_mean[["resist.D"]][[i]]=mean_univarriate_list[[i]][resist_to.D.IND$IND]
  univariate_mean[["resist.I"]][[i]]=mean_univarriate_list[[i]][resist_to.I.IND$IND]
}
for(i in names(univariate_similarity_list)){
  univariate_stability[[i]] = as.data.frame(univariate_stability[[i]])
  univariate_mean[[i]] = as.data.frame(univariate_mean[[i]])
}

names(univariate_stability)
names(univariate_stability$space.stab)

all(rownames(alpha_all)==rownames(map))

nrow(map)

####################### NPERMANOVA ###################
reps = 9999
permANOVA=NULL
for(i in names(data_SEM.list)){
  a=1-data_SEM.list[[i]]
  a[is.na(a)] = 1
  permANOVA[[i]]=adonis2(as.dist(a)~D*I*PR,data=map,permutations = reps)
}
permANOVA$Pl
permANOVA$SF
permANOVA$LF
permANOVA$Multifunc

####################### Mantel test ###################
mantel.re = list()
for(i in names(data_SEM.list)[1:3]){
  mantel.re[[paste0(i,".vs.","MultiFunc")]] = vegan::mantel(xdis = data_SEM.list[[i]], ydis = data_SEM.list[["Multifunc"]],permutations = reps, method = "sp")
}
for(i in names(data_SEM.list)[1:2]){
  for(j in names(data_SEM.list)[2:3]){
    if(i != j){
      print(paste0(i,".vs.",j))
      mantel.re[[paste0(i,".vs.",j)]] = vegan::mantel(xdis = data_SEM.list[[i]], ydis = data_SEM.list[[j]],permutations = reps,  method = "sp")
    }
  }
}
cbind(r = sapply(mantel.re,function(X){X$statistic}),
      p=sapply(mantel.re,function(X){X$signif}))
############# other codes not that dose not relate to the results in the current manuscript ################
########################### calculating CVs within group ###################
data_forcv=alpha_all[,uniq_functions]
group=paste(map$D,map$I,map$PR,sep=".")
Fun = function(x){if(sum(!is.na(x))>2){sd(x,na.rm=T)/mean(x+10^-16,na.rm=T)}else{NA}}
cvs=apply(data_forcv,2,FUN=function(x){group_calculate(x,Fun=Fun,group=group)})
plot_data=apply(map[c("D","I","PR")],2,FUN=function(x){group_calculate(x, Fun=function(x){x[1]},group=group)})
cvs=as.data.frame(cvs)
plot_data=as.data.frame(plot_data);
p=NULL
sem_cv=NULL
for(ii in names(cvs)){

  plot_data$x=cvs[,ii]
  p[[ii]]=ggplot(plot_data,aes(PR,x,color=as.factor(D),shape=as.factor(I)))+
    #facet_wrap(~paste(D,I))+
    geom_point(size=3)+
    scale_shape_manual(values = c(1,2))+
    scale_color_manual(values = c("green","blue","black"))+
    #geom_smooth(group)+
    labs(x=element_blank(), y = ii)+
    #scale_y_continuous(breaks = sort(c(-seq(bk.y,5,bk.y),seq(0,5,bk.y))))+
    theme_bw()+
    theme(text=element_text(size = 8, family = "sans", face = "plain"),axis.title.x = element_blank(),
          panel.grid.major.x = element_blank(),    panel.grid.minor.x =  element_blank(),  panel.grid.minor.y =  element_blank(),
          panel.background = element_rect(color = 'black', fill = 'transparent'),legend.position = "none")
  plot_data$x=stand(plot_data$x)
  sem_cv=rbind(sem_cv,plot_data)
}
pp <- plot_grid(p[[1]],p[[2]],p[[3]],p[[4]],
                p[[5]],p[[6]],p[[7]],p[[8]],
                p[[9]],p[[10]],p[[11]],p[[12]],
                p[[13]],p[[14]],
                hjust = 0,vjust = 0,
                label_fontfamily = "sans", label_size = 12,  label_x =c(0.15,0.15), label_y =1, ncol = 4, nrow=4)
#+ theme(plot.margin = unit(c(0.8,0,0,0), "cm"))
pp

ggplot(sem_cv,aes(D,x))+
  facet_wrap(~PR,nrow = 1)+
  geom_point(aes(x=D+0.35*(I-0.5),color=as.factor(D),shape=as.factor(I)),size=2)+
  scale_shape_manual(values = c(1,2))+
  scale_color_manual(values = c("black","#0092FF", "#FF0000"))+
  geom_smooth(method = "lm",se=T,linetype="solid",color="black")+
  stat_poly_eq(aes(label=paste(stat(rr.label), stat(p.value.label), sep = "*\", \"*")),
               formula = y ~ x,size=3.5, parse = T)+
  labs(x=element_blank(), y = ii)+
  #scale_y_continuous(breaks = sort(c(-seq(bk.y,5,bk.y),seq(0,5,bk.y))))+
  theme_bw()+
  theme(text=element_text(size = 8, family = "sans", face = "plain"),axis.title.x = element_blank(),
        panel.grid.major.x = element_blank(),    panel.grid.minor.x =  element_blank(),  panel.grid.minor.y =  element_blank(),
        panel.background = element_rect(color = 'black', fill = 'transparent'),legend.position = "none")

names(sem_cv)
library(lavaan)
library(semPlot)
fit=lavaan(add.variance('x~D+I+PR'),data = sem_cv,bootstrap=F)
summary(fit)
fitMeasures(fit,c("chisq","df","pvalue", "cfi", "rmsea","rmsea.pvalue"))

semPaths(fit,what = "est",whatLabels = "std.all",style="ram",residuals=T,nCharNodes=6,
         sizeMan=10,sizeMan2 = 5,sizeLat=20,sizeInt=10)

semPaths(fit,whatLabels = "std.all")

glm(x~D+I+PR,data = sem_cv)

###########################calculating CVs between Drought group ########################### 
library(matrixStats)
data_forcv=alpha_all[,uniq_functions]
group=data.frame(within=paste(map$D,map$I,map$PR,sep="."),bet.D=paste(map$I,map$PR,sep="."),bet.I=paste(map$D,map$PR,sep="."))

for(i in unique(group$bet.D)){
  data_forcv.bet.D=data_forcv[group$bet.D == i,]
}
data_forcv.bet.D=group_mean2(data_forcv$Pl.mass,group = group$within)
Fun = function(x){if(sum(!is.na(x))>2){sd(x,na.rm=T)/mean(x+10^-16,na.rm=T)}else{NA}}
cvs=apply(data_forcv,2,FUN=function(x){group_calculate(x,Fun=Fun,group=group$bet.I)})
plot_data=apply(map[c("D","I","PR")],2,FUN=function(x){group_calculate(x, Fun=function(x){x[1]},group=group$bet.I)})
cvs=as.data.frame(cvs)
plot_data=as.data.frame(plot_data);
p=NULL
sem_cv=NULL
for(ii in names(cvs)){
  
  plot_data$x=cvs[,ii]
  p[[ii]]=ggplot(plot_data,aes(PR,x,color=as.factor(D),shape=as.factor(I)))+
    #facet_wrap(~paste(D,I))+
    geom_point(size=3)+
    scale_shape_manual(values = c(1,2))+
    scale_color_manual(values = c("green","blue","black"))+
    #geom_smooth(group)+
    labs(x=element_blank(), y = ii)+
    #scale_y_continuous(breaks = sort(c(-seq(bk.y,5,bk.y),seq(0,5,bk.y))))+
    theme_bw()+
    theme(text=element_text(size = 8, family = "sans", face = "plain"),axis.title.x = element_blank(),
          panel.grid.major.x = element_blank(),    panel.grid.minor.x =  element_blank(),  panel.grid.minor.y =  element_blank(),
          panel.background = element_rect(color = 'black', fill = 'transparent'),legend.position = "none")
  plot_data$x=stand(plot_data$x)
  sem_cv=rbind(sem_cv,plot_data)
}
pp <- plot_grid(p[[1]],p[[2]],p[[3]],p[[4]],
                p[[5]],p[[6]],p[[7]],p[[8]],
                p[[9]],p[[10]],p[[11]],p[[12]],
                p[[13]],p[[14]],
                hjust = 0,vjust = 0,
                label_fontfamily = "sans", label_size = 12,  label_x =c(0.15,0.15), label_y =1, ncol = 4, nrow=4)
#+ theme(plot.margin = unit(c(0.8,0,0,0), "cm"))
pp

ggplot(sem_cv,aes(D,x))+
  facet_wrap(~PR,nrow = 1)+
  geom_point(aes(x=D+0.35*(I-0.5),color=as.factor(D),shape=as.factor(I)),size=2)+
  scale_shape_manual(values = c(1,2))+
  scale_color_manual(values = c("black","#0092FF", "#FF0000"))+
  geom_smooth(method = "lm",se=T,linetype="solid",color="black")+
  stat_poly_eq(aes(label=paste(stat(rr.label), stat(p.value.label), sep = "*\", \"*")),
               formula = y ~ x,size=3.5, parse = T)+
  labs(x=element_blank(), y = ii)+
  #scale_y_continuous(breaks = sort(c(-seq(bk.y,5,bk.y),seq(0,5,bk.y))))+
  theme_bw()+
  theme(text=element_text(size = 8, family = "sans", face = "plain"),axis.title.x = element_blank(),
        panel.grid.major.x = element_blank(),    panel.grid.minor.x =  element_blank(),  panel.grid.minor.y =  element_blank(),
        panel.background = element_rect(color = 'black', fill = 'transparent'),legend.position = "none")

names(sem_cv)
fit=lavaan(add.variance('x~D+PR'),data = sem_cv,bootstrap=TRUE)
summary(fit)
fitMeasures(fit,c("chisq","df","pvalue ", "cfi", "rmsea","rmsea.pvalue"))

semPaths(fit,what = "est",whatLabels = "std.all",style="ram",residuals=T,nCharNodes=6,
         sizeMan=10,sizeMan2 = 5,sizeLat=20,sizeInt=10)
semPaths(fit,whatLabels = "std.all")

gl=glm(x~D+I+PR,data = sem_cv)
anova(gl)

########################################## anova-univariables ##############
D=stand(map$D)
I=stand(map$I)
PR=stand(map$PR)
aov.RE=NULL
for(i in names(alpha_all)){
  md=lm(stand(alpha_all[[i]])~I*PR*D)
  aov.re=as.data.frame(anova(md))
  aov.re[]=round(aov.re[],digits = 4)
  aov.re$sig.sym=as.symbol.pvalue(aov.re$`Pr(>F)`)
  aov.re$RR=round(aov.re$`Sum Sq`/sum(aov.re$`Sum Sq`),4)
  aov.re$Coefficients=round(c(md$coefficients[rownames(aov.re)[rownames(aov.re)!="Residuals"]],md$coefficients[1]),3)
  aov.re$dependent=i;aov.re$factor=rownames(aov.re);
  aov.RE=rbind(aov.RE,aov.re)
}
aov.RE
#write.csv(aov.RE,"SEM/H1_4/aov.univariables.csv")

