# Please run step1--stability_calcuation.R
#setwd("E:/Fauna/Communications Biology/data")
library(vegan)
library(GUniFrac)
library(picante)
library(phangorn)
library(ggpmisc)
library(cowplot)
source("E:/Fauna/Communications Biology/data/R_sources/other_source.R")

taxa=read.csv("data_files/fauna_tax.csv",row.names = 1)
library(vegan)
order_fauna=F
same_PS=F
SQRT_biomass=F
std.Method="no"



abs.PR=0.5
map=read.csv("data_files/map.csv",row.names = 1)
nrow(map)
#map=map[colnames(tabs[["Pl"]]),]
#tabs$Pl=tabs$LF
#trees$Pl=trees$LF
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
  a=c(a,all(apply(tabs[[i]],2,sum) >0 ))
}
a

dim(tabs$Pl)
dim(tabs$SF)
dim(tabs$LF)

order_dist = NULL
for(i in c("LF","SF")){
  order_dist[[i]]=as.matrix(vegan::vegdist(t( rowsum(tabs[[i]],group = taxa[rownames(tabs[[i]]),]$Order) ),method = "bray"))
}
dim(order_dist[[i]])

i="SF"
tabulate.char(taxa[rownames(tabs[[i]]),]$Order,weight = rowSums(tabs[[i]]))


############################### PL #######################
plant_taxon = read.csv("data_files/plant_taxon.csv",row.names = 1,header = T)
tabi = "Pl"
TAB=tabs[[tabi]]

V=c(1:16,
    apply(combn(16,2),2,paste,collapse = "-"),
    apply(combn(16,3),2,paste,collapse = "-"),
    apply(combn(16,4),2,paste,collapse = "-"))
inv = resD = resI = Sim = dtf(NA,length(V),5,dimnames = list(V,c("r0","df","r1","part.r","part.r.p")))
nn=0
for (delete_subset in V){
  nn=nn + 1
  print(c(nn,delete_subset))
  SUBTAB=TAB
  #dim(TAB)
  rownames(SUBTAB)
  SUBTAB = SUBTAB[!rownames(SUBTAB) %in% rownames(SUBTAB)[unlist(strsplit(delete_subset,"-"))%>%as.numeric],]
  SUBTAB = SUBTAB[,colSums(SUBTAB)>0]
  #dim(SUBTAB)
  df = ncol(SUBTAB)
  mat_temp=matrix(0,nrow(map),nrow(map),dimnames = list(rownames(map),rownames(map)))
  
  unifracs=UniFrac(if(SQRT_biomass){sqrt(SUBTAB)}else{SUBTAB},trees[[tabi]],total.normalize = F)
  simi=1-as.matrix(unifracs$w_unifrac[,,1])
  mat_temp[]=0;mat_temp[rownames(simi),colnames(simi)]=simi[]
  mat_temp[!rownames(mat_temp)%in%colnames(SUBTAB),!rownames(mat_temp)%in%colnames(SUBTAB)]=NA
  
  j1 = !is.na(cov1 <-  mat_temp[invariability$IND])
  j2 = !is.na(cov2 <-  mat_temp[resist_to.I$IND])
  j3 = !is.na(cov3 <-  mat_temp[resist_to.D$IND])
  j4 = !is.na(cov4 <-  mat_temp[IND[judge.low.up]])
  
  r = cor(cov1,invariability$Multifunc,use = "pairwise.complete.obs")
  r2 = cor(invariability[[tabi]][j1],invariability$Multifunc[j1],use = "pairwise.complete.obs")
  inv[delete_subset,] = c(r=r,df = df, r2=r2, part.r = partial.corr(x = invariability[[tabi]], y = invariability$Multifunc, covariate = cov1)$par.cor)
  
  r = cor(cov2,resist_to.I$Multifunc,use = "pairwise.complete.obs")
  r2 = cor(resist_to.I[[tabi]][j2],resist_to.I$Multifunc[j2],use = "pairwise.complete.obs")
  resI[delete_subset,] = c(r=r,df = df, r2=r2, part.r = partial.corr(x = resist_to.I[[tabi]], y = resist_to.I$Multifunc, covariate = cov2)$par.cor)
  
  r = cor(cov3,resist_to.D$Multifunc,use = "pairwise.complete.obs")
  r2 = cor(resist_to.D[[tabi]][j3],resist_to.D$Multifunc[j3],use = "pairwise.complete.obs")
  resD[delete_subset,] = c(r=r,df = df, r2=r2, part.r = partial.corr(x = resist_to.D[[tabi]], y = resist_to.D$Multifunc, covariate = cov3)$par.cor)
  
  r = cor(cov4,stab$Multifunc,use = "pairwise.complete.obs")
  r2 = cor(stab[[tabi]][j4],stab$Multifunc[j4],use = "pairwise.complete.obs")
  Sim[delete_subset,] = c(r=r,df = df, r2=r2, part.r = partial.corr(x = stab[[tabi]], y = stab$Multifunc, covariate = cov4)$par.cor)
}
r = c(cor(invariability$Pl,invariability$Multifunc,use = "pairwise.complete.obs"),
      cor(resist_to.I$Pl,resist_to.I$Multifunc,use = "pairwise.complete.obs"),
      cor(resist_to.D$Pl,resist_to.D$Multifunc,use = "pairwise.complete.obs"),
      cor(stab$Pl,stab$Multifunc,use = "pairwise.complete.obs") )
PL_pair = list(inv=inv, resD=resD, resI=resI, Sim=Sim, id = V)
#save(PL_pair,file = "PL_pair.RData")


data1=NULL
for(i in 1:3){
  data1 = rbind(data1,cbind(PL_pair[[i]],PL_pair$id,names(PL_pair)[i]))
}
names(data1)[6:7] = c("id","aspect")

data1$n=apply(data1$id%>%strsplit2mat(split = "-"),1,function(x){length(unique(x[x!=""]))})
data1$sp = data1$id

all(rownames(data1) == data1$sp)
data1$dif = as.numeric(data1$r0) - as.numeric(data1$r1)
data1$sort = 1:nrow(data1)

for(i in names(data1)){
  if(all(!is.na(as.numeric(data1[[i]])))) data1[[i]] = as.numeric(data1[[i]])
}

for(i in 1:4){
  for(j in unique(data1$aspect))
  data1[data1$n == i & data1$aspect == j,] = sortby(data1[data1$n == i & data1$aspect == j,],"dif")
}


library(ggplot2)
p.PL <- ggplot(data1,aes(x = as.numeric(n), y = dif,color=as.factor(n)))+
  facet_wrap(~aspect)+
  geom_jitter(width = 0.2,height = 0,size= 0.1,shape = 1)+
  geom_boxplot(aes(group=as.factor(n)),fill=NA,color = "black",size = 0.2)+
  labs(y="sim change", x = "species deletion")+
  #scale_x_continuous(breaks = sort(c(-seq(bk.y,5,bk.y),seq(0,5,bk.y))))+
  theme_bw()+
  theme(text=element_text(size = 7, family = "sans", face = "plain"),
        panel.grid.major.y = element_blank(),    panel.grid.minor.y =  element_blank(),  panel.grid.minor.x =  element_blank(),
        panel.background = element_rect(color = 'black', fill = 'transparent'),legend.position = "none")
p.PL

#ggsave(filename = "Pl_stab_change_after_deletion.pdf",p,width = 160,height = 50,units = "mm")
  
realtive_importance = matrix(NA,length(rownames(TAB)), 4, dimnames = list(rownames(TAB),1:4))
#colnames(realtive_importance) = uni_class
xx = strsplit2mat(data1$id,split = "-") %>% as.data.frame()
for(i in 1:ncol(xx)) xx[[i]] %<>% as.numeric()

for(i in 1:4){
  j1 = data1$n == i
  for(ii in 1:length(rownames(TAB))){
    temp = data1[apply(xx==ii,1,sum,na.rm=T) > 0 & j1,]
    realtive_importance[ii,i] = mean(-temp$dif/temp$r1)
  }
}
realtive_importance%<>%as.data.frame()
realtive_importance$abund = apply(TAB,1,function(x){mean(x[x>0])})
realtive_importance$sd = apply(TAB,1,function(x){sd(x[x>0])})
realtive_importance[,c("abund","sd")] = realtive_importance[,c("abund","sd")]/sum(realtive_importance$abund)
sum(realtive_importance$abund)
realtive_importance$cv = realtive_importance$sd / realtive_importance$abund
realtive_importance$sp = rownames(realtive_importance)
realtive_importance_pl = realtive_importance
plot_df = NULL
for(i in 1:4){
  plot_df = rbind(plot_df,cbind(value = realtive_importance[,i],realtive_importance[,c(5:8)], delition = i))
}
ggplot(plot_df,aes(abund, value, color = as.factor(delition)))+
  geom_point()
ggplot(plot_df,aes(sd, value, color = as.factor(delition)))+
  geom_point()
ggplot(plot_df,aes(cv, value, color = as.factor(delition)))+
  geom_point()

plot_df$sp = plant_taxon[plot_df$sp,"Species"]
realtive_importance$sp=plant_taxon[realtive_importance$sp,"Species"]

sp1 <- realtive_importance$sp[1]
realtive_importance$sp[realtive_importance$sp == sp1] = paste0(rep("a",20),collapse = "")
plot_df$sp[plot_df$sp == sp1] = paste0(rep("a",20),collapse = "")

plot_df$sp = factor(plot_df$sp,levels = realtive_importance$sp[order(realtive_importance$abund,decreasing = T)])
realtive_importance$sp = factor(realtive_importance$sp,levels = realtive_importance$sp[order(realtive_importance$abund,decreasing = T)])



p.pl_importantSP = ggplot(realtive_importance,aes(abund, sp))+
  #geom_point()+
  geom_bar(aes(abund,sp),stat='identity', position='stack',fill="lightgray",color=NA,size=0.25,width = 0.75) +
  geom_errorbar(aes(xmin = abund -sd, xmax = abund +sd),width = 0.25,size=0.375) +
  geom_point(data = plot_df,aes(value,sp, color = as.factor(delition)),shape = 1)+
  theme_bw() +
  theme(text=element_text(size = 8, family = "sans", face = "plain"),axis.title.x = element_blank(),
        panel.grid.major.x = element_blank(),    panel.grid.minor.x =  element_blank(),  panel.grid.minor.y =  element_blank(),
        panel.background = element_rect(color = 'black', fill = 'transparent'),legend.position = "none")

data1$realtive_importance = -data1$dif/data1$r1
data1$abund = data1$sd = data1$cv = 0#sapply(X=1:nrow(data1),FUN = function(X,y=strsplit(data1$sp,split = "-"),...){sum(realtive_importance[as.numeric(y[[X]]),"abund"])})
y=strsplit(data1$sp,split = "-")
for(i in 1:nrow(data1)){
  data1$abund[i] = sum(realtive_importance[as.numeric(y[[i]]),"abund"])
  data1$sd[i] = sqrt(sum(realtive_importance[as.numeric(y[[i]]),"sd"]^2))
  data1$cv[i] = mean(realtive_importance[as.numeric(y[[i]]),"cv"])
}

p1.1=ggplot(data1[order(data1$n,decreasing = T),],aes(abund, realtive_importance, color = as.factor(n)))+
  facet_wrap(~aspect)+
  geom_point(aes(size = 5-n), size = 0.2, shape = 1, alpha = 1)+
  geom_smooth(method = "lm", se=F,size = 1) +
  theme(text=element_text(size = 7, family = "sans", face = "plain"),
        panel.grid.major.y = element_blank(),    panel.grid.minor.y =  element_blank(),  panel.grid.minor.x =  element_blank(),
        panel.background = element_rect(color = 'black', fill = 'transparent'),legend.position = "none")

p1.2=ggplot(data1[order(data1$n,decreasing = T),],aes(sd, realtive_importance, color = as.factor(n)))+
  facet_wrap(~aspect)+
  geom_point(aes(size = 5-n), size = 0.2, shape = 1, alpha = 1)+
  geom_smooth(method = "lm", se=F,size = 1) +
  theme(text=element_text(size = 7, family = "sans", face = "plain"),
        panel.grid.major.y = element_blank(),    panel.grid.minor.y =  element_blank(),  panel.grid.minor.x =  element_blank(),
        panel.background = element_rect(color = 'black', fill = 'transparent'),legend.position = "none")

p1.3=ggplot(data1[order(data1$n,decreasing = T),],aes(cv, realtive_importance, color = as.factor(n)))+
  facet_wrap(~aspect)+
  geom_point(aes(size = 5-n), size = 0.2, shape = 1, alpha = 1)+
  geom_smooth(method = "lm", se=F,size = 1) +
  theme(text=element_text(size = 7, family = "sans", face = "plain"),
        panel.grid.major.y = element_blank(),    panel.grid.minor.y =  element_blank(),  panel.grid.minor.x =  element_blank(),
        panel.background = element_rect(color = 'black', fill = 'transparent'),legend.position = "none")
data = data1[data1$n ==1,]
data$sp = plant_taxon[rownames(tabs[["Pl"]])[as.numeric(data$sp)],"Species"]
data$sp[data$sp == "Cichorium intybus"] = rep("a",20)%>%paste0(collapse = "")
data$sp = factor(data$sp,levels(realtive_importance$sp))

p1.4=ggplot(realtive_importance,aes(abund, sp))+
  #geom_point()+
  geom_errorbar(aes(xmin = abund/3, xmax = abund +sd),width = 0.25,size=0.375)+
  geom_bar(aes(abund,sp),stat='identity', position='stack',fill="lightgray",color=NA,size=0.25,width = 0.75) +
  #facet_wrap(~aspect) +
  geom_point(data=data,aes(realtive_importance,sp, color = as.factor(aspect)), shape =1, alpha = 1) +
  #scale_color_manual(values = c("black","#DC143C", "#CD853F"))
  #geom_smooth(method = "lm", se=F,size = 1) +
  theme(text=element_text(size = 7, family = "sans", face = "plain"),
        panel.grid.major.y = element_blank(),    panel.grid.minor.y =  element_blank(),  panel.grid.minor.x =  element_blank(),
        panel.background = element_rect(color = 'black', fill = 'transparent'),legend.position = "none")

p1.1
p1.2
p1.3
p1 = pp <- plot_grid(p1.1,p1.2,p1.3,
                     hjust = 0,vjust = 0,
                     label_fontfamily = "sans", label_size = 7,  label_x =c(0.15,0.15), label_y =1, ncol = 1, nrow=3)
p1
#+ theme(plot.margin = unit(c(0.8,0,0,0), "cm"))

#ggsave(filename = "Pl_p1.pdf",p1,width = 160,height = 150,units = "mm")
############################### LF #######################
tabi = "LF"
TAB=tabs[[tabi]]
tabulate.char(taxa[rownames(tabs[[tabi]]),]$Order,weight = rowSums(tabs[[tabi]]))

uni_class = unique(taxon <- taxa[rownames(TAB),"Order"])

V=c(1:length(uni_class),
    apply(combn(length(uni_class),2),2,paste,collapse = "-"),
    apply(combn(length(uni_class),3),2,paste,collapse = "-"),
    apply(combn(length(uni_class),4),2,paste,collapse = "-"))
inv = resD = resI = Sim = dtf(NA,length(V),5,dimnames = list(V,c("r0","df","r1","part.r","part.r.p")))
nn=0
for (delete_subset in V){
  nn=nn+1
  print(c(nn,delete_subset))
  SUBTAB=TAB
  dim(TAB)
  #rownames(SUBTAB)
  SUBTAB = SUBTAB[!taxon %in% uni_class[unlist(strsplit(delete_subset,"-"))%>%as.numeric],]
  SUBTAB = SUBTAB[,colSums(SUBTAB)>0]
  #dim(SUBTAB)
  df = ncol(SUBTAB)
  mat_temp=matrix(0,nrow(map),nrow(map),dimnames = list(rownames(map),rownames(map)))
  
  unifracs=UniFrac(if(SQRT_biomass){sqrt(SUBTAB)}else{SUBTAB},trees[[tabi]],total.normalize = F)
  simi=1-as.matrix(unifracs$w_unifrac[,,1])
  mat_temp[]=0;mat_temp[rownames(simi),colnames(simi)]=simi[]
  mat_temp[!rownames(mat_temp)%in%colnames(SUBTAB),!rownames(mat_temp)%in%colnames(SUBTAB)]=NA
  
  j1 = !is.na(cov1 <-  mat_temp[invariability$IND])
  j2 = !is.na(cov2 <-  mat_temp[resist_to.I$IND])
  j3 = !is.na(cov3 <-  mat_temp[resist_to.D$IND])
  j4 = !is.na(cov4 <-  mat_temp[IND[judge.low.up]])
  
  r = cor(cov1,invariability$Multifunc,use = "pairwise.complete.obs")
  r2 = cor(invariability[[tabi]][j1],invariability$Multifunc[j1],use = "pairwise.complete.obs")
  inv[delete_subset,] = c(r=r,df = df, r2=r2, part.r = partial.corr(x = invariability[[tabi]], y = invariability$Multifunc, covariate = cov1)$par.cor)
  
  r = cor(cov2,resist_to.I$Multifunc,use = "pairwise.complete.obs")
  r2 = cor(resist_to.I[[tabi]][j2],resist_to.I$Multifunc[j2],use = "pairwise.complete.obs")
  resI[delete_subset,] = c(r=r,df = df, r2=r2, part.r = partial.corr(x = resist_to.I[[tabi]], y = resist_to.I$Multifunc, covariate = cov2)$par.cor)
  
  r = cor(cov3,resist_to.D$Multifunc,use = "pairwise.complete.obs")
  r2 = cor(resist_to.D[[tabi]][j3],resist_to.D$Multifunc[j3],use = "pairwise.complete.obs")
  resD[delete_subset,] = c(r=r,df = df, r2=r2, part.r = partial.corr(x = resist_to.D[[tabi]], y = resist_to.D$Multifunc, covariate = cov3)$par.cor)
  
  r = cor(cov4,stab$Multifunc,use = "pairwise.complete.obs")
  r2 = cor(stab[[tabi]][j4],stab$Multifunc[j4],use = "pairwise.complete.obs")
  Sim[delete_subset,] = c(r=r,df = df, r2=r2, part.r = partial.corr(x = stab[[tabi]], y = stab$Multifunc, covariate = cov4)$par.cor)
  
}

#LF_pair=list(inv=inv,resD=resD,resI=resI,Sim=Sim,id = V)
#save(LF_pair,file = "LF_pair.RData")


data2=NULL
for(i in 1:3){
  data2 = rbind(data2,cbind(LF_pair[[i]],LF_pair$id,names(LF_pair)[i]))
}
names(data2)[6:7] = c("id","aspect")

data2$n=apply(data2$id%>%strsplit2mat(split = "-"),1,function(x){length(unique(x[x!=""]))})
data2$sp = data2$id

all(rownames(data2) == data2$sp)
data2$dif = as.numeric(data2$r0) - as.numeric(data2$r1)
data2$sort = 1:nrow(data2)
for(i in names(data2)){
  if(all(!is.na(as.numeric(data2[[i]])))) data2[[i]] = as.numeric(data2[[i]])
}

for(i in 1:4){
  for(j in unique(data2$aspect))
    data2[data2$n == i & data2$aspect == j,] = sortby(data2[data2$n == i & data2$aspect == j,],"dif")
}


library(ggplot2)
p.LF <- ggplot(data2,aes(x = as.numeric(n), y = dif,color=as.factor(n)))+
  facet_wrap(~aspect)+
  geom_jitter(width = 0.2,height = 0,size= 0.1,shape = 1)+
  geom_boxplot(aes(group=as.factor(n)),fill=NA,color = "black",size = 0.2)+
  labs(y="sim change", x = "species deletion")+
  #scale_x_continuous(breaks = sort(c(-seq(bk.y,5,bk.y),seq(0,5,bk.y))))+
  theme_bw()+
  theme(text=element_text(size = 7, family = "sans", face = "plain"),
        panel.grid.major.y = element_blank(),    panel.grid.minor.y =  element_blank(),  panel.grid.minor.x =  element_blank(),
        panel.background = element_rect(color = 'black', fill = 'transparent'),legend.position = "none")
p.LF
#'+  geom_smooth(method = "loess")'
#ggsave(filename = "LF_stab_change_after_deletion.pdf",p,width = 160,height = 50,units = "mm")

realtive_importance = matrix(NA,length(uni_class),4,dimnames = list(uni_class,1:4))
#colnames(realtive_importance) = uni_class
xx = strsplit2mat(data2$id,split = "-") %>% as.data.frame()
for(i in 1:ncol(xx)) xx[[i]] %<>% as.numeric()

for(i in 1:4){
  j1 = data2$n == i
  for(ii in 1:length(uni_class)){
    temp = data2[apply(xx==ii,1,sum,na.rm=T) > 0 & j1,]
    realtive_importance[ii,i] = mean(-temp$dif/temp$r1)
  }
}
order_tab = rowsum(TAB,group = taxa[rownames(TAB),]$Order)[uni_class,]
realtive_importance%<>%as.data.frame()
realtive_importance$abund = apply(order_tab,1,function(x){mean(x[x>0])})
realtive_importance$sd = apply(order_tab,1,function(x){sd(x[x>0])})
realtive_importance[,c("abund","sd")] = realtive_importance[,c("abund","sd")]/sum(realtive_importance$abund)
sum(realtive_importance$abund)
realtive_importance$cv = realtive_importance$sd / realtive_importance$abund
realtive_importance$sp = rownames(realtive_importance)
rownames(realtive_importance) == uni_class

plot_df = NULL
for(i in 1:4){
  plot_df = rbind(plot_df,cbind(value = realtive_importance[,i],realtive_importance[,c(5:8)], delition = i))
}
ggplot(plot_df,aes(abund, value, color = as.factor(delition)))+
  geom_point()
ggplot(plot_df,aes(sd, value, color = as.factor(delition)))+
  geom_point()
ggplot(plot_df,aes(cv, value, color = as.factor(delition)))+
  geom_point()


realtive_importance$sp[realtive_importance$sp == "Poduromorpha"] = paste0(rep("a",20),collapse = "")
plot_df$sp[plot_df$sp == "Poduromorpha"] = paste0(rep("a",20),collapse = "")

plot_df$sp = factor(plot_df$sp,levels = realtive_importance$sp[order(realtive_importance$abund,decreasing = T)])
realtive_importance$sp = factor(realtive_importance$sp,levels = realtive_importance$sp[order(realtive_importance$abund,decreasing = T)])
realtive_importance = realtive_importance[realtive_importance$sp %in% levels(realtive_importance$sp)[1:10],]
plot_df = plot_df[plot_df$sp %in% levels(realtive_importance$sp)[1:10],]

p.LF_importantSP = ggplot(realtive_importance,aes(abund, sp))+
  #geom_point()+
  geom_bar(aes(abund,sp),stat='identity', position='stack',fill="lightgray",color=NA,size=0.25,width = 0.75) +
  geom_errorbar(aes(xmin = abund -sd, xmax = abund +sd),width = 0.25,size=0.375) +
  geom_point(data = plot_df,aes(value,sp, color = as.factor(delition)),shape = 1)+
  theme_bw() +
  theme(text=element_text(size = 8, family = "sans", face = "plain"),axis.title.x = element_blank(),
        panel.grid.major.x = element_blank(),    panel.grid.minor.x =  element_blank(),  panel.grid.minor.y =  element_blank(),
        panel.background = element_rect(color = 'black', fill = 'transparent'),legend.position = "none")

data2$realtive_importance = -data2$dif/data2$r1
data2$abund = data2$sd = data2$cv = 0#sapply(X=1:nrow(data2),FUN = function(X,y=strsplit(data2$sp,split = "-"),...){sum(realtive_importance[as.numeric(y[[X]]),"abund"])})
y=strsplit(data2$sp,split = "-")
for(i in 1:nrow(data2)){
  data2$abund[i] = sum(realtive_importance[as.numeric(y[[i]]),"abund"])
  data2$sd[i] = sqrt(sum(realtive_importance[as.numeric(y[[i]]),"sd"]^2))
  data2$cv[i] = mean(realtive_importance[as.numeric(y[[i]]),"cv"])
}

p2.1 = ggplot(data2[order(data2$n,decreasing = T),],aes(abund, realtive_importance, color = as.factor(n)))+
  facet_wrap(~aspect)+
  geom_point(aes(size = 5-n), size = 0.2, shape = 1, alpha = 1)+
  geom_smooth(method = "lm", se=F,size = 1) +
  theme(text=element_text(size = 7, family = "sans", face = "plain"),
        panel.grid.major.y = element_blank(),    panel.grid.minor.y =  element_blank(),  panel.grid.minor.x =  element_blank(),
        panel.background = element_rect(color = 'black', fill = 'transparent'),legend.position = "none")

p2.2 = ggplot(data2[order(data2$n,decreasing = T),],aes(sd, realtive_importance, color = as.factor(n)))+
  facet_wrap(~aspect)+
  geom_point(aes(size = 5-n), size = 0.2, shape = 1, alpha = 1)+
  geom_smooth(method = "lm", se=F,size = 1) +
  theme(text=element_text(size = 7, family = "sans", face = "plain"),
        panel.grid.major.y = element_blank(),    panel.grid.minor.y =  element_blank(),  panel.grid.minor.x =  element_blank(),
        panel.background = element_rect(color = 'black', fill = 'transparent'),legend.position = "none")

p2.3 = ggplot(data2[order(data2$n,decreasing = T),],aes(cv, realtive_importance, color = as.factor(n)))+
  facet_wrap(~aspect)+
  geom_point(aes(size = 5-n), size = 0.2, shape = 1, alpha = 1)+
  geom_smooth(method = "lm", se=F,size = 1) +
  theme(text=element_text(size = 7, family = "sans", face = "plain"),
        panel.grid.major.y = element_blank(),    panel.grid.minor.y =  element_blank(),  panel.grid.minor.x =  element_blank(),
        panel.background = element_rect(color = 'black', fill = 'transparent'),legend.position = "none")

data = data2[data2$n ==1,]
data$sp = uni_class[as.numeric(data$sp)]
data$sp[data$sp == "Poduromorpha"] = rep("a",20)%>%paste0(collapse = "")
data$sp = factor(data$sp,levels(realtive_importance$sp))
data = data[data$sp %in% levels(realtive_importance$sp)[1:10],]

p2.4=ggplot(realtive_importance,aes(abund, sp))+
  #geom_point()+
  geom_errorbar(aes(xmin = abund/3, xmax = abund +sd),width = 0.25,size=0.375)+
  geom_bar(aes(abund,sp),stat='identity', position='stack',fill="lightgray",color=NA,size=0.25,width = 0.75) +
  #facet_wrap(~aspect) +
  geom_point(data=data,aes(realtive_importance,sp, color = as.factor(aspect)), shape =1, alpha = 1) +
  #scale_color_manual(values = c("black","#DC143C", "#CD853F"))
  #geom_smooth(method = "lm", se=F,size = 1) +
  theme(text=element_text(size = 7, family = "sans", face = "plain"),
        panel.grid.major.y = element_blank(),    panel.grid.minor.y =  element_blank(),  panel.grid.minor.x =  element_blank(),
        panel.background = element_rect(color = 'black', fill = 'transparent'),legend.position = "none")

p2 = pp <- plot_grid(p2.1,p2.2,p2.3,
                     hjust = 0,vjust = 0,
                     label_fontfamily = "sans", label_size = 7,  label_x =c(0.15,0.15), label_y =1, ncol = 1, nrow=3)
p2
#+ theme(plot.margin = unit(c(0.8,0,0,0), "cm"))

#ggsave(filename = "LF_p2.pdf",p2,width = 160,height = 150,units = "mm")
############################### SF #######################
tabi = "SF"
TAB=tabs[[tabi]]
tabulate.char(taxa[rownames(tabs[[tabi]]),]$Order,weight = rowSums(tabs[[tabi]]))

uni_class = unique(taxon <- taxa[rownames(TAB),"Order"])
uni_class = uni_class[!uni_class == "--"]
V=c(1:length(uni_class),
    apply(combn(length(uni_class),2),2,paste,collapse = "-"),
    apply(combn(length(uni_class),3),2,paste,collapse = "-"),
    apply(combn(length(uni_class),4),2,paste,collapse = "-"))
inv = resD = resI = Sim = dtf(NA,length(V),5,dimnames = list(V,c("r0","df","r1","part.r","part.r.p")))
nn=0
for (delete_subset in V){
  nn=nn+1
  print(c(nn,delete_subset))
  SUBTAB=TAB
  dim(TAB)
  #rownames(SUBTAB)
  SUBTAB = SUBTAB[!taxon %in% uni_class[unlist(strsplit(delete_subset,"-"))%>%as.numeric],]
  SUBTAB = SUBTAB[,colSums(SUBTAB)>0]
  #dim(SUBTAB)
  df = ncol(SUBTAB)
  mat_temp=matrix(0,nrow(map),nrow(map),dimnames = list(rownames(map),rownames(map)))
  
  unifracs=UniFrac(if(SQRT_biomass){sqrt(SUBTAB)}else{SUBTAB},trees[[tabi]],total.normalize = F)
  simi=1-as.matrix(unifracs$w_unifrac[,,1])
  mat_temp[]=0;mat_temp[rownames(simi),colnames(simi)]=simi[]
  mat_temp[!rownames(mat_temp)%in%colnames(SUBTAB),!rownames(mat_temp)%in%colnames(SUBTAB)]=NA
  
  j1 = !is.na(cov1 <-  mat_temp[invariability$IND])
  j2 = !is.na(cov2 <-  mat_temp[resist_to.I$IND])
  j3 = !is.na(cov3 <-  mat_temp[resist_to.D$IND])
  j4 = !is.na(cov4 <-  mat_temp[IND[judge.low.up]])
  
  r = cor(cov1,invariability$Multifunc,use = "pairwise.complete.obs")
  r2 = cor(invariability[[tabi]][j1],invariability$Multifunc[j1],use = "pairwise.complete.obs")
  inv[delete_subset,] = c(r=r,df = df, r2=r2, part.r = partial.corr(x = invariability[[tabi]], y = invariability$Multifunc, covariate = cov1)$par.cor)
  
  r = cor(cov2,resist_to.I$Multifunc,use = "pairwise.complete.obs")
  r2 = cor(resist_to.I[[tabi]][j2],resist_to.I$Multifunc[j2],use = "pairwise.complete.obs")
  resI[delete_subset,] = c(r=r,df = df, r2=r2, part.r = partial.corr(x = resist_to.I[[tabi]], y = resist_to.I$Multifunc, covariate = cov2)$par.cor)
  
  r = cor(cov3,resist_to.D$Multifunc,use = "pairwise.complete.obs")
  r2 = cor(resist_to.D[[tabi]][j3],resist_to.D$Multifunc[j3],use = "pairwise.complete.obs")
  resD[delete_subset,] = c(r=r,df = df, r2=r2, part.r = partial.corr(x = resist_to.D[[tabi]], y = resist_to.D$Multifunc, covariate = cov3)$par.cor)
  
  r = cor(cov4,stab$Multifunc,use = "pairwise.complete.obs")
  r2 = cor(stab[[tabi]][j4],stab$Multifunc[j4],use = "pairwise.complete.obs")
  Sim[delete_subset,] = c(r=r,df = df, r2=r2, part.r = partial.corr(x = stab[[tabi]], y = stab$Multifunc, covariate = cov4)$par.cor)
  
}

#SF_pair=list(inv=inv,resD=resD,resI=resI,Sim=Sim,id = V)
#save(SF_pair,file = "SF_pair.RData")


data3=NULL
for(i in 1:3){
  data3 = rbind(data3,cbind(SF_pair[[i]],SF_pair$id,names(SF_pair)[i]))
}
names(data3)[6:7] = c("id","aspect")

data3$n=apply(data3$id%>%strsplit2mat(split = "-"),1,function(x){length(unique(x[x!=""]))})
data3$sp = data3$id

all(rownames(data3) == data3$sp)
data3$dif = as.numeric(data3$r0) - as.numeric(data3$r1)
data3$sort = 1:nrow(data3)
for(i in names(data3)){
  if(all(!is.na(as.numeric(data3[[i]])))) data3[[i]] = as.numeric(data3[[i]])
}

for(i in 1:4){
  for(j in unique(data3$aspect))
    data3[data3$n == i & data3$aspect == j,] = sortby(data3[data3$n == i & data3$aspect == j,],"dif")
}


library(ggplot2)
p.SF <- ggplot(data3,aes(x = as.numeric(n), y = dif,color=as.factor(n)))+
  facet_wrap(~aspect)+
  geom_jitter(width = 0.2,height = 0,size= 0.1,shape = 1)+
  geom_boxplot(aes(group=as.factor(n)),fill=NA,color = "black",size = 0.2)+
  labs(y="sim change", x = "species deletion")+
  #scale_x_continuous(breaks = sort(c(-seq(bk.y,5,bk.y),seq(0,5,bk.y))))+
  theme_bw()+
  theme(text=element_text(size = 7, family = "sans", face = "plain"),
        panel.grid.major.y = element_blank(),    panel.grid.minor.y =  element_blank(),  panel.grid.minor.x =  element_blank(),
        panel.background = element_rect(color = 'black', fill = 'transparent'),legend.position = "none")
p.SF
#'+  geom_smooth(method = "loess")'
#ggsave(filename = "SF_stab_change_after_deletion.pdf",p,width = 160,height = 50,units = "mm")

realtive_importance = matrix(NA,length(uni_class),4,dimnames = list(uni_class,1:4))
#colnames(realtive_importance) = uni_class
xx = strsplit2mat(data3$id,split = "-") %>% as.data.frame()
for(i in 1:ncol(xx)) xx[[i]] %<>% as.numeric()

for(i in 1:4){
  j1 = data3$n == i
  for(ii in 1:length(uni_class)){
    temp = data3[apply(xx==ii,1,sum,na.rm=T) > 0 & j1,]
    realtive_importance[ii,i] = mean(-temp$dif/temp$r1)
  }
}

order_tab = rowsum(TAB,group = taxa[rownames(TAB),]$Order)[uni_class,]
realtive_importance%<>%as.data.frame()
realtive_importance$abund = apply(order_tab,1,function(x){mean(x[x>0])})
realtive_importance$sd = apply(order_tab,1,function(x){sd(x[x>0])})
realtive_importance[,c("abund","sd")] = realtive_importance[,c("abund","sd")]/sum(realtive_importance$abund)
sum(realtive_importance$abund)
realtive_importance$cv = realtive_importance$sd / realtive_importance$abund
realtive_importance$sp = rownames(realtive_importance)
rownames(realtive_importance) == uni_class

plot_df = NULL
for(i in 1:4){
  plot_df = rbind(plot_df,cbind(value = realtive_importance[,i],realtive_importance[,c(5:8)], delition = i))
}
ggplot(plot_df,aes(abund, value, color = as.factor(delition)))+
  geom_jitter(width = 0)
ggplot(plot_df,aes(sd, value, color = as.factor(delition)))+
  geom_jitter(width = 0)
ggplot(plot_df,aes(cv, value, color = as.factor(delition)))+
  geom_jitter(width = 0)


realtive_importance$sp[realtive_importance$sp == "Poduromorpha"] = paste0(rep("a",20),collapse = "")
plot_df$sp[plot_df$sp == "Poduromorpha"] = paste0(rep("a",20),collapse = "")

plot_df$sp = factor(plot_df$sp,levels = realtive_importance$sp[order(realtive_importance$abund,decreasing = T)])
realtive_importance$sp = factor(realtive_importance$sp,levels = realtive_importance$sp[order(realtive_importance$abund,decreasing = T)])
realtive_importance = realtive_importance[realtive_importance$sp %in% levels(realtive_importance$sp)[1:10],]
plot_df = plot_df[plot_df$sp %in% levels(realtive_importance$sp)[1:10],]

p.SF_importantSP = ggplot(realtive_importance,aes(abund, sp))+
  #geom_point()+
  geom_bar(aes(abund,sp),stat='identity', position='stack',fill="lightgray",color=NA,size=0.25,width = 0.75) +
  geom_errorbar(aes(xmin = abund -sd, xmax = abund +sd),width = 0.25,size=0.375) +
  geom_point(data = plot_df,aes(value,sp, color = as.factor(delition)),shape = 1)+
  theme_bw() +
  theme(text=element_text(size = 8, family = "sans", face = "plain"),axis.title.x = element_blank(),
        panel.grid.major.x = element_blank(),    panel.grid.minor.x =  element_blank(),  panel.grid.minor.y =  element_blank(),
        panel.background = element_rect(color = 'black', fill = 'transparent'),legend.position = "none")

data3$realtive_importance = -data3$dif/data3$r1
data3$abund = data3$sd = data3$cv = 0#sapply(X=1:nrow(data3),FUN = function(X,y=strsplit(data3$sp,split = "-"),...){sum(realtive_importance[as.numeric(y[[X]]),"abund"])})
y=strsplit(data3$sp,split = "-")
for(i in 1:nrow(data3)){
  data3$abund[i] = sum(realtive_importance[as.numeric(y[[i]]),"abund"])
  data3$sd[i] = sqrt(sum(realtive_importance[as.numeric(y[[i]]),"sd"]^2))
  data3$cv[i] = mean(realtive_importance[as.numeric(y[[i]]),"cv"])
}

p3.1 = ggplot(data3[order(data3$n,decreasing = T),],aes(abund, realtive_importance, color = as.factor(n)))+
  facet_wrap(~aspect)+
  geom_point(aes(size = 5-n), size = 0.2, shape = 1, alpha = 1)+
  geom_smooth(method = "lm", se=F,size = 1) +
  theme(text=element_text(size = 7, family = "sans", face = "plain"),
        panel.grid.major.y = element_blank(),    panel.grid.minor.y =  element_blank(),  panel.grid.minor.x =  element_blank(),
        panel.background = element_rect(color = 'black', fill = 'transparent'),legend.position = "none")

p3.2 = ggplot(data3[order(data3$n,decreasing = T),],aes(sd, realtive_importance, color = as.factor(n)))+
  facet_wrap(~aspect)+
  geom_point(aes(size = 5-n), size = 0.2, shape = 1, alpha = 1)+
  geom_smooth(method = "lm", se=F,size = 1) +
  theme(text=element_text(size = 7, family = "sans", face = "plain"),
        panel.grid.major.y = element_blank(),    panel.grid.minor.y =  element_blank(),  panel.grid.minor.x =  element_blank(),
        panel.background = element_rect(color = 'black', fill = 'transparent'),legend.position = "none")

p3.3 = ggplot(data3[order(data3$n,decreasing = T),],aes(cv, realtive_importance, color = as.factor(n)))+
  facet_wrap(~aspect)+
  geom_point(aes(size = 5-n), size = 0.2, shape = 1, alpha = 1)+
  geom_smooth(method = "lm", se=F,size = 1) +
  theme(text=element_text(size = 7, family = "sans", face = "plain"),
        panel.grid.major.y = element_blank(),    panel.grid.minor.y =  element_blank(),  panel.grid.minor.x =  element_blank(),
        panel.background = element_rect(color = 'black', fill = 'transparent'),legend.position = "none")

data = data3[data3$n ==1,]
data$sp = uni_class[as.numeric(data$sp)]
data$sp[data$sp == "Poduromorpha"] = rep("a",20)%>%paste0(collapse = "")
data$sp = factor(data$sp,levels(realtive_importance$sp))
data = data[data$sp %in% levels(realtive_importance$sp)[1:10],]

p3.4=ggplot(realtive_importance,aes(abund, sp))+
  #geom_point()+
  geom_errorbar(aes(xmin = abund/3, xmax = abund +sd),width = 0.25,size=0.375)+
  geom_bar(aes(abund,sp),stat='identity', position='stack',fill="lightgray",color=NA,size=0.25,width = 0.75) +
  #facet_wrap(~aspect) +
  geom_point(data=data,aes(realtive_importance,sp, color = as.factor(aspect)), shape =1, alpha = 1) +
  #scale_color_manual(values = c("black","#DC143C", "#CD853F"))
  #geom_smooth(method = "lm", se=F,size = 1) +
  theme(text=element_text(size = 7, family = "sans", face = "plain"),
        panel.grid.major.y = element_blank(),    panel.grid.minor.y =  element_blank(),  panel.grid.minor.x =  element_blank(),
        panel.background = element_rect(color = 'black', fill = 'transparent'),legend.position = "none")


p3 = pp <- plot_grid(p3.1,p3.3,
                     hjust = 0,vjust = 0,
                     label_fontfamily = "sans", label_size = 7,  label_x =c(0.15,0.15), label_y =1, ncol = 2, nrow=1)
p3
#+ theme(plot.margin = unit(c(0.8,0,0,0), "cm"))

#ggsave(filename = "SF_p3.pdf",p3,width = 160, height = 37,units = "mm")

pp <- plot_grid(p1.1,p1.3,
                p2.1,p2.3,
                p3.1,p3.3,
                     hjust = 0,vjust = 0,
                     label_fontfamily = "sans", label_size = 7,  label_x =c(0.15,0.15), label_y =1, ncol = 2, nrow=3)
pp
#+ theme(plot.margin = unit(c(0.8,0,0,0), "cm"))

#ggsave(filename = "Pl_LF_SF.pdf",pp,width = 160, height = 110,units = "mm")


pp1 <- plot_grid(p.PL,p.LF,p.SF,
                 hjust = 0,vjust = 0,
                 label_fontfamily = "sans", label_size = 7,  label_x =c(0.15,0.15), label_y =1, ncol = 1, nrow=3)
pp1
#+ theme(plot.margin = unit(c(0.8,0,0,0), "cm"))

#ggsave(filename = "Pl_LF_SF_1.pdf",pp1,width = 140, height = 160,units = "mm")


p.pl_importantSP
p.LF_importantSP
p.SF_importantSP
importantSP = plot_grid(p.pl_importantSP,p.LF_importantSP,p.SF_importantSP,
                        hjust = 0,vjust = 0,
                        label_fontfamily = "sans", label_size = 7,  label_x =c(0.15,0.15), label_y =1, ncol = 3, nrow=1)
importantSP

#ggsave(filename = "importantSP.pdf",importantSP,width = 250, height = 80,units = "mm")


importantSP1 = plot_grid(p1.4,p2.4,p3.4,
                        hjust = 0,vjust = 0,
                        label_fontfamily = "sans", label_size = 7,  label_x =c(0.15,0.15), label_y =1, ncol = 3, nrow=1)
importantSP1

#ggsave(filename = "importantSP1.pdf",importantSP1,width = 200, height = 80,units = "mm")


