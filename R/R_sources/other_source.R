UniFrac = function (OTU, tree, alpha = 1,total.normalize = FALSE, weighted = c("both","weighted","unweighted"),
                    membory.size=8,quiet=TRUE){
  #OTU = tabs$Pl.co
  #tree = trees$Pl.co
  #aaa = UniFrac(OTU,tree,total.normalize = F)
  #a1 = aaa$w_unifrac[,,1]
  #aaa = UniFrac(OTU,tree,total.normalize = F)
  #a1 = aaa$unifrac[,,1]
  #bbb = GUniFrac2(t(OTU),tree)
  #b1 = bbb$unifracs[,,"d_1"]
  #bbb = GUniFrac(t(OTU),tree)
  #b1 = bbb$unifracs[,,"d_1"]
  #bbb = UniFrac(OTU,tree,total.normalize = F)
  #b1 = bbb$w_unifrac[,,1]
  #a2 = a1 == b1
  #mantel(a1,b1)
  if(length(dim(OTU)) == 2){
    OTU = as.matrix(OTU)
    Multi = 1
    OTU = OTU[apply(OTU,1,sum)>0,apply(OTU,2,sum)>0]
    if(total.normalize){OTU = t(t(OTU)/apply(OTU,2,sum))}
    
  }else if(length(dim(OTU)) == 3){
    Multi=dim(OTU)[3]
    OTU=OTU[apply(OTU,1,sum)>0,apply(OTU,2,sum)>0,]
    if(total.normalize){OTU=OTU/rep(apply(OTU,c(2,3),sum),each=nrow(OTU))}
    
  }else{
    stop("Input format is wrong")
    
  }
  
  species_names=rownames(OTU)
  site_names=colnames(OTU)
  n.sites = ncol(OTU)
  if (is.null(site_names)) {
    colnames(OTU) = paste("comm", 1:n, sep = "_")
    site_names=colnames(OTU)
  }
  if (!ape::is.rooted(tree))
    stop("Rooted phylogenetic tree required!")
  
  if (sum(!(species_names %in% tree$tip.label)) != 0) {
    stop("The OTU table contains unknown OTUs! OTU names\n\t\t\t\t\tin the OTU table and the tree should match!")
  }
  absent = tree$tip.label[!(tree$tip.label %in% species_names)]
  if (length(absent) != 0) {
    tree = ape::drop.tip(tree, absent)
  }
  br.len = tree$edge.length
  if(all(br.len==0)){br.len[]=1}
  
  CUM=Edge_cum(OTU,tree)
  if(!all(rownames(CUM)==1:length(br.len))){stop("unknown error with tree or OTU")}
  
  result=calculate_unifrac_array(CUM, br.len, alpha, n.sites, site_names,membory.size = membory.size,quiet=quiet)
  
  return(result)
}

calculate_unifrac_array <- function(CUM, br.len,alpha, n.sites, site_names, membory.size=8,quiet=TRUE){
  Multi=if(length(dim(CUM))==2){1}else if(length(dim(CUM))==3){dim(CUM)[3]}else{stop("cum_edge wrong")}
  dist1 = matrix(NA,n.sites,n.sites)
  diag(dist1) = 0
  dist_mat_un = dist_mat = array(0,dim = c(n.sites,n.sites,Multi))
  j = which(dist1 == 0)
  w_unifrac = unifrac = array(NA, c(n.sites, n.sites, Multi), dimnames = list(site_names, site_names, as.character(1:Multi)))
  
  if(nrow(CUM)*n.sites*(n.sites-1)/2*Multi/100000000*32/membory.size > 1){
    warning("Membory size is limited. Calculate with cycle insted !")
    CUM = array(CUM,dim = c(nrow(CUM),(n.sites+1),Multi),
                dimnames = list(rownames(CUM),c(site_names[c(1:n.sites,1)]),as.factor(1:Multi)))
    dist1=matrix(NA,n.sites,n.sites)
    diag(dist1)=0
    dist_mat_un=dist_mat=array(0,dim=c(n.sites,ncol(CUM),Multi))
    j=which(dist1==0)
    
    for (i in 1:(n.sites-1)){
      if(!quiet){print(paste(date(),"----", site_names[i], "----",round(100-(n.sites-i)^2/n.sites^2*100,2),"% finished"))}
      cum1 <- CUM[,1:(n.sites-i+1),]
      cum2 <- CUM[,(i+1):(n.sites+1),]
      sum12 <- cum1+cum2
      diff = abs(cum1 - cum2)/(sum12+10^-32)
      if(alpha==1){w = br.len * sum12}else if(alpha==0){w = br.len}else if(alpha==0.5){w = br.len * sqrt(sum12)}else{w = br.len * sum12^alpha}
      dist_mat[rep(j[1:(n.sites-i+1)]+n.sites*i,Multi)+rep(((1:Multi)-1)*n.sites*ncol(CUM),each=n.sites-i+1)] = colSums(diff * w)/colSums(w)[]
      cum1_binary <- cum1>0
      cum2_binary <- cum2>0
      uni.occur = (cum1_binary + cum2_binary)==1
      w_binary <- br.len * (cum1_binary|cum2_binary)
      dist_mat_un[rep(j[1:(n.sites-i+1)]+n.sites*i,Multi)+rep(((1:Multi)-1)*n.sites*ncol(CUM),each=n.sites-i+1)] <- colSums(uni.occur * w_binary)/colSums(w_binary)
    }
    
    rm(list = c("cum1","cum2","w","sum12","diff","cum1_binary","cum2_binary","w_binary","uni.occur"))
    
    w_unifrac[]=dist_mat[1:n.sites,1:n.sites,]#+dist_mat[1:n,(1:n)+n]
    w_unifrac[]=w_unifrac + aperm.default(w_unifrac,perm = c(2,1,3))
    unifrac[]=dist_mat_un[1:n.sites,1:n.sites,]#+dist_mat_un[1:n,(1:n)+n]
    unifrac[]=unifrac+aperm.default(unifrac,perm = c(2,1,3))
    rm("dist_mat_un","dist_mat")
    gc()
    return(list(w_unifrac = w_unifrac,unifrac=unifrac))
  }else{
    cum1 = cum2 = array(NA,dim = c(nrow(CUM),n.sites*(n.sites-1)/2,Multi), dimnames = list(1:nrow(CUM), 1:(n.sites*(n.sites-1)/2), as.character(1:Multi)))
    ind = NULL
    for (i in 1:(n.sites-1)){
      x <- length(ind)
      ind = c(ind,j[1:(n.sites-i)]+n.sites*i)
      cum1[,(x+1):length(ind),] <- CUM[,1:(n.sites-i),]
      cum2[,(x+1):length(ind),] <- CUM[,(i+1):(n.sites),]
    }
    ind = rep(ind,Multi)+rep(((1:Multi)-1)*n.sites*n.sites,length(ind))
    sum12 = cum1 + cum2
    diff = abs(cum1 - cum2)/(sum12 + 10^-32)
    cum1_binary = cum1>0;rm(list=c("cum1"))
    cum2_binary = cum2>0;rm(list=c("cum2"))
    
    if(alpha==1){w = br.len * sum12}else if(alpha==0){w = br.len}else if(alpha==0.5){w = br.len * sqrt(sum12)}else{w = br.len * sum12^alpha}
    dist_mat[ind] = colSums(diff * w)/colSums(w)
    w_unifrac[,,1:Multi] = dist_mat+aperm.default(dist_mat,perm = c(2,1,3))
    
    rm("dist_mat","diff","w")
    gc()
    
    
    uni.occur = (cum1_binary + cum2_binary) == 1
    w_binary = br.len * (cum1_binary|cum2_binary)
    rm(list=c("cum1_binary","cum2_binary"))
    dist_mat_un[ind] = colSums(uni.occur * w_binary)/colSums(w_binary)
    unifrac[,,1:Multi] = dist_mat_un+aperm.default(dist_mat_un,perm = c(2,1,3))
    
    rm("dist_mat_un","uni.occur","w_binary")
    
    
    gc()
    return(list(w_unifrac = w_unifrac,unifrac = unifrac))
  }
  
  
}
Edge_cum <- function(OTU, tree){
  
  edge_cum_index.fun <- function(edge,ntip){
    f.edge.index <- function(i, edge, a.tip, c.len, ...){
      edge.id = rep(NA,c.len)
      cycle.node <- 0;
      node.loc <- edge$target == edge$source[edge$target == i]
      while (sum(node.loc)) {
        cycle.node <- cycle.node+1
        edge.id[cycle.node] <- edge$id[node.loc]
        node.loc <- edge$target == edge$source[node.loc]
      }
      return(if(cycle.node){edge.id[1:cycle.node]}else{NULL})
    }
    a.tip = data.frame(comm.id=edge$target[edge$target %in% 1:ntip],edge.id=edge$id[edge$target %in% 1:ntip])
    
    if(length(a.tip$comm.id)<nrow(edge)){
      re = sapply(1:nrow(a.tip), function(i,...) f.edge.index(i,edge,a.tip,max(round(2*sqrt(ntip))+1,50)))
      
      a.node=matrix(NA,length(unlist(re)),2);colnames(a.node)=c("comm.id","edge.id")
      k=0
      for(i in 1:length(re)){
        len=length(re[[i]])
        if(!is.null(re[[i]])){
          a.node[(1:len)+k,1]=i
          a.node[(1:len)+k,2]=re[[i]]
          k = k + len
        }
      }
      a.node <- a.node[1:k,]
    }else{
      a.node=c()
    }
    
    return(data.frame(comm.id=c(a.tip[,1],a.node[,1]),edge.id=c(a.tip[,2],a.node[,2])))
  }
  
  tip.label = tree$tip.label
  
  if(length(dim(OTU)) == 2){
    OTU = OTU[tip.label,]
    Multi=1
  }else if(length(dim(OTU)) == 3){
    Multi=dim(OTU)[3]
    OTU=OTU[tip.label,,]
  }
  
  br.len = tree$edge.length
  ntip = length(tip.label)
  edge = as.data.frame(tree$edge)
  nbr = nrow(edge)
  names(edge) = c("source","target");edge$id = 1:nbr
  
  edge_cum_index=edge_cum_index.fun(edge,ntip)
  
  CUM=array(NA,dim=c(length(unique(edge_cum_index$edge.id)),ncol(OTU),Multi),
            dimnames = list(sort(unique(edge_cum_index$edge.id)),
                            colnames(OTU),
                            as.character(1:Multi)))
  
  if(Multi>1){
    for(OTU_i in 1:Multi){
      CUM[,,OTU_i] = rowsum(OTU[,,OTU_i][edge_cum_index$comm.id,]
                            ,group = edge_cum_index$edge.id)
    }
  }else{
    CUM[,,1] = rowsum(OTU[edge_cum_index$comm.id,],
                      group = edge_cum_index$edge.id)
  }
  #rownames(CUM) = sort(unique(edge_cum_index$edge.id))
  colnames(CUM) = colnames(OTU)
  
  return(CUM)
}

missing.data.process <- function(x,method=c("rnorm","resample.self"),replace=TRUE){
  #par(mfrow=c(2,2))
  #x=c(rnorm(10000,mean = -5,sd = 1),rnorm(10000,mean =5,sd = 1))
  #hist(x,breaks = 100)
  #x[sample.int(20000,10000)]=NA
  #hist(x[!is.na(x)],breaks = 100)
  #hist(missing.data.process(x,method = "rnorm"),breaks = 100)
  #hist(missing.data.process(x,method = "re"),breaks = 100)
  x[is.na(x)]=NA
  na.logic=is.na(x)
  if(sum(na.logic)>0&sum(!na.logic)){
    if(method[1] == "rnorm"){
      rand.norm.set=rnorm(sum(na.logic)*10,mean(x,na.rm=T),sd(x,na.rm=T))
      j=rand.norm.set<=max(x[!na.logic])&rand.norm.set>=min(x[!na.logic])
      x[na.logic]=sample(rand.norm.set[j],sum(na.logic))
    }else{
      rand.norm.set=spline(x[!na.logic],n = sum(!na.logic)*(ceiling(sum(na.logic)/sum(!na.logic))+2))
      j=rand.norm.set$y<=max(x[!na.logic])&rand.norm.set$y>=min(x[!na.logic])
      x[na.logic]=sample(rand.norm.set$y[j],sum(na.logic),replace = replace)
    }
  }
  return(x)
}


row.mat=function(data=data,n.col=NULL){
  # function: Repeat a vector (or seqs) for 'n.col' times in a matrix's colums
  if(is.data.frame(data)){data=as.matrix(data)}
  data=if(length(data)==1){(as.matrix(1:data))}else{as.matrix(data)}
  row.name=if(length(rownames(data))==length(data)){rownames(data)}else{1:length(data)}
  n.col=if(is.null(n.col)){length(data)}else{n.col}
  row_mat=matrix(rep(data,times=n.col),nrow=length(data),ncol=n.col)
  rownames(row_mat)=if(!is.null(row.name)){row.name}else{1:length(data)}
  colnames(row_mat)=1:n.col
  return(row_mat)
}
group_calculate <- function(x = NULL, Fun = NULL, group = NULL, matrixout = TRUE, ...){
  
  if(is.null(Fun)){
    stop("please inplut a function or name of function")
  }
  
  if(is.character(Fun)){
    Fun = match.fun(FUN = Fun)
  }
  
  unique_g=unique(group)
  
  re = vector(mode = "list",length = length(unique_g))
  names(re) <- unique_g
  
  temp <- Fun(x[group == unique_g[1]], ...)
  
  for(i in unique_g){
    re[[i]] = Fun(x[group == i], ...)
  }
  
  if(matrixout){
    j <- unique(unlist(lapply(re, length)))
    if(length(j)==1){
      if(length(temp)==1){
        re <- unlist(re) ;names(re) <- unique_g
      }else{
        re <- dtf(unlist(re),nrow = length(temp),ncol = length(unique_g),dimnames = list(names(temp),unique_g))
      }
    }else if(length(j)>1){
      re <- strsplit2mat(re)%>%t%>%as.data.frame
    }else{
      return(NULL)
    }
  }
  
  return(re)
}
as.symbol.pvalue<-function(pvalue,
                           sig.level=c(0.001,0.01,0.05,0.1),
                           sig.symbol=c("***","**","*","."),
                           in.sig.symbol="ns",in.sig.out=TRUE,digits = 2,
                           align=c("none","left","right")){
  x=as.vector(as.numeric(pvalue[!is.na(pvalue)]))
  symbol.p=vector(length = length(pvalue));
  if(sum(is.na(pvalue))){symbol.p[is.na(pvalue)]="NA"}
  sym.p=vector(length = sum(!is.na(pvalue)))
  sig.symbol=sortby(sig.symbol,sig.level)
  sig.level=sort(sig.level)
  
  for(i in 1:length(sig.level)){
    logic.1=x<=sig.level[i]&x>=0
    if(sum(logic.1)){
      sym.p[logic.1]=sig.symbol[i]
      x[logic.1]=-1
    }
  }
  if(in.sig.out){
    sym.p[x>-1]=round(x[x>-1],digits)
  }else{
    sym.p[x>-1]=in.sig.symbol[1]
  }
  symbol.p[!is.na(pvalue)]=sym.p
  symbol.p=char.align(symbol.p, align =align[1], digits=digits )
  return(symbol.p)
}
sortby<-function(x,sort_vector=1, MARGIN = 2, decreasing = FALSE, na.last = TRUE){
  
  if (is.vector(x)){x=as.matrix(x)}
  
  if(length(sort_vector)==1){sort_vector=if(MARGIN == 2){x[,sort_vector]}else{x[sort_vector,]}}
  od <- order(sort_vector,decreasing=decreasing)
  if(all(dim(x)>1)){
    if (MARGIN==1){
      x <- x[,od]
    }else{
      x <- x[od,]
    }
  }else{
    if (MARGIN==1){
      x <- matrix(x[,od],ncol = ncol(x),dimnames = list(rownames(x),colnames(x)[od]))
    }else{
      x <- matrix(x[od,],nrow = nrow(x),dimnames = list(rownames(x)[od],colnames(x)))
    }
  }
  return(x)
}
char.align <- function(x,align=c("none","left","right"),digits=3){
  #(l=char.align(c(rep(NA,10),runif(100)/10),align = "left"))
  #(r=char.align(c(rep(NA,10),runif(100)/10),align = "right"))
  #paste(l,r,sep=" - ")
  align.X=vector(length = length(x));
  ISNA.X=is.na(x)
  if(sum(ISNA.X)){align.X[ISNA.X]="NA"}
  
  if(sum(!ISNA.X)){x=as.vector(x[!ISNA.X])}
  if(is.numeric(x)){x=round(x,digits)}
  if(!is.character(x)){x=as.character(x)}
  
  align.X[!ISNA.X]=x
  
  if(align[1]=="right"){
    temp=max(nchar(align.X))-nchar(align.X)
    align.char=apply(as.matrix(temp),1,FUN=function(x){return(if(x==0){""}else{paste0(rep(" ",x),collapse = "")})})
    align.X=paste0(align.char,align.X)
  }else if(align[1]=="left"){
    temp=max(nchar(align.X))-nchar(align.X)
    align.char=apply(as.matrix(temp),1,FUN=function(x){return(if(x==0){""}else{paste0(rep(" ",x),collapse = "")})})
    align.X=paste0(align.X,align.char)
  }
  return(align.X)
}
"dtf" <- function(data = NA, nrow = 1, ncol = 1, byrow = FALSE, dimnames = NULL){
  if (is.object(data) || !is.atomic(data))
    data <- as.vector(data)
  
  #if(is.null(rownames(a))) rownames(a) <- 1:nrow(a)
  #if(is.null(colnames(a))) colnames(a) <- 1:ncol(a)
  return(
    as.data.frame(
      .Internal(matrix(data, nrow, ncol, byrow, dimnames, missing(nrow),
                       missing(ncol)))
    )
  )
}
partial.corr=function(x,y,covariate=NULL,multiple_regress=T,interaction=T){
  #load("div_map.RData")
  #data=div_map
  #x=div_map$log_PR
  #y=div_map$biomass
  #z=data.frame(D=2-as.numeric(as.factor(as.character(div_map$D))),I=2-as.numeric(as.factor(as.character(div_map$I))))
  #p.cor=partial.corr(x,y,z)
  #p.cor$par.cor
  #plot(p.cor$resi.xy$resx,p.cor$resi.xy$resy)
  if(interaction){multiple_regress=T}
  z=covariate
  if(is.null(z)){
    rp=cor.test(x,y)
    return(c(r=as.numeric(rp$estimate),p=rp$p.value))
  }else{
    z=as.matrix(z)
    if(ncol(z)==1){loop=NULL}else{loop=2:ncol(z)}
    if(multiple_regress){
      if(is.null(loop)){
        x.lm=lm(x~z[])
        resx=x.lm$residuals
        y.lm=lm(y~z[])
        resy=y.lm$residuals
      }else{
        if(interaction){
          temp=loop_eval(loop=loop,start.express="lm(x~z[,1]",end.express=")",cycle.express=c("*z[,","cycle.i","]"))
          x.lm=eval(parse(text = temp))
          resx=x.lm$residuals
          temp=loop_eval(loop=loop,start.express="lm(y~z[,1]",end.express=")",cycle.express=c("*z[,","cycle.i","]"))
          y.lm=eval(parse(text = temp))
          resy=y.lm$residuals
        }else{
          temp=loop_eval(loop=loop,start.express="lm(x~z[,1]",end.express=")",cycle.express=c("+z[,","cycle.i","]"))
          x.lm=eval(parse(text = temp))
          resx=x.lm$residuals
          temp=loop_eval(loop=loop,start.express="lm(y~z[,1]",end.express=")",cycle.express=c("+z[,","cycle.i","]"))
          y.lm=eval(parse(text = temp))
          resy=y.lm$residuals
        }
      }
      
      
      
      
    }else{
      resx=x
      resy=y
      for(i in 1:ncol(z)){
        resx=lm(resx~z[,i])$residuals
        resy=lm(resy~z[,i])$residuals
      }
      resx=lm(resx~z[,i])$residuals
      resy=lm(resy~z[,i])$residuals
    }
    #plot(resx,resy)
    #plot(x,y)
    rp=cor.test(resx,resy)
  }
  return(list(par.cor=c(r=as.numeric(rp$estimate),p=rp$p.value),
              resi.xy=data.frame(resx=resx,resy=resy))
  )
}

tabulate.char<-function(char,sort=F,decreasing = F,partition=10000000,weight=NULL){
  #char: a vector of characters
  #sort: logical, sort the output or not
  if(!is.null(weight)){
    weight=as.matrix(weight)
  }else{
    weight=as.matrix(rep(1,length(char)))
  }
  weight <- as.data.frame(weight)
  if(is.null(colnames(weight))) colnames(weight) <- 1:ncol(weight)
  
  char=as.character(char)
  if(length(char)>partition){
    char.l <- as.list(1:floor(length(char)/partition))
    w <- as.list(1:floor(length(char)/partition))
    for(i in 1:length(char.l)){
      char.l[[i]] <- char[(1:partition)+partition*(i-1)]
      w[[i]] <- weight[(1:partition)+partition*(i-1),]
    }
    if(length(char)%%partition){
      char.l = c(char.l,list(char[(1:(length(char)%%partition))+partition*i]))
      w <- c(w,list(weight[(1:(length(char)%%partition))+partition*i,]))
    }
    for(i in 1:length(char.l)){
      char.l[[i]] <- tabulate.char(char = char.l[[i]],sort=sort, decreasing=decreasing,partition=partition,weight=w[[i]])
    }
    char_tabulate <- merge.list_mat(char.l, method = "sum")
    return(
      if(sort){
        sortby(char_tabulate,1,decreasing = decreasing)
      }else{
        char_tabulate
      }
    )
  }
  
  
  
  '  if(sum(is.na(char))){
    temp=matrix(sum(is.na(char)),1,1,dimnames = list("N.A",NULL))
    char=char[!is.na(char)]
    weight <- weight[!is.na(char),]
  }else{
    temp=NULL
  }'
  char[is.na(char)] <- "N.A"
  char[char==""]="..empty.."
  #unique_char=unique(char);
  char_tabulate <- rowsum(weight,group = char,reorder = FALSE)
  
  ######### 1
  #char_tabulate <- vector(mode = "numeric",length = length(unique_char));names(char_tabulate)=unique_char
  
  #for(i in char){
  #char_tabulate[i] <- char_tabulate[i]+1
  #}
  
  
  #if (!is.null(temp)) char_tabulate <- rbind(char_tabulate,temp)
  #rownames(char_tabulate) <- c(unique_char,names(temp));
  if(ncol(weight)>1){
    colnames(char_tabulate) = paste0("counts.",colnames(weight))
  }else{
    colnames(char_tabulate) = "counts"
  }
  
  
  ####### 2
  '  if(length(unique(char))>1){
    char_tabulate=matrix(0,length(unique_char),1,dimnames = list(unique_char,"number"))
    for(i in char){char_tabulate[i,] <- char_tabulate[i,]+1};
    #for(i in unique_char){char_tabulate[i,] <- sum(char==i)};
    if(sort){
      char_tabulate=sortby(char_tabulate,sort_vector = char_tabulate,decreasing = decreasing)
    }
    if(!is.null(temp)){
      char_tabulate=rbind(char_tabulate,N.A=temp[1])
    }

    colnames(char_tabulate)="number"

  }else{
    if(!is.null(temp)){
      char_tabulate=c(length(char),temp[1]);names(char_tabulate)=c(unique_char,"N.A")
    }else{
      char_tabulate=length(char);names(char_tabulate)=unique_char
    }
  }'
  
  return(
    if(sort & nrow(char_tabulate)>1){
      sortby(char_tabulate,1,decreasing = decreasing)
    }else{
      char_tabulate
    }
  )
  
  
}


"strsplit2mat" <- function(x,n.col=100,split=" ",...){
  if(is.character(x)) x <- strsplit(x,split)
  if(length(x)){
    if(length(x)==1){
      return(x[[1]])
    }else{
      max_col = max(unlist(lapply(x,function(i){length(i)})))
      if(is.null(n.col)) n.col = max_col
      n.col = min(max_col,n.col)
      
      Y <- lapply(x, FUN = function(y, null.char=null.char){
        leny <- length(y)
        if(leny >= n.col){
          return(y[1:n.col])
        }else{
          c(y,null.char[1:(n.col-leny)])
        }
      }, null.char = rep("",n.col))
      
      return(
        t(matrix(unlist(Y),nrow = n.col,dimnames = list(NULL,names(x))))
      )
      
    }
  }else{
    return(NULL)
  }
  
}















