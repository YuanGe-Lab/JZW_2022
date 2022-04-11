"all.pars" <- function(model) {x <- unique(strsplit(model,split = c("[~; (=:)(\\+)(\\=)(\n)]"))[[1]]);return(x[!is.na(x)&x!=""])}

"stepwise.model" <- function(fit,data=NULL,edge.pvalue.filter=0.1, MI.cut = 1,op=NULL,missing="listwise",change.ori.model.path=TRUE,man_path=NULL){
  library(magrittr)
  extracted.model=extract.sem.model(fit)
  Model=initial.model=extracted.model$model.nolable
  model.table=extracted.model$model.table

  #fit = lavaan(model = Model,data = data,missing=missing)
  I=fitMeasures(fit,c("chisq","df","pvalue", "cfi", "rmsea","rmsea.pvalue"))
  Stepwise.regression=matrix(0,1000,12);colnames(Stepwise.regression)=c("discard","chisq","df","pvalue","rmsea","Model","added","chisq1","df1","pvalue1","rmsea1","Model1")
  Stepwise.regression=as.data.frame(Stepwise.regression)
  Stepwise.regression[1,]=c(NA,I[c("chisq","df","pvalue","rmsea")],NA,NA,NA,NA,NA,NA,NA)
  n=1
  adding="noadding";re.add=FALSE
  re2=NULL
  for(cycle.k in 1:500){
    n=n+1
    Est <- parameterEstimates(fit, ci = FALSE, standardized = TRUE)

    judge_retain=Est$pvalue<edge.pvalue.filter&!is.na(Est$pvalue)

    sig.edge=Est[judge_retain,]
    sig.edge$label[sig.edge$label!=""] = paste(sig.edge$label[sig.edge$label!=""],sep="")

    all.eff.pars = sig.edge[sig.edge$op=="~",c(1,3)] %>%unlist %>% as.character()
    cov.pars =  sig.edge$op=="~~"#& sig.edge$lhs != sig.edge$rhs
    non.eff.covs = cov.pars&(!sig.edge$lhs%in%all.eff.pars | ! sig.edge$rhs%in%all.eff.pars)
    #which(non.eff.covs)
    sig.edge <- sig.edge[!non.eff.covs,]

    re=paste(sig.edge$lhs,sig.edge$op,sig.edge$rhs,sep="")
    Model=add.variance(paste0(re,collapse ="; "), add.path=man_path)

    if(sum(!judge_retain)){
      insig.edge=Est[!judge_retain,]
      insig.edge$label[insig.edge$label!=""]=paste(insig.edge$label[insig.edge$label!=""],sep="")
      re1=paste(insig.edge$lhs,insig.edge$op,insig.edge$rhs,sep="")
      if(sum(re1%in%re2)){
        re1=sort(unique(re1[re1%in%re2]))
        discard=paste0(re1,collapse ="; ")
      }else{
        discard=""
      }

    }else{
      discard=""
    }
    #print(paste("discarding : ",discard))
    fit = lavaan(model = Model,data = data,missing=missing)
    I=fitMeasures(fit,c("chisq","df","pvalue", "cfi", "rmsea","rmsea.pvalue"))
    Stepwise.regression$discard[n]=discard
    Stepwise.regression$Model[n]=Model
    Stepwise.regression[n,c("chisq", "df", "pvalue", "rmsea")]=I[c("chisq", "df", "pvalue", "rmsea")]

    print( Stepwise.regression[n,][,c("discard","chisq","df","pvalue","rmsea")])

    MI <- modificationIndices(fit,sort.=TRUE,minimum.value= (0.0001),op = NULL)

    rownames(MI)=MI$experession.nolable=paste(MI$lhs,MI$op,MI$rhs,sep="")
    temp=model.table[model.table$op=="~",]
    temp2=model.table[model.table$op=="~~",]

    MI=MI[!MI$experession.nolable%in%
            c(paste0(temp$rhs,"~~",temp$lhs),
              paste0(temp$lhs,"~~",temp$rhs),
              paste0(temp$rhs,temp$op,temp$lhs),
              paste0(temp2$lhs,"~",temp2$rhs),
              paste0(temp2$rhs,"~",temp2$lhs)
            ),]


    #sbuset.MI1$experession=model.table[rownames(sbuset.MI1),]$experession
    if(nrow(MI)>0){
      for(j in 1:length(MI$experession.nolable)){
        adding=MI$experession.nolable[j]
        if(grepl("~~",adding)){
          x <- extract.sem.model(fit)%>%.$model.nolable%>%strsplit(split = "[; ]")%>%unlist%>%.[!.==""]%>%
            .[!grepl("~~",.,fixed = T)]%>%strsplit(split = "[; ~]")%>%unlist%>%.[!.==""]%>%unique
          if(all(!strsplit(adding,"~~")%>%unlist%in%x)){adding="";next}
        }
        re.add=grepl(pattern = adding,x = discard )&Stepwise.regression$pvalue[n]<=Stepwise.regression$pvalue[n-1]
        re2=sort(unique(c(re,adding)))
        Model1=add.variance(paste0(re2,collapse ="; "),add.paths = man_path)
        if(!(Model1%in%unique(Stepwise.regression$Model1))){adding="";break}
      }

    }else{
      break
    }
    #print(paste("adding : " , adding))
    Stepwise.regression$Model1[n] = Model = Model1

    fit = lavaan(model = Model,data = data,missing=missing)
    I=fitMeasures(fit,c("chisq","df","pvalue", "cfi", "rmsea","rmsea.pvalue"))
    Stepwise.regression$added[n]=adding

    Stepwise.regression[n,c("chisq1", "df1", "pvalue1", "rmsea1")]=I[c("chisq", "df", "pvalue", "rmsea")]
    if(Model1%in%Stepwise.regression$Model1[1:(n-1)]){
      break
    }
    print( Stepwise.regression[n,][,c("added","chisq1","df1","pvalue1","rmsea1")])
  }

  Stepwise.regression=Stepwise.regression[1:n,]
  #if(max(Stepwise.regression$rmsea,na.rm = TRUE)>max(Stepwise.regression$rmsea1,na.rm = TRUE)){
  TEMP1=Stepwise.regression[,(1:6)+6];colnames(TEMP1)=colnames(Stepwise.regression)[1:6];rownames(TEMP1)=paste0("add",rownames(TEMP1))
  Stepwise.regression1=rbind(Stepwise.regression[,1:6],TEMP1)
  Stepwise.regression1=Stepwise.regression1[!is.na(Stepwise.regression1$Model),]

  p1_j=Stepwise.regression1$pvalue>=min(0.05,max(Stepwise.regression1$pvalue,na.rm = TRUE))&!is.na(Stepwise.regression1$pvalue)&
    Stepwise.regression1$rmsea<=max(0.08,min(Stepwise.regression1$rmsea,na.rm = TRUE));

  maxdf1=Stepwise.regression1$df==max(Stepwise.regression1$df[p1_j])
  Stepwise.regression1=Stepwise.regression1[maxdf1,]
  if(sum(maxdf1)>1){
    Stepwise.regression1=Stepwise.regression1[which(Stepwise.regression1$pvalue==max(Stepwise.regression1$pvalue))[1],]
  }


  best_Model=Stepwise.regression1$Model
  #}else{
  #best_Model=Stepwise.regression[2:n,]$Model1[which(round(Stepwise.regression[2:n,]$pvalue1,10)==round(max(Stepwise.regression[2:n,]$rmsea1,na.rm = TRUE),10))[1]]
  #}
  bestfit=lavaan(best_Model, data = data,missing=missing)
  all.p=all.pars(best_Model)
  x=extract.sem.model(bestfit)%>%.$model.nolable%>%strsplit(split = "[; ]")%>%unlist%>%.[!.==""]%>%
    .[!grepl("~~",.,fixed = T)]%>%strsplit(split = "[; ~]")%>%unlist%>%.[!.==""]%>%unique
  if(sum(!all.p%in%x)){
    best_Model <- delete.variables.in.model(bestfit,delete.pars = all.p[!all.p%in%x])
  }

  bestfit=lavaan(best_Model, data = data,missing=missing)
  bestfit.model.table=extract.sem.model(bestfit)
  addpaths=bestfit.model.table$model.table[!(bestfit.model.table$model.table$experession.nolable%in%model.table$experession.nolable),]
  return(list(bestfit=bestfit,best.fit.Model=best_Model,addpaths=addpaths,Stepwise.regression=Stepwise.regression,initial.model=list(initial.model=initial.model,model.table=model.table),bestfit.model.table=bestfit.model.table$model.table))
}

"delete.variables.in.model" <- function(fit=NULL,delete.pars=NULL,delete.paths=NULL,export.model.table=FALSE){
  #delete.paths=list(lhs=c("process"),rhs=c("SF","LF"),op="~")
  par.table <- partable(fit)
  par.table=par.table[par.table$op!=":=",]
  par.table$expression.nolable=with(data = par.table,paste0(lhs,op,rhs))

  no.lable=par.table$label==""
  if(sum(no.lable)){par.table$label[no.lable]=paste("par",1:sum(no.lable),sep="")}

  judge1=par.table$rhs%in%delete.pars | par.table$lhs%in%delete.pars

  if(!is.null(delete.paths)){
    delete.paths=as.data.frame(cbind(lhs=delete.paths$lhs,op=delete.paths$op,rhs=delete.paths$rhs))
    for(i in 1:nrow(delete.paths)){
      if(length(unique(delete.paths[i,c("lhs","rhs")]))==2){
        judge=judge1|(par.table$lhs==delete.paths[i,"lhs"]&par.table$rhs==delete.paths[i,"rhs"]&par.table$op==delete.paths[i,"op"])
      }
    }
  }else{
    judge=judge1
  }
  if(sum(judge)){
    print(paste0(sum(judge)," path",if(sum(judge)>1){"s were "}else{" was "},"deleted. < ",
                 paste(par.table$expression.nolable[judge],collapse = " & ")," >" ))
    par.table=par.table[!judge,]
  }else{
    warning("no paths deleted")
  }



  model=paste0(sort(paste(par.table$lhs,par.table$op,par.table$rhs,sep="")),collapse ="; ")
  if(export.model.table){
    model.table=par.table[,c("lhs","op","rhs","label")]
    model.table$experession=paste(par.table$lhs,par.table$op,par.table$label,"*",par.table$rhs,sep="")
    model.table$experession.nolable=paste(par.table$lhs,par.table$op,par.table$rhs,sep="")
    rownames(model.table)=model.table$experession.nolable
    return(list(model=model,model.table=model.table))
  }else{
    return(model=model)
  }
}

"group.effects.sem" <- function(fit,rh=NULL,lh=NULL,model=NULL,
                                group.rhs=NULL,
                                group.lhs=NULL,
                                via=NULL,exclude.paths=NULL,exclude.pars=NULL){


  equal.eff <- NULL
  #effects.model(rh=c("D","I","PR"),lh=c("Multi.stibility","M_Wellbeing"),edges,model=MD)
  par.tab_model=extract.sem.model(fit)
  par.tab=par.tab_model$model.table
  edges=par.tab[par.tab$op=="~",]
  edges0=edges
  if(!is.null(exclude.paths)){
    judge_exclude.paths=edges$experession.nolable%in%gsub(pattern = " ",replacement = "",x = exclude.paths)
  }else{
    judge_exclude.paths=FALSE
  }
  if(!is.null(exclude.pars)){
    judge_exclude.paths=judge_exclude.paths|(edges$lhs%in%exclude.pars|edges$rhs%in%exclude.pars)
  }
  if(sum(judge_exclude.paths)){
    print(c(paste(sum(judge_exclude.paths),if(sum(judge_exclude.paths)==1){"path <"}else{"paths <"},
                  paste(edges$experession.nolable[judge_exclude.paths],collapse = " & "),
                  if(sum(judge_exclude.paths)==1){"> was "}else{"> were "},"not included in the calculation of the direct, indirect or total effects")))

    edges=edges[!judge_exclude.paths,]
  }

  if(is.null(rh)){
    rh=rownames(sortby(tabulate.char(edges$rhs),sort_vector = 1,decreasing = T))
  }
  if(is.null(lh)){
    lh=rownames(sortby(tabulate.char(edges$lhs),1))
  }

  j1=rh%in%lh
  if(sum(j1)>0&sum(j1)<length(rh)){rh=c(rh[!j1],rh[j1])}
  j2=lh%in%rh
  if(sum(j2)>0&sum(j2)<length(lh)){lh=c(lh[j2],lh[!j2])}

  rh.id=NULL
  for(i in 1:length(rh)){
    rh.id[rh[i]]=i
  }
  lh.id=NULL
  for(j in 1:length(lh)){
    lh.id[lh[j]]=j
  }


  if(!is.null(group.rhs)){
    if(is.character(group.rhs)){
      group.rhs=gsub(" ","",group.rhs)
      group.rhs=strsplit(x = group.rhs,split = "\\+")
    }
    temp=NULL
    sorunique=NULL
    equal.eff=NULL
    for(i in names(group.rhs)){
      lh.rh=sort(unique(group.rhs[[i]]))
      if(length(lh.rh)>1){
        lh.rh=lh.rh[!lh.rh %in% edges[edges$lhs%in%lh.rh & edges$rhs%in%lh.rh,]$lhs]
        if(length(lh.rh)>0){
          if(!paste0(lh.rh,collapse = "")%in%sorunique){
            temp[[i]]=lh.rh
          }else{
            equal.eff[[i]]=paste0(lh.rh,collapse = "")
          }

        }
      }else if(length(lh.rh)==1){
        if(!paste0(lh.rh,collapse = "")%in%sorunique){
          temp[[i]]=lh.rh
        }else{
          equal.eff[[i]]=paste0(lh.rh,collapse = "")
        }
      }
      if(!is.null(lh.rh)){sorunique=unique(c(sorunique,paste0(lh.rh,collapse = "")))}
      lh.rh=NULL
    }
    group.rhs=temp
  }else{
    group.rhs=as.list(rh);
    names(group.rhs)=rh
  }
  if(!is.null(equal.eff)){
    for(i in names(equal.eff)){
      log.i=as.matrix(lapply(group.rhs,FUN = paste0,collapse=""))
      equal.eff[[i]]=rownames(log.i)[equal.eff[[i]]==log.i]
    }
  }

  if(!is.null(group.lhs)){
    if(is.character(group.lhs)){
      group.lhs=gsub(" ","",group.lhs)
      group.lhs=strsplit(x = group.lhs,split = "\\+")
    }
    temp=NULL
    sorunique=NULL
    for(i in names(group.lhs)){
      lh.rh=sort(unique(group.lhs[[i]]))
      if(length(lh.rh)>0 & !paste0(lh.rh,collapse = "")%in%sorunique){
        temp[[i]]=lh.rh
      }
      if(!is.null(lh.rh)){sorunique=unique(c(sorunique,paste0(lh.rh,collapse = "")))}
      lh.rh=NULL
    }
    group.lhs=temp
  }else{
    group.lhs=as.list(lh);
    names(group.lhs)=lh
  }

  if(!is.null(via)){
    #if(length(names(via))<length(via)){
    #names(via)=c(names(via),paste("via",(length(via)-length(names(via))+1):length(via)))
    #}
    BY=NULL
    if(is.null(names(via))){names(via)=gsub(" ","",via)}
    for(via.i in names(via)){
      BY[[via.i]]=unlist(strsplit(gsub(" ","",via[via.i]),split = "\\+"))
    }
    by1=BY[[1]];names(by1)=by1;by_pars=NULL
    for(tempi in by1){
      by_pars[[tempi]]=edges$label[edges$lhs==tempi]
    }
  }else{
    by_pars=NULL
  }

  group.effect_mat=matrix(NA,length(group.lhs),length(group.rhs),dimnames = list(names(group.lhs),names(group.rhs)))
  group.total.effect_mat=group.direct.effect_mat=group.indirect.effect_mat=group.effect_mat
  group.effect_id=group.effect_mat
  group.effect_id[]=""
  group.total.effect_id=group.direct.effect_id=group.indirect.effect_id=group.effect_id
  if(!is.null(via)){
    indirect.by.mat=indirect.by.id=NULL
    for(i in names(by_pars)){
      indirect.by.mat[[i]]=group.effect_mat;indirect.by.id[[i]]=group.effect_id
    }
    indirect.by.mat[["interac"]]=group.effect_mat;indirect.by.id[["interac"]]=group.effect_id
    indirect.by.mat[["others"]]=group.effect_mat;indirect.by.id[["others"]]=group.effect_id

  }

  group.element=NULL
  for(i in names(group.rhs)){
    for (j in names(group.lhs)){
      id.rh=rh.id[group.rhs[[i]]]
      id.lh=lh.id[group.lhs[[j]]]
      temp=paste0("total_rh",id.rh,"lh",rep(id.lh,each=length(id.rh)))
      group.element[[paste0(i,"_VS_",j)]]=temp
    }
  }

  if(is.null(model)){
    model=par.tab_model$model
  }

  mat=matrix(NA,length(rh),length(lh),dimnames = list(rh,lh))
  total.eff=direct.eff=indirect.eff=NULL
  for(i in 1:length(rh)){
    rh_i=rh[i]
    for(j in 1:length(lh)){
      lh_j=lh[j]
      effects=SEM.total.effect.expression(rh=rh_i,lh=lh_j,edges)

      total.eff[[paste("total_rh",i,"lh",j,sep="")]]=effects$total.effect.elem
      indirect.eff[[paste("indirect_rh",i,"lh",j,sep="")]]=effects$indirect.effect.elem
      direct.eff[[paste("direct_rh",i,"lh",j,sep="")]]=effects$direct.effect.elem
    }
  }

  ###########################################################
  all.factor=unique(c(rh,lh))
  direct.path.ids=indirect.path.ids=NULL
  #group.pars=c(group.rhs,group.lhs[names(group.lhs)[!(names(group.lhs)%in%names(group.rhs))]])
  if(!is.null(via)){
    for(g.rhsi in names(group.element)){
      rhilhi=unlist(strsplit(g.rhsi,split = "_VS_"));names(rhilhi)=c("rh","lh")
      indirect.factors=BY[[1]]
      indirect.path.ids[[g.rhsi]]=edges$label[ (edges$lhs%in%indirect.factors | edges$rhs%in%indirect.factors)]
      direct.path.ids[[g.rhsi]]=edges$label[!(edges$lhs%in%indirect.factors | edges$rhs%in%indirect.factors)]
      "for(g.rhsj in names(group.element)[!names(group.element)==g.rhsi]){
        indirect.path.ids[[g.rhsj]]=sort(unique(c(indirect.path.ids[[g.rhsj]],direct.path.ids[[g.rhsi ]])))
       }"
    }
  }else{
    for(g.rhsi in names(group.element)){
      rhilhi=unlist(strsplit(g.rhsi,split = "_VS_"));names(rhilhi)=c("rh","lh")

      direct.factors=c(unlist(group.rhs[rhilhi["rh"]]),unlist(group.lhs[rhilhi["lh"]]))
      if(length(direct.factors)>0){if(length(direct.factors)>length(unique(direct.factors))){warning(paste(g.rhsi,": repeated pars in rhs and lhs"))}}
      direct.factors=unique(direct.factors)
      direct.path.ids[[g.rhsi]]=edges$label[ edges$lhs%in%direct.factors & edges$rhs%in%direct.factors ]

      indirect.factors=all.factor[!all.factor%in%direct.factors]
      indirect.path.ids[[g.rhsi]]=edges$label[ (edges$lhs%in%indirect.factors | edges$rhs%in%indirect.factors)]

    }
  }
  # g.rhsi="Eco.stability_VS_Eco.stability"

  alleffect.temp.total=sapply(names(total.eff),function(iii,...) paste0(total.eff[[iii]],collapse = "+"))

  def.total=def.direct=def.indirect=NULL;def.indirect.by=NULL

  for(g.rhsi in names(group.element)){
    rhilhi=unlist(strsplit(g.rhsi,split = "_VS_"));names(rhilhi)=c("rh","lh")
    direct.group.total=indirect.group.total=NULL
    direct.group.total1=indirect.group.total1=NULL
    indirect.group.by=NULL

    for( ii in group.element[[g.rhsi]] ){
      if(!is.null(total.eff[[ii]])){
        direct.ind=gsub("total_","direct_",ii)
        inderect.ind=gsub("total_","indirect_",ii)
        direct.ii=direct.eff[[direct.ind]]
        indirect.ii=indirect.eff[[inderect.ind]]
        total.ii=total.eff[[ii]]
        if(!all(sort(unique(total.ii))==sort(unique(c(direct.ii,indirect.ii))))){stop("codes error, please debug!!!")}

        with_indirect.path.pars=vector(length = length(indirect.ii))
        with_direct.path.pars=with_indirect.path.pars

        if(!is.null(indirect.ii)&&length(indirect.ii)>0&&!is.na(indirect.ii)){
          for(tempi in 1:length(indirect.ii)){
            indirect.ii.temp=unlist(strsplit(indirect.ii[tempi],split = "\\*"))
            with_indirect.path.pars[tempi]=any(indirect.ii.temp%in%indirect.path.ids[[g.rhsi]])
            with_direct.path.pars[tempi]=all(indirect.ii.temp%in%direct.path.ids[[g.rhsi]])
          }
          if(sum(with_indirect.path.pars)){
            indirect.group.total[[ii]]=indirect.ii[with_indirect.path.pars]
          }

          if(sum(!with_direct.path.pars)){
            indirect.group.total1[[ii]]=indirect.ii[!with_direct.path.pars]
          }
        }


        if(sum(!with_indirect.path.pars)){
          direct.group.total[[ii]]=c(direct.ii,indirect.ii[!with_indirect.path.pars])
        }else if(!is.null(direct.ii)&&length(direct.ii)>0&&!is.na(direct.ii)){
          direct.group.total[[ii]]=direct.ii
        }

        if(sum(with_direct.path.pars)){
          direct.group.total1[[ii]]=c(direct.ii,indirect.ii[with_direct.path.pars])
        }else if(!is.null(direct.ii)&&length(direct.ii)>0&&!is.na(direct.ii)){
          direct.group.total1[[ii]]=direct.ii
        }


        #if(!is.null(indirect.group.total1)){print(".....")}
      }
    }

    if(!all(unlist(indirect.group.total)==unlist(indirect.group.total1))|!all(unlist(direct.group.total)==unlist(direct.group.total1))){
      stop("line 179 in 'group.effects: unknown error !")
    }

    if(!is.null(direct.group.total)){
      direct.effect.temp=sapply(names(direct.group.total),function(i,...) paste0(direct.group.total[[i]],collapse = "+"))
      group.ele.i=direct.effect.temp[names(direct.effect.temp)%in%group.element[[g.rhsi]]]
      group.ele.i=group.ele.i[!is.na(group.ele.i)&group.ele.i!=""]
      if(!is.null(group.ele.i)&length(group.ele.i)>0){
        group.direct.effect_mat[rhilhi["lh"],rhilhi["rh"]]=1
        group.direct.effect_id[rhilhi["lh"],rhilhi["rh"]]=paste0("direct_",g.rhsi)
        def.direct=c(def.direct,paste(paste0("direct_",g.rhsi),paste0(group.ele.i,collapse = "+"),sep=":="))
      }

    }

    if(!is.null(indirect.group.total)){
      if(!is.null(via)){
        name.ind=names(indirect.group.total)
        name.ind.i=0
        for(tempi in indirect.group.total){
          name.ind.i=name.ind.i+1
          name=name.ind[name.ind.i]
          for(tempk in 1:length(tempi)){
            temp1=strsplit(tempi[tempk],"\\*")[[1]]
            judge_in=vector(length = length(by_pars));names(judge_in)=names(by_pars)
            for(tempj in names(by_pars)){
              judge_in[tempj]=any(temp1[1:(length(temp1)-1)]%in%by_pars[[tempj]])
            }

            if(sum(judge_in)==1){
              indirect.group.by[[names(judge_in)[judge_in]]][[name]]=
                c(indirect.group.by[[names(judge_in)[judge_in]]][[name]],tempi[tempk])
            }else if(sum(judge_in)>1){
              indirect.group.by[["interac"]][[name]]=c(indirect.group.by[["interac"]][[name]],tempi[tempk])
            }else{
              indirect.group.by[["others"]][[name]]=c(indirect.group.by[["others"]][[name]],tempi[tempk])
            }
          }
        }

        for(tempi in names(indirect.group.by)) {
          indirect.group.by.temp=sapply(names(indirect.group.by[[tempi]]),function(i,...) paste0(indirect.group.by[[tempi]][[i]],collapse = "+"))
          group.ele.i=indirect.group.by.temp[names(indirect.group.by.temp)%in%group.element[[g.rhsi]]]
          group.ele.i=group.ele.i[!is.na(group.ele.i)&group.ele.i!=""]
          if(!is.null(group.ele.i)&length(group.ele.i)>0){
            indirect.by.mat[[tempi]][rhilhi["lh"],rhilhi["rh"]]=1
            indirect.by.id[[tempi]][rhilhi["lh"],rhilhi["rh"]]=paste0("by.",tempi,"_",g.rhsi)
            def.indirect.by[[tempi]]=c(def.indirect.by[[tempi]],paste(paste0("by.",tempi,"_",g.rhsi),paste0(group.ele.i,collapse = "+"),sep=":="))
          }
        }
      }

      indirect.effect.temp=sapply(names(indirect.group.total),function(i,...) paste0(indirect.group.total[[i]],collapse = "+"))
      group.ele.i=indirect.effect.temp[names(indirect.effect.temp)%in%group.element[[g.rhsi]]]
      group.ele.i=group.ele.i[!is.na(group.ele.i)&group.ele.i!=""]
      if(!is.null(group.ele.i)&length(group.ele.i)>0){
        group.indirect.effect_mat[rhilhi["lh"],rhilhi["rh"]]=1
        group.indirect.effect_id[rhilhi["lh"],rhilhi["rh"]]=paste0("indirect_",g.rhsi)
        def.indirect=c(def.indirect,paste(paste0("indirect_",g.rhsi),paste0(group.ele.i,collapse = "+"),sep=":="))
      }
    }

    total.effect.total=total.eff[names(total.eff)%in%group.element[[g.rhsi]]]
    if(!is.null(total.effect.total)){
      total.effect.temp=sapply(names(total.effect.total),function(i,...) paste0(total.effect.total[[i]],collapse = "+"))
      group.ele.i=total.effect.temp[names(total.effect.temp)%in%group.element[[g.rhsi]]]
      group.ele.i=group.ele.i[!is.na(group.ele.i)&group.ele.i!=""]
      if(!is.null(group.ele.i)&length(group.ele.i)>0){
        group.total.effect_mat[rhilhi["lh"],rhilhi["rh"]]=1
        group.total.effect_id[rhilhi["lh"],rhilhi["rh"]]=paste0("total_",g.rhsi)
        def.total=c(def.total,paste(paste0("total_",g.rhsi),paste0(group.ele.i,collapse = "+"),sep=":="))
      }
    }

  }
  group.total.effect_mat[is.na(group.total.effect_mat)]=0
  group.direct.effect_mat[is.na(group.direct.effect_mat)]=0
  group.indirect.effect_mat[is.na(group.indirect.effect_mat)]=0

  effects.group.total=paste0(def.total,collapse =";")
  effects.group.direct=paste0(def.direct,collapse =";")
  effects.group.indirect=paste0(def.indirect,collapse =";")

  effects=paste0(c(effects.group.total,effects.group.direct,effects.group.indirect),collapse = ";")


  effects.group.indirect.by=NULL;indirect.effects.by=NULL
  if(!is.null(via)){
    for(i in names(def.indirect.by)){
      indirect.by.mat[[i]][is.na(indirect.by.mat[[i]])]=0
      if(!is.null(def.indirect.by[[i]])){effects.group.indirect.by[[i]]=paste0(def.indirect.by[[i]],collapse =";")}
      indirect.effects.by[[i]]$mat=indirect.by.mat[[i]]
      indirect.effects.by[[i]]$id=indirect.by.id[[i]]
    }
    indirect.by.mat=indirect.by.mat[names(def.indirect.by)]
    indirect.by.id=indirect.by.id[names(def.indirect.by)]



    effects.by=paste0(effects.group.indirect.by,collapse = ";")
  }else{
    effects.by=NULL
  }




  model.with.effects=if(is.null(model)){NULL}else{paste(model,effects,effects.by,sep=";")}




  return(
    list(model=list(model.with.effects=model.with.effects, effects.defined=effects),
         total.effects=list(mat=group.total.effect_mat,id=group.total.effect_id),
         direct.effects=list(mat=group.direct.effect_mat,id=group.direct.effect_id),
         indirect.effects=list(mat=group.indirect.effect_mat,id=group.indirect.effect_id),
         indirect.effects.by=indirect.effects.by,
         equal.eff=equal.eff)

  )


}

"lavaan.path.stat" <- function(fit,data,effect.type=c("total.effects", "direct.effects", "indirect.effects"),
                               group.rhs=NULL,
                               group.lhs=NULL,
                               via=NULL,
                               rh=NULL,lh=NULL,standardized = TRUE,
                               bootstrap=FALSE,reps=1000,
                               return.boot=TRUE,missing="listwise",boot.depth=NULL,permutation=FALSE,
                               ind.in.dist=NULL,dist.list=NULL,boot.perm.simultaneous = FALSE,
                               strata=NULL,exclude.paths=NULL,exclude.pars=NULL){
  eff.type <- effect.type
  if(!is.null(ind.in.dist)&length(ind.in.dist)!=nrow(data)){stop("ind.in.dist input wrong")}
  MD.with.effects=group.effects.sem(fit,rh=NULL,lh=NULL,model=NULL,
                                    group.rhs=group.rhs,
                                    group.lhs=group.lhs,
                                    via=via,exclude.paths=exclude.paths,exclude.pars=exclude.pars)
  if(!is.null(MD.with.effects$indirect.effects.by)){
    indirect.by_ele=paste0("indirect.by_",names(MD.with.effects$indirect.effects.by))
    for(i in names(MD.with.effects$indirect.effects.by)){
      MD.with.effects[[paste0("indirect.by_",i)]]=MD.with.effects$indirect.effects.by[[i]]
    }
    eff.type=c(effect.type,indirect.by_ele)
  }
  fit = lavaan(model = MD.with.effects$model$model.with.effects,
               data = data,missing=missing, auto.var=TRUE, auto.fix.first=TRUE,auto.cov.lv.x=TRUE)
  Est <- parameterEstimates(fit, ci = FALSE, standardized = standardized,remove.system.eq = T,remove.nonfree=T,rsquare=T,fmi=F)
  effects.Est=Est[Est$op==":=",]
  rownames(effects.Est)=effects.Est$label

  if(permutation){boot.depth=NULL}
  if(!bootstrap){boot.perm.simultaneous=FALSE}

  effects=se.effects=p.effects=std_effects=se.std_effects=p.std_effects=judge.non.zero=NULL
  effects.summary=std.effects.summary=NULL
  if(!bootstrap){
    for(i in eff.type){
      temp.mat=MD.with.effects[[i]]$mat
      judge.non.zero[[i]]=which(temp.mat>0)
      temp.mat[]=NA
      effects[[i]]=se.effects[[i]]=p.effects[[i]]=temp.mat


      est=effects.Est[MD.with.effects[[i]]$id[judge.non.zero[[i]]],]

      effects[[i]][judge.non.zero[[i]]]=est$est
      se.effects[[i]][judge.non.zero[[i]]]=est$se
      p.effects[[i]][judge.non.zero[[i]]]=est$pvalue
      effects[[i]][is.na(effects[[i]])]=0

      effects.summary[[i]]=temp.mat
      summary.temp=paste(char.align(effects[[i]],align = "right",digits = 2),as.symbol.pvalue(p.effects[[i]],align = "left",digits = 2))
      effects.summary[[i]][]=summary.temp[]
      effects.summary[[i]] <- as.data.frame(effects.summary[[i]])


      if(standardized){
        std_effects[[i]]=se.std_effects[[i]]=p.std_effects[[i]]=temp.mat
        std_effects[[i]][judge.non.zero[[i]]]=est$std.all
        std_effects[[i]][is.na(std_effects[[i]])]=0

        std.effects.summary[[i]]=temp.mat
        summary.temp=paste(char.align(std_effects[[i]],align = "right",digits = 2),as.symbol.pvalue(p.effects[[i]],align = "left",digits = 2))
        std.effects.summary[[i]][]=summary.temp[]
        std.effects.summary[[i]] <- as.data.frame(std.effects.summary[[i]])
      }
    }


    if(standardized){
      return.result <- list(effects=effects,se.effects=se.effects,p.effects=p.effects,
                            std_effects=std_effects,
                            effects.summary=effects.summary,
                            std.effects.summary=std.effects.summary)
    }else{
      return.result <- list(effects=effects,se.effects=se.effects,p.effects=p.effects,
                            effects.summary=effects.summary)
    }
    # return.result <- list(effects=effects,se.effects=se.effects,p.effects=p.effects,std_effects=std_effects))
  }else{

    all.pars=unique(c(Est[!Est$op==":=",]$lhs,Est[!Est$op==":=",]$rhs))
    data.filter=data[,all.pars]
    effects_array=std_effects_array=NULL
    judge.non.zero=len=NULL

    for(i in eff.type){
      temp.mat=MD.with.effects[[i]]$mat;temp.mat[temp.mat==0]=NA
      effects_array[[i]]=array(NA,dim = c(nrow(temp.mat),ncol(temp.mat),reps))
      judge.non.zero[[i]]=which(MD.with.effects[[i]]$mat>0)
      len[[i]]=length(MD.with.effects[[i]]$mat)
      effects[[i]]=se.effects[[i]]=p.effects[[i]]=temp.mat

      effects.summary[[i]]=temp.mat;effects.summary[[i]][]=NA

      if(standardized){
        std_effects_array[[i]]=effects_array[[i]]
        std_effects[[i]]=se.std_effects[[i]]=p.std_effects[[i]]=temp.mat
        std.effects.summary[[i]]=effects.summary[[i]]
      }
    }
    re.samp.n=0
    for(reps.i in 1:reps){

      fit=NULL
      boot.data<-if(is.null(boot.depth)){
        data.filter[sample.int(nrow(data.filter),replace = T),]
      }else{
        data.filter[sample.int(n = nrow(data.filter),size = boot.depth,replace = T),]
      }
      #if(permutation){boot.data=apply(boot.data,MARGIN = 2,function(x) x[sample.int(length(x),replace = FALSE)])}
      fit = lavaan(model = MD.with.effects$model$model.with.effects,
                   data = boot.data, missing=missing, auto.var=TRUE, auto.fix.first=TRUE,auto.cov.lv.x=TRUE)
      while(is.null(fit)){
        re.samp.n=re.samp.n+1
        boot.data<-if(is.null(boot.depth)){
          data.filter[sample.int(nrow(data.filter),replace = T),]
        }else{
          data.filter[sample.int(n = nrow(data.filter),size = boot.depth,replace = T),]
        }
        #if(permutation){boot.data=apply(boot.data,MARGIN = 2,function(x) x[sample.int(length(x),replace = FALSE)])}
        fit = lavaan(model = MD.with.effects$model$model.with.effects,
                     data = boot.data, missing=missing, auto.var=TRUE, auto.fix.first=TRUE,auto.cov.lv.x=TRUE)
      }

      Est <- parameterEstimates(fit, ci = FALSE, standardized = standardized,remove.system.eq = T,remove.nonfree=T,rsquare=T,fmi=F)
      effects.Est=Est[Est$op==":=",]
      rownames(effects.Est)=effects.Est$label

      for(i in eff.type){
        est<-effects.Est[MD.with.effects[[i]]$id[judge.non.zero[[i]]],]

        effects_array[[i]][judge.non.zero[[i]] + (reps.i-1)*len[[i]]]=est$est
        if(standardized){
          std_effects_array[[i]][judge.non.zero[[i]]+(reps.i-1)*len[[i]]]=est$std.all
        }
      }
    }
    if(re.samp.n>0){warning(paste("Bootstraping waring:",re.samp.n,"moedles without solutions was recalculated"))}
    for(i in eff.type){
      effects[[i]][]=apply(X = effects_array[[i]],MARGIN = c(1,2),FUN = mean)
      se.effects[[i]][]=apply(X = effects_array[[i]],MARGIN = c(1,2),FUN = sd)
      p.effects[[i]][judge.non.zero[[i]]]=pt(-abs(effects[[i]][judge.non.zero[[i]]])/(se.effects[[i]][judge.non.zero[[i]]]+10^-26),df = reps-1)
      effects[[i]][is.na(effects[[i]])]=0

      summary.temp=paste(char.align(effects[[i]],align = "right",digits = 2),as.symbol.pvalue(p.effects[[i]],align = "left",digits = 2))
      effects.summary[[i]][]=summary.temp[]
      effects.summary[[i]] <- as.data.frame(effects.summary[[i]])

      if(standardized){
        std_effects[[i]][]=apply(X = std_effects_array[[i]],MARGIN = c(1,2),FUN = mean)
        se.std_effects[[i]][]=apply(X = std_effects_array[[i]],MARGIN = c(1,2),FUN = sd)
        p.std_effects[[i]][judge.non.zero[[i]]]=pt(-abs(std_effects[[i]][judge.non.zero[[i]]])/(se.std_effects[[i]][judge.non.zero[[i]]]+10^-26),df = reps-1)
        std_effects[[i]][is.na(std_effects[[i]])]=0

        summary.temp=paste(char.align(std_effects[[i]],align = "right",digits = 2),as.symbol.pvalue(p.std_effects[[i]],align = "left",digits = 2))
        std.effects.summary[[i]][]=summary.temp[]
        std.effects.summary[[i]] <- as.data.frame(std.effects.summary[[i]])
      }
    }

    if(standardized){
      if(return.boot){
        return.result <- list(effects=effects,se.effects=se.effects,p.effects=p.effects,
                              std_effects=std_effects,se.std_effects=se.std_effects,p.std_effects=p.std_effects,
                              effects.summary=effects.summary,
                              std.effects.summary=std.effects.summary,
                              std.effetc.boot=std_effects_array,
                              effects.boot=effects_array)
      }else{
        return.result <- list(effects=effects,se.effects=se.effects,p.effects=p.effects,
                              std_effects=std_effects,se.std_effects=se.std_effects,p.std_effects=p.std_effects,
                              effects.summary=effects.summary,
                              std.effects.summary=std.effects.summary)
      }

    }else{
      if(return.boot){
        return.result <- list(effects=effects,se.effects=se.effects,p.effects=p.effects,
                              effects.summary=effects.summary,
                              effects.boot=effects_array)
      }else{
        return.result <- list(effects=effects,se.effects=se.effects,p.effects=p.effects,
                              effects.summary=effects.summary)
      }

    }

  }

  if(permutation){
    all.pars=unique(c(Est[!Est$op==":=",]$lhs,Est[!Est$op==":=",]$rhs))
    perm.data=data[,all.pars]
    independ.pars=all.pars[!all.pars%in%Est[Est$op=="~",]$lhs]
    independ.df=data[,independ.pars]
    depend.pars=all.pars[!all.pars%in%independ.pars]
    depend.df=data[,depend.pars]
    perm.list=dist.list[depend.pars]


    temp.distmat=perm.list[[1]];temp.distmat[]<-0
    temp.distmat[ind.in.dist]=1
    subset.perm=sort(unique(c(which(apply(temp.distmat,1,sum)>0),which(apply(temp.distmat,2,sum)>0))))
    ind.in.dist.new=which(temp.distmat[subset.perm,subset.perm]>0)
    if(!all(temp.distmat[subset.perm,subset.perm][ind.in.dist.new]==temp.distmat[ind.in.dist])){
      stop("ind.in.dist parameter wrong!")
    }
    if(nrow(temp.distmat)>length(subset.perm)){
      ind.in.dist=ind.in.dist.new
      for(iii in 1:length(perm.list)){
        perm.list[[iii]]=perm.list[[iii]][subset.perm,subset.perm]
      }
    }

    effects_array=std_effects_array=NULL
    judge.non.zero=len=NULL

    for(i in eff.type){
      temp.mat=MD.with.effects[[i]]$mat;temp.mat[temp.mat==0]=NA
      effects_array[[i]]=array(NA,dim = c(nrow(temp.mat),ncol(temp.mat),reps))
      judge.non.zero[[i]]=which(MD.with.effects[[i]]$mat>0)
      len[[i]]=length(MD.with.effects[[i]]$mat)
      effects[[i]]=se.effects[[i]]=p.effects[[i]]=temp.mat

      effects.summary[[i]]=temp.mat;effects.summary[[i]][]=NA

      if(standardized){
        std_effects_array[[i]]=effects_array[[i]]
        std_effects[[i]]=se.std_effects[[i]]=p.std_effects[[i]]=temp.mat
        std.effects.summary[[i]]=effects.summary[[i]]
      }
    }

    #perm.mat=sapply(1:reps,function(iii,...) {sample.int(nrow(perm.list[[1]]))})
    perm.mat=t(vegan:::getPermuteMatrix(N = nrow(perm.list[[1]]),perm = min(reps*10,100000),strata = strata))
    perm.vec=vector(length = length(perm.list));names(perm.vec)=names(perm.list)
    len_perm=ncol(perm.mat)
    re.samp.n=0
    for(reps.i in 1:reps){
      print(reps.i)
      perm.vec[]<-sample.int(len_perm,length(perm.list))
      for(k in names(perm.list)){
        perm.data[[k]]=perm.list[[k]][perm.mat[,perm.vec[k]],perm.mat[,perm.vec[k]]][ind.in.dist]
      }
      if(standardized){
        for(tempi in depend.pars){
          perm.data[,tempi]<-(perm.data[,tempi]-mean(perm.data[,tempi],na.rm=T))/sd(perm.data[,tempi],na.rm=T)
        }
      }
      fit=NULL
      fit = lavaan(model = MD.with.effects$model$model.with.effects,
                   data = perm.data, missing=missing, auto.var=TRUE, auto.fix.first=TRUE,auto.cov.lv.x=TRUE)

      while(is.null(fit)){
        re.samp.n=re.samp.n+1
        perm.vec=sample.int(nrow(perm.list[[1]]))
        for(k in names(perm.list)){
          perm.data[[k]]=perm.list[[k]][perm.vec,perm.vec][ind.in.dist]
        }
        if(standardized){
          for(tempi in depend.pars){
            perm.data[,tempi]<-(perm.data[,tempi]-mean(perm.data[,tempi],na.rm=T))/sd(perm.data[,tempi],na.rm=T)
          }
        }
        fit = lavaan(model = MD.with.effects$model$model.with.effects,
                     data = perm.data, missing=missing, auto.var=TRUE, auto.fix.first=TRUE,auto.cov.lv.x=TRUE)
      }

      Est <- parameterEstimates(fit, ci = FALSE, standardized = standardized,remove.system.eq = T,remove.nonfree=T,rsquare=T,fmi=F)
      effects.Est=Est[Est$op==":=",]
      rownames(effects.Est)=effects.Est$label

      for(i in eff.type){
        est<-effects.Est[MD.with.effects[[i]]$id[judge.non.zero[[i]]],]

        effects_array[[i]][judge.non.zero[[i]] + (reps.i-1)*len[[i]]]=est$est
        if(standardized){
          std_effects_array[[i]][judge.non.zero[[i]]+(reps.i-1)*len[[i]]]=est$std.all
        }
      }
    }
    if(re.samp.n>0){warning(paste("Permutation waring:",re.samp.n,"moedles without solutions was recalculated"))}
    for(i in eff.type){
      if(boot.perm.simultaneous) {
        effects_array[[i]]=-effects_array[[i]]+return.result$effects.boot[[i]]
      }else{
        effects_array[[i]]=-effects_array[[i]]+rep(return.result$effects[[i]],dim(effects_array[[i]])[3])
      }
      effects[[i]][]=apply(X = effects_array[[i]],MARGIN = c(1,2),FUN = mean)
      se.effects[[i]][]=apply(X = effects_array[[i]],MARGIN = c(1,2),FUN = sd)
      p.effects[[i]][judge.non.zero[[i]]]=pt(-abs(effects[[i]][judge.non.zero[[i]]])/(se.effects[[i]][judge.non.zero[[i]]]+10^-26),df = reps-1)
      effects[[i]][is.na(effects[[i]])]=0

      summary.temp=paste(char.align(effects[[i]],align = "right",digits = 2),as.symbol.pvalue(p.effects[[i]],align = "left",digits = 2))
      effects.summary[[i]][]=summary.temp[]
      effects.summary[[i]] <- as.data.frame(effects.summary[[i]])

      if(standardized){
        if(boot.perm.simultaneous) {
          std_effects_array[[i]] = -std_effects_array[[i]] + return.result$std.effetc.boot[[i]]
        }else{
          std_effects_array[[i]] = -std_effects_array[[i]]+rep(return.result$effects[[i]],dim(std_effects_array[[i]])[3])
        }
        std_effects[[i]][]=apply(X = std_effects_array[[i]],MARGIN = c(1,2),FUN = mean)
        se.std_effects[[i]][]=apply(X = std_effects_array[[i]],MARGIN = c(1,2),FUN = sd)
        p.std_effects[[i]][judge.non.zero[[i]]]=pt(-abs(std_effects[[i]][judge.non.zero[[i]]])/(se.std_effects[[i]][judge.non.zero[[i]]]+10^-26),df = reps-1)
        std_effects[[i]][is.na(std_effects[[i]])]=0

        summary.temp=paste(char.align(std_effects[[i]],align = "right",digits = 2),as.symbol.pvalue(p.std_effects[[i]],align = "left",digits = 2))
        std.effects.summary[[i]][]=summary.temp[]
        std.effects.summary[[i]] <- as.data.frame(std.effects.summary[[i]])

      }
    }

    if(standardized){
      if(return.boot){
        return.result <- list(effects=effects,se.effects=se.effects,p.effects=p.effects,
                              std_effects=std_effects,se.std_effects=se.std_effects,p.std_effects=p.std_effects,
                              effects.summary=effects.summary,
                              std.effects.summary=std.effects.summary,
                              std.effetc.boot=std_effects_array,
                              effects.boot=effects_array)
      }else{
        return.result <- list(effects=effects,se.effects=se.effects,p.effects=p.effects,
                              std_effects=std_effects,se.std_effects=se.std_effects,p.std_effects=p.std_effects,
                              effects.summary=effects.summary,
                              std.effects.summary=std.effects.summary)
      }

    }else{
      if(return.boot){
        return.result <- list(effects=effects,se.effects=se.effects,p.effects=p.effects,
                              effects.summary=effects.summary,
                              effects.boot=effects_array)
      }else{
        return.result <- list(effects=effects,se.effects=se.effects,p.effects=p.effects,
                              effects.summary=effects.summary)
      }

    }
  }

  rown=rownames(return.result$effects$total.effects)
  coln=colnames(return.result$effects$total.effects)
  for(i in c("std.effetc.boot","effects.boot")){
    if(i %in% names(return.result)){
      for(j in 1:length(return.result[[i]])){
        rownames(return.result[[i]][[j]])=rown
        colnames(return.result[[i]][[j]])=coln
      }
    }
  }



  #########################################################################################################
  if(!is.null(MD.with.effects$equal.eff)){
    equal.eff=MD.with.effects$equal.eff
    warning(paste(c("The effects of <",names(equal.eff),"> is equal to the effects of <",unlist(equal.eff),">"),collapse = " "))

    temp.array=array(NA,dim = c(length(rown),length(coln)+length(equal.eff),reps),dimnames = list(rown,c(coln,names(equal.eff)),1:reps))
    temp.mat=temp.array[,,1]
    cycle=data=NULL
    loop.depth=1
    for(cycle_1 in 1:length(return.result)){
      data[[loop.depth]]=return.result[[cycle_1]]


      if(!(is.array(data[[loop.depth]])|is.matrix(data[[loop.depth]])|is.data.frame(data[[loop.depth]]))){
        loop.depth=loop.depth+1
        for(cycle_2 in 1:length(data[[loop.depth-1]])){
          data[[loop.depth]]=data[[loop.depth-1]][[cycle_2]]
          if(!(is.array(data[[loop.depth]])|is.matrix(data[[loop.depth]])|is.data.frame(data[[loop.depth]]))){
            loop.depth=loop.depth+1
            warning("debug ... loop.depth low")
          }else{
            if(length(dim(return.result[[cycle_1]][[cycle_2]]))==3){
              temp.array[,coln,]=return.result[[cycle_1]][[cycle_2]]
              temp.array[,names(equal.eff),]=return.result[[cycle_1]][[cycle_2]][,unlist(equal.eff),]
              return.result[[cycle_1]][[cycle_2]]=temp.array
            }else if(length(dim(return.result[[cycle_1]][[cycle_2]]))==2){
              temp.mat[,coln]=return.result[[cycle_1]][[cycle_2]]
              temp.mat[,names(equal.eff)]=return.result[[cycle_1]][[cycle_2]][,unlist(equal.eff)]
              return.result[[cycle_1]][[cycle_2]]=temp.mat
            }

          }

        }
        loop.depth=loop.depth-1

      }else{
        if(length(dim(return.result[[cycle_1]]))==3){
          temp.array[,coln,]=return.result[[cycle_1]]
          temp.array[,names(equal.eff),]=return.result[[cycle_1]][,unlist(equal.eff),]
          return.result[[cycle_1]]=temp.array
        }else if(length(dim(return.result[[cycle_1]]))==2){
          temp.mat[,coln]=return.result[[cycle_1]]
          temp.mat[,names(equal.eff)]=return.result[[cycle_1]][,unlist(equal.eff)]
          return.result[[cycle_1]]=temp.mat
        }
      }
    }
  }

  #gc()
  return(return.result)

}


"SEM.total.effect.expression" <- function(rh=NULL,lh=NULL,edges){
  #total.effect.expression(rh="D",lh="Multi.stibility",edges)
  fun1<-function(x,edges){
    "x=sortby(x,-(1:length(x)))";
    if(length(x)<2&length(unique(x))<length(x)){
      return(NULL)
    }else if(length(x)==2){
      return(edges[paste0(x[2:1],collapse = "~"),]$label)
    }else if(length(x)>2){
      temp=vector(length = length(x)-1)
      for(k in 2:length(x)){
        temp[k-1]=paste0(x[k:(k-1)],collapse ="~")
      }
      return(paste0(edges[temp,]$label,collapse="*"))
    }
  }

  if(sum(edges$lhs==lh)==0|sum(edges$rhs==rh)==0){
    finish=NULL
  }else{
    rh.temp=rh
    finish=NULL
    unfinish.member=rh
    n=0
    while (length(rh.temp)>0&n<(length(unique(edges$lhs)))){
      n=n+1
      unfinish.rh=NULL
      unfinish=NULL
      for(rhi in unique(rh.temp)){
        times.rhi=rh.temp==rhi

        lh.temp=edges$lhs[edges$rhs==rhi]
        if(length(lh.temp)>0){
          lh.temp=rep(lh.temp,each=sum(times.rhi))
          find.lh=lh.temp==lh
          #find.lh=rep(find.lh,times=sum(times.rhi))
          temp1=paste(unfinish.member[times.rhi],lh.temp,sep = "~")
          rh.new=if(sum(!find.lh)){lh.temp[!find.lh]}else{NULL}
          if(sum(find.lh)){
            finish=c(finish,temp1[find.lh]);
            temp1=temp1[!find.lh]
          }
          unfinish=c(unfinish,temp1)
          unfinish.rh=c(unfinish.rh,rh.new)
        }else{
          #print(rhi)
          next
        }

      }
      rh.temp=unfinish.rh
      unfinish.member=unfinish
    }

  }



  if(!is.null(finish)){
    finish.split=strsplit(x =finish,split = "~" )

    total.effect.elem=sapply(1:length(finish.split),function(i,...) fun1(finish.split[[i]],edges))

    direct_edge=paste(lh,rh,sep="~")==rownames(edges)
    direct.effect.elem=if(sum(direct_edge)){edges[direct_edge,]$label}else{NULL}
    indirect.effect.elem=
      if(is.null(direct.effect.elem)){
        total.effect.elem
      }else{
        if(length(total.effect.elem)>1){total.effect.elem[!total.effect.elem==direct.effect.elem]}else{NULL}
      }

  }else{
    total.effect.elem=direct.effect.elem=indirect.effect.elem=NULL
  }


  return(list(total.effect.elem=total.effect.elem,
              direct.effect.elem=direct.effect.elem,
              indirect.effect.elem=indirect.effect.elem
  ))
}


"add.variance" <- function(model=NULL,variance=NULL,covariance=NULL,add.paths=NULL){
  if(!is.null(model)){
    temp=strsplit(model,split = "[\n;]");
    md=NULL
    for(i in 1:length(temp)){
      md=c(md,temp[[i]])
    }
    md=gsub("[ ]*","",md)
    md=md[!grepl("^#",md)&!md==""]
    all.pars=unique(unlist(strsplit(md,"[~=:\\+]")));all.pars=all.pars[!all.pars==""]
    model.list=list(variance=NULL,covariance=NULL,effect=NULL,lv=NULL,def.v=NULL)
    if(sum(grepl("=~",md))){model.list$lv=md[grepl("=~",md)];md=md[!grepl("=~",md)]}

    if(sum(grepl(":=",md))){model.list$def.v=md[grepl(":=",md)];  md=md[!grepl(":=",md)]}

    if(sum(grepl("~~",md))){
      temp=md[grepl("~~",md)];
      md=md[!grepl("~~",md)]
      temp1=strsplit(temp,"~~")
      for(i in 1:length(temp1)){
        tempj=strsplit(temp1[[i]],"\\+")
        for(j1 in tempj[[1]]){
          for(j2 in tempj[[2]]){
            if(j1==j2){
              model.list$variance=c(model.list$variance,paste0(j1,"~~",j2))
            }else{
              model.list$covariance=c(model.list$covariance,paste0(j1,"~~",j2))
            }
          }
        }
      }
    }

    md=md[grepl("~",md)]
    temp1=strsplit(md,"~")
    if(length(temp1)>0){
      for(i in 1:length(temp1)){
        tempj=strsplit(temp1[[i]],"\\+")
        for(j1 in tempj[[1]]){
          for(j2 in tempj[[2]]){
            model.list$effect=c(model.list$effect,paste0(j1,"~",j2))
          }
        }
      }
    }

  }
  if(is.null(variance)){variance=all.pars}
  variance=unique(c(variance,all.pars))
  model.list$variance=unique(c(model.list$variance,paste0(variance,"~~",variance)))
  tab=t(as.data.frame(strsplit(x = unlist(model.list),split = "~[~]*")))

  if(length(covariance)>0){
    add.cov=NULL
    covariance=unique(covariance)
    for(j1 in 1:length(covariance)){
      for (j2 in j1:length(covariance)){
        if(j1==j2){
          model.list$variance=c(model.list$variance,paste0(covariance[j1],"~~",covariance[j2]))
        }else{
          add.cov=c(add.cov,paste0(covariance[j1],"~~",covariance[j2]))
        }
      }
    }
    add.cov=add.cov[!gsub("~","",add.cov)%in%c(paste0(tab[,1],tab[,2]),paste0(tab[,2],tab[,1]))]
    model.list$covariance=sort(unique(c(model.list$covariance,add.cov)))
    model.list$variance=unique(sort(model.list$variance))
  }

  for(i in 1:length(model.list)){
    if(!is.null(model.list[[i]])){model.list[[i]]=sort(unique(model.list[[i]]))}
  }
  if(length(model.list$covariance)){
    cov.name <- model.list$covariance%>%strsplit(.,"~~")%>%strsplit2mat
    if(length(cov.name==2)){
      cov.name <- rbind(cov.name,cov.name[2:1])%>%apply(.,1,paste0,collapse="~~")
    }else{
      cov.name <- rbind(cov.name,cov.name[,2:1])%>%apply(.,1,paste0,collapse="~~")
    }
  }else{cov.name <- NULL}

  model.u=unique(unlist(model.list));


  if(!is.null(add.paths)){
    add.paths <- add.variance(add.paths)%>%strsplit(.,";")%>%unlist
    add.paths <- add.paths[!add.paths%in%cov.name]
    #add.paths <- add.paths[!add.paths==""]
    model.u <- unique(c(model.u,add.paths[!add.paths%in%model.u]))
  }
  model=paste0(model.u[!model.u==""],collapse = ";")
  return(model)
}


"effects.sem" <- function(rh=NULL,lh=NULL,fit,model=NULL){
  #effects.model(rh=c("D","I","PR"),lh=c("Multi.stibility","M_Wellbeing"),edges,model=MD)
  par.tab_model=extract.sem.model(fit)
  par.tab=par.tab_model$model.table
  edges=par.tab[par.tab$op=="~",]

  if(is.null(rh)){
    if(length(unique(edges$lhs))==1){
      rh=unique(edges$rhs)
    }else{
      rh=rownames(sortby(tabulate.char(edges$rhs),sort_vector = 1,decreasing = T))
    }

  }
  if(is.null(lh)){
    if(length(unique(edges$lhs))==1){
      lh=unique(edges$lhs)
    }else{
      lh=rownames(sortby(tabulate.char(edges$lhs),1))
    }

  }

  j1=rh%in%lh
  if(sum(j1)>0&sum(j1)<length(rh)){rh=c(rh[!j1],rh[j1])}
  j2=lh%in%rh
  if(sum(j2)>0&sum(j2)<length(lh)){lh=c(lh[j2],lh[!j2])}




  if(is.null(model)){
    model=par.tab_model$model
  }

  mat=matrix(NA,length(rh),length(lh),dimnames = list(rh,lh))
  total.eff=matrix("",length(rh),length(lh),dimnames = list(rh,lh))
  direct.eff=indirect.eff=total.eff
  total.eff.id=indirect.eff.id=direct.eff.id=total.eff
  for(i in 1:length(rh)){
    #if(rh_i==rh[3]&lh_j==lh[1]){break}
    rh_i=rh[i]

    for(j in 1:length(lh)){
      lh_j=lh[j]

      effects=SEM.total.effect.expression(rh=rh_i,lh=lh_j,edges)
      total.eff.id[rh_i,lh_j]=paste("total_rh",i,"lh",j,sep="")
      total.eff[rh_i,lh_j]=if(!is.null(effects$total.effect.elem)){paste(total.eff.id[rh_i,lh_j],":=",paste0(effects$total.effect.elem,collapse ="+"),sep="")}else{""}

      indirect.eff.id[rh_i,lh_j]=paste("indirect_rh",i,"lh",j,sep="")
      indirect.eff[rh_i,lh_j]=if(!is.null(effects$indirect.effect.elem)){paste(indirect.eff.id[rh_i,lh_j],":=",paste0(effects$indirect.effect.elem,collapse ="+"),sep="")}else{""}

      direct.eff.id[rh_i,lh_j]=paste("direct_rh",i,"lh",j,sep="")
      direct.eff[rh_i,lh_j]=if(!is.null(effects$direct.effect.elem)){paste(direct.eff.id[rh_i,lh_j],":=",paste0(effects$direct.effect.elem,collapse ="+"),sep="")}else{""}
    }
  }
  all.effect=c(total.eff[],indirect.eff[],direct.eff[])
  all.effect=all.effect[!all.effect==""]
  effects=paste0(all.effect,collapse =";")

  model.with.effects=if(is.null(model)){NULL}else{paste(model,effects,sep=";")}

  mat[]=NA;mat[total.eff==""]=0;   mat[total.eff!=""]=1;   total.eff=mat
  mat[]=NA;mat[direct.eff==""]=0;  mat[direct.eff!=""]=1;  direct.eff=mat
  mat[]=NA;mat[indirect.eff==""]=0;mat[indirect.eff!=""]=1;indirect.eff=mat
  return(
    list(model=list(model.with.effects=model.with.effects, effects.defined=effects),
         total.effects=list(mat=t(total.eff),id=t(total.eff.id)),
         direct.effects=list(mat=t(direct.eff),id=t(direct.eff.id)),
         indirect.effects=list(mat=t(indirect.eff),id=t(indirect.eff.id))
    )
  )
}



"sem.effects.stat" <- function(fit,data,effect.type=c("total.effects", "direct.effects", "indirect.effects"),
                               rh=NULL,lh=NULL,standardized = TRUE,
                               bootstrap=FALSE,reps=1000,return.boot=TRUE,missing="listwise",boot.depth=NULL){
  MD.with.effects=effects.sem(rh,lh,fit,model=NULL)
  fit = lavaan(model = MD.with.effects$model$model.with.effects,
               data = data, missing=missing, auto.var=TRUE, auto.fix.first=TRUE,auto.cov.lv.x=TRUE)
  Est <- parameterEstimates(fit, ci = FALSE, standardized = standardized,remove.system.eq = T,remove.nonfree=T,rsquare=T,fmi=F)
  effects.Est=Est[Est$op==":=",]
  rownames(effects.Est)=effects.Est$label


  effects=se.effects=p.effects=std_effects=se.std_effects=p.std_effects=judge.non.zero=NULL
  effects.summary=std.effects.summary=NULL
  if(!bootstrap){
    for(i in effect.type){
      temp.mat=MD.with.effects[[i]]$mat
      judge.non.zero[[i]]=which(temp.mat>0)
      temp.mat[]=NA
      effects[[i]]=se.effects[[i]]=p.effects[[i]]=temp.mat


      est=effects.Est[MD.with.effects[[i]]$id[judge.non.zero[[i]]],]

      effects[[i]][judge.non.zero[[i]]]=est$est
      se.effects[[i]][judge.non.zero[[i]]]=est$se
      p.effects[[i]][judge.non.zero[[i]]]=est$pvalue
      effects[[i]][is.na(effects[[i]])]=0

      effects.summary[[i]]=temp.mat
      summary.temp=paste(char.align(effects[[i]],align = "right",digits = 2),as.symbol.pvalue(p.effects[[i]],align = "left",digits = 2))
      effects.summary[[i]][]=summary.temp[]
      effects.summary[[i]] <- as.data.frame(effects.summary[[i]])

      if(standardized){
        std_effects[[i]]=se.std_effects[[i]]=p.std_effects[[i]]=temp.mat
        std_effects[[i]][judge.non.zero[[i]]]=est$std.all
        std_effects[[i]][is.na(std_effects[[i]])]=0

        std.effects.summary[[i]]=temp.mat
        summary.temp=paste(char.align(std_effects[[i]],align = "right",digits = 2),as.symbol.pvalue(p.effects[[i]],align = "left",digits = 2))
        std.effects.summary[[i]][]=summary.temp[]
        std.effects.summary[[i]] <- as.data.frame(std.effects.summary[[i]])
      }
    }


    if(standardized){
      return(list(effects=effects,se.effects=se.effects,p.effects=p.effects,
                  std_effects=std_effects,
                  effects.summary=effects.summary,
                  std.effects.summary=std.effects.summary))
    }else{
      return(list(effects=effects,se.effects=se.effects,p.effects=p.effects,
                  effects.summary=effects.summary))
    }
    # return(list(effects=effects,se.effects=se.effects,p.effects=p.effects,std_effects=std_effects))
  }else{

    all.pars=unique(c(Est[!Est$op==":=",]$lhs,Est[!Est$op==":=",]$rhs))
    data.filter=data[,all.pars]
    effects_array=std_effects_array=NULL
    judge.non.zero=len=NULL

    for(i in effect.type){
      temp.mat=MD.with.effects[[i]]$mat;temp.mat[temp.mat==0]=NA
      effects_array[[i]]=array(NA,dim = c(nrow(temp.mat),ncol(temp.mat),reps))
      judge.non.zero[[i]]=which(MD.with.effects[[i]]$mat>0)
      len[[i]]=length(MD.with.effects[[i]]$mat)
      effects[[i]]=se.effects[[i]]=p.effects[[i]]=temp.mat

      effects.summary[[i]]=temp.mat;effects.summary[[i]][]=NA

      if(standardized){
        std_effects_array[[i]]=effects_array[[i]]
        std_effects[[i]]=se.std_effects[[i]]=p.std_effects[[i]]=temp.mat
        std.effects.summary[[i]]=effects.summary[[i]]
      }
    }

    for(reps.i in 1:reps){
      boot.data<-if(is.null(boot.depth)){
        data.filter[sample.int(nrow(data.filter),replace = T),]
      }else{
        data.filter[sample.int(n = nrow(data.filter),size = boot.depth,replace = T),]
      }
      fit = lavaan(model = MD.with.effects$model$model.with.effects,
                   data = boot.data, missing=missing, auto.var=TRUE, auto.fix.first=TRUE,auto.cov.lv.x=TRUE)
      Est <- parameterEstimates(fit, ci = FALSE, standardized = standardized,remove.system.eq = T,remove.nonfree=T,rsquare=T,fmi=F)
      effects.Est=Est[Est$op==":=",]
      rownames(effects.Est)=effects.Est$label

      for(i in effect.type){
        est<-effects.Est[MD.with.effects[[i]]$id[judge.non.zero[[i]]],]

        effects_array[[i]][judge.non.zero[[i]] + (reps.i-1)*len[[i]]]=est$est
        if(standardized){
          std_effects_array[[i]][judge.non.zero[[i]]+(reps.i-1)*len[[i]]]=est$std.all
        }
      }
    }

    for(i in effect.type){
      effects[[i]][]=apply(X = effects_array[[i]],MARGIN = c(1,2),FUN = mean)
      se.effects[[i]][]=apply(X = effects_array[[i]],MARGIN = c(1,2),FUN = sd)
      p.effects[[i]][judge.non.zero[[i]]]=pt(-abs(effects[[i]][judge.non.zero[[i]]])/(se.effects[[i]][judge.non.zero[[i]]]+10^-26),df = reps-1)
      effects[[i]][is.na(effects[[i]])]=0

      summary.temp=paste(char.align(effects[[i]],align = "right",digits = 2),as.symbol.pvalue(p.effects[[i]],align = "left",digits = 2))
      effects.summary[[i]][]=summary.temp[]
      effects.summary[[i]] <- as.data.frame(effects.summary[[i]])

      if(standardized){
        std_effects[[i]][]=apply(X = std_effects_array[[i]],MARGIN = c(1,2),FUN = mean)
        se.std_effects[[i]][]=apply(X = std_effects_array[[i]],MARGIN = c(1,2),FUN = sd)
        p.std_effects[[i]][judge.non.zero[[i]]]=pt(-abs(std_effects[[i]][judge.non.zero[[i]]])/(se.std_effects[[i]][judge.non.zero[[i]]]+10^-26),df = reps-1)
        std_effects[[i]][is.na(std_effects[[i]])]=0

        summary.temp=paste(char.align(std_effects[[i]],align = "right",digits = 2),as.symbol.pvalue(p.std_effects[[i]],align = "left",digits = 2))
        std.effects.summary[[i]][]=summary.temp[]
        std.effects.summary[[i]] <- as.data.frame(std.effects.summary[[i]])
      }
    }

    if(standardized){
      if(return.boot){
        return(list(effects=effects,se.effects=se.effects,p.effects=p.effects,
                    std_effects=std_effects,se.std_effects=se.std_effects,p.std_effects=p.std_effects,
                    effects.summary=effects.summary,
                    std.effects.summary=std.effects.summary,
                    std.effetc.boot=std_effects_array,
                    effects.boot=effects_array))
      }else{
        return(list(effects=effects,se.effects=se.effects,p.effects=p.effects,
                    std_effects=std_effects,se.std_effects=se.std_effects,p.std_effects=p.std_effects,
                    effects.summary=effects.summary,
                    std.effects.summary=std.effects.summary))
      }

    }else{
      if(return.boot){
        return(list(effects=effects,se.effects=se.effects,p.effects=p.effects,
                    effects.summary=effects.summary,
                    effects.boot=effects_array))
      }else{
        return(list(effects=effects,se.effects=se.effects,p.effects=p.effects,
                    effects.summary=effects.summary))
      }

    }

  }


}

"extract.sem.model" <- function(fit=NULL,export.model.table=TRUE){
  par.tab <- parTable(fit)
  no.lable=par.tab$label==""
  if(sum(no.lable)){par.tab$label[no.lable]=paste("par",1:sum(no.lable),sep="")}
  model.with.self_def=paste0(sort(paste(par.tab$lhs,par.tab$op,par.tab$label,"*",par.tab$rhs,sep="")),collapse ="; ")
  par.tab=par.tab[par.tab$op%in%c("~","~~","=~"),]
  model=paste0(sort(paste(par.tab$lhs,par.tab$op,par.tab$label,"*",par.tab$rhs,sep="")),collapse ="; ")
  model.nolable=paste0(sort(paste(par.tab$lhs,par.tab$op,par.tab$rhs,sep="")),collapse ="; ")
  if(export.model.table){

    model.table=par.tab[,c("lhs","op","rhs","label")]
    model.table$experession=paste(par.tab$lhs,par.tab$op,par.tab$label,"*",par.tab$rhs,sep="")
    model.table$experession.nolable=paste(par.tab$lhs,par.tab$op,par.tab$rhs,sep="")
    rownames(model.table)=model.table$experession.nolable
    return(list(model=model,model.with.self_def=model.with.self_def,model.table=model.table,model.nolable=model.nolable))
  }else{
    return(list(model=model,model.with.self_def=model.with.self_def,model.nolable=model.nolable))
  }
}


"effect.plot" <- function(EFF,sort_elem = NULL,mycolor=NULL,
                          sample.n.dots = 5, WEIGHT = FALSE, weight=NULL,INDEPEND.pars=NULL, DEPEND.pars=NULL,legend.position="none"){

  if(sum(names(EFF)%in%c("effects","se.effects","p.effects","std_effects","effects.summary"))==5){
    EFF=list(g0=EFF);
    if(!is.null(weight)) weight=list(g0=weight)
  }
  uniq_levels=INDEPEND=DEPEND=NULL
  for(ii in 1:length(EFF)){
    uniq_levels = unique(c(uniq_levels,names(EFF[[ii]][["effects"]])))
    INDEPEND = unique(c(INDEPEND,colnames(EFF[[ii]][["effects"]][["direct.effects"]])))
    DEPEND =  unique(c(DEPEND,rownames(EFF[[ii]][["effects"]][["direct.effects"]])))
  }
  if(!is.null(DEPEND.pars)){DEPEND<-DEPEND[DEPEND%in%DEPEND.pars]}
  if(!is.null(INDEPEND.pars)){INDEPEND<-INDEPEND[INDEPEND%in%INDEPEND.pars]}
  uniq_levels=uniq_levels[!uniq_levels%in%c("total.effects","indirect.effects","direct.effects")]
  uniq_levels=c("direct.effects",uniq_levels)
  if(!is.null(sort_elem)){
    sort_elem=sort_elem[sort_elem%in%uniq_levels]
    if(length(sort_elem)) uniq_levels[uniq_levels%in%sort_elem] = sort_elem
  }
  if(is.null(mycolor)){
    mycolor=rainbow(length(uniq_levels));names(mycolor)=uniq_levels
  }else{
    nm=names(mycolor);nm<-nm[nm%in%uniq_levels]
    if(length(nm)){
      if(length(nm) < length(uniq_levels)){
        mycolor <- c(mycolor[nm],rainbow(length(uniq_levels)-length(nm)))
        names(mycolor) <- c(nm,uniq_levels[!uniq_levels%in%nm])
      }
    }else if(length(mycolor)==length(uniq_levels)){
      names(mycolor) <- uniq_levels
    }
  }
  data.boot=databarplot=NULL
  lvls=c("total.effects","indirect.effects",uniq_levels)
  for(iii in names(EFF)){
    eff=EFF[[iii]]
    #eff$std.effects.summary
    #eff$std_effects

    rown=rownames(eff$std_effects$total.effects)
    coln=colnames(eff$std_effects$total.effects)



    for(independ in INDEPEND){
      for(depend in DEPEND){
        independ.i=which(coln==independ)
        depend.i=which(rown==depend)
        for(effect.i in lvls){

          temp0=eff$std.effetc.boot[[effect.i]][depend.i,independ.i,]
          if(!is.null(temp0)){
            tempbar=data.frame(effect.type=effect.i,
                               effects=eff$std_effects[[effect.i]][depend,independ],
                               p.value=eff$p.std_effects[[effect.i]][depend,independ],
                               se=eff$se.std_effects[[effect.i]][depend,independ],
                               ci=eff$se.std_effects[[effect.i]][depend,independ]*1.96,
                               group=paste0(depend,"~",independ),
                               depend=depend,
                               independ=independ,independ.SEMi=paste(independ,iii))
            temp=data.frame(effect.type=effect.i,
                            group=paste0(depend,"~",independ),
                            effects=eff$std.effetc.boot[[effect.i]][depend.i,independ.i,1:sample.n.dots],
                            group=paste0(depend,"~",independ),
                            depend=depend,
                            independ=independ,independ.SEMi=paste(independ,iii))

          }else{tempbar=temp=NULL}


          databarplot=rbind(databarplot,tempbar)


          data.boot=rbind(data.boot,temp)
        }
      }
    }
  }

  databarplot=databarplot[!is.na(databarplot$p.value),]
  data.boot=data.boot[!is.na(data.boot$effects),]
  if(WEIGHT){
    if(!all(databarplot$independ%in%names(weight[[iii]]))){stop(paste("weights of (",paste(INDEPEND,collapse = "/"),") must be supplied with 'WEIGHT=TRUE' argement"))}
    for(i in 1:nrow(databarplot)){
      databarplot[i,c("effects","se","ci")] = databarplot[i,c("effects","se","ci")]*weight[[iii]][databarplot$depend[i]]/weight[[iii]][databarplot$independ[i]]
    }

    for(i in 1:nrow(databarplot)){
      data.boot$effects[i] = data.boot$effects[i]*weight[[iii]][data.boot$depend[i]]/weight[[iii]][data.boot$independ[i]]
    }
  }


  #################################################################################################### barplot for effects on communities
  #mycolor=c("#FF6A6A","#FF00FF","#90EE90","#2E8B57","#00BFFF", "#6495ED","#EEC900","#CD950C")
  #mycolor=c("#FF00FF", "#FF6A6A", "#CD950C", "#EEC900", "gray")
  #mycolor
  plot_bar=NULL

  #temp.data=data.boot[data.boot$depend%in%DEPEND,];temp.data$color=temp.data$depend

  #dataforjiter=temp.data

  temp.data.bar = databarplot
  log.func = databarplot$depend%in%DEPEND
  data.bar.comm = databarplot[log.func,]

  data.bar.comm$depend = factor(data.bar.comm$depend,levels = DEPEND)
  #data.bar.comm$depend=paste(as.numeric(data.bar.comm$depend),data.bar.comm$depend)

  data.bar.comm$effect.type=factor(data.bar.comm$effect.type,levels = lvls[length(lvls):1])
  #temp.slop.p=slope.p[slope.p$group==i&slope.p$effect.type!="direct.effects",]
  for(j in INDEPEND){
    judge_total=data.bar.comm$independ==j&data.bar.comm$effect.type%in%c("total.effects")
    judge_indirect=data.bar.comm$independ==j&data.bar.comm$effect.type%in%c("indirect.effects")
    data1 = data.bar.comm[data.bar.comm$independ==j&!(judge_total|judge_indirect),]
    data_total = data.bar.comm[judge_total,]
    max.gap=max(c(0,data.bar.comm[judge_total,]$effects+data.bar.comm[judge_total,]$ci),na.rm=T)-min(c(data.bar.comm[judge_total,]$effects-data.bar.comm[judge_total,]$ci,0),na.rm=T)
    bk.y=max(c(round(max.gap/4,digits = 1),0.02), na.rm=T)


    plot_bar[[j]] = ggplot(data.bar.comm[data.bar.comm$independ==j,],aes(x = depend,y = effects))+
      facet_wrap(~paste(independ.SEMi),nrow=1)+
      #geom_smooth(data=temp.data,aes(x = PR,y = effects),method = "lm",se=F)+
      #geom_jitter(data = dataforjiter[dataforjiter$independ==j,],na.rm = T,
      #aes(x = depend,y = effects,col=color),height=0,width = 0.2,size=0.1,shape=16,alpha=0.8)+
      geom_bar(aes(x = depend,y = effects,fill=effect.type),stat='identity',position='stack',color=NA,size=0.25,alpha=1,width = 0.65,
               data = data1)+
      geom_bar(aes(x = depend,y = effects),stat='identity',position='stack',fill=NA,color="black",size=0.25,width = 0.65,
               data = data_total)+
      geom_errorbar(aes(ymin = effects-se, ymax = effects + se),width = 0.5,size=0.375,
                    data = data_total)+
      geom_text( aes(x = depend, y=effects + sign(effects)*(se+0.1*max.gap),
                     label=as.symbol.pvalue(p.value,sig.level=c(0.001,0.01,0.05))),size=4,
                 data = data_total)+
      scale_color_manual(values = mycolor[levels(data1$effect.type)[levels(data1$effect.type)%in%unique(data1$effect.type)]])+
      scale_fill_manual(values =  mycolor[levels(data1$effect.type)[levels(data1$effect.type)%in%unique(data1$effect.type)]])+
      labs(x=element_blank(), y = paste0(j," effects"))+
      scale_y_continuous(breaks = sort(c(-seq(bk.y,5,bk.y),seq(0,5,bk.y))))+
      theme_bw()+
      theme(text=element_text(size = 8, family = "sans", face = "plain"),axis.title.x = element_blank(),
            panel.grid.major.x = element_blank(),    panel.grid.minor.x =  element_blank(),  panel.grid.minor.y =  element_blank(),
            panel.background = element_rect(color = 'black', fill = 'transparent'),legend.position = legend.position)
  }

return(list(plot_bar=plot_bar,data.bar=data.bar.comm))
}




