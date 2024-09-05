###################################
### REQUIRED PACKAGES #############
###################################
.libPaths("~/data/volume_2")
library("ggplot2")
library("purrr")
library("matrixStats")
library("foreach")
library("ggpmisc")


###################################
### DEFINITIONS / CONSTANTS #######
###################################
# Metabolite markers selected to compare
MET65 <- tolower(c("Ala","Gln","His","Phe","Tyr","Ile","Leu","Val","Glc","Lac","Pyr","Cit","Ace",
                   "AcAce","bOHBut","Crea","Alb","Gp","XXL_VLDL_L","XL_VLDL_L","L_VLDL_L","M_VLDL_L",
                   "S_VLDL_L","XS_VLDL_L","IDL_L","L_LDL_L","M_LDL_L","S_LDL_L","XL_HDL_L","L_HDL_L",
                   "M_HDL_L","S_HDL_L","IDL_C","Serum_C","VLDL_C","LDL_C","HDL_C",
                   "VLDL_D","LDL_D","HDL_D","Serum_TG","TotPG","PC","SM","TotCho","ApoA1","ApoB",
                   "TotFA","DHA","LA","FAw3","FAw6","PUFA","MUFA","SFA","FAw3_FA","FAw6_FA","PUFA_FA",
                   "MUFA_FA","SFA_FA","UnsatDeg","xl_hdl_c","l_hdl_c","m_hdl_c","s_hdl_c"
))

# DEFINITIONS:
mort_betas <- data.frame(
  Abbreviation=c("pufa_fa","gp","glc","s_hdl_l","xxl_vldl_l","alb","phe","acace","ile","vldl_d",
                 "leu","val","his","lac"),
  Metabolite=c("Ratio of polyunsaturated fatty acids to total fatty acids (%)",
               "Glycoprotein acetyls",
               "Glucos",
               "Total lipids in small HD",
               "Total lipids in chylomicrons and extremely large VLDL",
               "Albumin",
               "Phenylalanine",
               "Acetoacetate",
               "Isoleucine",
               "Mean diameter for VLDL particles",
               "Leucine",
               "Valine",
               "Histidine",
               "Lactate"),
  Beta_value=c(-0.251264056,0.280179378,0.148392648,-0.141912413,-0.223685945,-0.113290967,
               0.124227231,0.078704926,0.205706859,-0.160137595,-0.195741112,-0.143698197,
               -0.069223326,0.062262328),
  stringsAsFactors = FALSE)

###########
#Functions#
###########
## Defining a function to report on missingness given counts:
print.miss.report <- function(dat,on_sample=TRUE,type="missingOrZero"){
  rep.miss <- function(counts,Nmeas,Nsamp,on_sample){
    A <- table(counts)
    B <- data.frame(as.numeric(names(A)),
                    as.numeric(names(A))/ifelse(on_sample,Nmeas,Nsamp)*100,
                    as.vector(A),
                    as.vector(A)/ifelse(on_sample,Nsamp,Nmeas)*100,stringsAsFactors=FALSE)
    if(on_sample){
      colnames(B) <- c("# missMeas","[%]","# samp","[%]")
    } else {
      colnames(B) <- c("# missSamp","[%]","# meas","[%]")
    }
    return(B)
  }
  COUNTS <- count.miss(dat,on_sample=on_sample)
  if(type=="missingOrZero"){
    res <- rep.miss(COUNTS[["miss_or_zero"]],Nmeas=COUNTS$Nmeas,Nsamp=COUNTS$Nsamp,on_sample=on_sample)
  }
  if(type=="missing"){
    res <- rep.miss(COUNTS[["miss"]],Nmeas=COUNTS$Nmeas,Nsamp=COUNTS$Nsamp,on_sample=on_sample)
  }
  if(type=="zero"){
    res <- rep.miss(COUNTS[["zero"]],Nmeas=COUNTS$Nmeas,Nsamp=COUNTS$Nsamp,on_sample=on_sample)
  }
  pander:::pander(res,justify = c('left',rep('right',3)),caption=paste0(type," in ",ifelse(on_sample,"samples","measurements")),digits=3)
}

## Defining a function for counting missing values:
count.miss <- function(dat,on_sample=TRUE){
  #IND <- grep("abnormal_macromolecule_a|low_glucose|low_glutamine_high_glutamate|low_protein_content|high_citrate|high_ethanol|high_lactate|high_pyruvate|serum_sample|unidentified_small_molecule_a|unidentified_small_molecule_b|unknown_acetylated_compound|isopropyl_alcohol|polysaccharides|aminocaproic_acid|fast",rownames(dat))
  #dat <- dat[-IND,]
  Dummi <- matrix(NA,ncol=ncol(dat),nrow=nrow(dat))
  Dummi[which(dat==0)] <- 0
  Dummi[which(dat!=0)] <- 1
  if(on_sample){
    COUNTS <- list(miss=colSums(is.na(Dummi)),zero=colSums(Dummi==0,na.rm=TRUE),
                   miss_or_zero=nrow(Dummi)-colSums(Dummi==1,na.rm=TRUE),Nmeas=nrow(Dummi),Nsamp=ncol(Dummi))
  } else {
    COUNTS <- list(miss=rowSums(is.na(Dummi)),zero=rowSums(Dummi==0,na.rm=TRUE),
                   miss_or_zero=ncol(Dummi)-rowSums(Dummi==1,na.rm=TRUE),Nmeas=nrow(Dummi),Nsamp=ncol(Dummi))
  }
  return(COUNTS)
}

### Defining a function for plotting a heatmap indcating missing & zero values:
plot.na.heatmap  <- function(dat, title="Missingness Plot"){
  Dummi <- matrix(NA,ncol=ncol(dat),nrow=nrow(dat))
  Dummi[which(!is.na(dat==0))] <- 1
  #Dummi[which(dat!=0)] <- 1
  layout(mat=matrix(c(1,2,3,4),ncol=2,byrow=TRUE),widths=c(4,1),heights=c(1,4))
  par(xaxs = "i")  # THIS IS EVIL MAGIC!
  par(yaxs = "i")
  par(mar=c(0,3,1,0))
  XCOUNT <- nrow(Dummi)-colSums(Dummi,na.rm=TRUE)
  XPERC  <- XCOUNT/nrow(Dummi)*100
  YCOUNT <- ncol(Dummi)-rowSums(Dummi,na.rm=TRUE)
  YPERC  <- YCOUNT/ncol(Dummi)*100
  if(max(XPERC)<10){
    ylim <- c(0,10)
  } else {
    ylim <- c(0,max(pretty(XPERC)))
  }
  par(mar=c(0.5,13,2.5,0))
  barplot(XPERC,axes=TRUE,main=title,col="lightblue",border="lightblue",
          las=2,cex.main=3,ylim=ylim,cex.axis=1.2,font.axis=1.5)
  par(mar=c(0,0,0,0))
  plot.new()
  plot.window(xlim=c(0,1),ylim=c(0,1))
  legend("topright",fill=c("white","grey30"),legend=c("missing","value"))
  par(mar=c(2,13,0,0))
  image(t(Dummi[nrow(Dummi):1,]),axes=F, col=c("grey30"),xlim=c(0,1),ylim=c(0,1))
  mtext(text=rownames(dat), side=2, line=0.3, at=c(0.985,seq(0.945,0.055, length=dim(dat)[1]-2),0.03),
        las=1, cex=0.8)
  mtext(text="samples",side=1,font=1.5,cex=1.5, line = 0.6)
  par(mar=c(2,0,0,0))
  if(max(YPERC)<10){
    xlim <- c(0,10)
  } else {
    xlim <- c(0,max(pretty(YPERC)))
  }
  barplot(rev(YPERC),horiz=TRUE,axes=TRUE,col="lightblue",border="lightblue",cex.axis=1,font.axis=2)
}

#Function to retrieve the same metabolites names  available in BBMRI-nl
find_BBMRI_names<-function(names, quiet=F){
  new_names <- names %>% purrr::map_chr(function(id) {
    # Look through the alternative_ids
    hits <-
      purrr::map_lgl(
        BBMRI_translator$alternative_names,
        ~ id %in% .
      )
    
    # If one unambiguous hit, return it.
    if (sum(hits) == 1L) {
      return(BBMRI_translator$BBMRI_names[hits])
      # If not found, give a warning and pass through the input.
    } else {
      if(!quiet){
        warning("Biomarker not found: ", id, call. = FALSE)
      }
      return(id)
    } 
  })
  n<-data.frame(uploaded=names,BBMRI_names=new_names)
  return(n)
}

#Functions to calculate the associatins between 2 matrices
cor.assoc <- function(dat,dat2,met,covID=NULL,method="pearson",quiet=FALSE){
    c<-t(sapply(met,function(x){
      cor(x=dat[,x],y=dat2[,x],method=method, use="pairwise.complete.obs")
    }))
  c<-data.frame(met=met,cor=as.numeric(c))
  rownames(c)<-met
  return(c)
}

report.dim<-function(x,header,trailing="50"){
  return(paste0(sprintf(paste0("%-",trailing,"s"),paste0("| ",header,": ")),sprintf("%4s",ncol(x))," metabolites x ",sprintf("%4s",nrow(x))," samples \n"))
}

subset.samples.miss<-function(x,Nmax=1,quiet=FALSE){
  MISS <- colSums(is.na(t(x)))
  x <- x[which(MISS<=Nmax),,drop=FALSE]
  if(!quiet){
    cat(report.dim(x,header=paste0("Pruning samples on missing values [Nmax>=",Nmax,"]")))
  }
  return(invisible(x))
}

subset.samples.zero<-function(x,Nmax=1,quiet=FALSE){
  ZERO <- colSums(t(x==0),na.rm=TRUE)
  x <- x[which(ZERO<=Nmax),,drop=FALSE]
  if(!quiet){
    cat(report.dim(x,header=paste0("Pruning samples on zero values [Nmax>=",Nmax,"]")))
  }
  return(invisible(x))
}

subset.samples.sd<-function(x,MEAN,SD,quiet=FALSE,d=5){
  MEAN <- MEAN[colnames(x)]
  SD <- SD[colnames(x)]
  Dummi <- x
  # Exclude persons being an outlier:
  outl_samp <- rownames(Dummi)[unique(which(((Dummi > t(replicate(nrow(Dummi),MEAN)) + d*t(replicate(nrow(Dummi),SD))) | (Dummi < t(replicate(nrow(Dummi),MEAN)) - d*t(replicate(nrow(Dummi),SD)))),arr.ind=TRUE)[,"row"])]
  sample_names <- setdiff(rownames(Dummi),outl_samp)
  x <- x[sample_names,,drop=FALSE]
  if(!quiet){
    cat(report.dim(x,header=paste0("Pruning samples on 5SD")))
  }
  return(invisible(x))
}



density_scatterplot<- function(x, p, xname, yname, title=xname, subtitle="", my.formula = y ~ x,  xlimits, ylimits) {
  
    pl<-ggplot(data.frame(outcome=x,predicted_outcome=p), 
               aes(x=outcome, y= predicted_outcome)) +
      xlab(xname) + ylab(yname)+
      labs(title = title,
           subtitle = subtitle)+
      geom_point(size=0.6, alpha = 0.5) +
      geom_density2d(aes(x=outcome, y= predicted_outcome), data.frame(outcome=x,predicted_outcome=p), 
                     size=0.7, colour="white")+
      geom_abline(intercept = 0, slope = 1,linetype="dotted", colour="red", size=0.8)+
      geom_smooth(method=lm, aes(fill=outcome), colour="deepskyblue",size=1.2)+
      coord_cartesian(xlim = xlimits, ylim = ylimits)+
      stat_poly_eq(formula = my.formula,size=5,
                   eq.with.lhs = "italic(hat(y))~`=`~",
                   aes(label = paste(..eq.label.., ..rr.label.., sep = "~`,`~")), 
                   parse = TRUE)+
      annotate("text", x=Inf, y=-Inf,hjust=1,vjust=-0.5,size=5,
               label = paste0("(R= ", as.character(round(cor(x, p, method = c("spearman"), use="pairwise.complete.obs"), digits= 3)),
                              ", med. error=", as.character(round(median(abs(x - p), na.rm=T), digits=3)),")"))+
      theme(plot.title = element_text(size = 16, face = "bold",hjust = 0.5),
            plot.subtitle = element_text(size = 12,hjust = 0.5),
            axis.text=element_text(size=12),
            axis.title=element_text(size=14,face="bold"))
  
  suppressMessages(suppressWarnings(plot(pl)))
  return(pl)
}


## Defining a function to prepare metabolite data for Joris' mortality score:
prep.mort_data <- function(dat,featID=c("pufa_fa","gp","glc","s_hdl_l","xxl_vldl_l","alb","phe","acace","ile","vldl_d",
                                        "leu","val","his","lac"),quiet=FALSE, rin=rin){
  if(!quiet){
    cat("|| Preparing data ... \n")
  }
  ## 1. Check for zeroes:
  to_fix <- names(which(colSums(dat[,featID]==0,na.rm=TRUE)>0))
  if(length(to_fix)>0){
    if(!quiet){
      cat(paste0("| Adding 1.0 to metabolites featuring zero's: '",paste(to_fix,collapse="', '"),"'\n"))
    }
    dat[,to_fix] <- dat[,to_fix] + 1
  } else {
    if(!quiet){
      cat(paste0("| No metabolites found featuring zero's \n"))
    }
  }
  ## 2. Scale:
  
  if(rin){
    dat[,featID] <- t(RIN(t(as.matrix(dat[,featID]))))
    if(!quiet){
      cat("| Perform RIN transform .. ")
    }
  }else{
    dat[,featID] <- scale(log(dat[,featID]),center=TRUE,scale=TRUE)
    if(!quiet){
      cat("| Perform log transform & scaling to zero mean and unity sd .. ")
    }
  }
  
  if(!quiet){
    cat("Done!\n")
  }
  return(dat)
}

## Defining a function to compute Joris' mortality score on metabolite data:
comp.mort_score <- function(dat,betas=mort_betas,quiet=FALSE, rin=FALSE){
  ## 1. Prepare data:
  if(!quiet){
    cat("=== Computing mortality score === \n")
  }
  prepped_dat <- prep.mort_data(dat,featID=betas$Abbreviation,quiet=quiet,rin=rin)
  ## 2. Compute:
  if(!quiet){
    cat("| Computing score .. \n")
  }
  mortScore <- as.vector(as.matrix(prepped_dat[,betas$Abbreviation]) %*% betas$Beta_value)
  if("LLnr" %in% colnames(dat)){
    names(mortScore) <- dat$LLnr
  }
  if(!quiet){
    cat("Done!\n\n")
  }
  return(mortScore)
}