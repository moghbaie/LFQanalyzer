library(tidyverse)
library(R6)
library(xlsx)
library(doParallel)
library(BayesFactor)
library(reshape2)
library(ggpubr)
library(ggrepel)
library(ggridges)


## self
## 01/30/2021
## Preparing data from MaxQuant or ...

## Data preparation consists of four stages:
##  1. Remove contaminants and reverse proteins
##  2. Log transformation
## t-test and fold change analysis ...


set.seed(123)

###########################################################################
### Proteomic template
###########################################################################
###########################################################################
### Proteomic template
###########################################################################
TemplateProtein <- R6Class("TemplateProtein",
                           list(
                             input.dir = NA,
                             input = NA,
                             intensity.type = NA,
                             variables = list(),
                             metadata = NA,
                             interactors = NA ,
                             category = NA,
                             df = NA,
                             df_log = NA,
                             df_imputed = list(),
                             significant_seeds = list(),
                             df_significant = list(),
                             df_significant_woI = list(),
                             df_significant_table = list()
                           )
)



TemplateProtein$set("public","importInput", function(input.dir){
  print("Importing Maxquant output")
  self$input.dir <- input.dir
  ## Import metadata
  self$metadata <- xlsx::read.xlsx(file.path(self$input.dir,"output","backup","targetfile.xlsx"),1)
  ## Import proteinGroup files
  self$input <-  read.delim(file.path(self$input.dir, 'proteinGroups.txt'), header = TRUE)
  
  invisible(self)
}
)

TemplateProtein$set("public","removeContaminant", function(intensity.type , variables){
  print("Removing contaminant proteins and reversed sequences")
  self[["intensity.type"]] <- intensity.type
  self[["variables"]] <- variables
  cols =colnames(self$input)
  dl = self[["input"]] 
  dl$Peptide.counts..unique. <-  apply(dl,1,function(x) strsplit(as.character(x["Peptide.counts..unique."]),";")[[1]][1])
  dl$Protein.ID <- apply(dl,1,function(x) strsplit(as.character(x["Protein.IDs"]),";")[[1]][1])
  dl$Gene.name <- apply(dl,1,function(x) strsplit(as.character(x["Gene.names"]),";")[[1]][1])
  dl$Gene.name[is.na(dl$Gene.name)] <- dl$Protein.ID[is.na(dl$Gene.name)]
  
  intensity.cols <- str_c(str_c(intensity.type, ".intensity."), self$variables)
  
  self$df <- dl %>% 
    filter(!(dl$Potential.contaminant=="+"|dl$Reverse=="+"|dl[["Q.value"]]>0.05))%>%
    select(c("Protein.ID","Gene.name","Peptide.counts..unique.","Fasta.headers",intensity.cols))
  
  invisible(self)
})



TemplateProtein$set("public","transformData", function(){
  print("Log transformation from intensities")
  self[["df_log"]] <- self[["df"]]
  
  self[["df_log"]][,grepl(self$intensity.type,colnames(self[["df_log"]]))] <- data.frame(apply( self[["df_log"]][,grepl(self$intensity.type,colnames(self[["df_log"]]))],2, function(x) ifelse(x>0,log2(as.numeric(x)),0)))
  
  invisible(self)
})


TemplateProtein$set("public","anovaAnalysis", function( imputation = "both",
                                                        adjustment = "both",
                                                        bf = "yes"){
  print("Performing normalization, imputation and t-test")
  
  grid <- self$metadata %>%  dplyr::select(Comparison) %>% unique() %>% .$Comparison
  df <- self$df_log
  
  colnames(df)
  
  for(g in grid){
    print(g)
    case <- self$metadata %>% filter(Comparison == g & Type=="case") %>% dplyr::select(Condition) %>% .$Condition
    control <- self$metadata %>% filter(Comparison == g & Type=="control") %>% dplyr::select(Condition) %>% .$Condition
    
    case_reps <- str_c(str_c(self$intensity.type,".intensity."), case)
    control_reps <- str_c(str_c(self$intensity.type,".intensity."), control)
    
    print(case_reps)
    print(control_reps)
    ## Removing records with zero intensities in all replicates
    dfnz <- df[apply(df%>% dplyr::select( case_reps, control_reps),1,sum)!=0,]
    
    dfi <- dfnz[,!grepl(str_c(self$intensity.type,".intensity."), colnames( dfnz))]
    
    cols2 <- c(colnames(dfi), case_reps, control_reps)
    dfi <- cbind(dfi, matrix(NA, ncol=length(c(case_reps, control_reps)), nrow= dim(dfi)[1]))
    colnames(dfi) <- cols2
    
    ### Calculating significant protein without imputed significant table
    
    dff <- dfnz[,!grepl("LFQ.intensity", colnames( dfnz))]
    dff[["avg.Case.log.intensity"]] <- apply(dfnz[,case_reps],1, function(x) ifelse(!is.na(mean(x[x!=0])), mean(x[x!=0]),0)) 
    dff[["avg.Control.log.intensity"]] <-apply(dfnz[,control_reps],1, function(x) ifelse(!is.na(mean(x[x!=0])), mean(x[x!=0]),0))
    dff[["avg.Case.num.reps"]] <- apply(dfnz[,case_reps],1, function(x) length(x[x!=0])) 
    dff[["avg.Control.num.reps"]] <- apply(dfnz[,control_reps],1, function(x) length(x[x!=0]))
    dff[["log2foldChange"]] <- as.numeric(dff[["avg.Case.log.intensity"]]) - as.numeric(dff[["avg.Control.log.intensity"]])
    dff["p.value"] <- NA
    
    tt <- foreach(k=1:dim(dff)[1], .combine=rbind)%do%{
      t.test(unlist(dfnz[k,case_reps]), unlist(dfnz[k,control_reps]))$p.value
    }
    
    dff[["p.value"]] <- data.frame(tt)$tt
    dff[["p.adjust.value"]] <- p.adjust(dff[["p.value"]] ,method="BH")
    dff[["Significant"]] <- ifelse((dff[['log2foldChange']] >1 &dff[["avg.Case.num.reps"]]>1)&dff[["p.adjust.value"]]<0.05,'Yes','No')
    
    self$df_significant_woI[[g]] <- dff
    
    if(imputation =="both"|imputation == "yes"){
      ### Calculating significant protein with imputation using 10 different seeds
      cl <- parallel::makeCluster(parallel::detectCores()-1)
      doParallel::registerDoParallel(cl)
      
      ddd <- foreach(seed=1:10, .combine='rbind', .packages = c("doParallel"), .export =c("impute", "calculate_stats_nonzeros","impute_all_zeros","min_col","count_zeros","na_zeros_impute","impute_partial_zeros")) %dopar% {
        set.seed(seed)
        dfi <- dfnz[,!grepl(paste0(self$intensity.type,".intensity"), colnames( dfnz))]
        cols2 <- c(colnames(dfi), case_reps, control_reps)
        dfi <- cbind(dfi, matrix(NA, ncol=length(c(case_reps, control_reps)), nrow= dim(dfi)[1]))
        colnames(dfi) <- cols2
        dfi[,case_reps] <- impute(dfnz[case_reps], amm = "2", pmm = "6")
        dfi[,control_reps] <- impute(dfnz[control_reps], amm = "2", pmm = "6")
        print("got here")
        return(dfi)
      }
      
      dfi <- ddd %>% group_by(Protein.ID, Gene.name, Peptide.counts..unique., Fasta.headers) %>%
        summarise(across(c(control_reps, case_reps), ~ mean(.x)))
      
      self$df_imputed[[g]] <- dfi
      
      dfx <- dfi[,!grepl(str_c(self$intensity.type,".intensity"), colnames( dfi))]
      dfx[["avg.Case.log.intensity"]] <- apply(dfi[,case_reps],1, function(x) ifelse(!is.na(mean(x[x!=0])), mean(x[x!=0]),0)) 
      dfx[["avg.Control.log.intensity"]] <- apply(dfi[,control_reps],1, function(x) ifelse(!is.na(mean(x[x!=0])), mean(x[x!=0]),0))
      dfx[["avg.Case.num.reps"]] <- dff[["avg.Case.num.reps"]][match(dfx$Protein.ID,dff$Protein.ID)]
      dfx[["avg.Control.num.reps"]] <- dff[["avg.Control.num.reps"]][match(dfx$Protein.ID,dff$Protein.ID)]
      dfx[["log2foldChange"]] <- as.numeric(dfx[["avg.Case.log.intensity"]]) - as.numeric(dfx[["avg.Control.log.intensity"]])
      
      dfx["p.value"] <- NA
      
      tt2 <- foreach(k=1:dim(dfi)[1], .combine=rbind)%do%{
        t.test(unlist(dfi[k,case_reps]),unlist(dfi[k,control_reps]))$p.value}
      
      dfx[["p.value"]] <- data.frame(tt2)$tt2
      dfx[["p.adjust.value"]] <- p.adjust(dfx[["p.value"]] ,method="BH")
      dfx[["Significant"]] <- ifelse(((dfx[['log2foldChange']] >1 &dfx[["avg.Case.num.reps"]]>1))&dfx[["p.adjust.value"]]<0.05,'Yes','No')
      
      
      if(adjustment == "both"|adjustment =="yes"){
        ### Unifying case and control median based on unimputed intensities
        dfi2 <- t(t(dfi[,c(case_reps,control_reps)])/(apply(dfnz[c(case_reps, control_reps)],2,function(x) median(x[x!=0]))/ min(apply(dfnz[c(case_reps, control_reps)],2,function(x) median(x[x!=0])))))
        dfi2 <- data.frame(dfi2)
        dfi2[["Protein.ID"]] <- dfi[["Protein.ID"]]
        dfi2[["Gene.name"]] <- dfi[["Gene.name"]]
        
        dfx[["log2foldChange.AC"]] <- apply(dfi2[,case_reps],1, function(x) ifelse(!is.na(mean(x[x!=0])), mean(x[x!=0]),0)) - apply(dfi2[,control_reps],1, function(x) ifelse(!is.na(mean(x[x!=0])), mean(x[x!=0]),0))
        dfx[["p.value.AC"]] <- NA
        
        tt3<- foreach(k=1:dim(dfi2)[1], .combine=rbind)%do%{
          t.test(unlist(dfi2[k,case_reps]),unlist(dfi2[k,control_reps]))$p.value}
        
        dfx[["p.value.AC"]] <- data.frame(tt3)$tt3
        dfx[["p.adjust.value.AC"]] <- p.adjust(dfx[["p.value.AC"]] ,method="BH")
        
        dfx[["Significant.AC"]] <- ifelse(((dfx[['log2foldChange.AC']] >1 &dfx[["avg.Case.num.reps"]]>1))&dfx[["p.adjust.value.AC"]]<0.05,'Yes','No')
        
        if(bf == "yes"){
          bfe <- foreach(k=1:dim(dfi2)[1], .combine=rbind, .packages = c("doParallel","BayesFactor"))%dopar%{
            dh <- data.frame(rbind(cbind(unlist(dfi2[k,case_reps]),rep("case",length(dfi2[k,case_reps]))),
                                   cbind(unlist(dfi2[k,control_reps]),rep("control",length(dfi2[k,control_reps])))))
            colnames(dh) <- c("intensity", "condition")
            dh$intensity <- as.numeric(dh$intensity)
            return(ttestBF(formula = intensity ~ condition, data = dh)@bayesFactor$bf)
          }
          
          dfx[["bf"]] <- data.frame(bfe)$bfe
          dfx[["SignificantB"]] <- ifelse(dfx$bf>3 & (dfx$log2foldChange.AC>1 & dfx$avg.Case.num.reps>1), "Yes","No")
        }
      }
      parallel::stopCluster(cl)
    }
    
    self$df_significant[[g]] <- dfx
  }
  
  invisible(self)
})



TemplateProtein$set("public","saveTable", function(){

  if(file.exists(file.path(self$input.dir,"output","table","significant.xlsx"))){
    file.remove(file.path(self$input.dir,"output","table","significant.xlsx"))}
  
  if(file.exists(file.path(self$input.dir,"output","table","unimputed.xlsx"))){
    file.remove(file.path(self$input.dir,"output","table","unimputed.xlsx"))}
  
  if(file.exists(file.path(self$input.dir,"output","table","imputed.xlsx"))){
    file.remove(file.path(self$input.dir,"output","table","imputed.xlsx"))}
  
  xlsx::write.xlsx(as.data.frame(self$df_log), file=file.path(self$input.dir,"output","table","unimputed.xlsx"),sheetName="original.log2.LFQ", row.names=FALSE, append=TRUE)
  
  for(name in names(self$df_imputed)){
    xlsx::write.xlsx(as.data.frame(self$df_imputed[[name]]), file=file.path(self$input.dir,"output","table","imputed.xlsx"),sheetName=name, row.names=FALSE, append=TRUE)
  }
  
  for(name in names(self$df_significant)){
    xlsx::write.xlsx(as.data.frame(self$df_significant[[name]]), file=file.path(self$input.dir,"output","table","significant.xlsx"),sheetName=name, row.names=FALSE, append=TRUE)
  }
  
  invisible(self)
})


TemplateProtein$set("public","visualize", function(){
  
  dt <- self$df_log
  dt1 <- melt(dt)
  dt1 <- dt1[!dt1$value==0,]
  dt2 <- dt1 %>%
    group_by(variable) %>%
    summarize(n=n())
  
  new.file <- self$metadata %>% filter(!is.na(Comparison)) %>% .$Condition
  
  file.variable <- str_c(str_c(self$intensity.type, ".intensity."), new.file)
  
  dt2[['Comparison']] <- self$metadata$Comparison[match(dt2$variable, str_c(str_c(self$intensity.type, ".intensity."),self$metadata$Condition))]
  
  dt2 <- dt2[!is.na(dt2$Comparison),]
  dt2 <- dt2[order(-dt2$n),]
  dt2$variable <- factor(dt2$variable, levels= dt2$variable[order(-dt2$n)])
  
  
  p1 <- ggplot(dt2, (aes(x = variable, y = n))) +
    geom_bar(aes(fill = n), stat = "identity") +
    scale_color_gradient(low="blue", high="green")+
    xlab("sample") +
    ylab("Number of records with non zero intensity") +
    ggtitle("Number of proteins with intensity > 0") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))+
    scale_fill_viridis_b()+
    facet_grid(.~ Comparison, scales='free')
  
  ggsave(file= file.path(self$input.dir,"output","img","intensity.barplot.pdf"), p1, width=5, height=5, dpi=200)
  ggsave(file= file.path(self$input.dir, "output","img","intensity.barplot.png"), p1, width=5, height=5, dpi=200)
  
  ### boxplot after imputation
  
  df <-  self$df_log
  df2 <- melt(df)
  df3 <- df2[!df2$value==0,]
  new.file <- self$metadata %>% filter( !is.na(Comparison)) %>% .$Condition
  
  df3$condition <- self$metadata$Comparison[match(df3$variable, str_c(str_c(self$intensity.type, ".intensity."),self$metadata$Condition))]
  
  df3 <- df3[!is.na(df3$condition),]
  df3[['Comparison']] <- self$metadata$Type[match(df3$variable, str_c(str_c(self$intensity.type, ".intensity."),self$metadata$Condition))]
  
  pbox <- ggpubr::ggboxplot(df3, 
                            x = "variable",
                            y = "value",
                            combine = TRUE,
                            color = "Comparison",
                            facet.by = c('Comparison',"condition"),
                            scales="free_x",
                            space="free_x",
                            palette = "jco",
                            ylab = "log LFQ.intensity", 
                            add = "jitter",                               # Add jittered points
                            add.params = list(size = 1, jitter = 0.2),  # Point size and the amount of jittering
                            label = "Gene.name",                # column containing point labels
                            label.select = list(top.up = 5, top.down = 0),# Select some labels to display
                            font.label = list(size = 11, face = "italic"), # label font
                            repel = TRUE  ,
                            max.overlaps= 35) + # Avoid overlapping
    rotate_x_text(angle=60)
  
  ggsave(file=file.path(self$input.dir,"output","img","boxplot.pdf"), pbox, width=6, height=7, dpi=200, limitsize = F)
  ggsave(file=file.path(self$input.dir,"output","img","boxplot.png"), pbox, width=6, height=7, dpi=200, limitsize = F)
  
  
  
  ### Pre and Post Imputation 
  
  
  for(cell in names(self$df_imputed)){
    cols <- colnames(self$df_imputed[[cell]])[grepl("LFQ.",colnames(self$df_imputed[[cell]]))]
    dz <- self$df_log[apply(self$df_log[,cols], 1,sum )>0, c("Gene.name",cols)]
    dw <- self$df_imputed[[cell]][,grepl("LFQ|Gene.name", colnames(self$df_imputed[[cell]]))]
    
    dzm <- melt(dz)
    dwm <- melt(dw)
    dzm[["category"]] <- "PreImputation"
    dwm[["category"]] <- "PostImputation"
    dm <- rbind(dzm,dwm)
    colnames(dm) <- c("Gene.name","variable","value", "category")
    dm$type <- self$metadata$Type[match(dm$variable, dm$variable[unlist(lapply(self$metadata%>% .$Condition, function(x) stringdist::amatch(x, dm$variable, maxDist=Inf)))])]
    
    
    pd <- ggplot(dm, aes(x = value, y = variable, color = category, point_color = category, fill = category)) +
      geom_density_ridges(
        jittered_points = TRUE, scale = .95, rel_min_height = .01,
        point_shape = "|", point_size = 3, size = 0.25,
        position = position_points_jitter(height = 0)
      ) +
      scale_y_discrete(expand = c(0, 0)) +
      scale_x_continuous(expand = c(0, 0), name = "log LFQ intensity") +
      scale_fill_manual(values = c("#D55E0050", "#0072B250"), labels = c("Post Imputation", "Pre Imputation")) +
      scale_color_manual(values = c("#D55E00", "#0072B2"), guide = "none") +
      scale_discrete_manual("point_color", values = c("#D55E00", "#0072B2"), guide = "none") +
      coord_cartesian(clip = "off") +
      guides(fill = guide_legend(
        override.aes = list(
          fill = c("#D55E00A0", "#0072B2A0"),
          color = NA, point_color = NA)
      )
      ) +
      ggtitle(paste("Density plot before and after imputation", cell)) +
      theme_ridges(center = TRUE)+
      theme(axis.text.y = element_text(size = 10))+
      facet_grid(type ~., scales="free", switch="y")
    
    
    ggsave(file= file.path(self$input.dir,"output","img",paste0(cell,".density.imputation.pdf")), pd, width=12, height=9, dpi=100)
    ggsave(file= file.path(self$input.dir,"output","img",paste0(cell,".density.imputation.png")), pd, width=12, height=9, dpi=100)
    
    
  }
  
  ## drawing volcano plot
  
  for(cell in names(self$df_significant)){
    print(cell)
    dx <- self$df_significant[[cell]]
    shade = data.frame(x1=c( 1), 
                       x2=c( Inf),
                       y1=c( -log10(0.05)), 
                       y2=c( Inf))
    fold_cutoff = 1
    pvalue_cutoff = 0.05
    
    
    
    
    dx2 <- subset(dx, (Significant=="Yes" |Significant.AC == "Yes"| SignificantB == "Yes"))
    
    p <- ggplot(dx) +
      theme_bw()+
      geom_rect(data=shade, 
                mapping=aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2), fill='grey95')+
      geom_point(data  = dx[!dx$Protein.ID %in% dx2$Protein.ID,], aes(x =log2foldChange, y = -log10(p.adjust.value)),colour = "gray", alpha = 0.5,size=1)+
      geom_point(data  = dx2, aes(x =log2foldChange, y = -log10(p.adjust.value)), colour = "red",alpha = 1,size=1)+
      geom_vline(xintercept = fold_cutoff, col = "blue")+
      geom_vline(xintercept = -fold_cutoff, col = "blue")+
      geom_hline(yintercept = -log10(pvalue_cutoff), col = "green2")+
      ggtitle(cell)+
      geom_text_repel(data  = dx2,
                      aes(x=log2foldChange, y=-log10(p.adjust.value),label=Gene.name),
                      segment.alpha =0.35,max.overlaps =25,
                      size = 2.5)
    
    ggsave(file=file.path(self$input.dir,"output","img", paste0(cell,"_volcanoplot.pdf")), p, width=10, height=8, dpi=200)
    ggsave(file=file.path(self$input.dir,"output","img", paste0(cell,"_volcanoplot.png")), p, width=10, height=8, dpi=200)
  }
  
  
  ## Comparing significant proteins
  
  
  for(cell in names(self$df_significant)){
    png(file=file.path(self$input.dir,"output","img",paste0(cell,".pvalue.distribution.png")), width=800, height=500)
    par(mfrow=c(1,2))
    plot(density(self$df_significant_woI[[cell]]$p.value), main=paste("Pre Imputation", cell), xlab="", cex=2,xlim=c(0,1))
    lines(density(self$df_significant_woI[[cell]]$p.adjust.value), main="",  add=T, col='red')
    
    plot(density(self$df_significant[[cell]]$p.value),  main=paste("Post Imputation",cell), xlab="", cex=2,xlim=c(0,1))
    lines(density(self$df_significant[[cell]]$p.adjust.value), main="", add=T, col='red')
    legend("top", inset=.01, title="",
           c("P.value","P.adjusted.value"), fill=c("black","red"), horiz=TRUE, cex=0.7)
    dev.off()
  }
  
  invisible(self)
})
