#!/usr/bin/env Rscript

#Script continuing from 1_DataProcessing_MLphylogenetics.sh, 2_SNPcalling.sh, 3_speedemon.sh. 
#This one is for running several species clustering analyses using the snp calling output and speedemon output.

#Using Clade1 of genus Acropora data (Scleractinia: Acroporidae)
#Author: Hanaka Mera   mailto:hanaka.mera[@]my.jcu.edu.au
    #Tutorials and papers followed (in this script):
    #Derkarabetian, S., Castillo, S., Koo, P.K., Ovchinnikov, S., Hedin, M., 2019. A demonstration of unsupervised machine learning in species delimitation. Molecular Phylogenetics and Evolution 139, 106562. https://doi.org/10.1016/j.ympev.2019.106562
    #https://github.com/shahanderkarabetian/uml_species_delim

#Recommended: open this file in R studio and toggle 'document outline' so you can see the headings

#Last edited: July 2023

###################################################################

## Calculate Fst ----
#install.packages("StAMPP") #https://www.biostars.org/p/309019/
library(StAMPP)
library(vcfR)
library(xlsx)

#Read vcf, convert to genlight, then convert to matrix
forstampp_refC99 <- vcfR2genlight(read.vcfR("data/refC99_snps_95.vcf"))
forstampp_refC99_matrix <- as.matrix(forstampp_refC99)
forstampp_refC99_samples <- row.names(forstampp_refC99_matrix)
pop.names <- metadata$pop[match(forstampp_refC99_samples,metadata$UCE_Tube_Label)]

#convert allele counts to frequency
forstampp_refC99_matrix = forstampp_refC99_matrix * (1/ploidy(forstampp_refC99)) #convert allele counts to frequency 
forstampp_refC99_matrix[is.na(forstampp_refC99_matrix)] = NaN 
format <- vector(length = length(forstampp_refC99_samples))
#format id for the genotype data
format[1:length(format)] = "freq"  
#convert to basic r data.frame suitable to stamppConvert 
stampp_refC99 <- as.data.frame(cbind(forstampp_refC99_samples, pop.names, ploidy(forstampp_refC99), format, forstampp_refC99_matrix)) 
geno_refC99 <- stamppConvert(stampp_refC99, 'r') 
#Calculate Fst
fst_refC99 <- stamppFst(geno_refC99,nboots=100,percent=95,nclusters=2)
fst_refC99$Fsts
#Calculate Nei's genetic distance (Nei 1972) between populations or individuals
neisd_refC99 <- stamppNeisD(geno_refC99,pop=TRUE)
colnames(neisd_refC99) <- rownames(neisd_refC99)
heatmap(neisd_refC99)


stamppAmova()
write.xlsx(x=fst_refC99$Fsts,file="export/FST.xlsx",sheetName="Fst")
write.xlsx(x=fst_refC99$Pvalues,file="export/FST.xlsx",sheetName="pvalue",append=TRUE)

###################################################################

## sNMF ----
library(LEA)
library(tidyverse)
library(pophelper)

#Below two functions were copied and pasted from somewhere, possible in a error discussion thread for pophelper
mergeQ2 <- function(qlist) {
  
  is.qlist(qlist)
  if(diff(range(as.integer(tabulateQ(qlist)$ind)))!=0) stop("mergeQ: Number of individuals differ between runs.")
  
  # Computes mean cell-wise across dataframes
  # x A list of numeric dataframes
  # 
  mergy <- function(x) {
    return(list(Reduce(`+`, x)/length(x)))
  }
  
  # if all runs have same K, merge as is
  if(diff(range(as.integer(tabulateQ(qlist)$k)))==0) {
    labels <- summariseQ(tabulateQ(qlist))$k
    x <- mergy(qlist)
    names(x) <- labels
  }else{
    # if runs have different K, split them and merge within sublists
    qlist <- sortQ2(qlist)
    labels <- summariseQ(tabulateQ(qlist,sorttable=FALSE))$k
    x <- unlist(lapply(splitQ(qlist),mergy),recursive=FALSE)
    names(x) <- labels
  }
  
  return(as.qlist(x))
}
sortQ2 <- function(qlist,by="k",decreasing=FALSE,debug=FALSE) {
  
  is.qlist(qlist)
  if(length(by)==0) stop("sortQ: Argument 'by' must not be length zero.")
  if(!is.character(by)) stop("sortQ: Argument 'by' must be a character.")
  if(!is.logical(decreasing)) stop("sortQ: Argument 'decreasing' must be a logical datatype.")
  
  fun1 <- function(x) as.matrix(unlist(attributes(x)))
  a <- lapply(qlist,fun1)
  if(debug) print(a)
  if(any(!sapply(a,function(x) any(grepl(paste0(by,collapse="|"),rownames(x)))))) {
    stop(paste0("One or more of the attributes provided in by (",by,") is missing in one or more runs. If 'ind' or 'k' is missing, use 'as.qlist()' to add them."))
  }
  
  # get df of attributes
  b <- as.data.frame(t(as.data.frame(lapply(a,function(x,y) x[y,],by),stringAsFactors=FALSE)),stringsAsFactors=FALSE)
  fun2 <- function(x) if(all(!is.na(as.numeric(as.character(x))))) {return(as.numeric(as.character(x)))}else{return(x)}
  b <- as.data.frame(sapply(b,fun2),stringAsFactors=FALSE)
  
  if(debug) {print(str(b)); print(b)}
  
  # order
  ord <- do.call(order,b[,by,drop=FALSE])
  if(decreasing) ord <- rev(ord)
  # sort qlist
  return(qlist[ord])
}

packageDescription("pophelper",fields='Version')
collist <- list("cols_1A"=c("#006666","#009999","#00CCCC","#00FFFF","#65FFFF","#B2FFFF","#E5FFFF","#FFE5CB","#FFCA99","#FFAD65","#FF8E33","#FF6E00","#CC5500","#993D00","#FFFFFF"),
                "cols_1B"=c("#762A83","#9970AB","#C2A5CF","#D9F0D3","#A6DBA0","#5AAE61","#FFFFFF"),
                "cols_1C"=c("#FDBEDE","#67001F","#A50021","#D8152F","#D6604D","#FFAC75","#FDDBC7","#FFD699","#FFF1BC","#D1E5F0","#BCF9FF","#99EAFF","#92C5DE","#75D3FF","#4393C3","#2166AC","#053061","#FFFFFF"))

#Convert .vcf to .geno
vcf2geno("data/refC99_snps_95.vcf")
#Run sNMF
sNMF_refC99 = snmf("data/refC99_snps_95.geno",K=2:50,entropy=T,alpha=100,
                        repetitions=5,project="new",iterations=1000,CPU=2)

## view snmf within the LEA package ----
#Plot cross-entropy criterion of all runs of the project
#"Smaller values of the cross-entropy criterion usually mean better runs. 
#... choose the value of K for which the cross-entropy curve exhibits a plateau"
summary(sNMF_refC99)
png("export/plots/refC99_sNMF_crossentropy.png",width=3,height=3,units="cm",res=200,pointsize=5)
par(mai=c(0.25,0.25,0.1,0.1))
plot(sNMF_refC99, col="blue",pch=16,cex=1,cex.axis=0.8,tck=-0.02,cex.lab=0.8,mgp=c(1.5,0.5,0))
dev.off()
# Best K?
summary(sNMF_refC99)$crossEntropy[2,][which.min(summary(sNMF_refC99)$crossEntropy[2,])]
# get the cross-entropy of all runs for K=10
ce_refC99 = cross.entropy(sNMF_refC99,K=10)
# select the run with the lowest cross-entropy for K=
best_refC99 = which.min(ce_refC99)
Q.matrix_refC99 <- round(Q(sNMF_refC99,K=10,run=best_refC99))
G.matrix_refC99 <- round(G(sNMF_refC99,K=10,run=best_refC99))
#Plot
#c(RColorBrewer::brewer.pal(n=5,'Blues'),RColorBrewer::brewer.pal(n=5,'Greens'))
barchart(sNMF_refC99, K=10, run=best_refC99, sort.by.Q=T,
         border=NA, space=0, col=c(collist$cols_1A,collist$cols_1B,collist$cols_1C),xlab="Individuals", ylab="Ancestry proportions") -> bp_refC99
#axis(1, at=1:length(bp$order),labels=bp$order, las=3, cex.axis=0.4)
#Rename labels
#system('bcftools query -l filt_filtered_snps_95.vcf > samples.list')
indivs_refC99 <- read.table('data/samples.list')
indivs_refC99 = indivs_refC99 %>% mutate(order=row_number())
match(indivs_refC99$order,bp_refC99$order)
bp_refC99$label <- indivs_refC99$V1[match(indivs_refC99$order,bp_refC99$order,nomatch=NULL)]
barchart(sNMF_refC99, K=7, run=best_refC99, sort.by.Q=T,
         border=NA, space=0, col=c(collist$cols_1A,collist$cols_1B,collist$cols_1C),xlab="Individuals", ylab="Ancestry proportions")
axis(1, at=1:length(bp_refC99$label),labels=bp_refC99$label, las=3, cex.axis=0.4)

## Export Q matrices from sNMF for pophelper ----
for (g in 2:40){
  for (i in 1:5){
    Q.matrix_refC99=round(Q(sNMF_refC99, K=g, run=i),3)
    write.table(Q.matrix_refC99, paste0('data/refC99_snps_95.snmf/K',g,'/run',i,'/K',g,'_run',i,'.txt'),col.names=F,row.names=F)
  }
}

## Plot with pophelper ----
sfiles_refC99 <- list.files(path='data/refC99_snps_95.snmf',full.names=T,recursive=T,pattern='txt')
slist_refC99 <- readQ(files=sfiles_refC99)
slist_refC99 <- lapply(slist_refC99,"rownames<-",indivs_refC99$V1)
tr_refC99 <- tabulateQ(slist_refC99)
sr_refC99 <- summariseQ(tr_refC99)
slist_refC99_aligned <- alignK(slist_refC99[c(181:185,1:5,136:140)]) #K=7,10,35
#slist_refC99_aligned_merged <- mergeQ(slist_refC99_aligned) #This used to work with Rv4.2 but spits out error with v4.3
slist_refC99_aligned_merged <- mergeQ2(slist_refC99_aligned)
plotq_refC99 <- plotQ(slist_refC99_aligned_merged,exportplot=F,returnplot=T,imgoutput='join',
                    basesize=8,sortind='all',sharedindlab=F,showindlab=T,useindlab=T,
                    showdiv=T,divsize=3) 
ggsave("export/refC99_sNMF.png",plot=plotq_refC99$plot[[1]],device="png",width=8,height=6,units="in")

###################################################################

## DAPC ----
library(adegenet)
library(tidyverse)
library(showtext)
font_add_google(name="Work Sans",family="work")
showtext_auto(TRUE)

### Read files ----
#system("/Users/bioinformatics/programs/plink/plink --vcf data/filt_filtered_snps_95.vcf --const-fid 0 --allow-extra-chr --recode 'structure' --out data/structure_snps_vcf95")
#system("mv data/structure_snps_vcf95.recode.strct_in data/structure_snps_vcf95.str")
data_refC99 <- read.structure("data/refC99_structure_snps_vcf95.str",n.ind=176, n.loc=1942, col.lab=1,col.pop=0,
                            col.others=NULL,row.marknames=NULL,onerowperind=TRUE,ask=FALSE)
  #refC99: genotypes=176, markers=1942, col for labels=1, pop=0, options=blank, 0, n
  #refcons: genotypes=176, markers=2019, col for labels=1, pop=0, options=blank, 0, n
metadata <- read.csv("../metadata_Cl1.csv") %>% 
  arrange(UCE_Tube_Label) %>% #slice(4:163) %>% filter(!(PSH=='Outgroup')) %>% 
  rename(pop=PSH) %>% 
  mutate(pop_Cnum=paste0(pop,"_",UCE_Tube_Label))
## Edit data with metadata so that the CLUSTER info is included
pop(data_refC99) <- metadata$subclade[match(indNames(data_refC99),metadata$UCE_Tube_Label)]
pop(data_refC99) <- metadata$pop[match(indNames(data_refC99),metadata$UCE_Tube_Label)]
indNames(data_refC99) <- metadata$pop_Cnum[match(indNames(data_refC99),metadata$UCE_Tube_Label)]  #to change individuals' labels

### PCA-DAPC ----
data_refC99_scaled <- scaleGen(data_refC99, center=FALSE, scale=FALSE, NA.method=c("zero"), nf)
# PCA, can adjust nf to include more components
pca_refC99 <- dudi.pca(data_refC99_scaled, center=TRUE, scale=TRUE, scannf=FALSE, nf=2)
plot(pca_refC99$li, xlab="Scaling Dimension 1", ylab="Scaling Dimension 2", main="PCA", pch=16)
s.label(pca_refC99$li, clabel=0.5, grid=0)

# DAPC (interactive, requires input. Determine optimal number of genetic clusters)
dev.off()
#Comment back in below to save the results as single .png after you figure out the # of PCs to retain and number of clusters
png("export/plots/refC99_PCA-DAPC_BIC.png",width=5,height=5,units="cm",res=200,pointsize=5)
#Change max.n.clust equal to number of pops
clusters_refC99 <- find.clusters(data_refC99, max.n.clust=40, n.iter=1e6, n.start=10)
#PCs to retain=150, number of clusters=5
dev.off()
#See quick assignment result
table(pop(data_refC99),clusters_refC99$grp)
#Run DAPC
results_refC99 <- dapc(data_refC99, grp=clusters_refC99$grp, perc.pca=NULL)
#PCs to retain=150, discriminant functions to retain=10
assignplot(results_refC99)
compoplot(results_refC99)
grp_k_refC99 <- nlevels(clusters_refC99$grp) 
scatter(results_refC99, scree.da=F, cell=0, cstar=0, solid=0.4, cex=2, clab=0.5)

# "PCA with DAPC informed clusters"
png("export/plots/refC99_PCA-DAPCclusters.png",width=5,height=5,units="cm",res=200,pointsize=7)
ggplot(data=pca_refC99$li,aes(y=Axis2,x=Axis1))+
  #geom_point(aes(colour=factor(results_refC99$grp,levels=c("1A","1B","1C"))),size=0.5)+
  geom_point(aes(colour=factor(results_refC99$grp)),size=0.5)+
  annotate(geom="text",x=0,y=72,label=paste0("K=",grp_k_refC99),size=2)+
  theme_bw()+ #scale_colour_manual(values=c("darkgreen","blue","red"))+
  labs(x="",y="",title="PCA with DAPC informed clusters")+
  theme(text=element_text(family="work",size=10),legend.title=element_blank(),
        legend.position='inside',legend.position.inside=c(0.2,0.7),
        legend.text=element_text(size=3,margin=ggplot2::margin(l=0.1,r=0.1,unit="mm")),
        legend.box.background=element_rect(colour='black',linewidth=0.1),
        legend.background=element_blank(),legend.margin=ggplot2::margin(0.1,0.1,0.1,0.1,"mm"),
        legend.key.size=unit(2,"mm"),panel.grid=element_line(linewidth=0.1),
        axis.text=element_text(size=4),axis.ticks=element_line(linewidth=0.05),
        axis.title=element_blank(),plot.title=element_text(size=5))
dev.off()
#Merge .pngs
system("montage export/plots/refC99_PCA-DAPC* -geometry +1+1 export/refC99_PCA-DAPC.png")
#Add to cluster result compilation
DAPCRFtSNEresults_refC99 <- as.data.frame(clusters_refC99$grp,nm="DAPC")

###################################################################

## RF ----
library(randomForest)
library(cluster)
library(PCDimension)
library(mclust)
library(factoextra)  # detach("package:adegenet",unload=TRUE) if throws ggplot error
library(MASS)

# Convert genind scaled data to factors for randomForest
data_refC99_conv <- as.data.frame(data_refC99_scaled)
data_refC99_conv[is.na(data_refC99_conv)] <- ""
data_refC99_conv[sapply(data_refC99_conv,is.integer)] <- lapply(data_refC99_conv[sapply(data_refC99_conv,is.integer)],as.factor)
data_refC99_conv[sapply(data_refC99_conv,is.character)] <- lapply(data_refC99_conv[sapply(data_refC99_conv,is.character)],as.factor)
nsamp_refC99 <- nrow(data_refC99_conv)
# Unsupervised random forest ~20min
rftest_refC99 <- randomForest(data_refC99_conv, ntree=5000)

### Classic MDS ---------------------------------------------------------
# cMDS with optimal number of components to retain using broken-stick
# may need to adjust number of dimensions if given error (originally nsamp-1)
# ALSO, Adjust plot margin first because it also spits out a 'figure margins too large' error
dev.off()
par(mar=c(0,0,0,0),mai=c(0,0,0,0))
pdf(file="export/plots/deletethis.pdf",width=100,height=100)  
cmdsplot1_refC99 <- randomForest::MDSplot(rftest_refC99, fac=results_refC99$grp, k=nsamp_refC99-40)
cmdsplot_bstick_refC99 <- bsDimension(cmdsplot1_refC99$eig)
cmdsplot2_refC99 <- randomForest::MDSplot(rftest_refC99, results_refC99$grp, cmdsplot_bstick_refC99)
dev.off()

#### 1. cMDS ----
# Use optimal k from DAPC 
png(file="export/plots/refC99_RF-cMDS-DAPCoptK.png",width=5,height=5,units="cm",res=200,pointsize=5)
par(mar=c(0.1,0.1,0.1,0.1),mai=c(0.2,0.2,0.2,0.1))
plot(cmdsplot2_refC99$points,main="RF cMDS: DAPC optimal K",xlab="",ylab="",
     col=results_refC99$grp, pch=16,cex=1.1,cex.axis=0.8,cex.main=1)
dev.off()

#### 2. PAM clustering on RF proximity scores ----
# Test mean silhouette width of k, select optimal k with highest silhouette width
optimumK <- data.frame(valueK="k=1")
for (i in 2:40) {
  RF_pam_clust_prox_refC99 <- pam(rftest_refC99$proximity, i) 
  output <- paste0("k=",i,": ",mean(silhouette(RF_pam_clust_prox_refC99)[, "sil_width"])) -> optimumK[i,]
  print(output)
}
optimumK %>% 
  mutate(K=as.numeric(str_extract(valueK,'k=([0-9]+): ([0-9].*)',group=1)),
         value=as.numeric(str_extract(valueK,'(k=.*) ([0-9].*)',group=2))) %>% 
  subset(value==max(value,na.rm=T)) %>% dplyr::select(K) -> optimumK
RF_pam_clust_prox_refC99 <- pam(rftest_refC99$proximity, optimumK)  
# cMDS with optimal k of DAPC and clusters via PAM
png(file="export/plots/refC99_RF-cMDS-PAM-prox.png",width=5,height=5,units="cm",res=200,pointsize=5)
par(mar=c(0.1,0.1,0.1,0.1),mai=c(0.2,0.2,0.2,0.1))
plot(cmdsplot2_refC99$points,main="RF cMDS: PAM clusters on RF proximity scores",xlab="",ylab="",#xlab="Scaling Dimension 1", ylab="Scaling Dimension 2",
     col=RF_pam_clust_prox_refC99$clustering, pch=16,cex=1.1,cex.axis=0.8,cex.main=1)
dev.off()
#add to cluster result compilation
DAPCRFtSNEresults_refC99 <- cbind(DAPCRFtSNEresults_refC99,RFc_PAMprox=RF_pam_clust_prox_refC99$clustering)

#### 3. PAM clustering on cMDS output  ----
# Test mean silhouette width of k, select optimal k with highest value
optimumK <- data.frame(valueK="k=1")
for (i in 2:40) {
  RF_pam_clust_cMDS_refC99 <- pam(cmdsplot1_refC99$points, i)
  output <- paste0("k=",i,": ",mean(silhouette(RF_pam_clust_cMDS_refC99)[, "sil_width"]))-> optimumK[i,]
  print(output)
}
optimumK %>% 
  mutate(K=as.numeric(str_extract(valueK,'k=([0-9]+): ([0-9].*)',group=1)),
         value=as.numeric(str_extract(valueK,'(k=.*) ([0-9].*)',group=2))) %>% 
  subset(value==max(value,na.rm=T)) %>% dplyr::select(K) -> optimumK
RF_pam_clust_cMDS_refC99 <- pam(cmdsplot1_refC99$points, optimumK)
png(file="export/plots/refC99_RF-cMDS-PAM-cMDS.png",width=5,height=5,units="cm",res=200,pointsize=5)
par(mar=c(0.1,0.1,0.1,0.1),mai=c(0.2,0.2,0.2,0.1))
plot(cmdsplot2_refC99$points,main="RF cMDS: PAM clusters on cMDS",xlab="",ylab="", # xlab="Scaling Dimension 1", ylab="Scaling Dimension 2", 
     col=RF_pam_clust_cMDS_refC99$clustering, pch=16,cex=1.1,cex.axis=0.8,cex.main=1)
dev.off()
#add to cluster result compilation
DAPCRFtSNEresults_refC99 <- cbind(DAPCRFtSNEresults_refC99, RFc_PAMcMDS=RF_pam_clust_cMDS_refC99$clustering)

#### 4. Determine optimal k from cMDS using gap statistic with PAM clusters ----
# adjust iter.max if it doesn't converge. Also can adjust k.max
cmds_nbclust_refC99 <- 
  fviz_nbclust(cmdsplot1_refC99$points,kmeans,nstart=25,method="gap_stat",nboot=1000,k.max=40,iter.max=50)+
  labs(subtitle="Gap statistic method")
cmds_nbclust_refC99   #K=15
cmds_nbclust_k_refC99 <- cmds_nbclust_refC99[["layers"]][[4]][["data"]][["xintercept"]]
# pam clustering with optimal k from gap statistic
cmds_nbclust_clust_refC99 <- pam(cmdsplot1_refC99$points, cmds_nbclust_k_refC99)
# cMDS with optimal k of RF via gap statistic and clusters via PAM 
png(file="export/plots/refC99_RF-cMDS-GapStatistic.png",width=5,height=5,units="cm",res=200,pointsize=5)
par(mar=c(0.1,0.1,0.1,0.1),mai=c(0.2,0.2,0.2,0.1))
plot(cmdsplot2_refC99$points,main="RF cMDS: Gap statistic optimal K and clusters",xlab="",ylab="", # xlab="Scaling Dimension 1", ylab="Scaling Dimension 2", 
     col=cmds_nbclust_clust_refC99$clustering, pch=16,cex=1.1,cex.axis=0.8,cex.main=1)
dev.off()
#add to cluster result compilation
DAPCRFtSNEresults_refC99 <- cbind(DAPCRFtSNEresults_refC99, RFc_GS=cmds_nbclust_clust_refC99$clustering)

#### 5. Determine optimal k from cMDS via hierarchical clustering with BIC ----
# adjust G option to reasonable potential cluster values, e.g. for up to 12 clusters, G=1:12
cmdsplot_clust_refC99 <- Mclust(cmdsplot2_refC99$points,G=1:40)
mclust_grps_cmdsplot_refC99 <- as.numeric(cmdsplot_clust_refC99$classification)
max(mclust_grps_cmdsplot_refC99)  #k=7
# cMDS with optimal k and clusters of RF via hierarchical clustering
png(file="export/plots/refC99_RF-cMDS-HierarchicalClustering.png",width=5,height=5,units="cm",res=200,pointsize=5)
par(mar=c(0.1,0.1,0.1,0.1),mai=c(0.2,0.2,0.2,0.1))
plot(cmdsplot2_refC99$points,main="RF cMDS: Hierarchical clustering optimal K and clusters",xlab="",ylab="", # xlab="Scaling Dimension 1", ylab="Scaling Dimension 2", 
     col=mclust_grps_cmdsplot_refC99, pch=16,cex=1.1,cex.axis=0.8,cex.main=1)
dev.off()
#add to cluster result compilation
DAPCRFtSNEresults_refC99 <- cbind(DAPCRFtSNEresults_refC99, RFc_HC=mclust_grps_cmdsplot_refC99)

### Isotonic MDS ---------------------------------------------------------
#In cMDS, you are given the matrix of distances between the objects with a goal to figure out the locations of these objects in 1D/2D/3D++space. 
#In non-metric MDS (isotonicMDS is one of them), you only know that objects 1 and 2 are more distant than objects 2 and 3, so you try to quantify that, on top of finding the dimensions and locations.
isomdsplot_refC99 <- MASS::isoMDS(1-rftest_refC99$proximity)

#### 1. isoMDS ----
# Use optimal k from DAPC 
png(file="export/plots/refC99_RF-iMDS-DAPCoptK.png",width=5,height=5,units="cm",res=200,pointsize=5)
par(mar=c(0.1,0.1,0.1,0.1),mai=c(0.2,0.2,0.2,0.1))
plot(isomdsplot_refC99$points,main="RF isoMDS: DAPC optimal K",xlab="",ylab="", # xlab="Scaling Dimension 1", ylab="Scaling Dimension 2", 
     col=results_refC99$grp, pch=16,cex=1.1,cex.axis=0.8,cex.main=1)
dev.off()

#### 2. PAM clustering on iMDS output ----
# Test mean silhouette width of k, select optimal k with highest value
optimumK <- data.frame(valueK="k=1")
for (i in 2:40) {
  RF_pam_clust_iso_refC99 <- pam(isomdsplot_refC99$points, i)
  output <- paste0("k=",i,": ",mean(silhouette(RF_pam_clust_iso_refC99)[, "sil_width"])) -> optimumK[i,]
  print(output)
}
optimumK %>% 
  mutate(K=as.numeric(str_extract(valueK,'k=([0-9]+): ([0-9].*)',group=1)),
         value=as.numeric(str_extract(valueK,'(k=.*) ([0-9].*)',group=2))) %>% 
  subset(value==max(value,na.rm=T)) %>% dplyr::select(K) -> optimumK
RF_pam_clust_iso_refC99 <- pam(isomdsplot_refC99$points, optimumK)
png(file="export/plots/refC99_RF-iMDS-PAM-iMDS.png",width=5,height=5,units="cm",res=200,pointsize=5)
par(mar=c(0.1,0.1,0.1,0.1),mai=c(0.2,0.2,0.2,0.1))
plot(isomdsplot_refC99$points,main="RF isoMDS: PAM clusters on isoMDS",xlab="",ylab="", # xlab="Scaling Dimension 1", ylab="Scaling Dimension 2", 
     col=RF_pam_clust_iso_refC99$clustering,pch=16,cex=1.1,cex.axis=0.8,cex.main=1)
dev.off()
#add to cluster result compilation
DAPCRFtSNEresults_refC99 <- cbind(DAPCRFtSNEresults_refC99, RFi_PAMiMDS=RF_pam_clust_iso_refC99$clustering)

#### 3. Determine optimal k using gap statistic ----
# can adjust k.max
isomds_nbclust_refC99 <-
  fviz_nbclust(isomdsplot_refC99$points,kmeans,nstart=25,method="gap_stat",nboot=1000,k.max=40)+
  labs(subtitle="Gap statistic method")
isomds_nbclust_refC99  #K=4
isomds_nbclust_k_refC99 <- isomds_nbclust_refC99[["layers"]][[4]][["data"]][["xintercept"]]
# pam clustering with optimal k from gap statistic
isomds_nbclust_clust2_refC99 <- pam(rftest_refC99$proximity, isomds_nbclust_k_refC99)
# isoMDS with optimal k of RF via gap statistic and clusters via PAM
png(file="export/plots/refC99_RF-iMDS-GapStatistic.png",width=5,height=5,units="cm",res=200,pointsize=5)
par(mar=c(0.1,0.1,0.1,0.1),mai=c(0.2,0.2,0.2,0.1))
plot(isomdsplot_refC99$points,main="RF isoMDS: Gap statistic optimal K and clusters",xlab="",ylab="",#xlab="Scaling Dimension 1", ylab="Scaling Dimension 2", 
     col=isomds_nbclust_clust2_refC99$clustering,pch=16,cex=1.1,cex.axis=0.8,cex.main=1)
dev.off()
#add to cluster result compilation
DAPCRFtSNEresults_refC99 <- cbind(DAPCRFtSNEresults_refC99, RFi_GS=isomds_nbclust_clust2_refC99$clustering)

#### 4. Determine optimal k of RF via hierarchical clustering with BIC ----
# adjust G option to reasonable potential cluster values, e.g. for up to 12 clusters, G=1:12
isomdsplot_clust_refC99 <- Mclust(isomdsplot_refC99$points,G=1:40)
mclust_grps_isomdsplot2_refC99 <- as.numeric(isomdsplot_clust_refC99$classification)
max(mclust_grps_isomdsplot2_refC99)  #5
# isoMDS with optimal k and clusters of RF via hierarchical clustering
png(file="export/plots/refC99_RF-iMDS-HierarchicalClustering.png",width=5,height=5,units="cm",res=200,pointsize=5)
par(mar=c(0.1,0.1,0.1,0.1),mai=c(0.2,0.2,0.2,0.1))
plot(isomdsplot_refC99$points, main="RF isoMDS: Hierarchical clustering optimal K and clusters", # xlab="Scaling Dimension 1", ylab="Scaling Dimension 2",
     col=mclust_grps_isomdsplot2_refC99,pch=16,cex=1.1,cex.axis=0.8,cex.main=1)
dev.off()
#add to cluster result compilation
DAPCRFtSNEresults_refC99 <- cbind(DAPCRFtSNEresults_refC99, RFi_HC=mclust_grps_isomdsplot2_refC99)
#Merge .pngs
system("montage export/plots/refC99_RF* -geometry +5+5 export/refC99_RF.png")

###################################################################

# t-SNE ----
# Libraries
library(tidyverse)
library(cluster)
library(factoextra) # detach("package:adegenet",unload=TRUE) if throws ggplot error
library(mclust)
library(tsne)

## Prepare plot labels and colours ----
colors_refC99 = rainbow(length(unique(data_refC99$pop)))
names(colors_refC99) = unique(data_refC99$pop)
ecb_refC99 = function(x,y){plot(x,t='n'); text(x, labels=data_refC99$pop, col=colors_refC99[data_refC99$pop])}

# t-SNE on principal components of scaled data 
#**TAKES TIME (~10min) but this is where you can see the clustering plot happen in real time
# adjust perplexity(originally 5), initial_dims(originally 5), max_iter(originally 5000)
dev.off()
tsne_refC99_p5 = tsne(pca_refC99$tab, epoch_callback=ecb_refC99, max_iter=15000, perplexity=5, initial_dims=5)

### 1. tSNE plot with DAPC groups ----
png(file="export/plots/refC99_tSNE-p5_DAPCoptK.png",width=5,height=5,units="cm",res=200,pointsize=5)
par(mar=c(0.1,0.1,0.1,0.1),mai=c(0.2,0.2,0.2,0.1))
plot(tsne_refC99_p5,main="tSNE p5: DAPC optimal K",xlab="",ylab="",
     col=results_refC99$grp,pch=16,cex=1.1,cex.axis=0.8,cex.main=1)
dev.off()

### 2. PAM clustering with optimal k from DAPC ----
#Test mean silhouette width of k, select optimal k with highest value
optimumK <- data.frame(valueK="k=1")
for (i in 2:40) {
  pam_clust2_tsne_refC99 <- pam(tsne_refC99_p5, i) 
  output <- paste0("k=",i,": ",mean(silhouette(pam_clust2_tsne_refC99)[, "sil_width"])) -> optimumK[i,]
  print(output)
}
optimumK %>% 
  mutate(K=as.numeric(str_extract(valueK,'k=([0-9]+): ([0-9].*)',group=1)),
         value=as.numeric(str_extract(valueK,'(k=.*) ([0-9].*)',group=2))) %>% 
  subset(value==max(value,na.rm=T)) %>% dplyr::select(K) -> optimumK
pam_clust2_tsne_refC99 <- pam(tsne_refC99_p5, optimumK)
png(file="export/plots/refC99_tSNE-p5_PAM.png",width=5,height=5,units="cm",res=200,pointsize=5)
par(mar=c(0.1,0.1,0.1,0.1),mai=c(0.2,0.2,0.2,0.1))
plot(tsne_refC99_p5,main="tSNE p5: PAM clusters",xlab="",ylab="",
     col=pam_clust2_tsne_refC99$clustering,pch=16,cex=1.1,cex.axis=0.8,cex.main=1)
dev.off()
#add to cluster result compilation
DAPCRFtSNEresults_refC99 <- cbind(DAPCRFtSNEresults_refC99, tSNE_PAM=pam_clust2_tsne_refC99$clustering)

### 3. Clustering for perplexity=5, determine optimal k using gap statistic ----
# can adjust k.max
tsne_refC99_p5_nbclust <- 
  fviz_nbclust(tsne_refC99_p5,kmeans,nstart=25,method="gap_stat",nboot=1000,k.max=40)+
  labs(subtitle="Gap statistic method")
tsne_refC99_p5_nbclust   #k=9
tsne_refC99_p5_nbclust_k <- tsne_refC99_p5_nbclust[["layers"]][[4]][["data"]][["xintercept"]]
# pam clustering with optimal k from gap statistic
tsne_refC99_p5_nbclust_clust <- pam(tsne_refC99_p5, tsne_refC99_p5_nbclust_k)
# t-SNE with optimal k of RF via gap statistic and clusters via PAM
png(file="export/plots/refC99_tSNE-p5_GapStatistic.png",width=5,height=5,units="cm",res=200,pointsize=5)
par(mar=c(0.1,0.1,0.1,0.1),mai=c(0.2,0.2,0.2,0.1))
plot(tsne_refC99_p5, main="tSNE p5: Gap statistic optimal K and clusters",xlab="",ylab="",
     col=tsne_refC99_p5_nbclust_clust$clustering,pch=16,cex=1.1,cex.axis=0.8,cex.main=1)
dev.off()
#add to cluster result compilation
DAPCRFtSNEresults_refC99 <- cbind(DAPCRFtSNEresults_refC99, tSNE_GS=tsne_refC99_p5_nbclust_clust$clustering)

### 4. Determine optimal k via hierarchical clustering with BIC ----
# adjust G option to reasonable potential cluster values, e.g. for up to 12 clusters, G=1:12
tsne_refC99_p5_clust <- Mclust(tsne_refC99_p5,G=1:40)
mclust_grps_tsne_refC99_p5 <- as.numeric(tsne_refC99_p5_clust$classification)
max(mclust_grps_tsne_refC99_p5)  #6
# t-SNE p5 with optimal k and clusters of RF via hierarchical clustering
png(file="export/plots/refC99_tSNE-p5_HierarchicalClustering.png",width=5,height=5,units="cm",res=200,pointsize=5)
par(mar=c(0.1,0.1,0.1,0.1),mai=c(0.2,0.2,0.2,0.1))
plot(tsne_refC99_p5,main="t-SNE p5: Hierarchical clustering optimal K and clusters",xlab="",ylab="",
     col=mclust_grps_tsne_refC99_p5,pch=16,cex=1.1,cex.axis=0.8,cex.main=1)
dev.off()
#add to cluster result compilation
DAPCRFtSNEresults_refC99 <- cbind(DAPCRFtSNEresults_refC99, tSNE_HC=mclust_grps_tsne_refC99_p5)

#Merge .pngs
system("montage export/plots/refC99_tSNE-p5* -geometry +5+5 export/refC99_tSNE-p5.png")
#Write out .csv
write.csv(DAPCRFtSNEresults_refC99,file="export/DAPC-RF-tSNEresults_refC99.csv",row.names=TRUE)
