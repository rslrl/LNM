library(psych)
library(corpcor)
library(parallel)
library(stringr)
library(data.table)
library(edgeR)
library(limma)
library(tidyverse)

#input and output file download address(Additional file 3): http://bclab.inha.ac.kr/LNM/

cancer='COAD'
# compute MIC between normal and tumor samples
e3=paste('/mnt/d/InCOB_ren/input/',cancer,'_normalized_count.csv',sep='')
Counts=fread(e3,header=T,check.names=F)
Counts=data.frame(Counts,check.names=F,stringsAsFactors = FALSE)
RNA_type=Counts[,1:2]
print(paste("#miRNAs: ", nrow(Counts[Counts$type=='miRNA',]),sep=""))
print(paste("#mRNAs: ", nrow(Counts[Counts$type=='mRNA',]),sep=""))
print(paste("#lncRNAs: ", nrow(Counts[Counts$type=='lncRNA',]),sep=""))
print(paste("#pseudogenes: ", nrow(Counts[Counts$type=='pseudogene',]),sep=""))

miRNA=Counts[Counts$type=='miRNA',]
RNA=Counts[Counts$type%in%c('mRNA','pseudogene','lncRNA'),]
RNA=RNA[,-which(colnames(RNA)=='type')]
miRNA=miRNA[,-which(colnames(miRNA)=='type')]
RNA=avereps(RNA[,-which(colnames(RNA)=='name')],ID=RNA$name)
miRNA=avereps(miRNA[,-which(colnames(miRNA)=='name')],ID=miRNA$name)
Counts=rbind(RNA,miRNA)
miRNA=t(miRNA)
RNA=t(RNA)

miRNA_name=colnames(miRNA)
mi_num=length(miRNA_name)
n1=length(which(str_sub(colnames(Counts),start = 14, end = 15)=='11'))
t1=length(which(str_sub(colnames(Counts),start = 14, end = 15)=='01'))

print(paste("#normal samples: ", n1,sep=""))
print(paste("#tumor samples: ", t1,sep=""))

e2=paste('/mnt/d/InCOB_ren/input/TCGA-',cancer,'-clinical.csv',sep='')
clinical=fread(e2)
clinical=clinical[,c('submitter_id','ajcc_pathologic_stage','ajcc_pathologic_m','ajcc_pathologic_t','ajcc_pathologic_n')]
colnames(clinical)[1]='id'

m_index=unique(clinical$ajcc_pathologic_m)
n_index=unique(clinical$ajcc_pathologic_n)
t_index=unique(clinical$ajcc_pathologic_t)
m_index=m_index[which(str_sub(m_index,1,2)!='MX')]
n_index=n_index[which(str_sub(n_index,1,2)!='NX')]
t_index=t_index[which(str_sub(t_index,1,2)!='TX')]

clinical=clinical[clinical$ajcc_pathologic_m%in%m_index,]
clinical=clinical[clinical$ajcc_pathologic_n%in%n_index,]
clinical=clinical[clinical$ajcc_pathologic_t%in%t_index,]

metastatic=clinical[which(str_sub(clinical$ajcc_pathologic_n,1,2)!='N0'),] 
metastatic=metastatic[which(str_sub(metastatic$ajcc_pathologic_t,1,2)!='T0'),] 
metastatic=metastatic[which(str_sub(metastatic$ajcc_pathologic_m,1,2)!='M1'),]

non_metastatic=clinical[which(str_sub(clinical$ajcc_pathologic_n,1,2)=='N0'),] 
non_metastatic=non_metastatic[which(str_sub(non_metastatic$ajcc_pathologic_t,1,2)!='T0'),] 
non_metastatic=non_metastatic[which(str_sub(non_metastatic$ajcc_pathologic_m,1,2)!='M1'),] 

print(paste("#lymph node metastasis: ", length(which(str_sub(colnames(Counts)[(n1+1):ncol(Counts)],1,12)%in%metastatic$id)),sep=""))
print(paste("#non-metastasis: ", length(which(str_sub(colnames(Counts)[(n1+1):ncol(Counts)],1,12)%in%non_metastatic$id)),sep=""))


func <- function(i) {
    n2=n1+1
    nt=n1+t1
    RNA_name=rep(colnames(RNA)[i],mi_num)
    x=as.matrix(RNA[1:n1,i])
    y=as.matrix(miRNA[1:n1,])
    mic_n1 <- as.numeric(cstats(x, y, alpha=0.6, C=15, est="mic_approx")[,3]) 

    n2=n1+1
    nt=n1+t1
    x=as.matrix(RNA[n2:nt,i])
    y=as.matrix(miRNA[n2:nt,])
    mic_t1 <- as.numeric(cstats(x, y, alpha=0.6, C=15, est="mic_approx")[,3])

    f1=data.frame(miRNA=miRNA_name,RNA=RNA_name,
                  mic_nor=mic_n1,
                  mic_tum=mic_t1,
                  stringsAsFactors = FALSE)
    return(f1)
}

i <- 1:ncol(RNA)
cl <- makeCluster(35) # Number of cpu cores
clusterExport(cl,c("RNA","miRNA",'n1','t1','mi_num','miRNA_name'))
clusterEvalQ(cl,library(minerva))
results <- parLapply(cl,i,func) 
res.df <- do.call('rbind',results) 
stopCluster(cl)

# Select miRNA_RNA pairs with MIC greater than 0.5 in normal samples or tumor samples
res.df=res.df[!duplicated(res.df),]
mic_nor_tum_0.5=res.df[abs(res.df$mic_nor)>0.5|abs(res.df$mic_tum)>0.5,]
print(paste("#features after MIC filtering : ", nrow(mic_nor_tum_0.5),sep=""))

# compute delta PCC for tumor samples
a1=c(1:n1)
func <- function(i) {
    col1=t1+1
    cox_all=matrix(nrow=1,ncol=col1)

    RNA1= as.character(mic_nor_tum_0.5[i,]$RNA)
    miRNA1= as.character(mic_nor_tum_0.5[i,]$miRNA)
    s=cbind(scale(Counts[RNA1,1:n1]),scale(Counts[miRNA1,1:n1]))
    s=s[abs(s[,1])<3&abs(s[,2])<3,]
    cor_s=cor(s)
    cox_all[1,1]=cor_s[2,1]

    for (i1 in seq(t1)){
        index=i1+1
        n1_1=n1+1
        a1[n1_1]=n1+i1
        
        s=cbind(scale(Counts[RNA1,a1]),scale(Counts[miRNA1,a1]))
        s1=s[abs(s[1:n1,1])<3&abs(s[1:n1,2])<3,]
        s1=rbind(s1,s[n1_1,])
        cor_s=cor(s1)
        
        cox_all[1,index]=cor_s[2,1]
    }
    return(cox_all)
}
i <- 1:dim(mic_nor_tum_0.5)[1]
cl <- makeCluster(35) 
clusterExport(cl,c("mic_nor_tum_0.5","Counts",'a1','n1','t1'))
results <- parLapply(cl,i,func) 
res.df <- do.call('rbind',results) 
stopCluster(cl)

res.df=data.frame(res.df)
n1_1=n1+1
ncol1=ncol(Counts)
colnames(res.df)=c('normal',colnames(Counts)[n1_1:ncol1])
res.df$pair=paste(paste(mic_nor_tum_0.5[,2], mic_nor_tum_0.5[,1], sep="_"))
res.df=data.frame(res.df,check.names = FALSE)

for (i in seq(t1)){
    i1=i+1
    res.df[,i1]=res.df[,i1]-res.df[,1]
}


#Divide samples into lymph node metastasis group and non-metastatic group based on clinical file.
sample=c('pair',
  colnames(res.df)[which(str_sub(colnames(res.df),1,12)%in%metastatic$id)],
  colnames(res.df)[which(str_sub(colnames(res.df),1,12)%in%non_metastatic$id)])
res.df=res.df[,sample]
res.df=res.df[!duplicated(res.df),]

m1=length(which(str_sub(colnames(res.df),1,12)%in%metastatic$id))+1
m1_1=m1+1
n_m1=ncol(res.df)
func <- function(i) {
    likes=as.numeric(res.df[i,2:m1])
    likes2=as.numeric(res.df[i,m1_1:n_m1])
    Pvalue=tryCatch({
       Pvalue=wilcox.test(likes,likes2)$p.value
    }, warning = function(w){
          return(2)
    }, 
       error = function(e) {
        return(3)
    }
    )
    return(Pvalue)
}
i <- 1:nrow(res.df)
cl <- makeCluster(35) 
clusterExport(cl,c("res.df",'m1','m1_1','n_m1'))
results <- parLapply(cl,i,func) 
res.df_1 <- do.call('cbind',results) 
stopCluster(cl)

# Select the pairs with p-value less than 0.01 in wilcox test
res.df$pvalue=as.vector(res.df_1)
res.df=res.df[res.df$pvalue<0.01,]
res.df=res.df[,-which(colnames(res.df)=='pvalue')]
print(paste("#features after Wilcox test : ", nrow(res.df),sep=""))

# Extract miRNA_RNA pairs to construct a ceRNA network
pairs=separate(res.df,pair, into= c('name',"miRNA"),sep= "_")
pairs=pairs[,1:2]
pairs=merge(pairs, RNA_type, all = FALSE,by = "name")
pairs=pairs[,c('miRNA','name','type')]
colnames(pairs)=c('miRNA','RNA','RNA_type')

e3=paste('/mnt/d/InCOB_ren/',cancer,'_ceRNA_network_Cytoscape.csv',sep='')
write.csv(pairs,e3,row.names=FALSE)


## Label the tumor samples
# label=c('pair',
#         rep(1,length(which(str_sub(colnames(res.df),1,12)%in%metastatic$id))),
#         rep(0,length(which(str_sub(colnames(res.df),1,12)%in%non_metastatic$id)))
#         )
# res.df=rbind(res.df,label)
