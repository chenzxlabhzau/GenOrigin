---
title: "Geneage"
author: "Sheng Qian"
date: "2019年8月7日"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE,cache = TRUE, highlight = TRUE)
```

##### Gene age analysis
```{r}
human <- read.csv("/home/chenyj/qiansheng/Homo_sapiens.csv",header = TRUE,sep = ",",stringsAsFactors = F)
human <- human[-1,]
ggplot(human,aes(gene_age))+geom_density()+Theme+
  xlab("Gene Age (myr)")+ ylab("Distribution")+
  labs(title="Human")+theme(plot.title = element_text(hjust = 0.5))
mouse <- read.csv("/home/chenyj/qiansheng/Mus_musculus.csv",header = TRUE,sep = ",",stringsAsFactors = F)
mouse <- mouse[-1,]
ggplot(mouse,aes(gene_age))+geom_density()+Theme+
  xlab("Gene Age (myr)")+ ylab("Distribution")+
  labs(title="Mouse")+theme(plot.title = element_text(hjust = 0.5))
#coord_cartesian(xlim = c(0,250))
```

```{bash}
#cd /home/qians/
#mkdir geneage
#cd geneage/
#cp /home/chenyj/qiansheng/*csv ./
for i in `ls *csv`; do head -2 $i| tail -n 1; done |sort| uniq -c #before delete the first two lines, make sure it 
sed -i '1d' *csv
sed -i '1d' *csv
cat *csv > all_species
```
```{r}
all <- read.csv("/home/qians/geneage/all_species",header = FALSE,sep = ",",stringsAsFactors = FALSE)
```

```{r}
ggplot(all,aes(V2))+geom_density()+Theme+
  xlab("Gene Age (myr)")+ ylab("Distribution")+
  labs(title="All species")+theme(plot.title = element_text(hjust = 0.5))
all$V3 <- 1
# for (i in 1:nrow(all)) {
#   if (all[i,2]>= 50) {
#     all[i,3] = 100
#   }else if (all[i,2]>= 100) {
#     all[i,3] = 150
#   }
# }
all[all$V2>=  0,3] = 5
all[all$V2>=  5,3] = 10
all[all$V2>= 10,3] = 15
all[all$V2>= 15,3] = 20
all[all$V2>= 20,3] = 25
all[all$V2>= 25,3] = 30
all[all$V2>= 30,3] = 35
all[all$V2>= 35,3] = 40
all[all$V2>= 40,3] = 45
all[all$V2>= 45,3] = 50
all[all$V2>= 50,3] = 100
all[all$V2>=100,3] = 150
all[all$V2>=150,3] = 200
all[all$V2>=200,3] = 250
all[all$V2>=250,3] = 300
all[all$V2>=300,3] = 350
all[all$V2>=350,3] = 400
all[all$V2>=400,3] = 450
all[all$V2>=450,3] = 500
all[all$V2>=500,3] = 550
all[all$V2>=550,3] = 600
all[all$V2>=600,3] = 650
all[all$V2>=650,3] = 700
all[all$V2>=700,3] = 750
all[all$V2>=750,3] = 800
all[all$V2>=800,3] = ">800"
age <- table(all$V3) %>% as.data.frame() 
age$Var1 <- factor(age$Var1,levels = c(sort(c(seq(5,50,by = 5),seq(100,800,by = 50))),">800"))
colnames(age) <- c("age","Number")
ggplot(age,aes(age,log(Number)))+geom_point()+Theme+xlab("Gene Age (myr)")+ylab("Gene Number (#log)")
write.table(age,file = "/home/qians/geneage/all_species_age",sep = '\t',row.names = FALSE,quote = FALSE)
```

```{r}
#Improved date method
human <- read.csv("/home/qians/geneage/Book2.csv",header = FALSE,sep = ",",stringsAsFactors = F)
ggplot(human,aes(V2))+geom_density()+Theme+
  xlab("Gene Age (myr)")+ ylab("Distribution")+
  labs(title="Human")+theme(plot.title = element_text(hjust = 0.5))
```
###### Conparison between GenOrigin and Gentree
```{r human gene}
h1 <- read.csv("/home/G20/geneAge/all_species_gene_age/Homo_sapiens.csv",header = FALSE,sep = ",",stringsAsFactors = F)
colnames(h1) <- c("ensembl_id","gene_age","gene_Interval")
h1$start <- lapply(h1$gene_Interval,function(x)unlist(strsplit(x,"-")[[1]][1]))
h1$end <- lapply(h1$gene_Interval,function(x)unlist(strsplit(x,"-")[[1]][2]))
hg <- read.csv("gentree_human_age.csv",header = T,sep = ",",stringsAsFactors = F)
h <- h1[,c(1,2,4,5)]
h[,-1] <- lapply(h[,-1],as.numeric)
h1$tree_age <- hg[match(h1$ensembl_id,hg$gene),"age"] 
h1 <- na.omit(h1[,-3])
h1 <- na.omit(h1)
h1$tmid <- lapply(h1$tree_age,function(x)unlist(strsplit(x,"/")[[1]][1])) %>% unlist() %>% as.numeric()
h1$trange <- lapply(h1$tree_age,function(x)unlist(strsplit(x,"/")[[1]][2])) %>% unlist()
h1$ts <- lapply(h1$trange,function(x)unlist(strsplit(x,"-")[[1]][1])) %>% unlist() %>% as.numeric()
h1$te <- lapply(h1$trange,function(x)unlist(strsplit(x,"-")[[1]][2])) %>% unlist() %>% as.numeric()
table(h1$tmid==h1$ts)
plot(h$tmid,h$gene_age,xlab="Gene Age in Gentree",ylab="Gene Age in GenOrigin")
#t <- factor(h$ts)
#levels(t) <- paste(rep("b",14),1:length(unique(t)),sep = "")
#h$group <- t
h1 <- h1[order(h1$ts),]
h1$trange <- factor(h1$trange,levels = unique(h1$trange)) #先排序 后指定level
cor.test(h1$V2,h1$tmid,method = "pearson")

ph <- ggplot2::ggplot(h,aes(trange,gene_age))+geom_boxplot(outlier.colour = "white")+
    xlab("Gene Age in Gentree")+ylab("Gene Age in GenOrigin")+
    theme(axis.text.x = element_text(angle = 45, hjust = 1))+
    ggtitle("Human gene")+ annotate("text",x=4,y=800,label="Pearson correlation = 0.61 \n P value < 2.2e-16",size=5)
```

```{r mouse gene}
m1 <- read.csv("/home/G20/geneAge/all_species_gene_age/Mus_musculus.csv",header = FALSE,sep = ",",stringsAsFactors = F)
colnames(m1) <- c("ensembl_id","gene_age","gene_Interval")
mg <- read.csv("gentree_mouse_age.csv",header = T,sep = ",",stringsAsFactors = F)
m1$start <- lapply(m1$gene_Interval,function(x)unlist(strsplit(x,"-")[[1]][1]))
m1$end <- lapply(m1$gene_Interval,function(x)unlist(strsplit(x,"-")[[1]][2]))
m1$tree_age <- mg[match(m1$ensembl_id,mg$gene),"age"]
m1 <- na.omit(m1[,-3])
m1$tmid <- lapply(m1$tree_age,function(x)unlist(strsplit(x,"/")[[1]][1])) %>% unlist() %>% as.numeric()
m1$trange <- lapply(m1$tree_age,function(x)unlist(strsplit(x,"/")[[1]][2])) %>% unlist()
m1$ts <- lapply(m1$trange,function(x)unlist(strsplit(x,"-")[[1]][1])) %>% unlist() %>% as.numeric()
m1$te <- lapply(m1$trange,function(x)unlist(strsplit(x,"-")[[1]][2])) %>% unlist() %>% as.numeric()
cor.test(m1$gene_age,m1$tmid,method = "pearson")
plot(m1$tmid,m1$gene_age,xlab="Gene Age in Gentree",ylab="Gene Age in GenOrigin")
m1 <- m1[order(m1$ts),]
m1$trange <- factor(m1$trange,levels = unique(m1$trange))
pm <- ggplot(m1,aes(trange,gene_age))+geom_boxplot(outlier.colour = "white")+
    xlab("Gene Age in Gentree")+ylab("Gene Age in GenOrigin")+
    theme(axis.text.x = element_text(angle = 45, hjust = 1))+
    ggtitle("Mouse gene")+ annotate("text",x=4,y=800,label="Pearson correlation = 0.72 \n P value < 2.2e-16",size=5)
```

```{r Fly}
f1 <- read.csv("/home/G20/geneAge/all_species_gene_age/Drosophila_melanogaster.csv",header = FALSE,sep = ",",stringsAsFactors = F)
colnames(f1) <- c("ensembl_id","gene_age","gene_Interval")
fg <- read.csv("gentree_fly_age.csv",header = T,sep = ",",stringsAsFactors = F)
f1$start <- lapply(f1$gene_Interval,function(x)unlist(strsplit(x,"-")[[1]][1])) %>% unlist() %>% as.numeric()
f1$end <- lapply(f1$gene_Interval,function(x)unlist(strsplit(x,"-")[[1]][2]))
f1$tree_age <- fg[match(f1$ensembl_id,fg$gene),"age"]
f1 <- na.omit(f1[,-3])
f1$tmid <- lapply(f1$tree_age,function(x)unlist(strsplit(x,"/")[[1]][1])) %>% unlist() %>% as.numeric()
f1$trange <- lapply(f1$tree_age,function(x)unlist(strsplit(x,"/")[[1]][2])) %>% unlist()
f1$ts <- lapply(f1$trange,function(x)unlist(strsplit(x,"-")[[1]][1])) %>% unlist() %>% as.numeric()
f1$te <- lapply(f1$trange,function(x)unlist(strsplit(x,"-")[[1]][2])) %>% unlist() %>% as.numeric()
cor.test(f1$gene_age,f1$tmid,method = "pearson")
f1 <- f1[order(f1$ts),]
f1$trange <- factor(f1$trange,levels = unique(f1$trange))
pf <- ggplot(f1,aes(trange,gene_age))+geom_boxplot(outlier.colour = "white")+
    xlab("Gene Age in Gentree")+ylab("Gene Age in GenOrigin")+
    theme(axis.text.x = element_text(angle = 45, hjust = 1))+ggtitle("Fly gene")
```

```{r Chicken}
c1 <- read.csv("/home/G20/geneAge/all_species_gene_age/Gallus_gallus.csv",header = FALSE,sep = ",",stringsAsFactors = F)
colnames(c1) <- c("ensembl_id","gene_age","gene_Interval")
cg <- read.csv("gentree_chicken_age.csv",header = T,sep = ",",stringsAsFactors = F)
c1$start <- lapply(c1$gene_Interval,function(x)unlist(strsplit(x,"-")[[1]][1])) %>% unlist() %>% as.numeric()
c1$end <- lapply(c1$gene_Interval,function(x)unlist(strsplit(x,"-")[[1]][2]))
c1$tree_age <- cg[match(c1$ensembl_id,cg$gene),"age"]
c1 <- na.omit(c1[,-3])
c1$tmid <- lapply(c1$tree_age,function(x)unlist(strsplit(x,"/")[[1]][1])) %>% unlist() %>% as.numeric()
c1$trange <- lapply(c1$tree_age,function(x)unlist(strsplit(x,"/")[[1]][2])) %>% unlist()
c1$ts <- lapply(c1$trange,function(x)unlist(strsplit(x,"-")[[1]][1])) %>% unlist() %>% as.numeric()
c1$te <- lapply(c1$trange,function(x)unlist(strsplit(x,"-")[[1]][2])) %>% unlist() %>% as.numeric()
cor.test(c1$gene_age,c1$tmid,method = "pearson")
c1 <- c1[order(c1$ts),]
c1$trange <- factor(c1$trange,levels = unique(c1$trange))
pc <- ggplot(c1,aes(trange,gene_age))+geom_boxplot(outlier.colour = "white")+
    xlab("Gene Age in Gentree")+ylab("Gene Age in GenOrigin")+
    theme(axis.text.x = element_text(angle = 45, hjust = 1))+ggtitle("Chicken gene")+
    annotate("text",x=2,y=1000,label="Pearson correlation = 0.33, P value < 2.2e-16",size=5)
```

```{r Rat}
r1 <- read.csv("/home/G20/geneAge/all_species_gene_age/Rattus_norvegicus.csv",header = FALSE,sep = ",",stringsAsFactors = F)
colnames(r1) <- c("ensembl_id","gene_age","gene_Interval")
rg <- read.csv("gentree_rat_age.csv",header = T,sep = ",",stringsAsFactors = F)
r1$start <- lapply(r1$gene_Interval,function(x)unlist(strsplit(x,"-")[[1]][1])) %>% unlist() %>% as.numeric()
r1$end <- lapply(r1$gene_Interval,function(x)unlist(strsplit(x,"-")[[1]][2]))
r1$tree_age <- rg[match(r1$ensembl_id,rg$gene),"age"]
r1 <- na.omit(r1[,-3])
r1$tmid <- lapply(r1$tree_age,function(x)unlist(strsplit(x,"/")[[1]][1])) %>% unlist() %>% as.numeric()
r1$trange <- lapply(r1$tree_age,function(x)unlist(strsplit(x,"/")[[1]][2])) %>% unlist()
r1$ts <- lapply(r1$trange,function(x)unlist(strsplit(x,"-")[[1]][1])) %>% unlist() %>% as.numeric()
r1$te <- lapply(r1$trange,function(x)unlist(strsplit(x,"-")[[1]][2])) %>% unlist() %>% as.numeric()
cor.test(r1$gene_age,r1$tmid,method = "pearson")
r1 <- r1[order(r1$ts),]
r1$trange <- factor(r1$trange,levels = unique(r1$trange))
pr <- ggplot(r1,aes(trange,gene_age))+geom_boxplot(outlier.colour = "white")+
    xlab("Gene Age in Gentree")+ylab("Gene Age in GenOrigin")+
    theme(axis.text.x = element_text(angle = 45, hjust = 1))+ggtitle("Rat gene")+
    annotate("text",x=2.5,y=850,label="Pearson correlation = 0.60 \n P value < 2.2e-16",size=5)
```

```{r Opposum}
o1 <- read.csv("/home/G20/geneAge/all_species_gene_age/Monodelphis_domestica.csv",header = FALSE,sep = ",",stringsAsFactors = F)
colnames(o1) <- c("ensembl_id","gene_age","gene_Interval")
og <- read.csv("gentree_opposum_age.csv",header = T,sep = ",",stringsAsFactors = F)
o1$start <- lapply(o1$gene_Interval,function(x)unlist(strsplit(x,"-")[[1]][1])) %>% unlist() %>% as.numeric()
o1$end <- lapply(o1$gene_Interval,function(x)unlist(strsplit(x,"-")[[1]][2]))
o1$tree_age <- og[match(o1$ensembl_id,og$gene),"age"]
o1 <- na.omit(o1[,-3])
o1$tmid <- lapply(o1$tree_age,function(x)unlist(strsplit(x,"/")[[1]][1])) %>% unlist() %>% as.numeric()
o1$trange <- lapply(o1$tree_age,function(x)unlist(strsplit(x,"/")[[1]][2])) %>% unlist()
o1$ts <- lapply(o1$trange,function(x)unlist(strsplit(x,"-")[[1]][1])) %>% unlist() %>% as.numeric()
o1$te <- lapply(o1$trange,function(x)unlist(strsplit(x,"-")[[1]][2])) %>% unlist() %>% as.numeric()
cor.test(o1$gene_age,o1$tmid,method = "pearson")
o1 <- o1[order(o1$ts),]
o1$trange <- factor(o1$trange,levels = unique(o1$trange))
po <- ggplot(o1,aes(trange,gene_age))+geom_boxplot(outlier.colour = "white")+
    xlab("Gene Age in Gentree")+ylab("Gene Age in GenOrigin")+
    theme(axis.text.x = element_text(angle = 45, hjust = 1))+ggtitle("Opposum gene")+
    annotate("text",x=2,y=1000,label="Pearson correlation = 0.53, P value < 2.2e-16",size=5)
```

```{r}
#human 
ha <- read.csv("/home/qians/geneage/annotation/human.csv",header = T,sep = "\t",stringsAsFactors = F)
h1$chr <- ha[match(h1$ensembl_id,ha$ensembl_gene_id),3]
h1$type <- ha[match(h1$ensembl_id,ha$ensembl_gene_id),2]
rm(ha)
#mouse
ma <- read.csv("/home/qians/geneage/annotation/mouse.csv",header = T,sep = "\t",stringsAsFactors = F)
m1$chr <- ma[match(m1$ensembl_id,ma$ensembl_gene_id),3]
m1$type <- ma[match(m1$ensembl_id,ma$ensembl_gene_id),2]
rm(ma)
#rat
ra <- read.csv("/home/qians/geneage/annotation/rat.csv",header = T,sep = "\t",stringsAsFactors = F)
r1$chr <- ra[match(r1$ensembl_id,ra$ensembl_gene_id),3]
r1$type <- ra[match(r1$ensembl_id,ra$ensembl_gene_id),2]
rm(ra)
#opposum
oa <- read.csv("/home/qians/geneage/annotation/opposum.csv",header = T,sep = "\t",stringsAsFactors = F)
o1$chr <- oa[match(o1$ensembl_id,oa$ensembl_gene_id),3]
o1$type <- oa[match(o1$ensembl_id,oa$ensembl_gene_id),2]
rm(oa)
#chicken
ca <- read.csv("/home/qians/geneage/annotation/chicken.csv",header = T,sep = "\t",stringsAsFactors = F)
c1$chr <- ca[match(c1$ensembl_id,ca$ensembl_gene_id),3]
c1$type <- ca[match(c1$ensembl_id,ca$ensembl_gene_id),2]
rm(ca)
```

```{r cowplot}
library(cowplot)
plot_grid(ph, pm, po, pr, pc, pf, labels = c("A", "B", "C", "D", "E", "F"), ncol = 2)
all <- rbind(h1,m1,o1,r1,c1)
all$group <- cut(all$tmid,breaks = seq(0,450,by = 50))
##all$group <- as.vector(all$group)
#all$group <- lapply(all$group,function(x)gsub("(","",x,fixed = T)) %>% unlist()
all$group <- lapply(all$group,function(x)gsub("(400,450]",">400",x,fixed = T)) %>% unlist()
all <- all[order(all$tmid),]
all$group <- factor(all$group,levels = unique(all$group))
cor.test(all$gene_age,all$tmid,method = "spearman")
ggplot(all,aes(group,gene_age))+geom_boxplot(outlier.colour = "white")+
  xlab("Gene Age in Gentree")+ylab("Gene Age in GenOrigin")+
  #theme(axis.text.x = element_text(angle = 45, hjust = 1))+Theme
  annotate("text",x=2,y=900,label="Spearman correlation = 0.57 \n P value < 2.2e-16",size=5)
```

```{r}
all$species <- lapply(all$ensembl_id,function(x)strsplit(x,"G0")[[1]][1]) %>% unlist()
table(all$chr=="X"|all$chr=="Z")/nrow(all)
young <- all[all$gene_age <50,]
nrow(young) / nrow(all)
table(young$chr=="X"|young$chr=="Z")/nrow(young)
table(all$chr=="X"|all$chr=="Z")/nrow(all)
table(all$chr=="X"|all$chr=="Z") %>% rbind(table(young$chr=="X"|young$chr=="Z")) %>% fisher.test()
table(young$type)
young[young$chr!="X"&young$chr!="Z","chr"]="A"
table(young$species,young$chr)
table(young$species)/table(all$species)
table(young[young$type!="protein_coding","species"],young[young$type!="protein_coding","chr"])
table(young[young$species=="ENSMOD","type"])
```

```{r}
#GO for young gene
##human
library(clusterProfiler)
library(org.Hs.eg.db,lib.loc = "/home/qians/R/x86_64-pc-linux-gnu-library/3.3")
bp <- enrichGO(young[young$species=="ENS",1], 
                OrgDb = org.Hs.eg.db, 
                keyType = "ENSEMBL", 
                ont = "BP",
                pvalueCutoff = 0.05
                )
bp_df <- as.data.frame(bp)
write.table(bp,file = "/home/G20/geneAge/all_species_gene_age/new/human_younggene_bp.csv",quote=FALSE,sep=",",row.names=FALSE)
barplot(bp, showCategory=100)
cc <- enrichGO(young[young$species=="ENS",1], 
                OrgDb = org.Hs.eg.db, 
                keyType = "ENSEMBL", 
                ont = "CC",
                pvalueCutoff = 0.05
                )
cc_df <- as.data.frame(cc)
write.table(cc_df,file = "/home/G20/geneAge/all_species_gene_age/new/human_younggene_cc.csv",quote=FALSE,sep=",",row.names=FALSE)
barplot(cc, showCategory=100)
mf <- enrichGO(young[young$species=="ENS",1], 
                OrgDb = org.Hs.eg.db, 
                keyType = "ENSEMBL", 
                ont = "MF",
                pvalueCutoff = 0.05
                )
mf_df <- as.data.frame(mf)
barplot(mf, showCategory=100)
write.table(mf_df,file = "/home/G20/geneAge/all_species_gene_age/new/human_younggene_mf.csv",quote=FALSE,sep=",",row.names=FALSE)

id <- toTable(org.Hs.egENSEMBL)
a <- young[young$species=="ENS",1] %>% as.data.frame()
a$id <- id[match(a$.,id$ensembl_id),1]
a <- na.omit(a)
kegg <- enrichKEGG(a$id, organism = 'hsa', keyType = 'kegg', 
                                  pvalueCutoff = 0.05,pAdjustMethod = 'BH', 
                                  minGSSize = 10, maxGSSize = 500,
                                  qvalueCutoff = 0.2,use_internal_data = FALSE)
kegg_df <- as.data.frame(kegg)
write.table(kegg_df,file = "/home/G20/geneAge/all_species_gene_age/new/human_younggene_kegg.csv",quote=FALSE,sep=",",row.names=FALSE)
```

##### RE-ANALYSIS
###### all species gene age distribution
```{r gene age distribution}
all <- read.csv("/home/qians/geneage/all_species",header = FALSE,sep = ",",stringsAsFactors = FALSE)
```
```{r}
ggplot(all,aes(V2))+geom_density()+Theme+
  xlab("Gene Age (myr)")+ ylab("Distribution")+
  labs(title="All species")+theme(plot.title = element_text(hjust = 0.5))
head(all)
table(all$V2)
all$V4 <- cut(all$V2,breaks = c(seq(5,50,5),seq(100,1000,50)))
table(all$V4)
age <- table(all$V4) %>% as.data.frame() 
#age$Var1 <- factor(age$Var1,levels = c(sort(c(seq(5,50,by = 5),seq(100,800,by = 50))),">800"))
colnames(age) <- c("age","Number")
age <- age[age$Number!=0,]
age$age <- as.vector(age$age)
age[23,1] <- ">800"
age$age <- factor(age$age,levels = age$age)
ggplot(age,aes(age,log(Number)))+geom_point()+Theme+xlab("Gene Age (myr)")+ylab("Gene Number (#log)")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
write.table(age,file = "/home/G20/geneAge/all_species_gene_age/new/all_species_age",sep = '\t',row.names = FALSE,quote = FALSE)
```

```{r gene age cumulative distribution}
#累加图
all <- read.csv("/home/G20/geneAge/all_species_gene_age/all_species",header = FALSE,sep = ",",stringsAsFactors = FALSE)
age <- table(all$V2) %>% as.data.frame()
head(age)
age$cum <- cumsum(age$Freq)
test <- as.data.frame(1:770)
colnames(test) <- "age"
test$number <- age[match(test$age,age$Var1),2]
test[is.na(test)] <- 0
test$age <- rev(test$age)
test$cum <- cumsum(test$number) 
ggplot(test,aes(age,log10(cum),group=1))+geom_point()+geom_line()+
  scale_x_discrete(name = "Age (myr)",breaks =seq(0,750,50))+ylab("Gene cumulative distribution (#log10)")
write.table(test,file = "/home/G20/geneAge/all_species_gene_age/new/all_species_age",sep = '\t',row.names = FALSE,quote = FALSE)
#ggplot(test,aes(.,cum,group=1))+geom_line()+geom_point()+
#  scale_x_continuous(limits=c(0, 800), breaks=c(200, 400))
#  scale_x_discrete(name = "Age (myr)",breaks =c(1,5,10,15,20,25,30,35,40,50,55,63,70,81,80,101,115,126,134,151,168,177,206,218,235,266,291,332,424,543,646,737,770))+theme(axis.text.x = element_text(angle = 45, hjust = 1))

```

```{r PMC2880965 gene age data}
DB1 <- read.csv("/home/G20/geneAge/all_species_gene_age/two_DB/PMC2880965_Homo_sapiens.csv",header = T,sep = ",",stringsAsFactors = FALSE)
DB1$origin <- h[match(DB1$ensembl_gene_id,h$ensembl_id),2]
DB1 <- DB1[DB1$gene_Interval!="1105-",]
DB1 <- DB1[order(DB1$gene_age),]
DB1$gene_Interval <- factor(DB1$gene_Interval,levels=unique(DB1$gene_Interval))
ggplot(DB1,aes(gene_Interval,origin))+geom_boxplot(outlier.colour = "white")+Theme+
  xlab("Gene Age in PMC2880965")+ylab("Gene Age in GenOrigin")+
  #theme(axis.text.x = element_text(angle = 45, hjust = 1))+Theme
  annotate("text",x=3,y=900,label="Spearman correlation = 0.38 \n P value < 2.2e-16",size=5)
```

```{r PMC4943184 gene age data}
#setwd("/home/G20/geneAge/all_species_gene_age/two_DB/")
#PMC4943184 <- list.files(pattern="*.csv",path = "/home/G20/geneAge/all_species_gene_age/two_DB/")
#PMC4943184 <- PMC4943184[-1]
#fcsv <-lapply(PMC4943184,read.csv)
#PMC4943184  <- data.frame(do.call(rbind,fcsv))
setwd("/home/G20/geneAge/all_species_gene_age/13species/")
origin <- list.files(pattern="*.csv",path = "/home/G20/geneAge/all_species_gene_age/13species/")
fcsv2 <- lapply(origin,read.csv)
origin <- data.frame(do.call(rbind,fcsv2))

Tdb <- read.csv("/home/G20/geneAge/all_species_gene_age/two_DB/othersource.csv",header = T,sep = ",",stringsAsFactors = FALSE)
PMC4943184 <- Tdb[Tdb$source=="PMC4943184",]
PMC4943184$origin <- origin[match(PMC4943184$ensembl_gene_id,origin$ensembl_gene_id),2]
PMC4943184 <- na.omit(PMC4943184[order(PMC4943184$gene_age),])
PMC4943184$gene_Interval <-  factor(PMC4943184$gene_Interval,levels = unique(PMC4943184$gene_Interval))
ggplot(PMC4943184,aes(gene_Interval,origin))+geom_boxplot(outlier.colour = "white")+Theme+
  xlab("Gene Age in PMC4943184")+ylab("Gene Age in GenOrigin")+
  annotate("text",x=6,y=150,label="Spearman correlation = 0.42 \n P value < 2.2e-16",size=5)
  
PMC6442393 <- Tdb[Tdb$source=="PMC6442393",]
PMC6442393$origin <-  origin[match(PMC6442393$ensembl_gene_id,origin$ensembl_gene_id),2]
PMC6442393 <- na.omit(PMC6442393[order(PMC6442393$gene_age),])
PMC6442393$gene_Interval <-  factor(PMC6442393$gene_Interval,levels = unique(PMC6442393$gene_Interval))
PMC6442393$group <- cut(PMC6442393$gene_age,breaks = seq(0,450,by = 50))
ggplot(PMC6442393,aes(group,origin))+geom_boxplot(outlier.colour = "white")+Theme+
    xlab("Gene Age in PMC6442393")+ylab("Gene Age in GenOrigin")+
    annotate("text",x=3,y=850,label="Spearman correlation = 0.31 \n P value < 2.2e-16",size=4)
```

```{r Gene age distribution}
all <- read.csv("/home/G20/geneAge/all_species_gene_age/Allspecies.csv",header = T,sep = ",",stringsAsFactors = FALSE)
Theme <-  theme_bw()+theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank(),axis.title = element_text(size = 15))
ggplot(all,aes(gene_age))+geom_density()+Theme+
  xlab("Gene Age (myr)")+ ylab("Distribution")+
  labs(title="All species")+theme(plot.title = element_text(hjust = 0.5))
# for (i in 1:nrow(all)) {
#   if (all[i,2]>= 50) {
#     all[i,3] = 100
#   }else if (all[i,2]>= 100) {
#     all[i,3] = 150
#   }
# }
all[all$gene_age>=  0,7] = 5
all[all$gene_age>=  5,7] = 10
all[all$gene_age>= 10,7] = 15
all[all$gene_age>= 15,7] = 20
all[all$gene_age>= 20,7] = 25
all[all$gene_age>= 25,7] = 30
all[all$gene_age>= 30,7] = 35
all[all$gene_age>= 35,7] = 40
all[all$gene_age>= 40,7] = 45
all[all$gene_age>= 45,7] = 50
all[all$gene_age>= 50,7] = 100
all[all$gene_age>=100,7] = 150
all[all$gene_age>=150,7] = 200
all[all$gene_age>=200,7] = 250
all[all$gene_age>=250,7] = 300
all[all$gene_age>=300,7] = 350
all[all$gene_age>=350,7] = 400
all[all$gene_age>=400,7] = 450
all[all$gene_age>=450,7] = 500
all[all$gene_age>=500,7] = 550
all[all$gene_age>=550,7] = 600
all[all$gene_age>=600,7] = 650
all[all$gene_age>=650,7] = 700
all[all$gene_age>=700,7] = 750
all[all$gene_age>=750,7] = 800
all[all$gene_age>=800,7] = ">800"
all$V8 <- cut(all$gene_age,breaks = c(seq(0,50,5),seq(100,800,50),1600))
age <- table(all$V8) %>% as.data.frame() 
#age$Var1 <- factor(age$Var1,levels = c(sort(c(seq(5,50,by = 5),seq(100,800,by = 50))),">800"))
colnames(age) <- c("age","Number")
ggplot(age,aes(age,log(Number)))+geom_point()+Theme+xlab("Gene Age (myr)")+ylab("Gene Number (#log)")
write.table(age,file = "/home/G20/geneAge/all_species_gene_age/all_species_age",sep = '\t',row.names = FALSE,quote = FALSE)
```

```{r correlation with other database}
other <- read.csv("/home/G20/geneAge/all_species_gene_age/two_DB/othersource.csv",header = T,sep = ",",stringsAsFactors = FALSE)
table(other$species,other$source)
other$origin <- all[match(other$ensembl_gene_id,all$ensembl_gene_id),2]
tapply(other[other$source=="PMC6442393",2],other[other$source=="PMC6442393",5],mean) #fly gene age in Gentree
tapply(other[other$source=="PMC4943184",2],other[other$source=="PMC4943184",5],mean)
#fly should be delete in Gentree (PMC6442393)
PMC6442393 <- other[other$source=="PMC6442393"& other$species!="Drosophila melanogaster",]
table(is.na(PMC6442393$origin))
PMC6442393 <- na.omit(PMC6442393)
PMC2880965 <- other[other$source=="PMC2880965",]
table(is.na(PMC2880965$origin))
cor.test(PMC2880965$origin,PMC2880965$gene_age,method = "spearman")
PMC4943184 <- other[other$source=="PMC4943184",]
table(is.na(PMC4943184$origin))
cor.test(PMC4943184$origin,PMC4943184$gene_age,method = "spearman")
all$PMC2880965 <- PMC2880965[match(all$ensembl_gene_id,PMC2880965$ensembl_gene_id),2]
all$PMC6442393 <- PMC6442393[match(all$ensembl_gene_id,PMC6442393$ensembl_gene_id),2]
all$PMC4943184 <- PMC4943184[match(all$ensembl_gene_id,PMC4943184$ensembl_gene_id),2]
cor.test(all$gene_age,all$PMC2880965,method = "spearman")
cor.test(all$gene_age,all$PMC6442393,method = "spearman")
cor.test(all$gene_age,all$PMC4943184,method = "spearman") 
row.names(all) <- all$ensembl_gene_id
a <- all[,c(2,9,10,11)]
library(corrplot,lib.loc = "/home/G4/R/x86_64-pc-linux-gnu-library/3.4")
b <- cor(a) %>% as.data.frame()
cor.test(a$gene_age,a$PMC2880965,method = "spearman")
cor.test(a$gene_age,a$PMC6442393,method = "spearman")
cor.test(a$gene_age,a$PMC4943184,method = "spearman")
cor.test(a$PMC2880965,a$PMC6442393,method = "spearman")
cor.test(a$PMC2880965,a$PMC4943184,method = "spearman")
b[1,4] = 0.4847178
b[1,3] = 0.5531828
b[1,2] = 0.3559807
b[3,4] = 0.354399
b[2,4] = 0.4156448
b[2,3] = 0.3214826
for (i in 1:3) {
  c[i] <- cor.test(human[,i],human[,i+1])
}
human <- all[all$common_name=="Human",c(2,9:11)]
M <- cor(human,method = "spearman")
row.names(M) <- c("GenOrigin","PMID20492640","PMID30862647","PMID27259914")
colnames(M) <- c("GenOrigin","PMID20492640","PMID30862647","PMID27259914")
corrplot(M, type = "upper", tl.pos = "d",cl.lim = c(0,1))
corrplot(M, add = TRUE, type = "lower",col = "black", method = "number",
         diag = FALSE, tl.pos = "n", cl.pos = "n")
```

```{r new gene analysis}
a <- all[,c(1,2,5,6,10)]
sum(is.na(a)) #NA的个数
sum(is.na(a[5]))
a <- na.omit(a) 
##cat /home/qians/geneage/annotation/*csv | sort| uniq > /home/G20/geneAge/all_species_gene_age/five_species_annotation
##for i in `ls /home/qians/geneage/annotation/`; do awk -F "," '{print $1}' /home/qians/geneage/annotation/$i | sort| uniq | wc -l; done
##for i in `ls /home/qians/geneage/annotation/`; do sort /home/qians/geneage/annotation/$i | uniq -c|wc -l; done
anno <- read.csv("/home/G20/geneAge/all_species_gene_age/five_species_annotation",header = T,sep = "\t",stringsAsFactors = F)
a$chr <- anno[match(a$ensembl_gene_id,anno$ensembl_gene_id),3]
a[a$chr!="X"&a$chr!="Z",6] = "A"
a$type <- anno[match(a$ensembl_gene_id,anno$ensembl_gene_id),2]
#result:
nrow(a)
cor.test(a$gene_age,a$PMC6442393,method = "spearman")
young <- a[a$gene_age <50,]
table(young$common_name)/table(a$common_name)
table(young$common_name)
nrow(young) / nrow(a)
nrow(young)

table(young$chr=="X"|young$chr=="Z")/nrow(young)
nrow(young[young$chr=="X"|young$chr=="Z",])
table(a$chr=="X"|a$chr=="Z")/nrow(a)
nrow(a[a$chr=="X"|a$chr=="Z",])
matrix(c(nrow(young[young$chr=="X"|young$chr=="Z",]),nrow(young),nrow(a[a$chr=="X"|a$chr=="Z",]),nrow(a)),nrow = 2) %>% fisher.test()
nrow(young[young$type!="protein_coding",])/nrow(young)
nrow(young[young$type!="protein_coding",])
nrow(a[a$type!="protein_coding",])/nrow(a)
nrow(a[a$type!="protein_coding",])
```

```{r}
#GO for young gene
##human
library(clusterProfiler)
library(org.Hs.eg.db,lib.loc = "/home/qians/R/x86_64-pc-linux-gnu-library/3.3")
table(young$common_name)
bp <- enrichGO(young[young$common_name=="Human",1], 
                OrgDb = org.Hs.eg.db, 
                keyType = "ENSEMBL", 
                ont = "BP",
                pvalueCutoff = 0.05
                )
bp_df <- as.data.frame(bp)
write.table(bp,file = "/home/G20/geneAge/all_species_gene_age/new/human_younggene_bp.csv",quote=FALSE,sep=",",row.names=FALSE)
barplot(bp, showCategory=100)
cc <- enrichGO(young[young$common_name=="Human",1], 
                OrgDb = org.Hs.eg.db, 
                keyType = "ENSEMBL", 
                ont = "CC",
                pvalueCutoff = 0.05
                )
cc_df <- as.data.frame(cc)
write.table(cc_df,file = "/home/G20/geneAge/all_species_gene_age/new/human_younggene_cc.csv",quote=FALSE,sep=",",row.names=FALSE)
barplot(cc, showCategory=100)
mf <- enrichGO(young[young$common_name=="Human",1], 
                OrgDb = org.Hs.eg.db, 
                keyType = "ENSEMBL", 
                ont = "MF",
                pvalueCutoff = 0.05
                )
mf_df <- as.data.frame(mf)
barplot(mf, showCategory=100)
write.table(mf_df,file = "/home/G20/geneAge/all_species_gene_age/new/human_younggene_mf.csv",quote=FALSE,sep=",",row.names=FALSE)

id <- toTable(org.Hs.egENSEMBL)
t <- young[young$common_name=="Human",1] %>% as.data.frame()
t$id <- id[match(t$.,id$ensembl_id),1]
t <- na.omit(t)
kegg <- enrichKEGG(t$id, organism = 'hsa', keyType = 'kegg', 
                                  pvalueCutoff = 0.05,pAdjustMethod = 'BH', 
                                  minGSSize = 10, maxGSSize = 500,
                                  qvalueCutoff = 0.2,use_internal_data = FALSE)
kegg_df <- as.data.frame(kegg)
barplot(kegg, showCategory=100)
write.table(kegg_df,file = "/home/G20/geneAge/all_species_gene_age/new/human_younggene_kegg.csv",quote=FALSE,sep=",",row.names=FALSE)
```
