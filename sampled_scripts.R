library('Seurat')
library('SeuratDisk')
library('dplyr')
library('patchwork')
library('scCATCH')
library('SingleR')
library('mindr')
library('ggplot2')

gene_short_name<-as.character(df_gene$gene_short_name)
gene_count_sampled@Dimnames[[1]]<-gene_short_name



Sampled_raw<-CreateSeuratObject(count=gene_count_sampled)

library(tidyverse)
data<-Sampled_raw@meta.data
data
data<-rownames_to_column(data,var="barcodes")
metadata_Sampled<-merge(data, df_cell, by.x='barcodes', by.y='sample')
metadata_Sampled<-column_to_rownames(metadata_Sampled, var="barcodes")
Sampled_raw<-AddMetaData(Sampled_raw, metadata_Sampled)
colnames(Sampled_raw@meta.data)
Sampled_processed<-subset(Sampled_raw, Main_cluster_name != 'NA')

Sampled_processed<-NormalizeData(Sampled_processed,normalization.method = "LogNormalize",scale.factor = 10000)
gc()
Sampled_processed<-FindVariableFeatures(Sampled_processed, selection.method = "vst", nfeatures = 40000)
gc()
Sampled_processed<-ScaleData(Sampled_processed)
gc()
Sampled_processed<-RunPCA(Sampled_processed, features=VariableFeatures(object=Sampled_processed))
gc()
print(Sampled_processed[["pca"]], dims = 1:15, nfeatures = 20)
pct <- Sampled_processed [["pca"]]@stdev / sum( Sampled_processed [["pca"]]@stdev)*100
# *100 
# Calculate cumulative percents for each PC
cumu <- cumsum(pct)
# Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
co1 <- which(cumu > 95 & pct < 5)[1]
co1
# Determine the difference between variation of PC and subsequent PC
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
# last point where change of % of variation is more than 0.1%.
co2
# Minimum of the two calculation
pcs <- min(co1, co2)
pcs
# Create a dataframe with values
plot_df <- data.frame(pct = pct,   cumu = cumu,   rank = 1:length(pct))
# Elbow plot to visualize 
dev.new()
ggplot(plot_df, aes(cumu, pct, label = rank, color = rank > pcs)) + 
  geom_text() + 
  geom_vline(xintercept = 95, color = "grey") + 
  geom_hline(yintercept = min(pct[pct > 5]), color = "grey") +
  theme_bw()
#???࣬????ŷ?Ͼ?????KNNͼ


Sampled_processed<-FindNeighbors(Sampled_processed, dims=1:50)


Sampled_processed<-RunUMAP(Sampled_processed, dims=1:50)


DimPlot(all_combine, reduction = "umap", label=T, label.size=5,raster=FALSE, group.by = "celltype")
DimPlot(Sampled_processed, reduction = "tsne", label=T, label.size=5,raster=FALSE, group.by = "Organ")
DimPlot(Sampled_processed, reduction = "umap", label=T, label.size=5,raster=FALSE, group.by = "Main_cluster_name")
DimPlot(Sampled_processed, reduction = "tsne", label=T, label.size=5,raster=FALSE, group.by = "Main_cluster_name")
saveRDS(Sampled_processed, file="Sampled_processed2.rds")

library('tidydr')
library('ggplot2')
pp4<-DotPlot(Sampled_processed, features = c("ENSG00000163399.11","ENSG00000147533.12","ENSG00000156599.6","ENSG00000140564.6","ENSG00000183762.8","ENSG00000141505.7","ENSG00000099250.13","ENSG00000184012.7","ENSG00000130234.6"),group.by = 'Development_day')+coord_flip()+
  theme_bw()+
  theme(panel.grid = element_blank(), axis.text.x=element_text(hjust = 1,vjust=0.5))+
  labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3))+
  scale_color_gradientn(values = seq(0,1,0.2),colours = c('#132651','#4D95C7','#FBC241','#EC581A'))
pp4

Placenta_processed$Main_cluster_name <- factor(Placenta_processed$Main_cluster_name,levels=c('AFP_ALB positive cells','Extravillous trophoblasts','IGFBP1_DKK1 positive cells','Lymphoid cells','Megakaryocytes','PAEP_MECOM positive cells','Myeloid cells','Smooth muscle cells','Syncytiotrophoblasts and villous cytotrophoblasts','Trophoblast giant cells','Stromal cells','Vascular endothelial cells'))

pp4<-DotPlot(Placenta_processed, features = c("ENSG00000163399.11","ENSG00000147533.12","ENSG00000156599.6","ENSG00000140564.6","ENSG00000183762.8","ENSG00000141505.7","ENSG00000099250.13","ENSG00000184012.7","ENSG00000130234.6"),group.by = 'Main_cluster_name')+coord_flip()+
  theme_bw()+
  theme(panel.grid = element_blank(), axis.text.x=element_text(hjust = 1,vjust=0.5))+
  labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3))+
  scale_color_gradientn(values = seq(0,1,0.2),colours = c('#132651','#4D95C7','#FBC241','#EC581A'))
pp4


'#132651','#4D95C7','#FBC241','#EC581A'

VlnPlot(Sampled_processed, features = c("ENSG00000130234.6","ENSG00000184012.7"),group.by="Main_cluster_name",
        stack=T, pt.size = 5,
        flip = T,
        add.noise = T)+#横纵轴不标记任何东西
  ggplot2::scale_fill_manual(values = c("#ff5f2e","#fcbe32"))
"#004e66","#e1eef6"#3C5488FF","#F39B7FFF","#8491B4FF"))
                                        
DimPlot(all_combine, reduction='umap',label.box = F,label = F, group.by = "Main_cluster_name",raster = FALSE)+theme_dr()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = 10,face = 'bold'),
        plot.title = element_text(hjust = 0.5,size = 16,face = 'bold'))+
  scale_fill_manual(values = c('white','white','white','white','white','white','white','white','white','white','white','white','white','white','white','white','white','white','white','white','white','white','white','white','white','white','white','white','white','white','white','white','white','white','white','white','white','white','white','white','white','white','white','white','white','white','white','white','white','white','white','white','white','white','white','white','white','white','white','white','white','white','white','white','white','white','white','white','white','white','white','white','white','white','white','white','white','white','white','white','white','white','white','white','white','white','white','white','white','white','white','white','white','white','white','white','white','white','white','white','white','white','white','white','white','white','white','white','white','white','white','white','white','white','white','white','white','white','white','white','white','white','white','white','white','white','white','white')) +
  scale_color_manual(values = c("#04111a","#041923","#05212c","#052835","#05303e","#063847","#06404f","#064858","#064f61","#07576a","#075f73","#086477","#096a7b","#096f7e","#0a7482","#0b7a86","#0c7f8a","#0d848e","#0d8991","#0e8f95","#0f9499","#1c9a9e","#29a1a2","#36a7a7","#43adab","#50b4b0","#5dbab5","#6ac0b9","#77c6be","#84cdc2","#91d3c7","#9ad4c4","#a3d4c1","#add5be","#b6d5bb","#bfd6b8","#c8d6b5","#d1d7b2","#dbd7af","#e4d8ac","#edd8a9","#edd298","#edcc88","#edc677","#edc066","#edbb56","#ecb545","#ecaf34","#eca923","#eca313","#ec9d02","#e99802","#e59202","#e28d02","#df8802","#dc8303","#d87d03","#d57803","#d27303","#ce6d03","#cb6803","#c96403","#c86003","#c65c03","#c55803","#c35503","#c15102","#c04d02","#be4902","#bd4502","#bb4102","#ba3e04","#b93b05","#b83707","#b73409","#b6310b","#b42e0c","#b32b0e","#b22710","#b12411","#b02113","#ae2115","#961C12","#821810","#6E140E","#5E110C")) +
  labs(title="Main_cluster")

DimPlot(all_combine, reduction='umap',label.box = F,label = F, group.by = "Main_cluster_name",raster = FALSE)+theme_dr()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = 10,face = 'bold'),
        plot.title = element_text(hjust = 0.5,size = 16,face = 'bold'))+
  scale_fill_manual(values = c('white','white','white','white','white','white','white','white','white','white','white','white','white','white','white','white','white','white','white','white','white','white','white','white','white','white','white','white','white','white','white','white','white','white','white','white','white','white','white','white','white','white','white','white','white','white','white','white','white','white','white','white','white','white','white','white','white','white','white','white','white','white','white','white','white','white','white','white','white','white','white','white','white','white','white','white','white','white','white','white','white','white','white','white','white','white','white','white','white','white','white','white','white','white','white','white','white','white','white','white','white','white','white','white','white','white','white','white','white','white','white','white','white','white','white','white','white','white','white','white','white','white','white','white','white','white','white','white')) +
  scale_color_manual(values = c("#D3D3D3","#D3D3D3","#D3D3D3","#D3D3D3","#D3D3D3","#D3D3D3","#D3D3D3","#D3D3D3","#D3D3D3","#D3D3D3","#D3D3D3","#D3D3D3","#D3D3D3","#D3D3D3","#D3D3D3","#D3D3D3","#D3D3D3","#D3D3D3","#D3D3D3","#D3D3D3","#D3D3D3","#D3D3D3","#D3D3D3","#D3D3D3","#D3D3D3","#D3D3D3","#D3D3D3","#D3D3D3","#D3D3D3","#D3D3D3","#D3D3D3","#D3D3D3","#D3D3D3","#D3D3D3","#D3D3D3","#D3D3D3","#D3D3D3","#D3D3D3","#D3D3D3","#D3D3D3","#D3D3D3","#D3D3D3","#D3D3D3","#D3D3D3","#D3D3D3","#D3D3D3","#D3D3D3","#D3D3D3","#D3D3D3","#D3D3D3","#D3D3D3","#D3D3D3","#D3D3D3","#D3D3D3","#D3D3D3","#D3D3D3","#D3D3D3","#D3D3D3","#D3D3D3","#D3D3D3","#D3D3D3","#D3D3D3","#D3D3D3","#D3D3D3","#D3D3D3","#D3D3D3","#D3D3D3","#b93b05","#D3D3D3","#D3D3D3","#D3D3D3","#D3D3D3","#D3D3D3","#D3D3D3","#D3D3D3","#b02113","#D3D3D3")) +
  labs(title="Sampled_Main_cluster")


DimPlot(all_combine, reduction='umap',label.box = F,label = F, group.by = "Organ",raster = FALSE)+theme_dr()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = 10,face = 'bold'),
        plot.title = element_text(hjust = 0.5,size = 16,face = 'bold'))+
  scale_fill_manual(values = c('white','white','white','white','white','white','white','white','white','white','white','white','white','white','white','white')) +
  scale_color_manual(values = c('#04111A',"#EDD8A9","#EC9D02","#075F73","#CB6803","#91D3C7","#BB4102","#0F9499","#3C5488FF","#9E2127","#8491B4FF","#F39B7FFF","#B02113","#FFC0CB","#bebada","#80b1d3")) +
  labs(title="Organs")

DimPlot(Sampled_processed, reduction='umap',label.box = F,label = F, group.by = "Organ",raster = FALSE)+theme_dr()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = 10,face = 'bold'),
        plot.title = element_text(hjust = 0.5,size = 16,face = 'bold'))+
  scale_fill_manual(values = c('white','white','white','white','white','white','white','white','white','white','white','white','white','white')) +
  scale_color_manual(values = c('#04111A',"#075F73","#91D3C7","#EDD8A9","#EC9D02","#CB6803","#BB4102","#B02113","#9E2127","#3C5488FF","#8491B4FF","#F39B7FFF","#FFC0CB","#bebada","#80b1d3")) +
  labs(title="Organs")

DimPlot(Sampled_processed, reduction='umap',label.box = T,label = T, group.by = "Organ",raster = FALSE)+theme_dr()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = 10,face = 'bold'),
        plot.title = element_text(hjust = 0.5,size = 16,face = 'bold'))+
    scale_color_manual(values = c('#04111A',"#EDD8A9","#EC9D02","#075F73","#CB6803","#91D3C7","#BB4102","#0F9499","#3C5488FF","#9E2127","#8491B4FF","#F39B7FFF")) +
  labs(title="Sampled_Development_day")

'#04111A',"#075F73","#0F9499","#91D3C7","#EDD8A9","#EC9D02","#CB6803","#9E2127"
DimPlot(Sampled_processed, reduction='umap',label.box = T,label = T, group.by = "Development_day",raster = FALSE)+theme_dr()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = 10,face = 'bold'),
        plot.title = element_text(hjust = 0.5,size = 16,face = 'bold'))+
  scale_fill_manual(values = c('white','white','white','white','white','white','white','white','white','white','white','white','white')) +
  scale_color_manual(values = c("#075F73","#EDD8A9","#EC9D02","#9E2127","#3C5488FF","#F39B7FFF")) +
  labs(title="Sampled_development_day")
FeaturePlot(Sampled_processed, features = c("ENSG00000130234.6","ENSG00000184012.7","ENSG00000099250.13","ENSG00000141505.7","ENSG00000183762.8","ENSG00000140564.6","ENSG00000156599.6","ENSG00000147533.12","ENSG00000163399.11"),cols = c("#D3D3D3" ,"#DE1F1F"),slot = "data",label.size = 6,pt.size = 0.5,raster=FALSE,max.cutoff = 1) 
FeaturePlot(Sampled_processed, features = c("ENSG00000099250.13","ENSG00000183762.8"),cols = c("#D3D3D3" ,"#DE1F1F"),slot = "data",label.size = 6,pt.size = 0.3,raster=FALSE,max.cutoff = 1) 
FeaturePlot(Sampled_processed, features = c("ENSG00000134853.7","ENSG00000113721.9","ENSG00000185551.8"),cols = c("#D3D3D3" ,"#DE1F1F"),slot = "data",label.size = 6,pt.size = 0.3,raster=FALSE,max.cutoff = 1) 


FeaturePlot(Sampled_processed, features = c("ENSG00000130234.6"),cols = c("#D3D3D3" ,"#DE1F1F"),slot = "data",label.size = 6,pt.size = 0.5,raster=FALSE,max.cutoff = 1) 
FeaturePlot(Sampled_processed, features = c("ENSG00000184012.7"),cols = c("#D3D3D3" ,"#DE1F1F"),slot = "data",label.size = 6,pt.size = 0.5,raster=FALSE,max.cutoff = 1) 
FeaturePlot(Sampled_processed, features = c("ENSG00000099250.13"),cols = c("#D3D3D3" ,"#DE1F1F"),slot = "data",label.size = 6,pt.size = 0.5,raster=FALSE,max.cutoff = 1) 
FeaturePlot(Sampled_processed, features = c("ENSG00000141505.7"),cols = c("#D3D3D3" ,"#DE1F1F"),slot = "data",label.size = 6,pt.size = 0.5,raster=FALSE,max.cutoff = 1) 
FeaturePlot(Sampled_processed, features = c("ENSG00000183762.8"),cols = c("#D3D3D3" ,"#DE1F1F"),slot = "data",label.size = 6,pt.size = 0.5,raster=FALSE,max.cutoff = 1) 
FeaturePlot(Sampled_processed, features = c("ENSG00000140564.6"),cols = c("#D3D3D3" ,"#DE1F1F"),slot = "data",label.size = 6,pt.size = 0.5,raster=FALSE,max.cutoff = 1) 
FeaturePlot(Sampled_processed, features = c("ENSG00000156599.6"),cols = c("#D3D3D3" ,"#DE1F1F"),slot = "data",label.size = 6,pt.size = 0.5,raster=FALSE,max.cutoff = 1) 
FeaturePlot(Sampled_processed, features = c("ENSG00000147533.12"),cols = c("#D3D3D3" ,"#DE1F1F"),slot = "data",label.size = 6,pt.size = 0.5,raster=FALSE,max.cutoff = 1) 
FeaturePlot(Sampled_processed, features = c("ENSG00000163399.11"),cols = c("#D3D3D3" ,"#DE1F1F"),slot = "data",label.size = 6,pt.size = 0.5,raster=FALSE,max.cutoff = 1) 


FeaturePlot(endothelial, features = c("ENSG00000130234.6","ENSG00000184012.7","ENSG00000099250.13"),cols = c("#D3D3D3" ,"#DE1F1F"),slot = "data",label.size = 6,pt.size = 0.5,max.cutoff = 2) 

FeaturePlot(Sampled_processed, features = c("ENSG00000261371.1","ENSG00000113721.9","ENSG00000110799.9","ENSG00000164736.5"),cols = c("lightgrey" ,"#DE1F1F"),slot = "data",label.size = 6,pt.size = 0.5,raster=FALSE,max.cutoff = 1) 
FeaturePlot(Sampled_processed, features = c("ENSG00000130234.6","ENSG00000184012.7","ENSG00000099250.13","ENSG00000141505.7","ENSG00000183762.8","ENSG00000140564.6","ENSG00000156599.6","ENSG00000147533.12","ENSG00000163399.11"),cols = c("#D3D3D3" ,"#DE1F1F"),slot = "data",label.size = 6,pt.size = 0.5,raster=FALSE) 
FeaturePlot(Sampled_processed, features = c("ENSG00000130234.6"),cols = c("#D3D3D3" ,"#DE1F1F"),slot = "data",label.size = 6,pt.size = 0.5,raster=FALSE) 


Sampled_processed<-subset(Sampled_processed, subset=Main_cluster_name %in% c('Sampled_processed cells')) 
Sampled_processed<-subset(Sampled_processed, subset=Organ %in% c('Heart','Intestine','Sampled','Lung','Pancreas','Placenta','Spleen')) 


saveRDS(Sampled_processed, file="Sampled_processed.rds")
library(tidydr)
library(ggsci)
DimPlot(Sampled_processed, reduction='umap',label.box = T,label = T, group.by = "Main_cluster_name",raster = FALSE)+theme_dr()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = 10,face = 'bold'),
        plot.title = element_text(hjust = 0.5,size = 16,face = 'bold'))+
   scale_color_npg()+labs(title="Sampled_Development_day")


#loom文件转换
devtools::install_github(repo = "hhoeflin/hdf5r")
devtools::install_github(repo = "mojaveazure/loomR", ref = "develop")
library(loomR)
/mnt/NAS/NAS2/Public/SeqShare/PMID33184181/zyy_analysis

Sampled_processed=connect('/mnt/NAS/NAS2/Public/SeqShare/PMID33184181/Processed_datasets/Human_RNA_processed.loom',mode = 'r+',skip.validate = TRUE)
Sampled_processed<-as.Seurat(Sampled_processed)
AllCount=Sampled_processed[["matrix"]][,] #提取.loom文件的matrix
Allcount<-Sampled_processed$matrix[,]
ds[['/matrix']]$methods()psoriasis=t(psoriasis) #发现那个matrix的gene和barcode颠倒了,换过来才符合seurat对象中的矩阵行为基因,列为barcode
dim(psoriasis)
[1] 19968 24234
gene=hc.2$row.attrs$Gene[] #提取基因名
barcode=hc.2$col.attrs$CellID[]# 提取barcode
length(gene)
[1] 19968
length(barcode)
[1] 24234
colnames(psoriasis)= barcode
row.names(psoriasis)= gene
dim(psoriasis)
[1] 19968 24234
####创建seurat对象####
psoriasis=CreateSeuratObject(counts = psoriasis,project = 'psoriasis',min.cells = 3, min.features = 200)
psoriasis

DimPlot(Vascular,reduction='umap', label=T, label.size=5)

###### merge data
Alldata<-merge(Cerebrum_processed,y=c(Cerebellum_processed, Spleen_processed, Placenta_processed, Pancreas_processed, Sampled_processed, Heart_processed, Liver_processed, Sampled_processed, Lung_processed))


save.image ("Seperate_samples") 

table(Sampled_processed$Main_cluster_name)



#vascular cells
Vascular<-subset(all_combine, subset=Main_cluster_name %in% c('Vascular endothelial cells')) 
Vascular<-RunPCA(Vascular, features=VariableFeatures(object=Vascular))
gc()
print(Vascular[["pca"]], dims = 1:15, nfeatures = 20)
pct <- Vascular [["pca"]]@stdev / sum( Vascular [["pca"]]@stdev)*100
# *100 
# Calculate cumulative percents for each PC
cumu <- cumsum(pct)
# Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
co1 <- which(cumu > 95 & pct < 5)[1]
co1
# Determine the difference between variation of PC and subsequent PC
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
# last point where change of % of variation is more than 0.1%.
co2
# Minimum of the two calculation
pcs <- min(co1, co2)
pcs
# Create a dataframe with values
plot_df <- data.frame(pct = pct,   cumu = cumu,   rank = 1:length(pct))
# Elbow plot to visualize 
dev.new()
ggplot(plot_df, aes(cumu, pct, label = rank, color = rank > pcs)) + 
  geom_text() + 
  geom_vline(xintercept = 95, color = "grey") + 
  geom_hline(yintercept = min(pct[pct > 5]), color = "grey") +
  theme_bw()
#???࣬????ŷ?Ͼ?????KNNͼ


Vascular<-FindNeighbors(Vascular, dims=1:17)
Vascular<- FindClusters(Vascular, resolution = 1)

Vascular<-RunUMAP(Vascular, dims=1:17)


DimPlot(Vascular, reduction = "umap", label=T, label.size=5,raster=FALSE)
DimPlot(Sampled_processed, reduction = "tsne", label=T, label.size=5,raster=FALSE, group.by = "Organ")
DimPlot(Sampled_processed, reduction = "umap", label=T, label.size=5,raster=FALSE, group.by = "Main_cluster_name")
DimPlot(Sampled_processed, reduction = "tsne", label=T, label.size=5,raster=FALSE, group.by = "Main_cluster_name")
saveRDS(Sampled_processed, file="Sampled_processed2.rds")


FeaturePlot(Vascular, features = c("PECAM1","CD34","TIE1","ANPEP","RGS5","PDGFRB"),cols = c("#D3D3D3" ,"#DE1F1F"),slot = "data",label.size = 6,pt.size = 0.5,max.cutoff = 2) 
FeaturePlot(Vascular, features = c("MYL9","TAGLN","SERPING1","LUM","ATP1A2","KCNJ8","CD248"),cols = c("#D3D3D3" ,"#DE1F1F"),slot = "data",label.size = 6,pt.size = 0.5,max.cutoff = 2) 

FeaturePlot(Vascular, features = c("ACE2","TMPRSS2","NRP1","ASGR1","KREMEN1","FURIN","ZDHHC5","GOLGA7","ATP1A1"),cols = c("#D3D3D3" ,"#DE1F1F"),slot = "data",label.size = 6,pt.size = 0.5,max.cutoff = 2) 
