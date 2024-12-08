#Figure2 Pannel A
Organ_positive<-subset(Sampled_processed,subset=Organ %in% c("Placenta","Pancreas","Lung","Intestine","Cerebrum"))
Organ_positive<-NormalizeData(Organ_positive,normalization.method = "LogNormalize",scale.factor = 10000)
Organ_positive<-FindVariableFeatures(Organ_positive, selection.method = "vst", nfeatures = 3000)
Organ_positive<-ScaleData(Organ_positive,features=rownames(Organ_positive))
Organ_positive<-RunPCA(Organ_positive, features=VariableFeatures(object=Organ_positive))
Organ_positive<-FindNeighbors(Organ_positive, dims=1:50)
Organ_positive<-RunUMAP(Organ_positive, dims=1:50)
DimPlot(Organ_positive, reduction = "umap", label=F, label.size=5,raster=FALSE, group.by = "Organ",cols = c('#132651','#4D95C7','#FBC241','#EC581A',"#91D3C7"))

#Figure2 Pannel B
DimPlot(Organ_positive, reduction = "umap", label=F, label.size=3,raster=FALSE, group.by = "Development_day",cols=c("#E64B35", "#8F8A8F", "#43B7CB", "#18A89F", "#0F8C87", "#306287", "#806E84", "#E7967F", "#BB9599", "#8494B4", "#8CB9BC" ,"#9FA99D", "#C93430", "#BE1E16", "#89543F", "#937A62", "#B09C85"))

#Figure2 Pannel C
receptors<-c("ENSG00000099250.13","ENSG00000172270.14","ENSG00000072274.8","ENSG00000106460.14","ENSG00000010610.5","ENSG00000167601.7","ENSG00000183762.8","ENSG00000130234.6","ENSG00000141505.7","ENSG00000104938.12","ENSG00000113249.8")
Organ_positive$Organ<-factor(Organ_positive$Organ,levels = c("Intestine","Lung","Pancreas","Placenta","Cerebrum"))
gene_cell_exp<-AverageExpression(Organ_positive, features=receptors, group.by='Organ',slot='data',layer = 'data')                                      
organ_receptors<-as.data.frame(gene_cell_exp[["RNA"]])
pheatmap(organ_receptors, color = c(colorRampPalette(colors = c('#176BA0','black','#EE9A3A'))),display_numbers = F,fontsize = 10, angle_col = 45,border_color = NA,cellwidth = 20,cellheight = 15,cluster_rows = F,cluster_cols = F)

#Figure2 Pannel D
interactors<-c("ENSG00000164733.16","ENSG00000163399.11","ENSG00000135047.10","ENSG00000073060.11","ENSG00000141458.8","ENSG00000156599.6","ENSG00000147533.12","ENSG00000140564.6","ENSG00000175426.6","ENSG00000184012.7")
gene_cell_exp<-AverageExpression(Organ_positive, features=interactors, group.by='Organ',slot='data',layer = 'data')                                      
organ_interactors<-as.data.frame(gene_cell_exp[["RNA"]])
pheatmap(organ_interactors, color = c(colorRampPalette(colors = c('#176BA0','black','#EE9A3A'))),display_numbers = F,fontsize = 10, angle_col = 45,border_color = NA,cellwidth = 20,cellheight = 15,cluster_rows = F,cluster_cols = F)

#Figure3 Pannel A
gene_cell_exp<-AverageExpression(Organ_positive, features=receptors, group.by='Main_cluster_name',slot='data',layer = 'data')                                      
celltype_receptors<-as.data.frame(gene_cell_exp[["RNA"]])
x = 650
fig = pheatmap(celltype_receptors, color = c(colorRampPalette(colors = c('#176BA0','black','#EE9A3A'), 
                                                       bias = 4)(x), colorRampPalette(colors = c('#EE9A3A', "red"), )(1000 - x)),
               display_numbers = F, fontsize = 10, angle_col = 45, border_color = NA,
               cellwidth = 20, cellheight = 15, cluster_rows = F,cluster_cols = T);fig
##Stack plot
cell.prop<-as.data.frame(prop.table(table(Positive2$Organ, Positive2$Main_cluster_name_ordered)))
colname(cell.prop)<-c("cluster","Main_cluster_name","proportion")
ggplot(cell.prop,aes(Var2,Freq,fill=Var1))+
  geom_bar(stat="identity",position="fill")+
  ggtitle("")+
  theme_bw()+
  guides(fill=guide_legend(title=NULL))+
  scale_fill_manual(values = c("Cerebrum"="#98d09d","Intestine"="#fbf398","Lung"="#dadada","Pancreas"="#f7a895","Placenta"="#9b8191"),
                    limits=c("Cerebrum","Intestine","Lung","Pancreas","Placenta"))+
  scale_y_continuous(
    labels = scales::percent_format())+
    theme(panel.background = element_blank(),
    axis.line = element_line(),
    legend.position = "bottom")+
  labs(x=NULL,y="Organ proportion (%)")+
  guides(fill=guide_legend(title = NULL,nrow = 1,byrow = FALSE))+
  theme(panel.grid =element_blank()) +
  theme(axis.text = element_text(size=12, colour = "black"))+
  theme(axis.title.y = element_text(size=12, colour = "black"))+
  theme(panel.border = element_rect(size = 1, linetype = "solid", colour = "black"))+
  theme(axis.text.x = element_text(angle = 45,hjust = 0.8, vjust = 0.6))


# Color Code
Covid_color = c("darkorange1", "deepskyblue3");names(Covid_color) = c("Covid-19-negative", "Covid-19-positive")

refine_cluster_color <- c("#f95b69", "#5452bd", "#E5D2DD", "#f056b4", "#2dddd3",
                          "#363c48", "#1fa7bc", "#af86f5")
names(refine_cluster_color) = c('Ex_Neuron', 'In_Neuron', 'IPC', 'RGC', 'OPC',
                                'Microglia', 'Pericytes', 'Endothelia')

# Figure4 Pannel A
Idents(Covid.sct.clean) = Covid.sct.clean$CellTypeNew
new_idents = c("Microglia", "Pericytes", "OPC", "Endothelia", "In_Neuron", "RGC", "IPC", "Ex_Neuron");names(new_idents) = c('Microglia', 'Peri', 'OPC', 'Endo', 'In_Neuron', 'RGC', 'IPC', 'Ex_Neuron')
Covid.sct.clean = RenameIdents(Covid.sct.clean, new_idents)

fig = DimPlot(Covid.sct.clean, cols = refine_cluster_color, shuffle = T, label = T, pt.size = 0.05) + 
  ggtitle("Major Type") +
  theme_void() + NoLegend() +
  fig.fix;fig

fig = DimPlot(Covid.sct.clean, cols = refine_cluster_color, shuffle = T, label = F, pt.size = 0.05) + 
  ggtitle("Major Type") +
  theme_void() + NoLegend() +
  fig.fix;fig

# Figure4 Pannel B
Seurat_test = subset(Covid.sct.clean, (`percent.mt` <= 5) & (`percent.rp` <= 25) & (nFeature_RNA <= 3200) )
Seurat_test$T_Type = factor(Seurat_test %>% Idents(), levels = c(
  "Ex_Neuron","IPC","RGC", "In_Neuron", "OPC", "Microglia", "Pericytes", "Endothelia"));
plotC <- table(Seurat_test$Covid_State, Seurat_test$T_Type) %>% melt();colnames(plotC) <- c("Sample", "CellType", "Number")
Seurat_test$refine_cluster = Seurat_test %>% Idents()

fig <- ggplot(data = plotC, aes(x = Sample, y = Number, fill = CellType)) +
  geom_bar(stat = "identity", width=0.8,aes(group=CellType),position="stack")+
  scale_fill_manual(values=refine_cluster_color) +
  # scale_fill_npg() +
  theme_bw()+
  theme(panel.grid =element_blank()) +
  labs(x="", y="Cell Number")+
  theme(axis.text = element_text(size=12, colour = "black"))+
  theme(axis.title.y = element_text(size=12, colour = "black"))+
  theme(panel.border = element_rect(size = 1, linetype = "solid", colour = "black"))+
  theme(axis.text.x = element_text(angle = 45,hjust = 0.8, vjust = 0.6)) +
  out_cover_sol + fig.fix;fig

fig <- ggplot(data = plotC, aes(x = Sample, y = Number, fill = CellType)) +
  geom_bar(stat = "identity", width=0.8,aes(group=CellType),position="fill")+
  scale_fill_manual(values=refine_cluster_color) +
  # scale_fill_npg() +
  theme_bw()+
  theme(panel.grid =element_blank()) +
  labs(x="",y="Cell Proportion")+
  scale_y_continuous(labels = percent)+ 
  theme(axis.text = element_text(size=12, colour = "black"))+
  theme(axis.title.y = element_text(size=12, colour = "black"))+
  theme(panel.border = element_rect(size = 1, linetype = "solid", colour = "black"))+
  theme(axis.text.x = element_text(angle = 45,hjust = 0.8, vjust = 0.6)) +
  out_cover_sol + fig.fix;fig

# Figure4 Pannel C
Idents(Seurat_test) <- 'Covid19' 
sce_markers <- FindMarkers(Seurat_test,ident.1 = 'negative',ident.2 = 'positive')
sce_markers <- sce_markers %>% tibble::rownames_to_column('gene')
sce_markers_sel <- subset(sce_markers,p_val_adj <= 0.05)
go_organism <- "org.Hs.eg.db"
kegg_organism <- 'hsa'
ids = bitr(sce_markers_sel$gene, fromType="SYMBOL", toType="ENTREZID", OrgDb=go_organism)
sce_markers_sel = merge(sce_markers_sel, ids, by.x='gene', by.y='SYMBOL')
sce_markers_sel$group <- factor(ifelse(sce_markers_sel$avg_log2FC < 0, -1, 1), levels = c('down-regulated', 'up-regulated'))
gcSample = split(sce_markers_sel$ENTREZID, sce_markers_sel$group)
KEGG <- compareCluster(gcSample, fun="enrichKEGG", organism = kegg_organism, pvalueCutoff=0.05)
write.xlsx(KEGG@results,file = '')

##Draw the picture
kegg_sel <- read.xlsx('') 
kegg_sel$GeneRatio[which(stringr::str_detect(kegg_sel$group,'^down'))] <- -1 * kegg_sel$GeneRatio
p <- ggplot(kegg_sel, aes(x = reorder(Description, GeneRatio), y = GeneRatio, fill = group)) +
  geom_bar(stat = "identity") + 
  labs(y = 'GeneRatio', fill = 'Group') +
  coord_flip() + 
  geom_text(aes(y = 0, label = Description, hjust = as.numeric(GeneRatio > 0))) +  
  theme(legend.position = c(0.25, 0.9)) + 
  guides(fill = guide_legend(ncol = 1)) +
  scale_fill_manual(values = c("#82D2BE", "#CD64B4"))  
ggsave2(filename = 'KEGG.pdf',plot = p)

# Figure4 Pannel D
Seurat_test$CellTypeNew_withID <- paste(Seurat_test$CellTypeNew,Seurat_test$Covid19,sep = '_') #需要师姐帮忙确定CellTypeNew和Covid19这两列
Idents(Seurat_test) <- 'CellTypeNew_withID'
p <- DotPlot(Seurat_test,cols = c('#281E50','#F0B932'),features = c('IL7','IL7R')) + 
  coord_flip()   

# FigureS3 Pannel A
fig = DimPlot(Covid.sct.clean, split.by = "Covid_State", shuffle = T, label = F, cols = Covid_color, pt.size = 0.05) + 
  ggtitle("Covid State") +
  # scale_color_npg() +
  theme_void() +
  fig.fix;fig

# FigureS3 Pannel B
df = Seurat_test@meta.data[, c('nCount_RNA', 'nFeature_RNA', paste0("percent.", c("mt", "rp", "ercc")), "Covid_State")]
df$color = Covid_color[df$Covid_State]

fig = ggplot(data = df,
             aes(x=Covid_State, y=`nFeature_RNA`, fill=Covid_State)) +
  geom_half_violin(side = "r", color=NA, alpha=0.55) +
  geom_half_boxplot(side = "r", errorbar.draw = T, width=0.15, linewidth=0.5, outlier.size = 0) +
  geom_half_point_panel(side = "l", shape=21, size=0.4, color="white") +
  scale_fill_manual(values = Covid_color) +
  # scale_y_continuous(limits = c(0, 100), expand = c(0, 0)) +
  # scale_x_discrete(labels = c('Artio','Carn','Prim','Rod','Others')) +
  labs(y="Number of genes per cell", x=NULL) +
  theme_classic() +
  theme(legend.position = "bottom",
        axis.title = element_text(size = 16, color = "black"),
        axis.text = element_text(size=13, color = "black"),
        axis.text.x = element_text(angle = 45,hjust = 0.8, vjust = 0.6)) +
  out_cover_sol +
  fig.fix + NoLegend();fig
fig = ggplot(data = df,
             aes(x=Covid_State, y=`nCount_RNA`, fill=Covid_State)) +
  geom_half_violin(side = "r", color=NA, alpha=0.55) +
  geom_half_boxplot(side = "r", errorbar.draw = T, width=0.15, linewidth=0.5, outlier.size = 0) +
  geom_half_point_panel(side = "l", shape=21, size=0.4, color="white") +
  scale_fill_manual(values = Covid_color) +
  # scale_y_continuous(limits = c(0, 100), expand = c(0, 0)) +
  # scale_x_discrete(labels = c('Artio','Carn','Prim','Rod','Others')) +
  labs(y="Number of UMIs per cell", x=NULL) +
  theme_classic() +
  theme(legend.position = "bottom",
        axis.title = element_text(size = 16, color = "black"),
        axis.text = element_text(size=13, color = "black"),
        axis.text.x = element_text(angle = 45,hjust = 0.8, vjust = 0.6)) +
  out_cover_sol +
  fig.fix + NoLegend();fig
fig = ggplot(data = df,
             aes(x=Covid_State, y=`percent.mt`, fill=Covid_State)) +
  geom_half_violin(side = "r", color=NA, alpha=0.55) +
  geom_half_boxplot(side = "r", errorbar.draw = T, width=0.15, linewidth=0.5, outlier.size = 0) +
  geom_half_point_panel(side = "l", shape=21, size=0.4, color="white") +
  scale_fill_manual(values = Covid_color) +
  scale_y_continuous(limits = c(0, 6), expand = c(0, 0), n.breaks = 6) +
  # scale_x_discrete(labels = c('Artio','Carn','Prim','Rod','Others')) +
  labs(y="Percentage of Mitochondrial Counts per cell (%)", x=NULL) +
  theme_classic() +
  theme(legend.position = "bottom",
        axis.title = element_text(size = 16, color = "black"),
        axis.text = element_text(size=13, color = "black"),
        axis.text.x = element_text(angle = 45,hjust = 0.8, vjust = 0.6)) +
  out_cover_sol +
  fig.fix + NoLegend();fig
fig = ggplot(data = df,
             aes(x=Covid_State, y=`percent.rp`, fill=Covid_State)) +
  geom_half_violin(side = "r", color=NA, alpha=0.55) +
  geom_half_boxplot(side = "r", errorbar.draw = T, width=0.15, linewidth=0.5, outlier.size = 0) +
  geom_half_point_panel(side = "l", shape=21, size=0.4, color="white") +
  scale_fill_manual(values = Covid_color) +
  scale_y_continuous(limits = c(0, 26), expand = c(0, 0)) +
  # scale_x_discrete(labels = c('Artio','Carn','Prim','Rod','Others')) +
  labs(y="Percentage of Ribosomal Counts per cell (%)", x=NULL) +
  theme_classic() +
  theme(legend.position = "bottom",
        axis.title = element_text(size = 16, color = "black"),
        axis.text = element_text(size=13, color = "black"),
        axis.text.x = element_text(angle = 45,hjust = 0.8, vjust = 0.6)) +
  out_cover_sol +
  fig.fix + NoLegend();fig

# FigureS3 Pannel C
DefaultAssay(Seurat_test) = "RNA"
Idents(Seurat_test) = factor(Seurat_test$refine_cluster, levels = c(
  "RGC", "IPC", "Ex_Neuron", "In_Neuron", "Microglia", "OPC", "Pericytes", "Endothelia"))
markers_test = FindAllMarkers(Seurat_test, only.pos = T, logfc.threshold = log2(1.2), return.thresh = 1e-5)
test.markers = markers_test;table(test.markers$cluster) 
test.markers_filtered = test.markers %>% group_by(cluster) %>%
  # slice_max(n = 10, order_by = avg_log2FC, with_ties = F) %>% as.data.frame()
  slice_min(n = 10, order_by = p_val_adj, with_ties = F) %>% as.data.frame()
Anno = c(# RGC
  "SOX2", "HES5", "HES1", "HOPX", 
  # "NES", "CLU", "RGS6", "MKI-67",
  # IPC
  "EGFR", "MKI67", "ASCL1", "EOMES",
  # "EGFR", "FABP7","DLX2",
  # Excited  GRIN1 GRIN2B BCL11B
  "SATB2", "NEUROD6", "NEUROD2", "SLC14A2",
  # interneuron  1 13 "NEUROD2",
  "GAD2", "CALB2", "SLC32A1", "DLX6-AS1",
  # Microglia
  "TMEM119", "P2RY12", "ABI3", "CX3CR1",
  # OPC PCDH15 EGFR
  "OLIG1", "OLIG2", "PDGFRA", "SOX10",
  # Pericyte
  "KCNJ8", "ANPEP", "RGS5", "PDGFRB",
  # Endothelium
  "MFSD2A", "TIE1", "CD34", "PECAM1");
gene_selected = c(union(test.markers_filtered$gene[1:10], Anno[1:4]),
                  union(test.markers_filtered$gene[11:20], Anno[5:8]),
                  union(test.markers_filtered$gene[21:30], Anno[9:12]),
                  union(test.markers_filtered$gene[31:40], Anno[13:16]),
                  union(test.markers_filtered$gene[41:50], Anno[17:20]),
                  union(test.markers_filtered$gene[51:60], Anno[21:24]),
                  union(test.markers_filtered$gene[61:70], Anno[25:28]),
                  union(test.markers_filtered$gene[61:70], Anno[29:32])
)
fig = AverageHeatmap(object = Seurat_test, assays = "RNA", slot = "data", 
                     markerGene = gene_selected, annoCol = TRUE, 
                     myanCol = refine_cluster_color[levels(Idents(Seurat_test))],
                     # width = 12, height = 20,
                     fontsize = 10,
                     showRowNames = F, markGenes = Anno);fig
