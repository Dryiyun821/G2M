
####Figure 2A######
if(T){
  library(survival)
  library(survminer)
  dev <- read_excel('G2M_Cli.csv')
  data.survdiff <- survdiff(Surv(dev$OS,dev$Status) ~dev$G2M)
  p.val = 1 - pchisq(data.survdiff$chisq, length(data.survdiff$n)- 1)
  HR = (data.survdiff$obs[2]/data.survdiff$exp[2])/(data.survdiff$obs[1]/data.survdiff$exp[1])
  up95 = exp(log(HR) + qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1]))
  low95 = exp(log(HR) - qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1]))
  fit<-survfit(Surv(OS,Status)~G2M,dev)
  
  ggsurvplot(
    fit,
    data = dev,
    
    # 生存曲线设置
    conf.int = TRUE,  # 是否显示置信区间
    conf.int.style = "step",  # 置信区间样式
    surv.median.line = "hv",  # 中位生存时间线
    
    # p值设置
    pval = TRUE,
    pval.size = 5,
    pval.coord = c(8, 0.4),  # p值位置 (可选)
    # 或者自定义p值文本
    # pval = paste0("p = ", round(p_value, 3)),
    
    # 坐标轴标签
    xlab = "Time (months)",
    ylab = "Overall Survival  Probability",
    break.time.by = 24,  # x轴刻度间隔
    
    # 风险表设置
    risk.table = TRUE,
    risk.table.title = "",
    risk.table.col = "strata",
    risk.table.height = 0.31,  # 风险表高度
    risk.table.y.text.col = TRUE,  # 风险表y轴文字颜色
    risk.table.y.text = FALSE,  # 是否显示风险表y轴文字
    
    # 图例设置
    legend.title = "G2M checkpoint",
    legend.labs = c("Low", "High"),  # 确保与fit中的分组顺序一致
    legend.position = c(0.8, 0.25),  # 图例位置 (x, y)
    
    # 配色
    palette = "jco",  # journal of clinical oncology配色
    # 主题
    ggtheme = theme_bw(),
    
    # 其他美化
    font.main = 14,
    font.x = 12,
    font.y = 12,
    font.tickslab = 10,
    censor = TRUE,  # 是否显示删失点
    censor.shape = "+",  # 删失点形状
    censor.size = 4,  # 删失点大小
    
    # 分面设置 (如果有多个变量)
    # facet.by = "another_variable",
    
    # 添加P值比较方法
    pval.method = FALSE  # 是否添加检验方法名称
  )
}

####Figure 2B###
if(T){
  ###ssGSEA分析###
  output<-read.csv('ssGSEA_G2M.csv')
  output$G2M <-factor(output$G2M,labels = c('Low','High'))
  #数据宽转长
  library(reshape2)
  library(ggplot2)
  library(ggpubr)
  library(dplyr)
  pdata_melt <- melt(output,
                     id.vars = c("ID",'G2M'),
                     variable.name ="Immune.cells",
                     value.name = "Immune.Score")
  
  attach(pdata_melt)
  # 输出到文件
  #write.csv(pdata_melt, "ssGSEA_input_TCGA.csv", row.names = F, quote = F)
  # 使用ggplo2画图
  library(ggsci)
  c <- ggplot(pdata_melt,
              aes(x=Immune.cells, y=Immune.Score, fill =MATH)) + #按类给边框着色
    geom_boxplot(notch = F, alpha = 0.95, 
                 outlier.shape = 16,
                 outlier.colour = "black", #outlier点用黑色
                 outlier.size = 0.65) +
    #自定义配色
    #scale_fill_manual(values= c("#E31A1C","#E7B800","#2E9FDF")) +
    #ggtitle("Gene signature score", "stratified by TME-cluster") + 
    theme_classic() +
    theme(axis.text.x = element_text(angle = 60, hjust = 1,size = 10), 
          axis.text.y = element_text(angle = 90, size = 12),
          axis.title.y = element_text(angle = 90, size = 15)) +
    theme(legend.position = "top")+
    scale_fill_jco()+
    labs(x = "Immune cells", y = "Enrichment score")+
    theme(plot.title = element_text(face = "bold",
                                    margin = margin(10, 0, 10, 0),
                                    size = 14))
  
  # 标注p value
  c + stat_compare_means(aes(label = paste0("p = ", ..p.format..)))
  # 标注*
  # `****` = 1e-04, `***` = 0.001, `**` = 0.01, `*` = 0.05, ns = 1
  p <- c + stat_compare_means(label = "p.signif")
  p
}

###Figure 2C####
if(T){
  library(ggplot2)
  library(ggridges)
  library(dplyr)
  
  # 输入数据：dev
  # 两列：Group (MPR/NMPR), G2M_checkpoint
  dev <- read.csv('')
  dev$Group <- factor(dev$Group, levels = c("MPR", "NMPR"))
  dev$y <- as.numeric(dev$Group)
  
  cols <- c("MPR" = "#0073C2FF", "NMPR" = "#EFC000FF")
  
  # Wilcoxon检验
  pval <- wilcox.test(G2M_checkpoint ~ Group, data = dev)$p.value
  sig <- ifelse(pval < 0.001, "***",
                ifelse(pval < 0.01, "**",
                       ifelse(pval < 0.05, "*", "ns")))
  
  gg <- ggplot(dev, aes(fill = Group, color = Group)) +
    
    # 山峦密度图
    geom_density_ridges(
      aes(x = G2M_checkpoint, y = y + 0.15, group = Group),
      scale = 0.35,
      alpha = 1,
      color = NA
    ) +
    
    # 箱线图
    geom_boxplot(
      aes(x = G2M_checkpoint, y = y, group = Group),
      orientation = "y",
      width = 0.055,
      outlier.shape = NA,
      color = "black",
      linewidth = 0.5
    ) +
    
    # 散点
    geom_point(
      aes(x = G2M_checkpoint, y = y - 0.12),
      position = position_jitter(width = 0, height = 0.035),
      size = 0.9,
      alpha = 0.9
    ) +
    
    # 显著性括号
    annotate("segment", x = 2.55, xend = 2.55,
             y = 1, yend = 2, linewidth = 0.5) +
    annotate("segment", x = 2.49, xend = 2.55,
             y = 1, yend = 1, linewidth = 0.5) +
    annotate("segment", x = 2.49, xend = 2.55,
             y = 2, yend = 2, linewidth = 0.5) +
    annotate("text", x = 2.65, y = 1.5,
             label = sig, size = 4) +
    
    scale_fill_manual(values = cols) +
    scale_color_manual(values = cols) +
    scale_y_continuous(
      breaks = c(1, 2),
      labels = c("MPR", "NMPR")
    ) +
    coord_cartesian(xlim = c(-1.4, 2.8),
                    ylim = c(0.4, 2.65)) +
    
    labs(x = "G2M checkpoint", y = NULL) +
    
    theme_classic(base_size = 14) +
    theme(
      legend.position = "none",
      axis.text.y = element_text(
        face = "bold.italic", color = "black"
      ),
      axis.text.x = element_text(color = "black"),
      axis.title.x = element_text(
        face = "bold.italic", size = 16
      ),
      axis.ticks.y = element_blank(),
      axis.line = element_line(linewidth = 0.6)
    )
  
  gg
}

###Figure 2D###
if(T){
  library(qs)
  library(Seurat)
  sub_sce <- qread('')
  sub_sce$celltype <- factor(sub_sce$celltype,levels=c( "Epithelial cells",
                                              "T/NK cells",
                                              "B cells",
                                              "Plasma cells",
                                              "Proliferative cells",
                                              "Myeloid cells",
                                              "Mast cells",
                                              "Neutrophil cells",
                                              "Fibroblast cells",
                                              "Endothelial cells"))
  
  sub_sce$Cluster <-sub_sce$celltype
  table(sub_sce$Cluster)
  g.colSet1 <- c("#00468BFF", "#ED0000FF", "#42B540FF", "#0099B4FF", "#925E9FFF",
                 "#FDAF91FF", "#AD002AFF", "#ADB6B6FF", "#1B1919FF", "#F39B7FFF",
                 "#8491B4FF", "#91D1C2FF", "#DC0000FF", "#7E6148FF", "#B09C85FF",
                 "#3C5488FF", "#E64B35FF", "#00A087FF", "#4DBBD5FF", "#6C5483FF")
  sub_sce<-AddMetaData(sub_sce,sub_sce@reductions$umap@cell.embeddings,col.name = colnames(sub_sce@reductions$umap@cell.embeddings))
  class_avg <- sub_sce@meta.data %>%
    group_by(Cluster) %>%
    dplyr::summarise(
      umap_1 = median(umap_1),
      umap_2 = median(umap_2)
    )
  
  sub_sce$Cluster <- sub_sce$celltype
  ggplot(sub_sce@meta.data ,aes(x=umap_1,y=umap_2,color=Cluster))+
    # geom_point(aes(color=SubCluster),size=0.1) +
    geom_point(aes(color=Cluster),size = 0.65,shape = 16, stroke = 0) +
    scale_color_manual(breaks = c(levels(sub_sce$Cluster)),
                       labels= c( "Epithelial cells",
                                  "T/NK cells",
                                  "B cells",
                                  "Plasma cells",
                                  "Proliferative cells",
                                  "Myeloid cells",
                                  "Mast cells",
                                  "Neutrophil cells",
                                  "Fibroblast cells",
                                  "Endothelial cells"                       ),
                       values = g.colSet1)+
    ggtitle("Cluster")+
    geom_text_repel(
      aes(label = Cluster), 
      data = class_avg,
      color = "black",
      size = 3,
      box.padding = 0.5,  # 标签与点之间的间距
      point.padding = 0.5,  # 标签与其他标签之间的间距
      max.overlaps = Inf,  # 允许所有标签显示
      segment.color = "grey50",  # 连接线颜色
      segment.size = 0.3 ) + # 连接线粗细
    theme_classic()+
    #theme(text=element_text(family="Arial",size=10)) +
    theme(plot.title = element_text(size = 12,hjust = 0.5),
          ##axis.title = element_text(size = 10),
          panel.border = element_rect(fill=NA,color="grey",size=0.5),
          axis.text = element_blank(),
          axis.title = element_blank(),
          axis.line = element_blank(),
          axis.ticks = element_blank(),
          legend.position = "right",
          legend.text = element_text(size = 9),
          strip.text  = element_text(size = 10),
          strip.background = element_blank())+ 
    guides(colour = guide_legend(override.aes = list(size=3)))
  
  sub_sce$RECIST <- factor(sub_sce$RECIST,levels = c('PreNMPR','PreMPR','PostNMPR','PostMPR'))
  
  
  
  ggplot(sub_sce@meta.data ,aes(x=umap_1,y=umap_2,color=SubCluster))+
    geom_point(aes(color=SubCluster),size=0.5) +
    # geom_point(aes(color=SubCluster),size = 0.65,shape = 16, stroke = 0) +
    scale_color_manual(breaks = c(levels(sub_sce$SubCluster)),
                       labels= c( "Epithelial cells",
                                  "T/NK cells",
                                  "B cells",
                                  "Plasma cells",
                                  "Proliferative cells",
                                  "Myeloid cells",
                                  "Mast cells",
                                  "Neutrophil cells",
                                  "Fibroblast cells",
                                  "Endothelial cells" 
                       ),
                       values = g.colSet1)+
    ggtitle("Cluster")+
    # geom_text_repel(
    #   aes(label = SubCluster), 
    #   data = class_avg,
    #   color = "black",
    #   size = 3,
    #   box.padding = 0.5,  # 标签与点之间的间距
    #   point.padding = 0.5,  # 标签与其他标签之间的间距
    #   max.overlaps = Inf,  # 允许所有标签显示
    #   segment.color = "grey50",  # 连接线颜色
    #   segment.size = 0.3 ) + # 连接线粗细
    theme_classic()+
    #theme(text=element_text(family="Arial",size=10)) +
    theme(plot.title = element_text(size = 12,hjust = 0.5),
          ##axis.title = element_text(size = 10),
          panel.border = element_rect(fill=NA,color="grey",size=0.5),
          axis.text = element_blank(),
          axis.title = element_blank(),
          axis.line = element_blank(),
          axis.ticks = element_blank(),
          legend.position = "none",
          legend.text = element_text(size = 9),
          strip.text  = element_text(size = 10),
          strip.background = element_blank())+ 
    guides(colour = guide_legend(override.aes = list(size=3)))+ 
    facet_wrap(~RECIST,ncol=2)
}

###Figure 2E###
if(T){
###sub_sce 为Seurat文件###
  library(Seurat)
  library(ggplot2)
  
  # 1. Marker基因
  markers <- c(
    "EPCAM", "KRT19", "KRT8",           # Epithelial
    "CD3D", "CD3E", "NKG7",             # T/NK
    "MS4A1", "CD19", "CD79A",           # B
    "MZB1", "JCHAIN", "IGKC",           # Plasma
    "MKI67", "STMN1", "PCNA",           # Proliferative
    "LYZ", "CD68", "APOC1",             # Myeloid
    "CPA3", "TPSAB1", "MS4A2",          # Mast
    "FCGR3B", "CSF3R", "S100A8",        # Neutrophil
    "COL1A1", "DCN", "LUM",             # Fibroblast
    "CLDN5", "PECAM1", "VWF"           # Endothelial
  )
  
  # 2. 设置细胞类型顺序（从下到上）
  cell_order <- c(
    "Epithelial cells",
    "T/NK cells",
    "B cells",
    "Plasma cells",
    "Proliferative cells",
    "Myeloid cells",
    "Mast cells",
    "Neutrophil cells",
    "Fibroblast cells",
    "Endothelial cells"
  )
  
  sub_sce$celltype <- factor(sub_sce$celltype, levels = cell_order)
  
  DefaultAssay(sub_sce) <- "RNA"
  
  # 3. 绘制DotPlot
  p <- DotPlot(
    sub_sce,
    features = markers,
    group.by = "celltype",
    dot.scale = 9,
    scale = TRUE
  ) +
    scale_color_gradientn(
      colors = c("#E8EEB2", "#4BAEB5", "#172653"),
      name = "Average"
    ) +
    labs(x = NULL, y = NULL, size = "Percent") +
    guides(
      size = guide_legend(order = 1),
      color = guide_colorbar(order = 2)
    ) +
    theme_bw(base_size = 12) +
    theme(
      panel.grid.major = element_line(
        color = "#E5E5E5", linewidth = 0.4
      ),
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(
        angle = 50, hjust = 1, vjust = 1,
        color = "#444444", size = 10
      ),
      axis.text.y = element_text(
        color = "#444444", size = 11
      ),
      legend.position = "right",
      panel.border = element_rect(
        color = "#BDBDBD", linewidth = 0.6
      )
    )
  
  p
}

###Figure 2F###
if(T){
  sub_sce$RECIST <- factor(sub_sce$RECIST,levels = c('PreNMPR','PreMPR','PostNMPR','PostMPR'))
  sub_sce$SubCluster <- sub_sce$celltype
  Cellratio <- prop.table(table(sub_sce$SubCluster, sub_sce$RECIST), margin = 2) %>% 
    as.data.frame() %>% 
    set_names(c("celltype", "group", "ratio"))
  
  ## 2. 计算各组细胞比例
  library(tidyverse)
  library(ggalluvial)
  library(scales)
  ## 图例从上到下的顺序
  ## 对应堆叠柱中从上到下的顺序
  group_order <-c('PreNMPR','PreMPR','PostNMPR','PostMPR')
  celltype_order <- c(
    "Epithelial cells",
    "T/NK cells",
    "B cells",
    "Plasma cells",
    "Proliferative cells",
    "Myeloid cells",
    "Mast cells",
    "Neutrophil cells",
    "Fibroblast cells",
    "Endothelial cells"
  )
  
  
  ## 检查每组比例之和是否为1
  Cellratio %>%
    group_by(group) %>%
    summarise(
      total_ratio = sum(ratio),
      .groups = "drop"
    )
  colour <- c(
    "Epithelial cells"     = "#00468BFF",
    "T/NK cells"          = "#ED0000FF",
    "B cells"             = "#42B540FF",
    "Plasma cells"        = "#0099B4FF",
    "Proliferative cells" = "#925E9FFF",
    "Myeloid cells"       = "#FDAF91FF",
    "Mast cells"          = "#AD002AFF",
    "Neutrophil cells"    = "#ADB6B6FF",
    "Fibroblast cells"    = "#1B1919FF",
    "Endothelial cells"   = "#F39B7FFF"
  )
  
  p <- ggplot(
    Cellratio,
    aes(
      x = group,
      y = ratio,
      stratum = celltype,
      alluvium = celltype,
      fill = celltype
    )
  ) +
    
    ## 组间流线
    geom_flow(
      stat = "flow",
      width = 0.60,
      alpha = 0.28,
      knot.pos = 0.35,
      curve_type = "sine",
      color = "white",
      linewidth = 0.55,
      aes.bind = "alluvia",
      lode.guidance = "frontback",
      show.legend = FALSE
    ) +
    
    ## 每组堆叠柱
    geom_stratum(
      width = 0.60,
      color = "white",
      linewidth = 0.55,
      alpha = 1
    ) +
    
    scale_fill_manual(
      values = colour,
      breaks = celltype_order,
      limits = celltype_order,
      drop = FALSE
    ) +
    
    scale_y_continuous(
      limits = c(0, 1),
      breaks = seq(0, 1, by = 0.25),
      labels = number_format(
        accuracy = 0.01
      ),
      expand = expansion(
        mult = c(0, 0)
      )
    ) +
    
    scale_x_discrete(
      limits = group_order,
      expand = expansion(
        mult = c(0.07, 0.07)
      )
    ) +
    
    labs(
      x = NULL,
      y = "Percentage (%)",
      fill = NULL
    ) +
    
    coord_cartesian(
      clip = "off"
    ) +
    
    theme_minimal(
      base_size = 15,
      base_family = "Arial"
    ) +
    
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      
      panel.grid.major.y = element_line(
        color = "#D9D9D9",
        linewidth = 0.55
      ),
      
      axis.line = element_blank(),
      
      axis.ticks = element_blank(),
      
      axis.title.y = element_text(
        size = 16,
        color = "black",
        margin = margin(r = 12)
      ),
      
      axis.text.y = element_text(
        size = 13,
        color = "black"
      ),
      
      axis.text.x = element_text(
        size = 14,
        face = "bold",
        color = "black",
        margin = margin(t = 7)
      ),
      
      legend.position = "right",
      
      legend.text = element_text(
        size = 13,
        color = "black"
      ),
      
      legend.key = element_rect(
        fill = NA,
        color = NA
      ),
      
      legend.key.width = unit(
        0.55,
        "cm"
      ),
      
      legend.key.height = unit(
        0.55,
        "cm"
      ),
      
      legend.spacing.y = unit(
        0.02,
        "cm"
      ),
      
      plot.margin = margin(
        t = 10,
        r = 15,
        b = 10,
        l = 10
      )
    ) +
    
    guides(
      fill = guide_legend(
        title = NULL,
        reverse = FALSE,
        override.aes = list(
          alpha = 1,
          color = NA
        )
      )
    )
  
  p
}

####Figure 2G####
if(T){
  library(Seurat)
  library(dplyr)
  library(pheatmap)
  group_order <- c(
    "PreNMPR",
    "PreMPR",
    "PostNMPR",
    "PostMPR"
  )
  
  sub_sce$RECIST <- factor(
    sub_sce$RECIST,
    levels = group_order
  )
  observed_matrix <- table(
    sub_sce$RECIST,
    sub_sce$Cluster
  )
  
  observed_matrix <- as.matrix(
    observed_matrix
  )
  
  # 确保分组顺序正确
  observed_matrix <- observed_matrix[
    group_order,
    ,
    drop = FALSE
  ]
  
  observed_matrix
  # ============================================================
  # 4. 计算随机分布情况下的期望细胞数
  #
  # expected_ij =
  # 第i组细胞总数 × 第j类细胞总数 / 全部细胞总数
  # ============================================================
  
  expected_matrix <- outer(
    rowSums(observed_matrix),
    colSums(observed_matrix)
  ) / sum(observed_matrix)
  
  rownames(expected_matrix) <- rownames(
    observed_matrix
  )
  
  colnames(expected_matrix) <- colnames(
    observed_matrix
  )
  
  
  # ============================================================
  # 5. 计算Observed/Expected ratio
  # ============================================================
  
  roe_matrix <- observed_matrix /
    expected_matrix
  
  # 转置：
  # 行为细胞亚群，列为RECIST分组
  heatdata <- t(
    roe_matrix
  )
  
  heatdata
  # ============================================================
  # 6. 设置颜色范围
  # ============================================================
  
  color_min <- 0
  color_max <- 2.5
  
  # 用于着色的矩阵
  heatdata_color <- heatdata
  
  heatdata_color[
    heatdata_color < color_min
  ] <- color_min
  
  heatdata_color[
    heatdata_color > color_max
  ] <- color_max
  
  # 用于显示真实数值的字符矩阵
  heatdata_number <- matrix(
    sprintf(
      "%.2f",
      heatdata
    ),
    nrow = nrow(heatdata),
    ncol = ncol(heatdata),
    dimnames = dimnames(heatdata)
  )
  
  heatmap_colors <- colorRampPalette(
    c(
      "#0044BB",
      "#92C5DE",
      "#F7F7F7",
      "#F4A582",
      "#C00038"
    )
  )(100)
  
  heatmap_breaks <- seq(
    color_min,
    color_max,
    length.out = 101
  )
  # ============================================================
  # 7. 绘制热图
  # ============================================================
  
  p_ROE <- pheatmap::pheatmap(
    mat = heatdata_color,
    
    color = heatmap_colors,
    breaks = heatmap_breaks,
    
    # 在Pre和Post之间添加间隔
    gaps_col = 2,
    
    # 显示原始真实数值
    display_numbers = heatdata_number,
    number_color = "black",
    fontsize_number = 11,
    
    # 单元格设置
    border_color = NA,
    cellwidth = 68,
    cellheight = 48,
    
    # 标题与字体
    main = "Ro/e",
    fontsize = 15,
    fontsize_row = 14,
    fontsize_col = 13,
    
    # 参考图中列名为竖直方向
    angle_col = "90",
    
    show_rownames = TRUE,
    show_colnames = TRUE,
    
    # 聚类设置
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    
    clustering_distance_rows = "euclidean",
    clustering_method = "complete",
    
    treeheight_row = 60,
    treeheight_col = 0,
    
    # 图例
    legend = TRUE,
    legend_breaks = c(
      0,
      0.5,
      1,
      1.5,
      2,
      2.5
    ),
    
    legend_labels = c(
      "0",
      "0.5",
      "1",
      "1.5",
      "2",
      "2.5"
    )
  )
  
  p_ROE
}

####Figure 2H###
if(T){
  #####-----------------------========================------------------------------
  Enrich_type<-"Hallmark"   ##GO_BP,GO_CC,GO_MF,Reactome,KEGG,
  # Type="All"
  Select_cluster_list<- cluster_list
  total_table<-data.frame()
  enriched_pathway<-c()
  for (cluster in Select_cluster_list){
    path<-paste("",cluster,"/",Enrich_type,"_enrich.xls",sep="")
    enrich_table<-read.table(path,sep="\t",header = T,row.names = 1)
    enrich_table$Cluster<-cluster
    
    if (nrow(enrich_table)>5){
      enrich_table<-enrich_table[which(enrich_table$p.adjust<0.05),]
      enrich_table<-head(enrich_table,5)
    }
    
    enriched_pathway<-c(enriched_pathway,enrich_table$Description)
    total_table<-rbind(total_table,enrich_table)
    
  }
  
  total_table$Description<-as.character(total_table$Description)
  total_table$Description<-factor(total_table$Description,levels=rev(unique(total_table$Description)))
  total_table$Cluster<-factor(total_table$Cluster,levels=Select_cluster_list)
  
  S1<- ggplot(total_table, aes(x= Cluster, y=Description, size=Count, color=p.adjust, group=Cluster)) + 
    geom_point(alpha = 0.8) + 
    theme_classic()+
    scale_color_gradient(low = "#FF7F0E",  high = "#1F77B4", space = "Lab")+
    scale_size(range = c(2, 6))+
    ggtitle(paste0(Enrich_type," in CTL"))+
    theme_bw()+
    theme(plot.title = element_text(size = 10,hjust = 0.5),
          axis.title.x = element_text(size = 10),
          axis.title.y = element_text(size = 10),
          axis.text.x = element_text(size = 10,angle = 45,hjust = 1,vjust = 1),
          axis.text.y = element_text(size = 10),
          legend.title = element_text(size = 10),
          
          legend.text = element_text(size = 10)) 
  
  S1 
}

###Figure 2I####
if(T){
library(Seurat)
library(msigdbr)
library(dplyr)
library(ggplot2)
library(ggrepel)

# 1. 获取HALLMARK_G2M_CHECKPOINT基因集
hallmark <- msigdbr(
  species = "Homo sapiens",
  collection = "H"
)

genes <- hallmark %>%
  filter(gs_name == "HALLMARK_G2M_CHECKPOINT") %>%
  pull(gene_symbol) %>%
  unique()

# 2. 计算G2M checkpoint评分
DefaultAssay(sub_sce) <- "RNA"

# Seurat v5多layer数据处理
if (inherits(sub_sce[["RNA"]], "Assay5") &&
    length(Layers(sub_sce[["RNA"]], search = "counts")) > 1) {
  sub_sce <- JoinLayers(sub_sce, assay = "RNA")
}

sub_sce <- NormalizeData(sub_sce, verbose = FALSE)

genes <- intersect(genes, rownames(sub_sce))
stopifnot(length(genes) >= 10)

set.seed(123)
sub_sce <- AddModuleScore(
  sub_sce,
  features = list(genes),
  name = "G2M_checkpoint"
)

sub_sce$G2M_checkpoint <- sub_sce$G2M_checkpoint1

# 3. 提取UMAP坐标及细胞类型
df <- as.data.frame(Embeddings(sub_sce, "umap"))
colnames(df)[1:2] <- c("UMAP_1", "UMAP_2")

df$celltype <- sub_sce@meta.data[rownames(df), "celltype"]

centers <- df %>%
  group_by(celltype) %>%
  summarise(
    UMAP_1 = median(UMAP_1),
    UMAP_2 = median(UMAP_2),
    .groups = "drop"
  )

# 4. 绘图
p <- FeaturePlot(
  sub_sce,
  features = "G2M_checkpoint",
  reduction = "umap",
  cols = c("#E8EEB2", "#4BAEB5", "#172653"),
  pt.size = 0.1,
  order = TRUE
) +
  geom_text_repel(
    data = centers,
    aes(x = UMAP_1, y = UMAP_2, label = celltype),
    inherit.aes = FALSE,
    size = 3.5,
    color = "black",
    max.overlaps = Inf
  ) +
  labs(color = "G2M checkpoint") +
  theme_classic(base_size = 14) +
  theme(
    axis.title = element_text(size = 14),
    axis.line = element_line(linewidth = 0.7),
    legend.title = element_text(size = 11)
  )

p
}

###Figure S1##
if(T){
  library(survival)
  library(survminer)
  dev <- read_excel('G2M_GSE30219.csv')
  data.survdiff <- survdiff(Surv(dev$OS,dev$Status) ~dev$G2M)
  p.val = 1 - pchisq(data.survdiff$chisq, length(data.survdiff$n)- 1)
  HR = (data.survdiff$obs[2]/data.survdiff$exp[2])/(data.survdiff$obs[1]/data.survdiff$exp[1])
  up95 = exp(log(HR) + qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1]))
  low95 = exp(log(HR) - qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1]))
  fit<-survfit(Surv(OS,Status)~G2M,dev)
  
  ggsurvplot(
    fit,
    data = dev,
    
    # 生存曲线设置
    conf.int = TRUE,  # 是否显示置信区间
    conf.int.style = "step",  # 置信区间样式
    surv.median.line = "hv",  # 中位生存时间线
    
    # p值设置
    pval = TRUE,
    pval.size = 5,
    pval.coord = c(8, 0.4),  # p值位置 (可选)
    # 或者自定义p值文本
    # pval = paste0("p = ", round(p_value, 3)),
    
    # 坐标轴标签
    xlab = "Time (months)",
    ylab = "Overall Survival  Probability",
    break.time.by = 24,  # x轴刻度间隔
    
    # 风险表设置
    risk.table = TRUE,
    risk.table.title = "",
    risk.table.col = "strata",
    risk.table.height = 0.31,  # 风险表高度
    risk.table.y.text.col = TRUE,  # 风险表y轴文字颜色
    risk.table.y.text = FALSE,  # 是否显示风险表y轴文字
    
    # 图例设置
    legend.title = "G2M checkpoint",
    legend.labs = c("Low", "High"),  # 确保与fit中的分组顺序一致
    legend.position = c(0.8, 0.25),  # 图例位置 (x, y)
    
    # 配色
    palette = "jco",  # journal of clinical oncology配色
    # 主题
    ggtheme = theme_bw(),
    
    # 其他美化
    font.main = 14,
    font.x = 12,
    font.y = 12,
    font.tickslab = 10,
    censor = TRUE,  # 是否显示删失点
    censor.shape = "+",  # 删失点形状
    censor.size = 4,  # 删失点大小
    
    # 分面设置 (如果有多个变量)
    # facet.by = "another_variable",
    
    # 添加P值比较方法
    pval.method = FALSE  # 是否添加检验方法名称
  )
  
  ###森林图
  library(forestmodel) # 森林图
  library(ggplot2) # 使用到里面的数据集“msleep”
  library(dplyr) # 数据处理
  library(survival) # 生存分析
  dev <- read.csv('foreatplot_GSE30219.csv')
  dev$Age <- as.numeric(dev$Age)
  dev$Gender <- factor(dev$Gender,levels=c('Female','Male'))
  dev$fT <- factor(dev$fT,levels=c('T1','T2','T3','T4'))
  dev$fN <- factor(dev$fN,levels=c('N0','N1','N2','N3'))
  dev$G2M.checkpoint <- factor(dev$G2M.checkpoint,levels=c('Low','High'))
  # 拟合生存分析模型
  surv_model <- coxph(Surv(OS, Status) ~ Age + Gender + fT+fN+G2M.checkpoint, data = dev)
  # 制作森林图
  forest_model(surv_model)
  forest_model(
    surv_model,
    format_options = forest_model_format_options(colour = "#0073C2FF", # 修改颜色
                                                 shape = 15, # 改变形状
                                                 text_size = 3.8, # 字体文本的大小
                                                 point_size = 6, # 森林图中“方框”的大小
                                                 banded = TRUE), # 如“FALSE”，可去掉森林图背景条带
    factor_separate_line = TRUE # 如为分类变量，需换行；比如Male，Female与gender不在同一行。
  )
}

###Figure S2##
if(T){
  library(Seurat)
  library(ggplot2)
  library(patchwork)
  library(Nebulosa)
  library(ggrepel)
  
  DefaultAssay(sub_sce) <- "RNA"
  
  #Figure S2A-C#
  sub_sce[["percent.mito"]] <- PercentageFeatureSet(
    sub_sce, pattern = "^MT-"
  ) / 100
  
  features <- c("nCount_RNA", "nFeature_RNA", "percent.mito")
  
  p_qc <- lapply(features, function(x) {
    VlnPlot(sub_sce, features = x, group.by = "Sample",
            pt.size = 0) +
      labs(x = "Sample ID", y = NULL, title = x) +
      theme_classic() +
      theme(
        axis.text.x = element_text(angle = 90, size = 6),
        legend.position = "none"
      )
  })
  
  # 细胞类型标签
  centers <- aggregate(
    Embeddings(sub_sce, "umap")[, 1:2],
    by = list(celltype = sub_sce$celltype),
    FUN = median
  )
  colnames(centers)[2:3] <- c("x", "y")
  
  labels <- geom_text_repel(
    data = centers,
    aes(x, y, label = celltype),
    inherit.aes = FALSE,
    size = 2.5,
    max.overlaps = Inf
  )
  
  # D：RECIST分组UMAP
  cols <- c(
    "PreNMPR"  = "#00468BFF",
    "PreMPR"   = "#ED0000FF",
    "PostNMPR" = "#42B540FF",
    "PostMPR"  = "#0099B4FF"
  )
  ###Figure S2D###
  pD <- DimPlot(
    sub_sce, group.by = "RECIST",
    cols = cols, pt.size = 0.1
  ) + labels + ggtitle("All Cluster") + theme_void()
  
  ###Figure S2E###
  # E：样本UMAP
  pE <- DimPlot(
    sub_sce, group.by = "Sample",
    pt.size = 0.1
  ) + labels + ggtitle("All Cluster") + theme_void()
  
  ###Figure S2F###
  # F：Marker基因表达密度图
  markers <- c(
    "EPCAM", "CD3D", "NKG7",
    "CD79A", "JCHAIN", "MKI67",
    "TOP2A", "C1QB", "CPA3",
    "CSF3R", "COL1A1", "VWF"
  )
  
  pF <- plot_density(
    sub_sce,
    features = markers,
    reduction = "umap",
    combine = TRUE
  )
  
  # 拼接A-F
  p <- (
    (p_qc[[1]] | p_qc[[2]]) /
      (p_qc[[3]] | pD) /
      (pE | pF)
  ) + plot_annotation(tag_levels = "A")
  
  p
}

###Figure 3A####
if(T){
  library(qs)
  library(Seurat)
  sub_sce <- qread('Epi.qs')
  sub_sce$celltype <- factor(sub_sce$celltype,levels=c(   "Epi_C0_SFTPA1",
                                                          "Epi_C1_KRT17",
                                                          "Epi_C2_TPX2",
                                                          "Epi_C3_AKR1C1",
                                                          "Epi_C4_IDO1",
                                                          "Epi_C5_SCGB1A1",
                                                          "Epi_C6_PIP",
                                                          "Epi_C7_MMP10",
                                                          "Epi_C8_AGER"))
  
  sub_sce$Cluster <-sub_sce$celltype
  table(sub_sce$Cluster)
  g.colSet1 <- c("#00468BFF", "#ED0000FF", "#42B540FF", "#0099B4FF", "#925E9FFF",
                 "#FDAF91FF", "#AD002AFF", "#ADB6B6FF", "#1B1919FF", "#F39B7FFF",
                 "#8491B4FF", "#91D1C2FF", "#DC0000FF", "#7E6148FF", "#B09C85FF",
                 "#3C5488FF", "#E64B35FF", "#00A087FF", "#4DBBD5FF", "#6C5483FF")
  sub_sce<-AddMetaData(sub_sce,sub_sce@reductions$umap@cell.embeddings,col.name = colnames(sub_sce@reductions$umap@cell.embeddings))
  class_avg <- sub_sce@meta.data %>%
    group_by(Cluster) %>%
    dplyr::summarise(
      umap_1 = median(umap_1),
      umap_2 = median(umap_2)
    )
  
  sub_sce$Cluster <- sub_sce$celltype
  ggplot(sub_sce@meta.data ,aes(x=umap_1,y=umap_2,color=Cluster))+
    # geom_point(aes(color=SubCluster),size=0.1) +
    geom_point(aes(color=Cluster),size = 0.65,shape = 16, stroke = 0) +
    scale_color_manual(breaks = c(levels(sub_sce$Cluster)),
                       labels= c(   "Epi_C0_SFTPA1",
                                    "Epi_C1_KRT17",
                                    "Epi_C2_TPX2",
                                    "Epi_C3_AKR1C1",
                                    "Epi_C4_IDO1",
                                    "Epi_C5_SCGB1A1",
                                    "Epi_C6_PIP",
                                    "Epi_C7_MMP10",
                                    "Epi_C8_AGER"                      ),
                       values = g.colSet1)+
    ggtitle("Cluster")+
    geom_text_repel(
      aes(label = Cluster), 
      data = class_avg,
      color = "black",
      size = 3,
      box.padding = 0.5,  # 标签与点之间的间距
      point.padding = 0.5,  # 标签与其他标签之间的间距
      max.overlaps = Inf,  # 允许所有标签显示
      segment.color = "grey50",  # 连接线颜色
      segment.size = 0.3 ) + # 连接线粗细
    theme_classic()+
    #theme(text=element_text(family="Arial",size=10)) +
    theme(plot.title = element_text(size = 12,hjust = 0.5),
          ##axis.title = element_text(size = 10),
          panel.border = element_rect(fill=NA,color="grey",size=0.5),
          axis.text = element_blank(),
          axis.title = element_blank(),
          axis.line = element_blank(),
          axis.ticks = element_blank(),
          legend.position = "right",
          legend.text = element_text(size = 9),
          strip.text  = element_text(size = 10),
          strip.background = element_blank())+ 
    guides(colour = guide_legend(override.aes = list(size=3)))
  
  sub_sce$RECIST <- factor(sub_sce$RECIST,levels = c('PreNMPR','PreMPR','PostNMPR','PostMPR'))
  
  
  
  ggplot(sub_sce@meta.data ,aes(x=umap_1,y=umap_2,color=SubCluster))+
    geom_point(aes(color=SubCluster),size=0.5) +
    # geom_point(aes(color=SubCluster),size = 0.65,shape = 16, stroke = 0) +
    scale_color_manual(breaks = c(levels(sub_sce$SubCluster)),
                       labels= c(  "Epi_C0_SFTPA1",
                                   "Epi_C1_KRT17",
                                   "Epi_C2_TPX2",
                                   "Epi_C3_AKR1C1",
                                   "Epi_C4_IDO1",
                                   "Epi_C5_SCGB1A1",
                                   "Epi_C6_PIP",
                                   "Epi_C7_MMP10",
                                   "Epi_C8_AGER"),
                       values = g.colSet1)+
    ggtitle("Epi Cluster")+
    # geom_text_repel(
    #   aes(label = SubCluster), 
    #   data = class_avg,
    #   color = "black",
    #   size = 3,
    #   box.padding = 0.5,  # 标签与点之间的间距
    #   point.padding = 0.5,  # 标签与其他标签之间的间距
    #   max.overlaps = Inf,  # 允许所有标签显示
    #   segment.color = "grey50",  # 连接线颜色
    #   segment.size = 0.3 ) + # 连接线粗细
    theme_classic()+
    #theme(text=element_text(family="Arial",size=10)) +
    theme(plot.title = element_text(size = 12,hjust = 0.5),
          ##axis.title = element_text(size = 10),
          panel.border = element_rect(fill=NA,color="grey",size=0.5),
          axis.text = element_blank(),
          axis.title = element_blank(),
          axis.line = element_blank(),
          axis.ticks = element_blank(),
          legend.position = "none",
          legend.text = element_text(size = 9),
          strip.text  = element_text(size = 10),
          strip.background = element_blank())+ 
    guides(colour = guide_legend(override.aes = list(size=3)))+ 
    facet_wrap(~RECIST,ncol=2)
}

####Figure 3B###
if(T){
  library(Seurat)
  library(ggplot2)
  library(patchwork)
  library(ggnewscale)
  library(MASS)
  
  # 1. 上皮亚群及Marker
  genes <- c(
    "SFTPA1", "KRT17", "TPX2",
    "AKR1C1", "IDO1", "SCGB1A1",
    "PIP", "MMP10", "AGER"
  )
  
  clusters <- paste0("Epi_C", 0:8, "_", genes)
  
  # 2. 提取UMAP及表达矩阵
  df <- as.data.frame(Embeddings(sub_sce, "umap"))
  colnames(df)[1:2] <- c("x", "y")
  df$celltype <- sub_sce@meta.data[rownames(df), "celltype"]
  df <- cbind(df, FetchData(sub_sce, vars = genes,
                            cells = rownames(df)))
  
  # 3. 计算细胞空间密度
  den <- MASS::kde2d(df$x, df$y, n = 150)
  ix <- pmax(1, pmin(150, findInterval(df$x, den$x)))
  iy <- pmax(1, pmin(150, findInterval(df$y, den$y)))
  df$density <- den$z[cbind(ix, iy)]
  
  # 4. 绘制9个亚群
  plots <- lapply(seq_along(genes), function(i) {
    
    d <- df[df$celltype == clusters[i], ]
    d$expr <- d[[genes[i]]]
    d$expr <- pmin(d$expr / max(quantile(d$expr, 0.99), 1e-6), 1)
    
    ggplot(df, aes(x, y)) +
      geom_point(aes(color = density), size = 0.08) +
      scale_color_gradientn(
        colors = c("#101820", "#722442", "#E94D71"),
        name = "Spot density"
      ) +
      ggnewscale::new_scale_color() +
      geom_point(
        data = d[d$expr > 0, ],
        aes(color = expr),
        size = 0.12
      ) +
      scale_color_gradientn(
        colors = c("#555849", "#A6A52B", "#FFE600", "#FFFFFF"),
        limits = c(0, 1),
        name = "Gene expression"
      ) +
      ggtitle(clusters[i]) +
      theme_void() +
      theme(
        panel.background = element_rect(fill = "#101820", color = NA),
        plot.background = element_rect(fill = "#101820", color = "grey40"),
        plot.title = element_text(
          color = "white", face = "bold.italic",
          hjust = 1, size = 12
        )
      )
  })
  
  # 5. 拼接并保存
  p <- wrap_plots(plots, ncol = 3, guides = "collect") &
    theme(legend.position = "bottom")
  
  p
}


####Figure 3C##
if(T){
  library(Seurat)
  library(ggplot2)
  
  # 1. 上皮细胞亚群
  cell_order <- c(
    "Epi_C0_SFTPA1", "Epi_C1_KRT17", "Epi_C2_TPX2",
    "Epi_C3_AKR1C1", "Epi_C4_IDO1", "Epi_C5_SCGB1A1",
    "Epi_C6_PIP", "Epi_C7_MMP10", "Epi_C8_AGER"
  )
  
  # 2. Marker基因
  markers <- c(
    "SFTPA1", "SFTPB", "NAPSA",
    "KRT17", "KRT5", "TP63",
    "TPX2", "MKI67", "TOP2A",
    "AKR1C1", "AKR1C2", "AKR1C3",
    "IDO1", "IFI27", "ISG15", "CXCL11",
    "SCGB1A1", "SCGB3A1", "CYP2F1",
    "PIP", "LTF", "PRR4",
    "MMP10", "IL33", "MMP1",
    "AGER", "CAV1", "CAV2"
  )
  
  # 3. 绘图
  Epi$celltype <- factor(Epi$celltype, levels = cell_order)
  
  p <- DotPlot(
    Epi,
    features = markers,
    group.by = "celltype",
    dot.scale = 7,
    scale = TRUE
  ) +
    scale_color_gradientn(
      colors = c("#E8EEB2", "#4BAEB5", "#172653"),
      name = "Average"
    ) +
    labs(x = NULL, y = NULL, size = "Percent") +
    theme_bw(base_size = 12) +
    theme(
      panel.grid.major = element_line(color = "#EAEAEA"),
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(
        angle = 50, hjust = 1, vjust = 1, size = 9
      ),
      axis.text.y = element_text(size = 10),
      legend.position = "right"
    )
  
  p
}

###Figure 3D##
if(T){
  Epi$RECIST <- factor(Epi$RECIST,levels = c('PreNMPR','PreMPR','PostNMPR','PostMPR'))
  Epi$SubCluster <- Epi$celltype
  Cellratio <- prop.table(table(Epi$SubCluster, Epi$RECIST), margin = 2) %>% 
    as.data.frame() %>% 
    set_names(c("celltype", "group", "ratio"))
  
  ## 2. 计算各组细胞比例
  library(tidyverse)
  library(ggalluvial)
  library(scales)
  ## 图例从上到下的顺序
  ## 对应堆叠柱中从上到下的顺序
  group_order <-c('PreNMPR','PreMPR','PostNMPR','PostMPR')
  celltype_order <- c(
    "Epi_C0_SFTPA1",
    "Epi_C1_KRT17",
    "Epi_C2_TPX2",
    "Epi_C3_AKR1C1",
    "Epi_C4_IDO1",
    "Epi_C5_SCGB1A1",
    "Epi_C6_PIP",
    "Epi_C7_MMP10",
    "Epi_C8_AGER"
  )
  
  
  ## 检查每组比例之和是否为1
  Cellratio %>%
    group_by(group) %>%
    summarise(
      total_ratio = sum(ratio),
      .groups = "drop"
    )
  colour <- c(
    "Epi_C0_SFTPA1"  = "#00468BFF",
    "Epi_C1_KRT17"   = "#ED0000FF",
    "Epi_C2_TPX2"    = "#42B540FF",
    "Epi_C3_AKR1C1"  = "#0099B4FF",
    "Epi_C4_IDO1"    = "#925E9FFF",
    "Epi_C5_SCGB1A1" = "#FDAF91FF",
    "Epi_C6_PIP"     = "#AD002AFF",
    "Epi_C7_MMP10"   = "#ADB6B6FF",
    "Epi_C8_AGER"    = "#F39B7FFF"
  )
  
  p <- ggplot(
    Cellratio,
    aes(
      x = group,
      y = ratio,
      stratum = celltype,
      alluvium = celltype,
      fill = celltype
    )
  ) +
    
    ## 组间流线
    geom_flow(
      stat = "flow",
      width = 0.60,
      alpha = 0.28,
      knot.pos = 0.35,
      curve_type = "sine",
      color = "white",
      linewidth = 0.55,
      aes.bind = "alluvia",
      lode.guidance = "frontback",
      show.legend = FALSE
    ) +
    
    ## 每组堆叠柱
    geom_stratum(
      width = 0.60,
      color = "white",
      linewidth = 0.55,
      alpha = 1
    ) +
    
    scale_fill_manual(
      values = colour,
      breaks = celltype_order,
      limits = celltype_order,
      drop = FALSE
    ) +
    
    scale_y_continuous(
      limits = c(0, 1),
      breaks = seq(0, 1, by = 0.25),
      labels = number_format(
        accuracy = 0.01
      ),
      expand = expansion(
        mult = c(0, 0)
      )
    ) +
    
    scale_x_discrete(
      limits = group_order,
      expand = expansion(
        mult = c(0.07, 0.07)
      )
    ) +
    
    labs(
      x = NULL,
      y = "Percentage (%)",
      fill = NULL
    ) +
    
    coord_cartesian(
      clip = "off"
    ) +
    
    theme_minimal(
      base_size = 15,
      base_family = "Arial"
    ) +
    
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      
      panel.grid.major.y = element_line(
        color = "#D9D9D9",
        linewidth = 0.55
      ),
      
      axis.line = element_blank(),
      
      axis.ticks = element_blank(),
      
      axis.title.y = element_text(
        size = 16,
        color = "black",
        margin = margin(r = 12)
      ),
      
      axis.text.y = element_text(
        size = 13,
        color = "black"
      ),
      
      axis.text.x = element_text(
        size = 14,
        face = "bold",
        color = "black",
        margin = margin(t = 7)
      ),
      
      legend.position = "right",
      
      legend.text = element_text(
        size = 13,
        color = "black"
      ),
      
      legend.key = element_rect(
        fill = NA,
        color = NA
      ),
      
      legend.key.width = unit(
        0.55,
        "cm"
      ),
      
      legend.key.height = unit(
        0.55,
        "cm"
      ),
      
      legend.spacing.y = unit(
        0.02,
        "cm"
      ),
      
      plot.margin = margin(
        t = 10,
        r = 15,
        b = 10,
        l = 10
      )
    ) +
    
    guides(
      fill = guide_legend(
        title = NULL,
        reverse = FALSE,
        override.aes = list(
          alpha = 1,
          color = NA
        )
      )
    )
  
  p
  
  
}

###Figure 3E###
if(T){
  library(Seurat)
  library(dplyr)
  library(pheatmap)
  group_order <- c(
    "PreNMPR",
    "PreMPR",
    "PostNMPR",
    "PostMPR"
  )
  
  Epi$RECIST <- factor(
    Epi$RECIST,
    levels = group_order
  )
  observed_matrix <- table(
    Epi$RECIST,
    Epi$Cluster
  )
  
  observed_matrix <- as.matrix(
    observed_matrix
  )
  
  # 确保分组顺序正确
  observed_matrix <- observed_matrix[
    group_order,
    ,
    drop = FALSE
  ]
  
  observed_matrix
  # ============================================================
  # 4. 计算随机分布情况下的期望细胞数
  #
  # expected_ij =
  # 第i组细胞总数 × 第j类细胞总数 / 全部细胞总数
  # ============================================================
  
  expected_matrix <- outer(
    rowSums(observed_matrix),
    colSums(observed_matrix)
  ) / sum(observed_matrix)
  
  rownames(expected_matrix) <- rownames(
    observed_matrix
  )
  
  colnames(expected_matrix) <- colnames(
    observed_matrix
  )
  
  
  # ============================================================
  # 5. 计算Observed/Expected ratio
  # ============================================================
  
  roe_matrix <- observed_matrix /
    expected_matrix
  
  # 转置：
  # 行为细胞亚群，列为RECIST分组
  heatdata <- t(
    roe_matrix
  )
  
  heatdata
  # ============================================================
  # 6. 设置颜色范围
  # ============================================================
  
  color_min <- 0
  color_max <- 2.5
  
  # 用于着色的矩阵
  heatdata_color <- heatdata
  
  heatdata_color[
    heatdata_color < color_min
  ] <- color_min
  
  heatdata_color[
    heatdata_color > color_max
  ] <- color_max
  
  # 用于显示真实数值的字符矩阵
  heatdata_number <- matrix(
    sprintf(
      "%.2f",
      heatdata
    ),
    nrow = nrow(heatdata),
    ncol = ncol(heatdata),
    dimnames = dimnames(heatdata)
  )
  
  heatmap_colors <- colorRampPalette(
    c(
      "#0044BB",
      "#92C5DE",
      "#F7F7F7",
      "#F4A582",
      "#C00038"
    )
  )(100)
  
  heatmap_breaks <- seq(
    color_min,
    color_max,
    length.out = 101
  )
  # ============================================================
  # 7. 绘制热图
  # ============================================================
  
  p_ROE <- pheatmap::pheatmap(
    mat = heatdata_color,
    
    color = heatmap_colors,
    breaks = heatmap_breaks,
    
    # 在Pre和Post之间添加间隔
    gaps_col = 2,
    
    # 显示原始真实数值
    display_numbers = heatdata_number,
    number_color = "black",
    fontsize_number = 11,
    
    # 单元格设置
    border_color = NA,
    cellwidth = 68,
    cellheight = 48,
    
    # 标题与字体
    main = "Ro/e",
    fontsize = 15,
    fontsize_row = 14,
    fontsize_col = 13,
    
    # 参考图中列名为竖直方向
    angle_col = "90",
    
    show_rownames = TRUE,
    show_colnames = TRUE,
    
    # 聚类设置
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    
    clustering_distance_rows = "euclidean",
    clustering_method = "complete",
    
    treeheight_row = 60,
    treeheight_col = 0,
    
    # 图例
    legend = TRUE,
    legend_breaks = c(
      0,
      0.5,
      1,
      1.5,
      2,
      2.5
    ),
    
    legend_labels = c(
      "0",
      "0.5",
      "1",
      "1.5",
      "2",
      "2.5"
    )
  )
  
  p_ROE
  
}

###Figure 3I###
if(T){
  library(Seurat)
  library(infercnv)
  library(AnnoProbe)
  library(ggplot2)
  library(ggpubr)
  library(dplyr)
  
  # 1. 读取数据
  combined1 <- qs::qread()
  
  epi <- sort(grep("^Epi_C", unique(combined1$SubCluster), value = TRUE))
  ref <- c("CD4T_Tfh_CXCL13", "NK_FCGR3A")
  
  sub_sce <- subset(
    combined1,
    subset = SubCluster %in% c(epi, ref)
  )
  
  DefaultAssay(sub_sce) <- "RNA"
  sub_sce <- JoinLayers(sub_sce, assay = "RNA")
  
  # 2. 提取原始Counts
  counts <- GetAssayData(
    sub_sce, assay = "RNA", layer = "counts"
  )
  
  # 3. 基因染色体位置
  geneInfor <- annoGene(
    rownames(counts), "SYMBOL", "human"
  )
  
  geneInfor <- geneInfor[, c("SYMBOL", "chr", "start", "end")]
  geneInfor <- geneInfor[
    !duplicated(geneInfor$SYMBOL) &
      geneInfor$SYMBOL %in% rownames(counts), 
  ]
  
  chr_order <- paste0("chr", 1:22)
  geneInfor <- geneInfor[geneInfor$chr %in% chr_order, ]
  
  geneInfor <- geneInfor[
    order(match(geneInfor$chr, chr_order), geneInfor$start),
  ]
  
  counts <- counts[geneInfor$SYMBOL, , drop = FALSE]
  
  rownames(geneInfor) <- geneInfor$SYMBOL
  geneInfor <- geneInfor[, c("chr", "start", "end")]
  
  # 4. 细胞注释
  anno <- data.frame(
    SubCluster = as.character(sub_sce$SubCluster),
    row.names = Cells(sub_sce)
  )
  
  # 5. 创建inferCNV对象
  infercnv_obj <- CreateInfercnvObject(
    raw_counts_matrix = counts,
    annotations_file = anno,
    gene_order_file = geneInfor,
    ref_group_names = ref
  )
  
  # 6. 运行inferCNV
  infercnv_obj <- infercnv::run(
    infercnv_obj,
    cutoff = 0.1,
    out_dir = "output/Epithelial_infercnv",
    cluster_by_groups = TRUE,
    denoise = TRUE,
    HMM = FALSE,
    num_threads = 4,
    write_expr_matrix = TRUE
  )
  
  # 7. 保存结果
  qs::qsave(
    infercnv_obj,
    "output/Epithelial_infercnv/infercnv_obj.qs"
  )
  # 读取inferCNV结果
  infercnv_obj <- qs::qread(
    "output/Epithelial_infercnv/infercnv_obj.qs"
  )
  
  # 计算CNV偏移评分
  expr <- infercnv_obj@expr.data
  
  cnv_score <- colMeans((expr - 1)^2)
  
  # 合并Seurat元数据
  meta <- sub_sce@meta.data[names(cnv_score), , drop = FALSE]
  
  meta$cnv_score <- cnv_score[rownames(meta)]
  meta$SubCluster <- factor(
    meta$SubCluster,
    levels = c(epi, ref)
  )
  # 配色
  cols <- c(
    setNames(scales::hue_pal()(length(epi)), epi),
    "NK_C1_NCAM1" = "#ADB6B6FF",
    "NK_C2_ZFP36" = "#1B1919FF"
  )
  
  p <- ggplot(
    meta,
    aes(x = SubCluster, y = cnv_score, fill = SubCluster)
  ) +
    geom_boxplot(
      width = 0.7,
      outlier.shape = NA,
      linewidth = 0.5
    ) +
    scale_fill_manual(values = g.colSet1) +
    labs(x = NULL, y = "inferCNV score") +
    theme_classic(base_size = 13) +
    theme(
      axis.text.x = element_text(
        angle = 45, hjust = 1, color = "black"
      ),
      axis.text.y = element_text(color = "black"),
      legend.position = "none"
    )
  
  p
}

###Figure 3J###
if(T){
  library(Seurat)
  library(msigdbr)
  library(dplyr)
  library(ComplexHeatmap)
  library(circlize)
  
  # 1. Hallmark基因集
  pathways <- c(
    "E2F targets"    = "HALLMARK_E2F_TARGETS",
    "G2M checkpoint" = "HALLMARK_G2M_CHECKPOINT",
    "DNA repair"     = "HALLMARK_DNA_REPAIR",
    "MYC targets V1" = "HALLMARK_MYC_TARGETS_V1",
    "MYC targets V2" = "HALLMARK_MYC_TARGETS_V2"
  )
  
  hallmark <- msigdbr(species = "Homo sapiens", collection = "H")
  
  genesets <- lapply(pathways, function(x)
    intersect(hallmark$gene_symbol[hallmark$gs_name == x],
              rownames(sub_sce)))
  
  # 2. 计算评分
  DefaultAssay(sub_sce) <- "RNA"
  set.seed(123)
  sub_sce <- AddModuleScore(
    sub_sce, features = genesets, name = "Hallmark_"
  )
  
  # 3. 亚群平均评分
  cell_order <- c(
    "Epi_C0_SFTPA1", "Epi_C1_KRT17", "Epi_C2_TPX2",
    "Epi_C3_AKR1C1", "Epi_C4_IDO1", "Epi_C5_SCGB1A1",
    "Epi_C6_PIP", "Epi_C7_MMP10", "Epi_C8_AGER"
  )
  
  df <- sub_sce@meta.data %>%
    filter(celltype %in% cell_order) %>%
    group_by(celltype) %>%
    summarise(across(paste0("Hallmark_", 1:5), mean), .groups = "drop")
  
  mat <- t(as.matrix(df[match(cell_order, df$celltype), -1]))
  rownames(mat) <- names(pathways)
  colnames(mat) <- cell_order
  
  # 4. Z-score及综合增殖评分
  mat <- t(scale(t(mat)))
  mat <- rbind(mat, "Residual proliferation\ncomposite" = colMeans(mat))
  
  # 5. 热图
  ha <- HeatmapAnnotation(
    celltype = ifelse(colnames(mat) == "Epi_C2_TPX2", "TPX2", "Other"),
    col = list(celltype = c(TPX2 = "#E53935", Other = "#D9D9D9"))
  )
  
  p <- Heatmap(
    mat,
    name = "Z-score",
    col = colorRamp2(c(-3, 0, 3),
                     c("#2166AC", "white", "#D7191C")),
    top_annotation = ha,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    row_names_side = "right",
    column_names_rot = 90,
    rect_gp = gpar(col = "white", lwd = 1.5),
    cell_fun = function(j, i, x, y, width, height, fill)
      grid.text(sprintf("%.2f", mat[i, j]), x, y,
                gp = gpar(fontsize = 9))
  )
  
  draw(p)
  
  
}

###Figure 3K###
if(T){
  library(Seurat)
  library(progeny)
  library(tidyverse)
  library(pheatmap)
  
  # 1. 提取上皮细胞
  seu <- sub_sce
  
  clusters_use <- c(
    "Epi_C0_SFTPA1", "Epi_C1_KRT17", "Epi_C2_TPX2",
    "Epi_C3_AKR1C1", "Epi_C4_IDO1", "Epi_C5_SCGB1A1",
    "Epi_C6_PIP", "Epi_C7_MMP10", "Epi_C8_AGER"
  )
  
  seu <- subset(seu, subset = celltype %in% clusters_use)
  DefaultAssay(seu) <- "RNA"
  
  # 2. PROGENy通路评分
  seu <- progeny::progeny(
    seu,
    scale = FALSE,
    organism = "Human",
    top = 500,
    perm = 1,
    return_assay = TRUE
  )
  
  # 3. 提取通路活性
  scores <- as.matrix(
    GetAssayData(seu, assay = "progeny", layer = "data")
  )
  
  pathways_show <- c(
    "VEGF", "PI3K", "EGFR", "MAPK", "JAK-STAT",
    "TGFb", "Hypoxia", "Estrogen", "WNT",
    "NFkB", "TNFa", "p53", "Trail"
  )
  
  # 4. 计算亚群平均通路活性
  df <- as.data.frame(t(scores)) %>%
    rownames_to_column("cell") %>%
    mutate(celltype = seu$celltype[match(cell, colnames(seu))]) %>%
    group_by(celltype) %>%
    summarise(across(where(is.numeric), mean), .groups = "drop")
  
  mat <- t(as.matrix(df[match(clusters_use, df$celltype), -1]))
  colnames(mat) <- clusters_use
  mat <- mat[intersect(pathways_show, rownames(mat)), , drop = FALSE]
  
  # 5. 行Z-score
  mat <- t(scale(t(mat)))
  mat[!is.finite(mat)] <- 0
  mat <- pmax(pmin(mat, 1.5), -1.5)
  
  # 6. 顶部注释
  anno <- data.frame(
    celltype = ifelse(
      colnames(mat) == "Epi_C2_TPX2",
      "Epi_C2_TPX2", "Other epithelial clusters"
    ),
    row.names = colnames(mat)
  )
  
  anno_colors <- list(
    celltype = c(
      "Epi_C2_TPX2" = "#E53935",
      "Other epithelial clusters" = "#D9D9D9"
    )
  )
  
  # 7. 绘制热图
  p <- pheatmap(
    mat,
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    color = colorRampPalette(
      c("#2166AC", "white", "#B2182B")
    )(100),
    breaks = seq(-1.5, 1.5, length.out = 101),
    annotation_col = anno,
    annotation_colors = anno_colors,
    border_color = "white",
    display_numbers = TRUE,
    number_format = "%.2f",
    fontsize_number = 8,
    angle_col = 90,
    cellwidth = 35,
    cellheight = 25
  )
}

###Figure 3L##
if(T){
  library(scTenifoldKnk)
  library(data.table)
  library(dplyr)
  library(Seurat)
  library(tidyverse)
  library(patchwork)
  library(ggplot2)
  library(ggrepel)
  library(remotes)
  library(qs)
  result <- readRDS('result_TPX2.rds')
  
  # ============================================================
  # 第五部分：结果处理与导出
  # ============================================================
  # 设置显著性阈值
  significant_threshold <- 0.05
  # 提取差异调控基因 ----
  differential_genes <- result$diffRegulation %>%
    # 计算log2倍数变化
    mutate(log_fold_change = log2(FC)) %>%
    # 排除目标基因本身
    filter(gene != target_gene) %>%
    # 筛选显著差异基因
    filter(p.value < significant_threshold)
  # 转换数值列类型（确保数据类型正确）
  numeric_columns <- 2:7
  differential_genes[, numeric_columns] <- sapply(differential_genes[, numeric_columns], as.numeric)
  write.csv(differential_genes,'TPX2/Differantial_genes.csv')
  # 保存结果到Excel文件
  output_filename <- paste0('1_', target_gene, "_虚拟扰动结果.xlsx")
  write.xlsx(differential_genes, output_filename)
  # 预览结果
  head(differential_genes)
  
  # 1. 条形图：显示前20个扰动最显著的基因 ----
  top_n_genes <- 20
  # 选择绝对值最大的前20个基因
  differential_genes_top <- differential_genes %>%
    mutate(abs_logfc = abs(log_fold_change)) %>%
    arrange(desc(abs_logfc)) %>%
    slice_head(n = top_n_genes)
  # 创建条形图
  ggplot(differential_genes_top, 
         aes(x = reorder(gene, log_fold_change), 
             y = log_fold_change)) +  # 把fill移出aes()
    geom_bar(stat = 'identity', 
             alpha = 0.8,
             fill = '#1F77B4FF') +  # 直接在geom_bar中指定颜色
    # 添加数值标签
    geom_text(aes(label = sprintf("%.2f", log_fold_change)),
              hjust = ifelse(differential_genes_top$log_fold_change > 0, -0.2, 1.2),
              size = 3,
              color = "black") +
    coord_flip() +
    labs(title = paste("Top", top_n_genes, "Differentially Regulated Genes"),
         subtitle = paste("KO gene:", target_gene),
         x = "Gene", 
         y = "log2(FC)") +
    theme_classic() +
    theme(
      axis.text = element_text(size = 10, color = "black"),
      axis.title = element_text(size = 12, face = "bold"),
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 11, hjust = 0.5, color = "gray50"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "none"
    ) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))  # 为标签留出空间
  
  
}

###Figure 3M###
if(T){
  library(ggplot2)
  library(dplyr)
  library(patchwork)
  
  dev$ONTOLOGY <- factor(dev$ONTOLOGY, levels = c("BP", "CC", "MF"))
  
  # 绘图函数
  plot_go <- function(type, xvar, col, top = FALSE, show_x = TRUE) {
    
    df <- dev %>%
      filter(ONTOLOGY == type) %>%
      slice_max(.data[[xvar]], n = 5, with_ties = FALSE) %>%
      arrange(desc(.data[[xvar]]))
    
    df$Description <- factor(
      df$Description, levels = rev(df$Description)
    )
    
    ggplot(df, aes(x = .data[[xvar]], y = Description)) +
      geom_col(fill = col, width = 0.75) +
      facet_grid(ONTOLOGY ~ .) +
      scale_x_continuous(
        position = if (top) "top" else "bottom",
        expand = expansion(mult = c(0, 0.08))
      ) +
      labs(
        x = if (top) "log2(FC)" else if (show_x) "Gene Count" else NULL,
        y = NULL
      ) +
      theme_bw(base_size = 12) +
      theme(
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "#D9D9D9"),
        strip.text.y.right = element_text(face = "bold", angle = 90),
        axis.title.x = element_text(face = "bold"),
        axis.text.y = element_text(size = 11),
        axis.text.x = if (show_x || top)
          element_text() else element_blank()
      )
  }
  
  # 三类GO通路
  p1 <- plot_go("BP", "log2FC", "#08549A", top = TRUE)
  p2 <- plot_go("CC", "Count",  "#F00000", show_x = FALSE)
  p3 <- plot_go("MF", "Count",  "#39B54A")
  
  # 拼接
  p <- (p1 / p2 / p3) +
    plot_layout(heights = c(1, 1, 1))
  
  p
  
  ggsave("GO_BP_CC_MF.pdf", p, width = 11, height = 8)
  
}

###Figure 4###
if(T){
  library(tidyverse)
  library(ggplot2)
  library(ggpubr)
  library(gghalves)
  library(survival)
  library(survminer)
  dev_expr <- read.csv('')
  cols <- c(
    "MPR" = "#0073C2FF",
    "NMPR" = "#EFC000FF"
  )
  
  dev_expr$Response <- factor(dev_expr$Response,
                              levels = c("MPR", "NMPR"))
  
  dev_pair$Response <- factor(dev_pair$Response,
                              levels = c("MPR", "NMPR"))
  
  dev$Response <- factor(dev$Response,
                         levels = c("MPR", "NMPR"))

  ###Figure 4A#####
  pA <- ggplot(dev_expr, aes(Response, TPX2, fill = Response)) +
    geom_half_violin(side = "r", width = 0.8) +
    geom_boxplot(width = 0.13, fill = "white",
                 outlier.shape = NA) +
    geom_jitter(aes(color = Response),
                width = 0.08, size = 1) +
    stat_compare_means(
      comparisons = list(c("MPR", "NMPR")),
      method = "wilcox.test",
      label = "p.signif"
    ) +
    scale_fill_manual(values = cols) +
    scale_color_manual(values = cols) +
    labs(x = NULL, y = "TPX2 Expression") +
    theme_classic() +
    theme(legend.position = "none")
  
  pA
  
  ####Figure 4C####
  
  dfC <- dev_pair %>%
    pivot_longer(
      cols = c(Pre, Post),
      names_to = "Time",
      values_to = "Score"
    ) %>%
    mutate(Time = factor(Time, levels = c("Pre", "Post")))
  
  pC <- ggplot(dfC, aes(Time, Score, group = Patient)) +
    geom_line(color = "grey70", linewidth = 0.5) +
    geom_point(aes(color = Response), size = 2.5) +
    facet_wrap(~Response, nrow = 1) +
    stat_compare_means(
      aes(group = Time),
      method = "wilcox.test",
      paired = TRUE,
      label = "p.format"
    ) +
    scale_color_manual(values = cols) +
    labs(x = NULL, y = "TPX2 IHC Score") +
    theme_classic() +
    theme(legend.position = "none")
  
  pC
  ###Figure 4D###
  dev_pair$Delta <- dev_pair$Post - dev_pair$Pre
  
  pD <- ggplot(dev_pair,
               aes(Response, Delta, fill = Response)) +
    geom_boxplot(width = 0.6, outlier.shape = NA) +
    geom_jitter(width = 0.08, size = 1.5) +
    geom_hline(yintercept = 0, linetype = 2,
               color = "grey60") +
    stat_compare_means(
      comparisons = list(c("MPR", "NMPR")),
      method = "wilcox.test",
      label = "p.format"
    ) +
    scale_fill_manual(values = cols) +
    labs(x = NULL,
         y = "TPX2 IHC Score (post-treatment - pretreatment)") +
    theme_classic() +
    theme(legend.position = "none")
  
  pD
  ####Figure 4E-F###
  dev$Delta <- dev$Post - dev$Pre
  
  dev$DeltaGroup <- factor(
    ifelse(dev$Delta > 0, "High", "Low"),
    levels = c("Low", "High")
  )
  
  df_surv <- dev %>% filter(!is.na(DeltaGroup))
  
  km_cols <- c("#0073C2FF", "#EFC000FF")
  km_plot <- function(data, time, event, group, title) {
    
    fit <- survfit(
      as.formula(paste0("Surv(", time, ",", event, ") ~ ", group)),
      data = data
    )
    
    ggsurvplot(
      fit, data = data,
      conf.int = TRUE,
      pval = TRUE,
      risk.table = TRUE,
      risk.table.height = 0.25,
      risk.table.y.text = FALSE,
      palette = km_cols,
      legend.title = title,
      legend.labs = c("Low", "High"),
      xlim = c(0, 72),
      break.time.by = 24,
      xlab = "Time (months)",
      ylab = "Survival Probability",
      ggtheme = theme_classic()
    )
  }
  ###Figure 4E###
  pE <- km_plot(
    df_surv, "OS_time", "OS_event",
    "DeltaGroup", "ΔTPX2 IHC Score"
  )
  
  print(pE)
  
  
  ####4F：Disease-free survival###
  pF <- km_plot(
    df_surv, "DFS_time", "DFS_event",
    "DeltaGroup", "ΔTPX2 IHC Score"
  )
  print(pF)
  
  ####Figure 4G###
  dev$PostGroup <- factor(
    paste0("Post", dev$Response),
    levels = c("PostMPR", "PostNMPR")
  )
  
  post_cols <- c(
    "PostMPR" = "#0073C2FF",
    "PostNMPR" = "#EFC000FF"
  )
  
  pG <- ggplot(dev, aes(PostGroup, Post, fill = PostGroup)) +
    geom_boxplot(width = 0.6, outlier.shape = NA) +
    geom_jitter(width = 0.12, size = 1.4) +
    stat_compare_means(
      comparisons = list(c("PostMPR", "PostNMPR")),
      method = "wilcox.test",
      label = "p.signif"
    ) +
    scale_fill_manual(values = post_cols) +
    labs(x = NULL, y = "TPX2 IHC Score") +
    theme_classic() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
  
  pG
  
  ####Figure I：TPX2 IHC与残余活肿瘤比例相关性###
  cor_res <- cor.test(
    dev$Post, dev$Residual,
    method = "spearman", exact = FALSE
  )
  
  pI <- ggplot(dev, aes(Post, Residual)) +
    geom_smooth(
      method = "loess", se = TRUE,
      color = "grey40", fill = "grey80"
    ) +
    geom_point(aes(color = Response), size = 2) +
    scale_color_manual(values = cols) +
    annotate(
      "text", x = Inf, y = Inf,
      label = sprintf("Spearman rho = %.2f\nP = %.3g",
                      cor_res$estimate, cor_res$p.value),
      hjust = 1.1, vjust = 1.5, size = 3.5
    ) +
    labs(
      x = "Post-treatment surgical specimen TPX2 IHC Score",
      y = "Residual viable tumor (%)",
      color = NULL
    ) +
    theme_classic()
  
  pI
  ###Figure 4J-K
  dev$SurgicalGroup <- factor(
    ifelse(dev$Post > 1, "High", "Low"),
    levels = c("Low", "High")
  )
  pJ <- km_plot(
    dev, "OS_time", "OS_event",
    "SurgicalGroup", "Surgical resection TPX2 IHC Score"
  )
  
  print(pJ)
  pK <- km_plot(
    dev, "DFS_time", "DFS_event",
    "SurgicalGroup", "Surgical resection TPX2 IHC Score"
  )
  
  print(pK)
}

####Figure 5A###
if(T){
  library(qs)
  library(Seurat)
  sub_sce <- qread('Prolif.qs')
  sub_sce$celltype <- factor(sub_sce$celltype,levels=c('Prolif_CD4_FOXP3','Prolif_CD8_GZMH','Prolif_Myeloid_MARCO','Prolif_B_TCL1A'))
  
  sub_sce$Cluster <-sub_sce$celltype
  table(sub_sce$Cluster)
  g.colSet1 <- c("#00468BFF", "#ED0000FF", "#42B540FF", "#0099B4FF", "#925E9FFF",
                 "#FDAF91FF", "#AD002AFF", "#ADB6B6FF", "#1B1919FF", "#F39B7FFF",
                 "#8491B4FF", "#91D1C2FF", "#DC0000FF", "#7E6148FF", "#B09C85FF",
                 "#3C5488FF", "#E64B35FF", "#00A087FF", "#4DBBD5FF", "#6C5483FF")
  sub_sce<-AddMetaData(sub_sce,sub_sce@reductions$umap@cell.embeddings,col.name = colnames(sub_sce@reductions$umap@cell.embeddings))
  class_avg <- sub_sce@meta.data %>%
    group_by(Cluster) %>%
    dplyr::summarise(
      umap_1 = median(umap_1),
      umap_2 = median(umap_2)
    )
  
  sub_sce$Cluster <- sub_sce$celltype
  ggplot(sub_sce@meta.data ,aes(x=umap_1,y=umap_2,color=Cluster))+
    # geom_point(aes(color=SubCluster),size=0.1) +
    geom_point(aes(color=Cluster),size = 0.65,shape = 16, stroke = 0) +
    scale_color_manual(breaks = c(levels(sub_sce$Cluster)),
                       labels= c(  'Prolif_CD4_FOXP3','Prolif_CD8_GZMH','Prolif_Myeloid_MARCO','Prolif_B_TCL1A'                    ),
                       values = g.colSet1)+
    ggtitle("Cluster")+
    geom_text_repel(
      aes(label = Cluster), 
      data = class_avg,
      color = "black",
      size = 3,
      box.padding = 0.5,  # 标签与点之间的间距
      point.padding = 0.5,  # 标签与其他标签之间的间距
      max.overlaps = Inf,  # 允许所有标签显示
      segment.color = "grey50",  # 连接线颜色
      segment.size = 0.3 ) + # 连接线粗细
    theme_classic()+
    #theme(text=element_text(family="Arial",size=10)) +
    theme(plot.title = element_text(size = 12,hjust = 0.5),
          ##axis.title = element_text(size = 10),
          panel.border = element_rect(fill=NA,color="grey",size=0.5),
          axis.text = element_blank(),
          axis.title = element_blank(),
          axis.line = element_blank(),
          axis.ticks = element_blank(),
          legend.position = "right",
          legend.text = element_text(size = 9),
          strip.text  = element_text(size = 10),
          strip.background = element_blank())+ 
    guides(colour = guide_legend(override.aes = list(size=3)))
    sub_sce$SubCluster <- sub_sce$celltype
    ggplot(sub_sce@meta.data ,aes(x=umap_1,y=umap_2,color=SubCluster))+
      geom_point(aes(color=SubCluster),size=0.5) +
      # geom_point(aes(color=SubCluster),size = 0.65,shape = 16, stroke = 0) +
      scale_color_manual(breaks = c(levels(sub_sce$SubCluster)),
                         labels= c('Prolif_CD4_FOXP3','Prolif_CD8_GZMH','Prolif_Myeloid_MARCO','Prolif_B_TCL1A'
                         ),
                         values = g.colSet1)+
      ggtitle("")+
      # geom_text_repel(
      #   aes(label = SubCluster), 
      #   data = class_avg,
      #   color = "black",
      #   size = 3,
      #   box.padding = 0.5,  # 标签与点之间的间距
      #   point.padding = 0.5,  # 标签与其他标签之间的间距
      #   max.overlaps = Inf,  # 允许所有标签显示
      #   segment.color = "grey50",  # 连接线颜色
      #   segment.size = 0.3 ) + # 连接线粗细
      theme_classic()+
      #theme(text=element_text(family="Arial",size=10)) +
      theme(plot.title = element_text(size = 12,hjust = 0.5),
            ##axis.title = element_text(size = 10),
            panel.border = element_rect(fill=NA,color="grey",size=0.5),
            axis.text = element_blank(),
            axis.title = element_blank(),
            axis.line = element_blank(),
            axis.ticks = element_blank(),
            legend.position = "none",
            legend.text = element_text(size = 9),
            strip.text  = element_text(size = 10),
            strip.background = element_blank())+ 
      guides(colour = guide_legend(override.aes = list(size=3)))+ 
      facet_wrap(~RECIST,ncol=2) }

###Figure 5B###
if(T){
  library(Seurat)
  library(Nebulosa)
  library(ggplot2)
  library(patchwork)
  
  seu <- Prolif
  DefaultAssay(seu) <- "RNA"
  
  # 1. 基因组合
  gene_sets <- list(
    "Prolif_CD4_FOXP3"     = c("MKI67", "CD4", "FOXP3"),
    "Prolif_CD8_GZMH"      = c("MKI67", "CD8A", "GZMH"),
    "Prolif_Myeloid_MARCO" = c("MKI67", "CD163", "MARCO"),
    "Prolif_B_TCL1A"       = c("MKI67", "MS4A1", "TCL1A")
  )
  
  # 2. 配色
  cols <- c(
    "#352A46", "#46327E", "#365C8D", "#277F8E",
    "#1FA187", "#4AC16D", "#9FDA3A", "#FDE725"
  )
  
  # 3. 白色虚线椭圆
  um <- Embeddings(seu, "umap")[, 1:2]
  xr <- range(um[, 1]); yr <- range(um[, 2])
  t <- seq(0, 2*pi, length.out = 200)
  a <- -15*pi/180
  
  rx <- 0.075*diff(xr); ry <- 0.105*diff(yr)
  x0 <- xr[1] + 0.17*diff(xr)
  y0 <- yr[1] + 0.35*diff(yr)
  
  ellipse <- data.frame(
    x = x0 + rx*cos(t)*cos(a) - ry*sin(t)*sin(a),
    y = y0 + rx*cos(t)*sin(a) + ry*sin(t)*cos(a)
  )
  
  # 4. 联合密度绘图函数
  plot_joint <- function(genes, title) {
    
    p <- Nebulosa::plot_density(
      seu,
      features = genes,
      reduction = "umap",
      joint = TRUE,
      combine = FALSE,
      method = "ks",
      size = 0.35,
      raster = TRUE
    )
    
    p[[length(p)]] +
      scale_color_gradientn(colors = cols, name = "Joint density") +
      geom_path(
        data = ellipse, aes(x, y),
        inherit.aes = FALSE,
        color = "white", linetype = "dashed",
        linewidth = 0.6
      ) +
      labs(title = title) +
      theme_void() +
      theme(
        plot.title = element_text(
          color = "#FF3030",
          face = "bold.italic", size = 12
        )
      )
  }
  
  # 5. 绘制并拼接
  plots <- Map(plot_joint, gene_sets, names(gene_sets))
  
  p <- wrap_plots(plots, ncol = 2)
  p
}

###Figure 5C###
###if(T)
{
###请进行细胞注释
library(qs)
library(Seurat)
sub_sce <- qread('Prolif.qs')
sub_sce$SubCluster <- sub_sce$celltype
cluster_order <- c(
  "Prolif_CD4_FOXP3",
  "Prolif_CD8_GZMH",
  "Prolif_Myeloid_MARCO",
  "Prolif_B_TCL1A"
)

## 建议用于注释展示的marker（按分群顺序排列）
marker_use <- c(
  ## Prolif_CD4_FOXP3
  "FOXP3", "TNFRSF4", "CCR8", "CD3D", "CD2",
  ## Prolif_CD8_GZMH
  "CD8A", "GZMB", "NKG7", "PRF1", "GZMH",
  ## Prolif_Myeloid_MARCO
  "AIF1", "FCER1G", "C1QA", "C1QB", "C1QC", "MS4A6A", "MARCO", "MSR1", "CD163",
  ## Prolif_B_TCL1A
  "MS4A1", "CD79A", "VPREB3", "TCL1A", "MEF2C", "LRMP",
  ## common proliferative markers
  "MKI67", "TOP2A", "TYMS"
)

## 仅保留对象中实际存在的基因
marker_use <- marker_use[marker_use %in% rownames(scobj2)]

##----------------------------------------------------------
## 3. 生成DotPlot原始数据
## cols只是占位，最终颜色用ggplot重绘
##----------------------------------------------------------
dp <- DotPlot(
  object   = sub_sce,
  features = marker_use,
  group.by = "SubCluster",
  assay    = DefaultAssay(scobj2),
  dot.scale = 8
)$data

##----------------------------------------------------------
## 4. 数据整理
## 为了和示例图更像：
## - y轴从上到下：Epi_C8 -> Epi_C0
## - x轴按marker_use顺序
## - 对avg.exp.scaled做轻度截断，避免极端值把颜色拉爆
##----------------------------------------------------------
dp <- dp %>%
  mutate(
    features.plot = factor(features.plot, levels = marker_use),
    id = factor(id, levels = rev(cluster_order)),
    avg.exp.scaled = pmax(pmin(avg.exp.scaled, 2.5), -1.0),
    pct.exp = pmax(pct.exp, 0)
  )

##----------------------------------------------------------
## 5. 作图
## 更贴近示例图的关键点：
## - 蓝-青绿-浅黄配色
## - 灰白背景+浅灰网格
## - 较大的点
## - 旋转x轴标签
##----------------------------------------------------------
p <- ggplot(dp, aes(x = features.plot, y = id)) +
  geom_point(aes(size = pct.exp, color = avg.exp.scaled)) +
  
  scale_size_continuous(
    name   = "Percent",
    range  = c(0.2, 10),
    breaks = c(0, 25, 50, 75),
    limits = c(0, 75)
  ) +
  
  scale_color_gradient2(
    name     = "Average",
    low      = "#D8DF8A",   # 浅黄绿
    mid      = "#4DB4B5",   # 青绿色
    high     = "#1E2D6E",   # 深蓝
    midpoint = 0,
    breaks   = c(-1, 0, 1, 2)
  ) +
  
  labs(x = NULL, y = NULL) +
  
  theme_bw(base_size = 16) +
  theme(
    panel.background   = element_rect(fill = "white", color = NA),
    plot.background    = element_rect(fill = "white", color = NA),
    
    panel.grid.major.x = element_line(color = "#D0D0D0", linewidth = 0.35),
    panel.grid.major.y = element_line(color = "#D0D0D0", linewidth = 0.35),
    panel.grid.minor   = element_blank(),
    
    axis.text.x = element_text(
      angle = 50,
      hjust = 1,
      vjust = 1,
      size  = 13,
      color = "black"
    ),
    axis.text.y = element_text(
      size  = 15,
      color = "black"
    ),
    axis.ticks = element_blank(),
    
    legend.title = element_text(size = 17, color = "black"),
    legend.text  = element_text(size = 13, color = "black"),
    
    legend.key.height = unit(0.95, "cm"),
    legend.key.width  = unit(0.95, "cm"),
    
    plot.margin = margin(10, 18, 10, 10)
  ) +
  
  guides(
    size = guide_legend(
      order = 1,
      override.aes = list(color = "black")
    ),
    color = guide_colorbar(
      order = 2
    )
  )
}  

####Figure 5D###
if(T){
  sub_sce$RECIST <- factor(sub_sce$RECIST,levels = c('PreNMPR','PreMPR','PostNMPR','PostMPR'))
  
  sub_sce$SubCluster <-  sub_sce$celltype
  Cellratio <- prop.table(table( sub_sce$SubCluster,  sub_sce$RECIST), margin = 2) %>% 
    as.data.frame() %>% 
    set_names(c("celltype", "group", "ratio"))
  
  ## 2. 计算各组细胞比例
  library(tidyverse)
  library(ggalluvial)
  library(scales)
  ## 图例从上到下的顺序
  ## 对应堆叠柱中从上到下的顺序
  celltype_order <- c(
    "Prolif_CD4_FOXP3",
    "Prolif_CD8_GZMH",
    "Prolif_Myeloid_MARCO",
    "Prolif_B_TCL1A"
  )
  
  Cellratio <- Cellratio %>%
    mutate(
      group = factor(
        group,
        levels = group_order
      ),
      celltype = factor(
        celltype,
        levels = celltype_order
      )
    ) %>%
    arrange(
      group,
      celltype
    )
  
  
  ## 检查每组比例之和是否为1
  Cellratio %>%
    group_by(group) %>%
    summarise(
      total_ratio = sum(ratio),
      .groups = "drop"
    )
  colour <- c(
    "Prolif_CD4_FOXP3"   = "#00468B",
    "Prolif_CD8_GZMH"         = "#ED0000",
    "Prolif_Myeloid_MARCO"            = "#42B540",
    "Prolif_B_TCL1A"       = "#0099B4"
  )
  
  p <- ggplot(
    Cellratio,
    aes(
      x = group,
      y = ratio,
      stratum = celltype,
      alluvium = celltype,
      fill = celltype
    )
  ) +
    
    ## 组间流线
    geom_flow(
      stat = "flow",
      width = 0.60,
      alpha = 0.28,
      knot.pos = 0.35,
      curve_type = "sine",
      color = "white",
      linewidth = 0.55,
      aes.bind = "alluvia",
      lode.guidance = "frontback",
      show.legend = FALSE
    ) +
    
    ## 每组堆叠柱
    geom_stratum(
      width = 0.60,
      color = "white",
      linewidth = 0.55,
      alpha = 1
    ) +
    
    scale_fill_manual(
      values = colour,
      breaks = celltype_order,
      limits = celltype_order,
      drop = FALSE
    ) +
    
    scale_y_continuous(
      limits = c(0, 1),
      breaks = seq(0, 1, by = 0.25),
      labels = number_format(
        accuracy = 0.01
      ),
      expand = expansion(
        mult = c(0, 0)
      )
    ) +
    
    scale_x_discrete(
      limits = group_order,
      expand = expansion(
        mult = c(0.07, 0.07)
      )
    ) +
    
    labs(
      x = NULL,
      y = "Percentage (%)",
      fill = NULL
    ) +
    
    coord_cartesian(
      clip = "off"
    ) +
    
    theme_minimal(
      base_size = 15,
      base_family = "Arial"
    ) +
    
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      
      panel.grid.major.y = element_line(
        color = "#D9D9D9",
        linewidth = 0.55
      ),
      
      axis.line = element_blank(),
      
      axis.ticks = element_blank(),
      
      axis.title.y = element_text(
        size = 16,
        color = "black",
        margin = margin(r = 12)
      ),
      
      axis.text.y = element_text(
        size = 13,
        color = "black"
      ),
      
      axis.text.x = element_text(
        size = 14,
        face = "bold",
        color = "black",
        margin = margin(t = 7)
      ),
      
      legend.position = "right",
      
      legend.text = element_text(
        size = 13,
        color = "black"
      ),
      
      legend.key = element_rect(
        fill = NA,
        color = NA
      ),
      
      legend.key.width = unit(
        0.55,
        "cm"
      ),
      
      legend.key.height = unit(
        0.55,
        "cm"
      ),
      
      legend.spacing.y = unit(
        0.02,
        "cm"
      ),
      
      plot.margin = margin(
        t = 10,
        r = 15,
        b = 10,
        l = 10
      )
    ) +
    
    guides(
      fill = guide_legend(
        title = NULL,
        reverse = FALSE,
        override.aes = list(
          alpha = 1,
          color = NA
        )
      )
    )
}

###Figure 5E###
if(T){
  library(Seurat)
  library(dplyr)
  library(pheatmap)
  group_order <- c(
    "PreNMPR",
    "PreMPR",
    "PostNMPR",
    "PostMPR"
  )
  
  Prolif$RECIST <- factor(
    Prolif$RECIST,
    levels = group_order
  )
  observed_matrix <- table(
    Prolif$RECIST,
    Prolif$Cluster
  )
  
  observed_matrix <- as.matrix(
    observed_matrix
  )
  
  # 确保分组顺序正确
  observed_matrix <- observed_matrix[
    group_order,
    ,
    drop = FALSE
  ]
  
  observed_matrix
  # ============================================================
  # 4. 计算随机分布情况下的期望细胞数
  #
  # expected_ij =
  # 第i组细胞总数 × 第j类细胞总数 / 全部细胞总数
  # ============================================================
  
  expected_matrix <- outer(
    rowSums(observed_matrix),
    colSums(observed_matrix)
  ) / sum(observed_matrix)
  
  rownames(expected_matrix) <- rownames(
    observed_matrix
  )
  
  colnames(expected_matrix) <- colnames(
    observed_matrix
  )
  
  
  # ============================================================
  # 5. 计算Observed/Expected ratio
  # ============================================================
  
  roe_matrix <- observed_matrix /
    expected_matrix
  
  # 转置：
  # 行为细胞亚群，列为RECIST分组
  heatdata <- t(
    roe_matrix
  )
  
  heatdata
  # ============================================================
  # 6. 设置颜色范围
  # ============================================================
  
  color_min <- 0
  color_max <- 2.5
  
  # 用于着色的矩阵
  heatdata_color <- heatdata
  
  heatdata_color[
    heatdata_color < color_min
  ] <- color_min
  
  heatdata_color[
    heatdata_color > color_max
  ] <- color_max
  
  # 用于显示真实数值的字符矩阵
  heatdata_number <- matrix(
    sprintf(
      "%.2f",
      heatdata
    ),
    nrow = nrow(heatdata),
    ncol = ncol(heatdata),
    dimnames = dimnames(heatdata)
  )
  
  heatmap_colors <- colorRampPalette(
    c(
      "#0044BB",
      "#92C5DE",
      "#F7F7F7",
      "#F4A582",
      "#C00038"
    )
  )(100)
  
  heatmap_breaks <- seq(
    color_min,
    color_max,
    length.out = 101
  )
  # ============================================================
  # 7. 绘制热图
  # ============================================================
  
  p_ROE <- pheatmap::pheatmap(
    mat = heatdata_color,
    
    color = heatmap_colors,
    breaks = heatmap_breaks,
    
    # 在Pre和Post之间添加间隔
    gaps_col = 2,
    
    # 显示原始真实数值
    display_numbers = heatdata_number,
    number_color = "black",
    fontsize_number = 11,
    
    # 单元格设置
    border_color = NA,
    cellwidth = 68,
    cellheight = 48,
    
    # 标题与字体
    main = "Ro/e",
    fontsize = 15,
    fontsize_row = 14,
    fontsize_col = 13,
    
    # 参考图中列名为竖直方向
    angle_col = "90",
    
    show_rownames = TRUE,
    show_colnames = TRUE,
    
    # 聚类设置
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    
    clustering_distance_rows = "euclidean",
    clustering_method = "complete",
    
    treeheight_row = 60,
    treeheight_col = 0,
    
    # 图例
    legend = TRUE,
    legend_breaks = c(
      0,
      0.5,
      1,
      1.5,
      2,
      2.5
    ),
    
    legend_labels = c(
      "0",
      "0.5",
      "1",
      "1.5",
      "2",
      "2.5"
    )
  )
  
  p_ROE
}

###Figure 5J###
if(T){
  library(Seurat)
  library(progeny)
  library(tidyverse)
  library(pheatmap)
  
  # 1. 读取数据及设置亚群
  Prolif <- qs::qread("Prolif.qs")
  
  clusters_use <- c(
    "Prolif_CD4_FOXP3",
    "Prolif_CD8_GZMH",
    "Prolif_Myeloid_MARCO",
    "Prolif_B_TCL1A"
  )
  
  seu <- subset(Prolif, subset = celltype %in% clusters_use)
  seu$SubCluster <- factor(seu$celltype, levels = clusters_use)
  DefaultAssay(seu) <- "RNA"
  
  # 2. PROGENy通路评分
  seu <- progeny::progeny(
    seu,
    scale = FALSE,
    organism = "Human",
    top = 500,
    perm = 1,
    return_assay = TRUE
  )
  
  seu <- ScaleData(
    seu,
    assay = "progeny",
    features = rownames(seu[["progeny"]]),
    verbose = FALSE
  )
  
  # 3. 提取评分矩阵
  scores <- as.matrix(
    GetAssayData(seu, assay = "progeny", layer = "scale.data")
  )
  
  # 4. 各亚群平均通路活性
  groups <- as.character(seu@meta.data[colnames(scores), "SubCluster"])
  
  mat <- sapply(clusters_use, function(x) {
    rowMeans(scores[, groups == x, drop = FALSE])
  })
  
  # 5. 选择通路
  pathways_show <- c(
    "VEGF", "PI3K", "EGFR", "MAPK",
    "JAK-STAT", "TGFb", "Hypoxia",
    "Estrogen", "WNT", "NFkB",
    "TNFa", "p53", "Trail"
  )
  
  mat <- mat[intersect(pathways_show, rownames(mat)), , drop = FALSE]
  
  # 保存原始平均评分
  write.csv(mat, "Prolif_PROGENy_scores.csv")
  
  # 6. 行Z-score标准化
  mat_z <- t(scale(t(mat)))
  mat_z[!is.finite(mat_z)] <- 0
  mat_z <- pmax(pmin(mat_z, 1.5), -1.5)
  
  # 7. 亚群注释
  cell_colors <- c(
    "Prolif_CD4_FOXP3"      = "#00468B",
    "Prolif_CD8_GZMH"       = "#ED0000",
    "Prolif_Myeloid_MARCO"  = "#42B540",
    "Prolif_B_TCL1A"        = "#0099B4"
  )
  
  annotation_col <- data.frame(
    Celltype = factor(colnames(mat_z), levels = clusters_use),
    row.names = colnames(mat_z)
  )
  
  annotation_colors <- list(Celltype = cell_colors)
  
  # 8. 绘制热图
  p <- pheatmap(
    mat_z,
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    color = colorRampPalette(
      c("#2166AC", "#FFFFFF", "#B2182B")
    )(100),
    breaks = seq(-1.5, 1.5, length.out = 101),
    annotation_col = annotation_col,
    annotation_colors = annotation_colors,
    annotation_names_col = FALSE,
    border_color = "white",
    fontsize = 11,
    fontsize_row = 11,
    fontsize_col = 10,
    angle_col = 45,
    cellwidth = 42,
    cellheight = 25,
    treeheight_row = 40,
    treeheight_col = 0,
    display_numbers = TRUE,
    number_format = "%.2f",
    fontsize_number = 8
  )
}

###Figure 5K###
if(T){
  library(Seurat)
  library(CytoTRACE2)
  library(tidyverse)
  library(ggpubr)
  library(rstatix)
  
  # 1. CytoTRACE2评分
  scobj <- Prolif
  DefaultAssay(scobj) <- "RNA"
  
  scobj <- cytotrace2(
    scobj,
    is_seurat = TRUE,
    slot_type = "counts",
    species = "human",
    seed = 666666
  )
  
  qs::qsave(scobj, "Prolif_cytotrace2.qs")
  
  # 2. UMAP评分图
  p1 <- FeaturePlot(
    scobj,
    features = "CytoTRACE2_Relative",
    reduction = "umap",
    cols = c("#E8EEB2", "#4BAEB5", "#172653"),
    order = TRUE,
    pt.size = 0.1
  ) + theme_classic()
  
  ggsave("CytoTRACE_UMAP.pdf", p1, width = 5, height = 5)
  
  # 3. 设置亚群及颜色
  clusters <- c(
    "Prolif_CD4_FOXP3",
    "Prolif_CD8_GZMH",
    "Prolif_Myeloid_MARCO",
    "Prolif_B_TCL1A"
  )
  
  cols <- setNames(
    c("#00468B", "#ED0000", "#42B540", "#0099B4"),
    clusters
  )
  
  df <- scobj@meta.data %>%
    filter(celltype %in% clusters) %>%
    mutate(celltype = factor(celltype, levels = clusters))
  
  # 4. Wilcoxon检验：以MARCO亚群为参考
  stat <- df %>%
    wilcox_test(
      CytoTRACE2_Relative ~ celltype,
      ref.group = "Prolif_Myeloid_MARCO",
      exact = FALSE
    ) %>%
    adjust_pvalue(method = "BH") %>%
    add_significance(p.col = "p.adj") %>%
    mutate(
      y.position = max(df$CytoTRACE2_Relative, na.rm = TRUE) +
        0.06 * row_number()
    )
  
  # 5. 箱线图
  p2 <- ggboxplot(
    df,
    x = "celltype",
    y = "CytoTRACE2_Relative",
    color = "celltype",
    palette = cols,
    outlier.shape = NA,
    width = 0.65
  ) +
    stat_pvalue_manual(
      stat,
      label = "p.adj.signif",
      hide.ns = TRUE,
      tip.length = 0.005
    ) +
    labs(x = NULL, y = "CytoTRACE2") +
    theme_classic(base_size = 13) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
}

###Figure 5I###
if(T){
  Prolif <- qread('Prolif.qs')
  Idents(Prolif) <- Prolif$celltype
  sce.markers <- FindAllMarkers(Prolif, only.pos = TRUE, min.pct = 0.2, return.thresh = 0.01)
  head(sce.markers)
  
  # 查看top10
  top10 <- sce.markers %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 1) %>%
    slice_head(n = 10) %>%
    ungroup()
  # 基因转成ENTREZID
  Symbol <- mapIds(get("org.Hs.eg.db"), keys = sce.markers$gene, keytype = "SYMBOL", column="ENTREZID")
  head(Symbol)
  ids <- bitr(sce.markers$gene,fromType = "SYMBOL",toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
  head(ids)
  # 合并ENTREZID到obj.markers中
  data <- merge(sce.markers, ids, by.x="gene", by.y="SYMBOL")
  head(data)
  gcSample <- split(data$ENTREZID, data$cluster)
  gcSample
  gcSample <- lapply(
    gcSample,
    function(x) {
      unique(as.character(x))
    }
  )
  
  names(gcSample)
  hallmark_raw <- msigdbr( species = "Homo sapiens", category = "H" )
  head(hallmark_raw)
  ## 构建Hallmark通路—ENTREZID对应关系
  library(dplyr)
  ## 自动识别新旧版本中的Entrez ID列
  entrez_col <- intersect(
    c("ncbi_gene", "entrez_gene"),
    colnames(hallmark_raw)
  )
  if (length(entrez_col) == 0) {
    stop(
      "hallmark_raw中未找到ncbi_gene或entrez_gene列。",
      "请运行colnames(hallmark_raw)检查实际列名。"
    )
  }
  
  entrez_col <- entrez_col[1]
  
  ## 构建Hallmark通路—Entrez ID对应表
  hallmark_t2g <- hallmark_raw %>%
    transmute(
      term = gs_name,
      gene = as.character(.data[[entrez_col]])
    ) %>%
    filter(
      !is.na(gene),
      gene != ""
    ) %>%
    distinct()
  xx_hallmark <- compareCluster(
    gcSample,
    fun = "enricher",
    TERM2GENE = hallmark_t2g,
    pvalueCutoff = 1,
    qvalueCutoff = 1,
    pAdjustMethod = "BH",
    minGSSize = 10,
    maxGSSize = 500
  )
  
  ## 提取富集结果
  res <- as.data.frame(
    xx_hallmark@compareClusterResult
  )
  if (nrow(res) > 0) {
    
    all_entrez <- unique(
      unlist(
        strsplit(
          as.character(res$geneID),
          split = "/",
          fixed = TRUE
        )
      )
    )
    
    id_map <- AnnotationDbi::select(
      org.Hs.eg.db,
      keys = all_entrez,
      keytype = "ENTREZID",
      columns = "SYMBOL"
    ) %>%
      filter(
        !is.na(SYMBOL)
      ) %>%
      distinct(
        ENTREZID,
        .keep_all = TRUE
      )
    
    entrez_to_symbol <- setNames(
      id_map$SYMBOL,
      id_map$ENTREZID
    )
    
    res$geneID <- vapply(
      strsplit(
        as.character(res$geneID),
        split = "/",
        fixed = TRUE
      ),
      function(ids) {
        
        gene_symbols <- unname(
          entrez_to_symbol[ids]
        )
        
        gene_symbols <- unique(
          gene_symbols[
            !is.na(gene_symbols) &
              gene_symbols != ""
          ]
        )
        
        paste(
          gene_symbols,
          collapse = "/"
        )
      },
      FUN.VALUE = character(1)
    )
  }

  enrich <- res %>%
    filter(
      !is.na(pvalue),
      !is.na(p.adjust)
    ) %>%
    group_by(Cluster) %>%
    slice_min(
      order_by = p.adjust,
      n = 5,
      with_ties = FALSE
    ) %>%
    ungroup()
  dt <- enrich %>%
    mutate(
      ## 避免p值极小时-log10出现Inf
      plot_pvalue = pmax(
        pvalue,
        .Machine$double.xmin
      ),
      
      ## 简化Hallmark通路名称
      Description = Description %>%
        str_remove("^HALLMARK_") %>%
        str_replace_all("_", " ") %>%
        str_to_upper()
    ) %>%
    arrange(
      Cluster,
      pvalue
    )
  
  dt$Description <- factor(
    dt$Description,
    levels = rev(
      unique(dt$Description)
    )
  )
  colors <- c(
    "Prolif_CD4_FOXP3" = "#00468B",
    "Prolif_CD8_GZMH" = "#ED0000",
    "Prolif_Myeloid_MARCO" = "#42B540",
    "Prolif_B_TCL1A" = "#0099B4"
  )
  p <- ggplot(dt) +
    geom_link(
      aes(
        x = 0,
        y = Description,
        xend = -log10(plot_pvalue),
        yend = Description,
        alpha = after_stat(index),
        color = Cluster,
        size = after_stat(index)
      ),
      n = 500,
      show.legend = TRUE
    )
  
  p1 <- p +
    geom_point(
      aes(
        x = -log10(plot_pvalue),
        y = Description
      ),
      color = "black",
      fill = "white",
      size = 6,
      shape = 21
    ) +
    theme_classic() +
    theme(
      panel.grid = element_blank(),
      strip.text = element_text(
        face = "bold.italic"
      ),
      axis.line = element_line(
        color = "black",
        linewidth = 0.6
      ),
      axis.text = element_text(
        face = "bold",
        color = "black"
      ),
      axis.title = element_text(
        size = 13
      )
    ) +
    xlab("-Log10 P value") +
    ylab("") +
    scale_color_manual(
      values = colors
    )
  
  p1
  
  ## 按细胞亚群分面
  p2 <- p1 +
    facet_wrap(
      ~Cluster,
      scales = "free",
      ncol = 2
    )
  
  p2
}

###Figure 5L###
if(T){
  options(stringsAsFactors = F)
  library(Seurat)
  library(tidyverse)
  library(data.table)
  library(ggrepel)
  # philentropy包用于计算JS散度，关于这个package的详细介绍在此：http://blog.fens.me/r-entropy/
  library(philentropy) # install.packages("philentropy")
  source("regulon_specificity.R")
  
  #### 计算RSS矩阵(Regulon Specificity Score) ####
  
  seu <- qs::qread("03-2.ifnb_pbmc.seurat.aucell.qs")
  ###提取治疗后的，比较不同的转录因子情况
  scobj <- seu
  table(scobj$celltype)
  DefaultAssay(scobj) <- "AUCell"
  
  # rows are cells and columns are features (regulons)
  rasMat <- t(scobj[["AUCell"]]@data)
  dim(rasMat)
  
  ## 2.2 计算cell type indicate matrix
  ctMat <- calIndMat(scobj$celltype)
  dim(ctMat)
  
  ## 2.3 计算RSS matrix
  rssMat <- calRSSMat(rasMat, ctMat)
  dim(rssMat)
  
  
  ## 2.4 可视化
  library(patchwork)
  ## Overall
  p1 <- PlotRegulonRank(rssMat, "Prolif_CD4_FOXP3",topn = 6)
  p2 <- PlotRegulonRank(rssMat, "Prolif_CD8_GZMH",topn = 6)
  p3 <- PlotRegulonRank(rssMat, "Prolif_Myeloid_MARCO",topn = 6)
  p4 <- PlotRegulonRank(rssMat, "Prolif_B_TCL1A",topn = 6)
}

####Figure S3A####
if(T){
  library(qs)
  library(Seurat)
  sub_sce <- qread('T.qs')
  sub_sce$celltype <- factor(sub_sce$celltype,levels=c("CD4T_Naive_TCF7",
  "CD4T_Tfh_CXCL13",
  "CD4T_Treg_FOXP3",
  "CD8T_Tem_GZMK",
  "CD8T_Tex_LAG3",
  "NK_FCGR3A",
  "T_HSPA1A",
  "T_IFN",
  "γδT"))
  
  sub_sce$Cluster <-sub_sce$celltype
  table(sub_sce$Cluster)
  g.colSet1 <- c("#00468BFF", "#ED0000FF", "#42B540FF", "#0099B4FF", "#925E9FFF",
                 "#FDAF91FF", "#AD002AFF", "#ADB6B6FF", "#1B1919FF", "#F39B7FFF",
                 "#8491B4FF", "#91D1C2FF", "#DC0000FF", "#7E6148FF", "#B09C85FF",
                 "#3C5488FF", "#E64B35FF", "#00A087FF", "#4DBBD5FF", "#6C5483FF")
  sub_sce<-AddMetaData(sub_sce,sub_sce@reductions$umap@cell.embeddings,col.name = colnames(sub_sce@reductions$umap@cell.embeddings))
  class_avg <- sub_sce@meta.data %>%
    group_by(Cluster) %>%
    dplyr::summarise(
      umap_1 = median(umap_1),
      umap_2 = median(umap_2)
    )
  
  sub_sce$Cluster <- sub_sce$celltype
  ggplot(sub_sce@meta.data ,aes(x=umap_1,y=umap_2,color=Cluster))+
    # geom_point(aes(color=SubCluster),size=0.1) +
    geom_point(aes(color=Cluster),size = 0.65,shape = 16, stroke = 0) +
    scale_color_manual(breaks = c(levels(sub_sce$Cluster)),
                       labels= c( "CD4T_Naive_TCF7",
                                  "CD4T_Tfh_CXCL13",
                                  "CD4T_Treg_FOXP3",
                                  "CD8T_Tem_GZMK",
                                  "CD8T_Tex_LAG3",
                                  "NK_FCGR3A",
                                  "T_HSPA1A",
                                  "T_IFN",
                                  "γδT"),
                       values = g.colSet1)+
    ggtitle("TNK Cluster")+
    geom_text_repel(
      aes(label = Cluster), 
      data = class_avg,
      color = "black",
      size = 3,
      box.padding = 0.5,  # 标签与点之间的间距
      point.padding = 0.5,  # 标签与其他标签之间的间距
      max.overlaps = Inf,  # 允许所有标签显示
      segment.color = "grey50",  # 连接线颜色
      segment.size = 0.3 ) + # 连接线粗细
    theme_classic()+
    #theme(text=element_text(family="Arial",size=10)) +
    theme(plot.title = element_text(size = 12,hjust = 0.5),
          ##axis.title = element_text(size = 10),
          panel.border = element_rect(fill=NA,color="grey",size=0.5),
          axis.text = element_blank(),
          axis.title = element_blank(),
          axis.line = element_blank(),
          axis.ticks = element_blank(),
          legend.position = "right",
          legend.text = element_text(size = 9),
          strip.text  = element_text(size = 10),
          strip.background = element_blank())+ 
    guides(colour = guide_legend(override.aes = list(size=3)))
  sub_sce$SubCluster <- sub_sce$celltype
  ggplot(sub_sce@meta.data ,aes(x=umap_1,y=umap_2,color=SubCluster))+
    geom_point(aes(color=SubCluster),size=0.5) +
    # geom_point(aes(color=SubCluster),size = 0.65,shape = 16, stroke = 0) +
    scale_color_manual(breaks = c(levels(sub_sce$SubCluster)),
                       labels= c('Prolif_CD4_FOXP3','Prolif_CD8_GZMH','Prolif_Myeloid_MARCO','Prolif_B_TCL1A'
                       ),
                       values = g.colSet1)+
    ggtitle("")+
    # geom_text_repel(
    #   aes(label = SubCluster), 
    #   data = class_avg,
    #   color = "black",
    #   size = 3,
    #   box.padding = 0.5,  # 标签与点之间的间距
    #   point.padding = 0.5,  # 标签与其他标签之间的间距
    #   max.overlaps = Inf,  # 允许所有标签显示
    #   segment.color = "grey50",  # 连接线颜色
    #   segment.size = 0.3 ) + # 连接线粗细
    theme_classic()+
    #theme(text=element_text(family="Arial",size=10)) +
    theme(plot.title = element_text(size = 12,hjust = 0.5),
          ##axis.title = element_text(size = 10),
          panel.border = element_rect(fill=NA,color="grey",size=0.5),
          axis.text = element_blank(),
          axis.title = element_blank(),
          axis.line = element_blank(),
          axis.ticks = element_blank(),
          legend.position = "none",
          legend.text = element_text(size = 9),
          strip.text  = element_text(size = 10),
          strip.background = element_blank())+ 
    guides(colour = guide_legend(override.aes = list(size=3)))+ 
    facet_wrap(~RECIST,ncol=2) }

###Figure S3B
if(T){
  library(Seurat)
  library(ggplot2)
  
  # 1. Seurat对象
  seu <- sub_sce
  DefaultAssay(seu) <- "RNA"
  
  # 2. 细胞亚群（从下到上）
  cell_order <- c(
    "CD4T_Naive_TCF7",
    "CD4T_Tfh_CXCL13",
    "CD4T_Treg_FOXP3",
    "CD8T_Tem_GZMK",
    "CD8T_Tex_LAG3",
    "NK_FCGR3A",
    "T_HSPA1A",
    "T_IFN",
    "γδT"
  )
  
  # 3. Marker基因
  markers <- c(
    "IL7R", "TCF7", "CCR7", "LEF1", "SELL",
    "CXCL13", "PDCD1", "ICOS", "MAF", "TOX2",
    "FOXP3", "IL2RA", "CTLA4", "TIGIT", "IKZF2",
    "GZMK", "GZMA", "CCL5",
    "LAG3", "HAVCR2", "TOX", "ENTPD1",
    "FCGR3A", "KLRD1", "GNLY", "NKG7", "FGFBP2",
    "HSPA1A", "HSPA1B", "DNAJB1", "JUN", "FOS",
    "ISG15", "MX1", "IFIT1", "IFIT3", "RSAD2",
    "TRDC", "TRGC1", "TRGC2", "GZMH", "PRF1"
  )
  
  # 4. 绘图
  seu$celltype <- factor(seu$celltype, levels = cell_order)
  
  p <- DotPlot(
    seu,
    features = markers,
    group.by = "celltype",
    dot.scale = 8,
    scale = TRUE,
    col.min = -1,
    col.max = 2
  ) +
    scale_color_gradientn(
      colors = c("#E8F2B2", "#7BC7B5", "#3195BF", "#172653"),
      values = scales::rescale(c(-1, 0, 1, 2)),
      limits = c(-1, 2),
      name = "Average"
    ) +
    scale_size_area(
      max_size = 8,
      breaks = c(0, 25, 50, 75),
      name = "Percent"
    ) +
    labs(x = NULL, y = NULL) +
    theme_bw(base_size = 12) +
    theme(
      panel.grid.major = element_line(color = "#E5E5E5"),
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(
        angle = 45, hjust = 1, color = "black", size = 10
      ),
      axis.text.y = element_text(color = "black", size = 12),
      legend.position = "right"
    )
  
  p
}

###Figure S4A###
if(T){
  Epi <- qread('Epi.qs')
  Idents(Epi) <- Epi$celltype
  sce.markers <- FindAllMarkers(Epi, only.pos = TRUE, min.pct = 0.2, return.thresh = 0.01)
  head(sce.markers)
  
  # 查看top10
  top10 <- sce.markers %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 1) %>%
    slice_head(n = 10) %>%
    ungroup()
  
  
  # 基因转成ENTREZID
  Symbol <- mapIds(get("org.Hs.eg.db"), keys = sce.markers$gene, keytype = "SYMBOL", column="ENTREZID")
  head(Symbol)
  ids <- bitr(sce.markers$gene,fromType = "SYMBOL",toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
  head(ids)
  
  # 合并ENTREZID到obj.markers中
  data <- merge(sce.markers, ids, by.x="gene", by.y="SYMBOL")
  head(data)
  
  gcSample <- split(data$ENTREZID, data$cluster)
  gcSample
  xx_go <- compareCluster(gcSample, fun="enrichGO", OrgDb="org.Hs.eg.db", ont="BP", pvalueCutoff=1, qvalueCutoff=1)
  res <- xx_go@compareClusterResult
  res <- xx_go@compareClusterResult
  head(res)
  
  ## 将富集结果中的 ENTREZID 重新转为 SYMBOL
  for (i in 1:dim(res)[1]) {
    arr = unlist(strsplit(as.character(res[i,"geneID"]), split="/"))
    gene_names = paste(unique(names(Symbol[Symbol %in% arr])), collapse="/")
    res[i,"geneID"] = gene_names
  }
  head(res)
  
  ## 通路筛选Top5
  enrich <- res %>% 
    group_by(Cluster) %>% 
    top_n(n = 5, wt = -pvalue) 
  
  dt <- enrich
  dt <- dt[order(dt$pvalue,dt$Cluster, decreasing = F), ]
  dt$Description <- factor(dt$Description, levels = unique(dt$Description))
  colnames(dt)
  library(ggforce)
  head(dt)
  table(dt$Cluster)
  colors <- c(
    "Epi_C0_SFTPA1"  = "#00468BFF",
    "Epi_C1_KRT17"   = "#ED0000FF",
    "Epi_C2_TPX2"    = "#42B540FF",
    "Epi_C3_AKR1C1"  = "#0099B4FF",
    "Epi_C4_IDO1"    = "#925E9FFF",
    "Epi_C5_SCGB1A1" = "#FDAF91FF",
    "Epi_C6_PIP"     = "#AD002AFF",
    "Epi_C7_MMP10"   = "#ADB6B6FF",
    "Epi_C8_AGER"    = "#F39B7FFF"
  )
  
  p <- ggplot(dt) +
    geom_link(aes(x = 0, y = Description,
                  xend = -log10(pvalue), yend = Description,
                  alpha = after_stat(index),
                  color = Cluster,
                  size = after_stat(index)), 
              n = 500, show.legend = T) 
  p
  
  
  
  p1 <- p +
    geom_point(aes(x = -log10(pvalue),y = Description), color = "black", fill = "white",size = 6,shape = 21) +
    geom_text(aes(x = -log10(pvalue), y = Description), label=dt$Count, size=3, nudge_x=0.05) + 
    theme_classic() +
    theme(panel.grid = element_blank(),
          strip.text = element_text(face = "bold.italic"),
          #axis.text = element_text(color = "black"),
          axis.line = element_line(color = "black", size = 0.6), # 加粗x轴和y轴的线条
          axis.text = element_text(face = "bold"), # 加粗x轴和y轴的标签
          axis.title = element_text( size = 13)    # 加粗x轴和y轴的标题
    ) +
    xlab("-Log10 Pvalue") + ylab("") + 
    scale_color_manual(values = colors)
  p1
  p2 <- p1 + 
    facet_wrap(~Cluster,scales = "free",ncol = 2) 
  p2
}

###Figure S4B###
if(T){
  library(Seurat)
  library(CytoTRACE2)
  library(ggplot2)
  library(ggrepel)
  library(dplyr)
  
  # 1. CytoTRACE2评分
  scobj <- Epi
  DefaultAssay(scobj) <- "RNA"
  
  scobj <- cytotrace2(
    scobj,
    is_seurat = TRUE,
    slot_type = "counts",
    species = "human",
    seed = 666666
  )
  
  qs::qsave(scobj, "Epi_cytotrace2.qs")
  
  # 2. 提取细胞亚群标签坐标
  umap <- as.data.frame(Embeddings(scobj, "umap"))
  colnames(umap)[1:2] <- c("UMAP_1", "UMAP_2")
  umap$celltype <- scobj$celltype[rownames(umap)]
  
  labels <- umap %>%
    group_by(celltype) %>%
    summarise(
      UMAP_1 = median(UMAP_1),
      UMAP_2 = median(UMAP_2),
      .groups = "drop"
    )
  
  # 3. 设置配色
  cols <- c(
    "#081D58", "#253494", "#225EA8",
    "#1D91C0", "#41B6C4", "#41AB5D",
    "#A1D76A", "#FFE54C", "#F46D43",
    "#D73027", "#7F0000"
  )
  
  # 4. UMAP绘图
  p <- FeaturePlot(
    scobj,
    features = "CytoTRACE2_Relative",
    reduction = "umap",
    order = TRUE,
    pt.size = 0.5,
    raster = TRUE,
    combine = FALSE
  )[[1]] +
    scale_color_gradientn(
      colors = cols,
      limits = c(0, 1),
      name = "CytoTRACE2\nrelative score"
    ) +
    geom_label_repel(
      data = labels,
      aes(UMAP_1, UMAP_2, label = celltype),
      inherit.aes = FALSE,
      fill = scales::alpha("white", 0.88),
      color = "black",
      fontface = "bold",
      size = 3,
      max.overlaps = Inf
    ) +
    labs(
      title = "CytoTRACE2 relative score",
      x = "UMAP 1", y = "UMAP 2"
    ) +
    coord_fixed() +
    theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      axis.text = element_blank(),
      axis.ticks = element_blank()
    )
}

###Figure S4C+D###
if(T){
  ###Python代码
  library(Seurat)
  library(Matrix)
  
  DefaultAssay(Epi) <- "RNA"
  Epi <- JoinLayers(Epi, assay = "RNA")
  
  dir.create("Epi_sctour", showWarnings = FALSE)
  
  mat <- GetAssayData(Epi, layer = "counts")
  
  writeMM(mat, "Epi_sctour/counts.mtx")
  write.table(rownames(mat), "Epi_sctour/genes.tsv",
              row.names = FALSE, col.names = FALSE, quote = FALSE)
  write.table(colnames(mat), "Epi_sctour/cells.tsv",
              row.names = FALSE, col.names = FALSE, quote = FALSE)
  
  um <- Embeddings(Epi, "umap")[colnames(mat), 1:2]
  
  write.csv(
    data.frame(
      cell = colnames(mat),
      celltype = Epi$celltype[colnames(mat)],
      UMAP_1 = um[, 1],
      UMAP_2 = um[, 2]
    ),
    "Epi_sctour/meta.csv",
    row.names = FALSE
  )
  ###scTour完整分析（Python）###
  import scanpy as sc
  import sctour as sct
  import pandas as pd
  import numpy as np
  import matplotlib.pyplot as plt
  from scipy.io import mmread
  
  # 1. 读取Seurat导出的数据
  path = "Epi_sctour/"
  
  genes = pd.read_csv(path + "genes.tsv", header=None)[0].astype(str)
  cells = pd.read_csv(path + "cells.tsv", header=None)[0].astype(str)
  meta = pd.read_csv(path + "meta.csv").set_index("cell")
  
  adata = sc.AnnData(
    X=mmread(path + "counts.mtx").T.tocsr().astype(np.float32),
    obs=meta.loc[cells].copy(),
    var=pd.DataFrame(index=genes)
  )
  adata.var_names_make_unique()
  adata.obsm["X_umap"] = adata.obs[["UMAP_1", "UMAP_2"]].to_numpy()
  
  # 2. 数据预处理
  sc.pp.calculate_qc_metrics(
    adata, percent_top=None, log1p=False, inplace=True
  )
  sc.pp.filter_genes(adata, min_cells=20)
  sc.pp.highly_variable_genes(
    adata, flavor="seurat_v3",
    n_top_genes=2000, subset=True
  )
  
  # 3. scTour模型训练
  tnode = sct.train.Trainer(
    adata,
    loss_mode="nb",
    random_state=123
  )
  tnode.train()
  
  # 4. 拟时序和向量场
  adata.obs["ptime"] = tnode.get_time()
  
  adata.obsm["X_TNODE"] = tnode.get_latentsp(
    alpha_z=0.5,
    alpha_predz=0.5
  )[0]
  
  adata.obsm["X_VF"] = tnode.get_vector_field(
    adata.obs["ptime"].values,
    adata.obsm["X_TNODE"]
  )
  
  # 5. 设置细胞亚群配色
  cell_colors = {
    "Epi_C0_SFTPA1": "#00468B",
    "Epi_C1_KRT17": "#ED0000",
    "Epi_C2_TPX2": "#42B540",
    "Epi_C3_AKR1C1": "#0099B4",
    "Epi_C4_IDO1": "#925E9F",
    "Epi_C5_SCGB1A1": "#FDAF91",
    "Epi_C6_PIP": "#AD002A",
    "Epi_C7_MMP10": "#ADB6B6",
    "Epi_C8_AGER": "#1B1919"
  }
  
  adata.obs["celltype"] = pd.Categorical(
    adata.obs["celltype"],
    categories=list(cell_colors),
    ordered=True
  )
  adata.uns["celltype_colors"] = list(cell_colors.values())
  
  # 6. Figure C + D
  fig, axes = plt.subplots(1, 2, figsize=(14, 5))
  
  # C: Pseudotime
  sc.pl.umap(
    adata,
    color="ptime",
    cmap="viridis",
    size=8,
    frameon=False,
    title="C  Pseudotime",
    ax=axes[0],
    show=False
  )
  
  # D: Cell-state vector field
  sct.vf.plot_vector_field(
    adata,
    zs_key="X_TNODE",
    vf_key="X_VF",
    use_rep_neigh="X_TNODE",
    t_key="ptime",
    color="celltype",
    size=12,
    alpha=0.5,
    stream=True,
    stream_color="black",
    stream_density=1.5,
    stream_linewidth=0.7,
    legend_loc="right margin",
    frameon=False,
    ax=axes[1],
    show=False
  )
  
  axes[1].set_title("D  Cell-state vector field")
  
  plt.tight_layout()
  plt.savefig("Epi_scTour_C_D.pdf", bbox_inches="tight")
  plt.show()
  
  # 7. 保存结果
  adata.write_h5ad("Epi_scTour_result.h5ad")
  adata.obs[["celltype", "ptime"]].to_csv("Epi_scTour_ptime.csv")
  tnode.save_model(".", "Epi_scTour_model")
}

###Figure S4E+F###
if(T){
  library(Seurat)
  library(monocle3)
  library(ggplot2)
  library(patchwork)
  library(qs)
  
  # 1. Seurat转换为Monocle3
  seu <- Epi
  seu$SubCluster <- seu$celltype
  DefaultAssay(seu) <- "RNA"
  seu <- JoinLayers(seu, assay = "RNA")
  
  counts <- GetAssayData(seu, layer = "counts")
  
  cds <- new_cell_data_set(
    counts,
    cell_metadata = seu@meta.data[colnames(counts), ],
    gene_metadata = data.frame(
      gene_short_name = rownames(counts),
      row.names = rownames(counts)
    )
  )
  
  cds <- preprocess_cds(cds, num_dim = 30)
  cds <- reduce_dimension(cds)
  
  # 使用Seurat原始UMAP
  SingleCellExperiment::reducedDim(cds, "UMAP") <-
    Embeddings(seu, "umap")[colnames(cds), ]
  
  # 2. 构建拟时序轨迹
  cds <- cluster_cells(cds, resolution = 0.01, k = 20)
  
  cds <- learn_graph(
    cds,
    use_partition = FALSE,
    learn_graph_control = list(
      minimal_branch_len = 50,
      geodesic_distance_ratio = 0.4
    )
  )
  
  # 3. 设置起始细胞
  xy <- SingleCellExperiment::reducedDim(cds, "UMAP")
  root_cells <- colnames(cds)[cds$SubCluster == "Epi_C0_SFTPA1"]
  root_xy <- xy[root_cells, , drop = FALSE]
  
  root <- root_cells[which.min(
    rowSums(sweep(root_xy, 2, apply(root_xy, 2, median))^2)
  )]
  
  cds <- order_cells(cds, root_cells = root)
  
  qs::qsave(cds, "Epi_monocle3.qs")
  
  # 4. 设置配色及分组顺序
  cell_colors <- c(
    "Epi_C0_SFTPA1"  = "#00468BFF",
    "Epi_C1_KRT17"   = "#ED0000FF",
    "Epi_C2_TPX2"    = "#42B540FF",
    "Epi_C3_AKR1C1"  = "#0099B4FF",
    "Epi_C4_IDO1"    = "#925E9FFF",
    "Epi_C5_SCGB1A1" = "#FDAF91FF",
    "Epi_C6_PIP"     = "#AD002AFF",
    "Epi_C7_MMP10"   = "#ADB6B6FF",
    "Epi_C8_AGER"    = "#E64B35FF"
  )
  
  cds$SubCluster <- factor(cds$SubCluster, levels = names(cell_colors))
  cds$RECIST <- factor(
    cds$RECIST,
    levels = c("PostMPR", "PostNMPR", "PreMPR", "PreNMPR")
  )
  
  # 5. 通用绘图主题
  theme_trajectory <- theme_classic() +
    theme(
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.line = element_blank(),
      panel.border = element_rect(
        fill = NA, color = "grey80"
      ),
      plot.title = element_text(hjust = 0.5)
    )
  
  # 6. Figure E：拟时序轨迹
  pE <- plot_cells(
    cds,
    color_cells_by = "pseudotime",
    show_trajectory_graph = TRUE,
    label_cell_groups = FALSE,
    label_branch_points = TRUE,
    label_leaves = TRUE,
    cell_size = 0.5,
    trajectory_graph_segment_size = 0.8
  ) +
    scale_color_gradientn(
      colors = c(
        "#E9F0BA", "#BDDDAE", "#57B4B6",
        "#2682B1", "#21579B", "#1A264C"
      ),
      name = "Pseudotime"
    ) +
    labs(title = "Cell trajectory of epithelial cells") +
    theme_trajectory
  
  # 7. Figure F：全部亚群轨迹
  pF1 <- plot_cells(
    cds,
    color_cells_by = "SubCluster",
    show_trajectory_graph = TRUE,
    label_cell_groups = FALSE,
    label_branch_points = FALSE,
    label_leaves = FALSE,
    cell_size = 0.4
  ) +
    scale_color_manual(values = cell_colors) +
    labs(title = "Cell trajectory by subcluster") +
    theme_trajectory +
    theme(legend.position = "none")
  
  # 8. Figure F：RECIST分组轨迹
  pF2 <- plot_cells(
    cds,
    color_cells_by = "SubCluster",
    show_trajectory_graph = TRUE,
    label_cell_groups = FALSE,
    label_branch_points = FALSE,
    label_leaves = FALSE,
    cell_size = 0.3
  ) +
    scale_color_manual(values = cell_colors) +
    facet_wrap(~RECIST, ncol = 2) +
    theme_trajectory
  
  # 9. 拼接并保存
  pF <- pF1 | pF2
  p <- pE | pF
}

###Figure S4G###
if(T){
  library(Seurat)
  library(tidyverse)
  library(philentropy)
  library(patchwork)
  library(ggrepel)
  
  source("regulon_specificity.R")
  
  # 1. Seurat对象
  scobj <- Epi
  DefaultAssay(scobj) <- "AUCell"
  
  # 2. 提取AUCell活性矩阵
  rasMat <- t(as.matrix(
    GetAssayData(scobj, assay = "AUCell", layer = "data")
  ))
  
  # 3. 计算RSS
  ctMat <- calIndMat(scobj$celltype[rownames(rasMat)])
  rssMat <- calRSSMat(rasMat, ctMat)
  
  # 4. 设置亚群
  clusters_use <- c(
    "Epi_C0_SFTPA1",
    "Epi_C1_KRT17",
    "Epi_C2_TPX2",
    "Epi_C3_AKR1C1",
    "Epi_C4_IDO1",
    "Epi_C5_SCGB1A1",
    "Epi_C6_PIP",
    "Epi_C7_MMP10",
    "Epi_C8_AGER"
  )
  
  # 5. 绘制每个亚群Top 6 Regulons
  plots <- lapply(clusters_use, function(x) {
    PlotRegulonRank(rssMat, x, topn = 6) +
      ggtitle(x) +
      theme_classic(base_size = 11) +
      theme(
        plot.title = element_text(hjust = 0.5, size = 12),
        axis.title = element_text(size = 10)
      )
  })
  
  # 6. 拼接3×3图
  p <- wrap_plots(plots, ncol = 3)
  p
}

###Figure S4H#
if(T){
  library(Seurat)
  library(tidyverse)
  library(patchwork)
  
  # 1. 提取AUCell评分
  seu <- Epi
  DefaultAssay(seu) <- "AUCell"
  
  TFs <- c("MXD3(+)", "E2F7(+)", "E2F1(+)")
  
  mat <- GetAssayData(seu, assay = "AUCell", layer = "data")
  
  df <- as.data.frame(t(as.matrix(mat[TFs, , drop = FALSE])))
  df$celltype <- seu@meta.data[rownames(df), "celltype"]
  
  # 2. 亚群顺序及配色
  cell_colors <- c(
    "Epi_C0_SFTPA1"  = "#00468BFF",
    "Epi_C1_KRT17"   = "#ED0000FF",
    "Epi_C2_TPX2"    = "#42B540FF",
    "Epi_C3_AKR1C1"  = "#0099B4FF",
    "Epi_C4_IDO1"    = "#925E9FFF",
    "Epi_C5_SCGB1A1" = "#FDAF91FF",
    "Epi_C6_PIP"     = "#AD002AFF",
    "Epi_C7_MMP10"   = "#ADB6B6FF",
    "Epi_C8_AGER"    = "#1B1919FF"
  )
  
  df$celltype <- factor(df$celltype, levels = names(cell_colors))
  
  # 3. 绘图函数
  plot_tf <- function(gene) {
    ggplot(df, aes(x = celltype, y = .data[[gene]], fill = celltype)) +
      geom_violin(
        trim = TRUE,
        width = 0.85,
        linewidth = 0.35,
        color = "grey35"
      ) +
      scale_fill_manual(values = cell_colors) +
      labs(title = gene, x = NULL, y = "TF activity score") +
      theme_classic(base_size = 12) +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none"
      )
  }
  
  # 4. 纵向拼接
  plots <- lapply(TFs, plot_tf)
  
  p <- wrap_plots(plots, ncol = 1)
  p
}

###igure S4I###
if(T){
  library(Seurat)
  library(monocle3)
  library(ggplot2)
  library(patchwork)
  
  # 1. 提取Monocle3拟时序
  Epi$Pseudotime <- pseudotime(cds)[Cells(Epi)]
  
  # 2. Marker基因
  genes <- c("SFTPA1", "KRT17", "TPX2", "AKR1C1")
  
  # 3. 分组及颜色
  Epi$Type <- factor(
    Epi$RECIST,
    levels = c("PreNMPR", "PreMPR", "PostNMPR", "PostMPR")
  )
  
  cols <- c(
    "PreNMPR"  = "#00468BFF",
    "PreMPR"   = "#ED0000FF",
    "PostNMPR" = "#42B540FF",
    "PostMPR"  = "#0099B4FF"
  )
  
  # 4. 提取表达数据
  plotdata <- FetchData(
    Epi,
    vars = c(genes, "Pseudotime", "Type")
  )
  
  plotdata <- plotdata[is.finite(plotdata$Pseudotime), ]
  
  # 5. 绘图函数
  plot_gene <- function(gene) {
    ggplot(plotdata, aes(
      x = Pseudotime,
      y = .data[[gene]],
      color = Type
    )) +
      geom_smooth(
        method = "loess",
        se = TRUE,
        linewidth = 1
      ) +
      scale_color_manual(values = cols) +
      labs(
        title = gene,
        x = "Pseudotime",
        y = "Expression",
        color = NULL
      ) +
      theme_gray(base_size = 12) +
      theme(
        panel.background = element_rect(fill = "gray95"),
        panel.grid = element_line(color = "white"),
        panel.border = element_rect(fill = NA, color = "black"),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.position = "bottom"
      )
  }
  
  # 6. 拼接及保存
  plots <- lapply(genes, plot_gene)
  
  p <- wrap_plots(plots, ncol = 2, guides = "collect") &
    theme(legend.position = "bottom")
}

###Figure S4J###
if(T){
  library(Seurat)
  library(tidyverse)
  library(ggpubr)
  library(rstatix)
  library(patchwork)
  
  # 1. 数据及分组
  seu <- Epi
  DefaultAssay(seu) <- "AUCell"
  
  TFs <- c("MXD3(+)", "E2F7(+)", "E2F1(+)")
  
  clusters <- c(
    "Epi_C0_SFTPA1", "Epi_C1_KRT17", "Epi_C2_TPX2",
    "Epi_C3_AKR1C1", "Epi_C4_IDO1", "Epi_C5_SCGB1A1",
    "Epi_C6_PIP", "Epi_C7_MMP10", "Epi_C8_AGER"
  )
  
  cols <- c(
    "PreNMPR"  = "#00468BFF",
    "PreMPR"   = "#ED0000FF",
    "PostNMPR" = "#42B540FF",
    "PostMPR"  = "#0099B4FF"
  )
  
  # 2. 提取AUCell评分
  mat <- GetAssayData(seu, assay = "AUCell", layer = "data")
  
  df <- as.data.frame(t(as.matrix(mat[TFs, , drop = FALSE]))) %>%
    rownames_to_column("Cell") %>%
    pivot_longer(-Cell, names_to = "TF", values_to = "Score") %>%
    mutate(
      celltype = seu$celltype[match(Cell, Cells(seu))],
      RECIST = seu$RECIST[match(Cell, Cells(seu))],
      celltype = factor(celltype, levels = clusters),
      RECIST = factor(RECIST, levels = names(cols)),
      TF = factor(TF, levels = TFs)
    ) %>%
    filter(!is.na(celltype), !is.na(RECIST))
  
  # 3. Wilcoxon检验（同一亚群内比较）
  stat <- df %>%
    group_by(TF, celltype) %>%
    wilcox_test(
      Score ~ RECIST,
      comparisons = list(
        c("PreNMPR", "PreMPR"),
        c("PostNMPR", "PostMPR")
      )
    ) %>%
    ungroup() %>%
    mutate(
      p.adj = p.adjust(p, method = "BH"),
      label = case_when(
        p.adj < 0.0001 ~ "****",
        p.adj < 0.001 ~ "***",
        p.adj < 0.01 ~ "**",
        p.adj < 0.05 ~ "*",
        TRUE ~ ""
      )
    )
  
  # 4. 显著性标注位置
  pos <- df %>%
    group_by(TF, celltype) %>%
    summarise(y = max(Score), .groups = "drop")
  
  stat <- left_join(stat, pos, by = c("TF", "celltype")) %>%
    mutate(
      x = as.numeric(celltype) +
        ifelse(group1 == "PreNMPR", -0.2, 0.2),
      y = y * 1.10
    )
  
  # 5. 绘制Figure J
  p <- ggplot(df, aes(celltype, Score, fill = RECIST)) +
    geom_boxplot(
      position = position_dodge(0.8),
      width = 0.7,
      outlier.shape = NA,
      linewidth = 0.3
    ) +
    geom_text(
      data = filter(stat, label != ""),
      aes(x = x, y = y, label = label),
      inherit.aes = FALSE,
      size = 3
    ) +
    facet_wrap(~TF, ncol = 1, scales = "free_y") +
    scale_fill_manual(values = cols) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.2))) +
    labs(x = NULL, y = "TF Scores", fill = "RECIST") +
    theme_classic(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      strip.background = element_blank(),
      strip.text = element_text(size = 13),
      legend.position = "right"
    )
  
  p
}

###Figure 5###
if(T){
  library(tidyverse)
  library(ggpubr)
  library(survival)
  library(survminer)
  library(lmerTest)
  library(emmeans)
  library(patchwork)
  
  # 分组及配色
  dev$Response <- factor(dev$Response, levels = c("MPR", "NMPR"))
  dev$PreGroup <- factor(dev$PreGroup, levels = c("Low", "High"))
  
  cols <- c("MPR" = "#0073C2FF", "NMPR" = "#EFC000FF")
  km_cols <- c("Low" = "#0073C2FF", "High" = "#EFC000FF")
  
  # A：治疗前TPX2 IHC评分
  pA <- ggplot(dev, aes(Response, Pre, fill = Response)) +
    geom_boxplot(width = 0.55, outlier.shape = NA) +
    geom_jitter(width = 0.08, size = 1.5) +
    stat_compare_means(
      comparisons = list(c("MPR", "NMPR")),
      method = "wilcox.test", label = "p.signif"
    ) +
    scale_fill_manual(values = cols) +
    labs(x = NULL, y = "TPX2 IHC Score") +
    theme_classic() +
    theme(legend.position = "none")
  
  # B：TPX2与残余活肿瘤比例
  pB <- ggplot(dev, aes(Pre, Residual)) +
    geom_point(shape = 21, fill = "white", size = 2) +
    geom_smooth(method = "loess", se = TRUE,
                color = "#00468B", fill = "#DCEBF5") +
    stat_cor(method = "spearman", label.x.npc = "right",
             label.y.npc = "top", size = 3.5) +
    labs(x = "Pretreatment TPX2 IHC score",
         y = "Residual viable tumor (%)") +
    theme_classic()
  
  # C-D：生存分析
  km_plot <- function(time, event, ylab) {
    
    fit <- survfit(
      as.formula(paste0("Surv(", time, ",", event, ") ~ PreGroup")),
      data = dev
    )
    
    ggsurvplot(
      fit, data = dev,
      palette = km_cols,
      conf.int = TRUE,
      pval = TRUE,
      risk.table = TRUE,
      risk.table.height = 0.25,
      risk.table.y.text = FALSE,
      legend.title = "Treatment-naive TPX2 IHC Score",
      legend.labs = c("Low", "High"),
      xlim = c(0, 72),
      break.time.by = 24,
      xlab = "Time (months)",
      ylab = ylab,
      ggtheme = theme_bw()
    )
  }
  
  pC <- km_plot("OS_time", "OS_event",
                "Overall Survival Probability")
  
  pD <- km_plot("DFS_time", "DFS_event",
                "Disease Free Survival Probability")
  
  # E：治疗前后TPX2纵向变化（混合效应模型）
  dfE <- dev %>%
    select(Patient, Response, Pre, Post) %>%
    pivot_longer(
      cols = c(Pre, Post),
      names_to = "Time",
      values_to = "Score"
    ) %>%
    mutate(Time = factor(Time, levels = c("Pre", "Post")))
  
  fitE <- lmer(
    Score ~ Time * Response + (1 | Patient),
    data = dfE
  )
  
  # 时间×疗效组交互作用
  p_int <- anova(fitE)["Time:Response", "Pr(>F)"]
  
  # 估计边际均值及95%CI
  emm <- as.data.frame(emmeans(fitE, ~ Time * Response))
  
  pE <- ggplot(
    emm,
    aes(Time, emmean, group = Response, color = Response)
  ) +
    geom_line(linewidth = 1.1) +
    geom_point(size = 3) +
    geom_errorbar(
      aes(ymin = lower.CL, ymax = upper.CL),
      width = 0.08, linewidth = 0.7
    ) +
    annotate(
      "text", x = 1.5, y = Inf,
      label = paste0(
        "Time × response: P = ",
        format.pval(p_int, digits = 2, eps = 0.001)
      ),
      vjust = 1.5, size = 4
    ) +
    scale_color_manual(values = cols) +
    scale_x_discrete(labels = c(
      "Pretreatment\nbiopsy",
      "Post-treatment\nresection"
    )) +
    labs(x = NULL, y = "Estimated TPX2 score", color = NULL) +
    theme_classic(base_size = 12) +
    theme(legend.position = "top")
  
  # 显示图像
  pA
  pB
  print(pC)
  print(pD)
  pE
  
  
  
}

####Figure 7###
if(T){
  library(Seurat)
  library(CellChat)
  library(tidyverse)
  library(patchwork)
  
  groups <- c("PreMPR", "PreNMPR")
  
  types <- c(
    "Epi_C2_TPX2",
    "Prolif_Myeloid_MARCO",
    "Prolif_CD8_GZMH",
    "Prolif_CD4_FOXP3",
    "Prolif_B_TCL1A"
  )
  
  seu <- subset(
    combined1,
    subset = RECIST %in% groups & SubCluster %in% types
  )
  
  DefaultAssay(seu) <- "RNA"
  seu <- JoinLayers(seu, assay = "RNA")
  seu <- NormalizeData(seu, verbose = FALSE)
  
  run_cc <- function(g) {
    obj <- subset(seu, subset = RECIST == g)
    
    cc <- createCellChat(obj, group.by = "SubCluster", assay = "RNA")
    cc@DB <- CellChatDB.human
    cc <- subsetData(cc)
    cc <- identifyOverExpressedGenes(cc)
    cc <- identifyOverExpressedInteractions(cc)
    cc <- computeCommunProb(cc, type = "truncatedMean", trim = 0.1)
    cc <- filterCommunication(cc, min.cells = 10)
    cc <- computeCommunProbPathway(cc)
    cc <- aggregateNet(cc)
    cc
  }
  
  cc <- setNames(lapply(groups, run_cc), groups)
  qs::qsave(cc, "CellChat_PreMPR_PreNMPR.qs")
  ###Figure 7A###
  epi <- "Epi_C2_TPX2"
  marco <- "Prolif_Myeloid_MARCO"
  targets <- setdiff(types, c(epi, marco))
  
  pairs <- c(
    paste(marco, c("Prolif_CD8_GZMH",
                   "Prolif_CD4_FOXP3",
                   "Prolif_B_TCL1A", epi), sep = " -> "),
    paste(epi, c(marco, targets), sep = " -> ")
  )
  
  pathways <- c(
    "ANNEXIN", "CCL", "CD45", "CD46", "CD99",
    "CDH1", "CPA3", "FN1", "GALECTIN",
    "ICAM", "JAM", "MHC-I", "MIF", "MK",
    "PECAM1", "Prostaglandin", "SPP1", "TNF"
  )
  
  netA <- purrr::imap_dfr(cc, function(obj, g) {
    x <- reshape2::melt(obj@netP$prob, value.name = "prob")
    colnames(x) <- c("source", "target", "pathway", "prob")
    x %>% mutate(group = g, pair = paste(source, target, sep = " -> "))
  })
  
  dfA <- expand_grid(
    group = groups, pair = pairs, pathway = pathways
  ) %>%
    left_join(netA, by = c("group", "pair", "pathway")) %>%
    mutate(prob = replace_na(prob, 0)) %>%
    select(group, pair, pathway, prob) %>%
    pivot_wider(names_from = group, values_from = prob) %>%
    mutate(
      diff = PreNMPR - PreMPR,
      pair = factor(pair, levels = rev(pairs)),
      pathway = factor(pathway, levels = pathways)
    )
  
  pA <- ggplot(dfA, aes(pathway, pair, fill = diff)) +
    geom_tile(color = "white", linewidth = 0.5) +
    scale_fill_gradient2(
      low = "#2166AC", mid = "white", high = "#D66060",
      midpoint = 0, name = "PreNMPR - PreMPR"
    ) +
    labs(x = NULL, y = NULL) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid = element_blank()
    )
  
  pA
  ##Figure 7B-C
  lr_names <- c(
    "MIF_CD74_CXCR4"  = "MIF-(CD74+CXCR4)",
    "MIF_CD74_CD44"   = "MIF-(CD74+CD44)",
    "SPP1_CD44"       = "SPP1-CD44",
    "SPP1_ITGA4_ITGB1" = "SPP1-(ITGA4+ITGB1)"
  )
  
  netLR <- purrr::imap_dfr(cc, function(obj, g) {
    subsetCommunication(obj) %>% mutate(group = g)
  }) %>%
    mutate(
      LR = paste(ligand, receptor, sep = "_"),
      LR = factor(LR, levels = names(lr_names),
                  labels = unname(lr_names)),
      significance = ifelse(pval <= 0.01, "P<=0.01", "0.01<P<0.05")
    ) %>%
    filter(!is.na(LR))
  
  plot_lr <- function(sender, receivers, title) {
    
    d <- netLR %>%
      filter(source == sender, target %in% receivers) %>%
      mutate(
        pair = factor(
          paste(target, group, sep = "\n"),
          levels = unlist(lapply(groups, function(g)
            paste(receivers, g, sep = "\n")))
        )
      )
    
    ggplot(d, aes(pair, LR)) +
      geom_point(aes(size = significance, color = prob)) +
      scale_size_manual(
        values = c("0.01<P<0.05" = 2.5, "P<=0.01" = 4.5),
        name = "Significance"
      ) +
      scale_color_gradientn(
        colors = c("#DDE68D", "#60BBC0", "#243D80"),
        name = "Communication\nprobability"
      ) +
      labs(title = title, x = NULL, y = NULL) +
      theme_bw(base_size = 10) +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold"),
        axis.text.x = element_text(angle = 90, hjust = 1),
        panel.grid = element_line(color = "grey92")
      )
  }
  
  pB <- plot_lr(
    epi, setdiff(types, epi),
    "Epi_C2_TPX2-derived signaling"
  )
  
  pC <- plot_lr(
    marco, setdiff(types, marco),
    "Prolif_Myeloid_MARCO-derived signaling"
  )
  
  pB
  pC
  ####igure D–E：配受体表达气泡图###
  lr_ref <- tibble(
    LR = unname(lr_names),
    ligand = c("MIF", "MIF", "SPP1", "SPP1"),
    receptor = c("CD74_CXCR4", "CD74_CD44",
                 "CD44", "ITGA4_ITGB1")
  )
  
  calc_lr <- function(obj, sender, receiver, ligand, receptor) {
    
    x <- obj@data
    ids <- as.character(obj@idents)
    s <- colnames(x)[ids == sender]
    t <- colnames(x)[ids == receiver]
    r <- strsplit(receptor, "_", fixed = TRUE)[[1]]
    
    if (!all(c(ligand, r) %in% rownames(x)) ||
        !length(s) || !length(t)) {
      return(tibble(score = NA_real_, pct = NA_real_))
    }
    
    L <- as.numeric(x[ligand, s])
    R <- as.matrix(x[r, t, drop = FALSE])
    
    receptor_mean <- exp(mean(log(rowMeans(R) + 1e-8)))
    
    tibble(
      score = sqrt(mean(L) * receptor_mean),
      pct = 100 * sqrt(
        mean(L > 0) * mean(colSums(R > 0) == length(r))
      )
    )
  }
  
  plot_expression <- function(sender, receivers, lr_use, title) {
    
    d <- expand_grid(
      group = groups, target = receivers, LR = lr_use
    ) %>%
      left_join(lr_ref, by = "LR") %>%
      rowwise() %>%
      mutate(result = list(calc_lr(
        cc[[group]], sender, target, ligand, receptor
      ))) %>%
      unnest(result) %>%
      ungroup() %>%
      mutate(
        pair = factor(
          paste(sender, "->", target),
          levels = paste(sender, "->", receivers)
        ),
        LR = factor(LR, levels = unname(lr_names))
      )
    
    ggplot(d, aes(pair, LR)) +
      geom_point(aes(size = pct, color = score)) +
      facet_grid(~group) +
      scale_color_gradientn(
        colors = c("#DDE68D", "#60BBC0", "#172653"),
        name = "LR expression\nscore"
      ) +
      scale_size_area(
        max_size = 7,
        limits = c(0, 100),
        breaks = c(0, 25, 50, 75, 100),
        name = "Combined\npositive %"
      ) +
      labs(title = title, x = NULL, y = NULL) +
      theme_bw(base_size = 10) +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold"),
        axis.text.x = element_text(angle = 50, hjust = 1),
        panel.grid = element_blank(),
        strip.background = element_rect(fill = "grey93")
      )
  }
  
  pD <- plot_expression(
    epi, setdiff(types, epi), unname(lr_names),
    "Epi_C2_TPX2-derived signaling before treatment"
  )
  
  pE <- plot_expression(
    marco, setdiff(types, marco), unname(lr_names)[3:4],
    "Prolif_Myeloid_MARCO-derived signaling before treatment"
  )
  
  pD
  pE
}