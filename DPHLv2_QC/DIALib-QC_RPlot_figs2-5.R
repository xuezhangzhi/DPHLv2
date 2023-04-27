# Reference: DIALib-QC_RPlot.pl

# --------- reset environment --------------------
rm(list = ls())
libs <- c('RColorBrewer', 'ggplot2', 'scales', 'ggpubr')
#sapply(libs, require, character.only=TRUE)
table(sapply(libs, require, character.only=TRUE))
rm(libs)
ModifiedPeptideSequence
# ------- function: plotting ------------------------
plotting <- function(df, pdfName, rlt_dir = './', allColor = c('#A2565B', '#6E8E84', '#1A476F', '#E37E00', '#90353A'),
                     allFill = c('#C79A9D', '#B7C7C2', '#8DA3B7', '#F1BE80', '#C89A9D'))
{
  # DIA-LibQC (For FragPipe-Easypqp library format)
  df_mz <- df %>% dplyr::distinct(ModifiedPeptideSequence, PrecursorCharge, PrecursorMz)
  # Plot A: precursor m/z distributaion
  p1 <- ggplot(df_mz, aes(x = PrecursorMz))+
    geom_histogram(aes(y = (..count..)/sum(..count..)),
                   color = allColor[1],
                   fill = allFill[1])+
    #geom_text(aes(y = ((..count..)/sum(..count..)), label = scales::percent((..count..)/sum(..count..))), stat = 'count', vjust = -10)+
    #xlim(350, 1300)+
    scale_y_continuous(labels = scales::percent)+
    labs(title = 'Precusor m/z', x = "Precusor m/z", y = "Frequency in percent (%)")+
    theme_bw()+
    theme(legend.text = element_text(size = 18, color = "black"),legend.position = 'top',
          legend.title = element_text(size = 18, color="black") ,
          panel.grid.major =element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"))+
    theme(panel.grid =element_blank())+
    theme(axis.text = element_text(size = 15,color = "black"))+
    theme(axis.text.x = element_text(size = 15,color = "black", angle = 0))+
    theme(plot.title = element_text(size = 20, face = "bold"))+
    theme(plot.subtitle=element_text(size = 20, hjust = 0, color="black"))+
    theme(axis.title.x=element_text(size = 18, hjust = 0.5, color="black"))+
    theme(axis.title.y=element_text(size = 18, hjust = 0.5, color="black"))
  # Plot B: precursor charge distributaion
  p2 <- ggplot(df_mz, aes(x = PrecursorCharge))+
    geom_bar(aes(y = (..count..)/sum(..count..)),
             width = 0.6,
             color = allColor[2],
             fill = allFill[2],
    )+
    geom_text(aes(y = ((..count..)/sum(..count..)), label = scales::percent((..count..)/sum(..count..))), stat = 'count', size = 6, hjust = 0.5, vjust = -0.1, position = "stack", color = 'black')+
    scale_y_continuous(labels = scales::percent)+
    labs(title = 'Precursor charge', x = "Precursor charge", y = "Precursor in percent (%)")+
    theme_bw()+
    theme(legend.text = element_text(size = 18, color = "black"),legend.position = 'top',
          legend.title = element_text(size = 18, color="black") ,
          panel.grid.major =element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"))+
    theme(panel.grid =element_blank())+
    theme(axis.text = element_text(size = 15,color = "black"))+
    theme(axis.text.x = element_text(size = 15,color = "black", angle = 0))+
    theme(plot.title = element_text(size = 20, face = "bold"))+
    theme(plot.subtitle=element_text(size = 20, hjust = 0, color="black"))+
    theme(axis.title.x=element_text(size = 18, hjust = 0.5, color="black"))+
    theme(axis.title.y=element_text(size = 18, hjust = 0.5, color="black"))
  
  
  df_rt <- df %>% dplyr::distinct(ModifiedPeptideSequence, PrecursorCharge, NormalizedRetentionTime)
  df_rt <- df_rt[df_rt$PrecursorCharge %in% c(2, 3), ]
  tb_tmp <- table(df_rt$ModifiedPeptideSequence)
  pairs <- names(tb_tmp[tb_tmp == 2])
  df_rt <- df_rt[df_rt$ModifiedPeptideSequence %in% pairs, ]
  df_rt <- reshape2::dcast(df_rt, ModifiedPeptideSequence~PrecursorCharge, value.var = 'NormalizedRetentionTime')
  colnames(df_rt) <- c('modSeq', 'RT2value', 'RT3value')
  #df_rt %<>% dplyr::filter(RT2value <= 250) %>% dplyr::filter(RT3value <= 250)
  # Plot C: +2/+3 RT linear regression
  correlation <- sprintf("%1f", cor(df_rt$RT2value, df_rt$RT3value, method = "pearson"))
  N <- nrow(df_rt)
  lbl = paste0('n=', N, ', R2=', correlation)
  p3 <- ggplot(df_rt, aes(x = RT2value, y = RT3value))+
    geom_point(shape = 1, size = 1, color = allFill[3])+
    geom_smooth( method = "lm", formula = 'y~x', color = 'black', alpha = 0, size = .3 )+
    labs(title = 'Library +2/+3 pair iRT correlation', x = "+2 iRT", y = "+3 iRT")+
    annotate("text", x = min(df_rt$RT2value), y = max(df_rt$RT3value) * 1.05, label = lbl, hjust = 0, vjust = 0, size = 5)+
    theme_bw()+
    theme(legend.text = element_text(size = 18, color = "black"),legend.position = 'top',
          legend.title = element_text(size = 18, color="black") ,
          panel.grid.major =element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"))+
    theme(panel.grid =element_blank())+
    theme(axis.text = element_text(size = 15,color = "black"))+
    theme(axis.text.x = element_text(size = 15,color = "black", angle = 0))+
    theme(plot.title = element_text(size = 20, face = "bold"))+
    theme(plot.subtitle=element_text(size = 20, hjust = 0, color="black"))+
    theme(axis.title.x=element_text(size = 18, hjust = 0.5, color="black"))+
    theme(axis.title.y=element_text(size = 18, hjust = 0.5, color="black"))
  
  
  df_len <- df %>% dplyr::distinct(PeptideSequence)
  df_len$PeptideLength <- nchar(df_len$PeptideSequence)
  # Plot D: Peptide length distribution
  p4 <- ggplot(df_len, aes(x = PeptideLength))+
    geom_bar(aes(y = (..count..)/sum(..count..)),
             width = 1,
             color = allColor[4],
             fill = allFill[4])+
    scale_y_continuous(labels = scales::percent)+
    labs(title = 'Peptide length', x = 'Peptide length', y = 'Peptide in percent (%)')+
    theme_bw()+
    theme(legend.text = element_text(size = 18, color = "black"),legend.position = 'top',
          legend.title = element_text(size = 18, color="black") ,
          panel.grid.major =element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"))+
    theme(panel.grid =element_blank())+
    theme(axis.text = element_text(size = 15,color = "black"))+
    theme(axis.text.x = element_text(size = 15,color = "black", angle = 0))+
    theme(plot.title = element_text(size = 20, face = "bold"))+
    theme(plot.subtitle=element_text(size = 20, hjust = 0, color="black"))+
    theme(axis.title.x=element_text(size = 18, hjust = 0.5, color="black"))+
    theme(axis.title.y=element_text(size = 18, hjust = 0.5, color="black"))
  
  
  df_mod <- df %>% dplyr::distinct(ModifiedPeptideSequence)
  # df_mod['[+42]'] <- unlist(lapply(df_mod$ModifiedPeptideSequence, function(e){ return(grepl('(UniMod:1)', e)) }))
  df_mod['[+57]'] <- unlist(lapply(df_mod$ModifiedPeptideSequence, function(e){ return(grepl('(UniMod:4)', e)) }))
  df_mod['[+16]'] <- unlist(lapply(df_mod$ModifiedPeptideSequence, function(e){ return(grepl('(UniMod:35)', e)) }))
  vec_mod <- unlist(lapply(df_mod[, -1], sum))
  df_mod <- data.frame(mod_type = names(vec_mod), freq = vec_mod)
  # Plot E: Modification counts
  p5 <- ggplot(df_mod, aes(mod_type, freq))+
    geom_col(color = allColor[5],
             fill = allFill[5],
             width = 0.5)+
    geom_text(aes(label = freq),
              size = 6, hjust = 0.5, vjust = 0, position = "stack", color = 'black')+
    labs(title = 'Modifications', x = 'Modification type', y = 'Number of modifications')+
    theme_bw()+
    theme(legend.text = element_text(size = 18, color = "black"),legend.position = 'top',
          legend.title = element_text(size = 18, color="black") ,
          panel.grid.major =element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(color = "black"))+
    theme(panel.grid =element_blank())+
    theme(axis.text = element_text(size = 15,color = "black"))+
    theme(axis.text.x = element_text(size = 15,color = "black", angle = 0))+
    theme(plot.title = element_text(size = 20, face = "bold"))+
    theme(plot.subtitle=element_text(size = 20, hjust = 0, color="black"))+
    theme(axis.title.x=element_text(size = 18, hjust = 0.5, color="black"))+
    theme(axis.title.y=element_text(size = 18, hjust = 0.5, color="black"))
  
  
  df_peppro <- df %>% dplyr::distinct(ProteinId, ModifiedPeptideSequence)
  df_peppro <- data.frame(table(table(df_peppro$ProteinId))); colnames(df_peppro) <- c('pep_num', 'count')
  if (sum(as.numeric(as.character(df_peppro$pep_num)) >= 8))
  {
    df_peppro$pep_num <- as.numeric(as.character(df_peppro$pep_num))
    over8 <- sum(df_peppro[df_peppro$pep_num >= 8, 'count'])
    df_peppro <- df_peppro[df_peppro$pep_num < 8, ]
    df_peppro <- rbind(df_peppro, c('>8', over8))
    df_peppro$pep_num <- factor(df_peppro$pep_num, levels = c('1', '2', '3', '4', '5', '6', '7', '>8'), ordered = T)
    df_peppro$count <- as.numeric(df_peppro$count)
    #df_peppro[1, 'sum'] <- df_peppro$count[1]
    #for (i in 2:nrow(df_peppro)){
    #  df_peppro[i, 'sum'] <- df_peppro[i-1, 'sum'] + df_peppro$count[i]
    #}
    #df_peppro <- df_peppro[which(df_peppro$sum < (0.95 * sum(df_peppro$count))), ]
  }
  # Plot F: Peptides per protein
  p6 <- ggplot(df_peppro, aes(pep_num, count))+
    geom_col(color = allColor[6], fill = allFill[6])+
    geom_text(aes(label = count), size = 6, hjust = 0.5, vjust = 0, position = "stack", color = 'black')+
    labs(title = 'Peptides per protein', x = 'Peptides per protein', y = 'Number of protiens')+
    scale_x_discrete("Peptides per protein",
                     #limits = 
    )+
    theme_bw()+
    theme(legend.text = element_text(size = 18, color = "black"),legend.position = 'top',
          legend.title = element_text(size = 18, color="black") ,
          panel.grid.major =element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"))+
    theme(panel.grid =element_blank())+
    theme(axis.text = element_text(size = 15,color = "black"))+
    theme(axis.text.x = element_text(size = 15,color = "black", angle = 0))+
    theme(plot.title = element_text(size = 20, face = "bold"))+
    theme(plot.subtitle=element_text(size = 20, hjust = 0, color="black"))+
    theme(axis.title.x=element_text(size = 18, hjust = 0.5, color="black"))+
    theme(axis.title.y=element_text(size = 18, hjust = 0.5, color="black"))
  
  
  df$Precursor <- paste0(df$ModifiedPeptideSequence, '_', df$'PrecursorCharge')
  df_frpr <- df %>% dplyr::distinct(Precursor, Annotation)
  df_frpr <- data.frame(table(table(df_frpr$Precursor))); colnames(df_frpr) <- c('fr_num', 'count')
  if(sum(as.numeric(as.character(df_frpr$fr_num)) >= 6))
  {
    df_frpr$fr_num <- as.numeric(as.character(df_frpr$fr_num))
    over6 <- sum(df_frpr[df_frpr$fr_num >= 6, 'count'])
    df_frpr <- df_frpr[df_frpr$fr_num < 6, ]
    df_frpr <- rbind(df_frpr, c('>6', over6))
    df_frpr$fr_num <- factor(df_frpr$fr_num, levels = c('1', '2', '3', '4', '5', '>6'), ordered = T)
    df_frpr$count <- as.numeric(df_frpr$count)
  }
  df_frpr$freq <- df_frpr$count / sum(df_frpr$count)
  # Plot G: Fragments per precursor
  p7 <- ggplot(df_frpr, aes(fr_num, freq))+
    geom_col(width = 0.5, color = allColor[7], fill = allFill[7])+
    geom_text(aes(label = paste0(sprintf('%.3f', freq * 100), '%')), size = 6, hjust = 0.5, vjust = 0, position = "stack", color = 'black')+
    labs(title = 'Fragments per precursor', x = 'Fragments per precursor', y = 'Frequency in percent (%)')+
    scale_x_discrete("Fragments per precursor",
                     #limits = 
    )+
    scale_y_continuous(labels = scales::percent)+
    theme_bw()+
    theme(legend.text = element_text(size = 18, color = "black"),legend.position = 'top',
          legend.title = element_text(size = 18, color="black") ,
          panel.grid.major =element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"))+
    theme(panel.grid =element_blank())+
    theme(axis.text = element_text(size = 15,color = "black"))+
    theme(axis.text.x = element_text(size = 15,color = "black", angle = 0))+
    theme(plot.title = element_text(size = 20, face = "bold"))+
    theme(plot.subtitle=element_text(size = 20, hjust = 0, color="black"))+
    theme(axis.title.x=element_text(size = 18, hjust = 0.5, color="black"))+
    theme(axis.title.y=element_text(size = 18, hjust = 0.5, color="black"))
  
  
  df_frtype <- df %>% dplyr::distinct(Precursor, Annotation, FragmentType)
  df_frtype <- data.frame(table(df_frtype$FragmentType)); colnames(df_frtype) <- c('fr_type', 'count')
  df_frtype$freq <- df_frtype$count / sum(df_frtype$count)
  # Plot H: Fragment ion type distribution
  p8 <- ggplot(df_frtype, aes(fr_type, freq))+ theme_bw()+
    geom_col(color = allColor[8], fill = allFill[8], width = 0.3)+
    geom_text(aes(label = paste0(sprintf('%.2f', freq * 100), '%')),
              size = 6, hjust = 0.5, vjust = 0, position = "stack", color = 'black')+
    labs(title = 'Fragment ion', x = 'Fragment ion type', y = 'Frequency in percent (%)')+
    scale_y_continuous(labels = scales::percent)+
    theme_bw()+
    theme(legend.text = element_text(size = 18, color = "black"),legend.position = 'top',
          legend.title = element_text(size = 18, color="black") ,
          panel.grid.major =element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"))+
    theme(panel.grid =element_blank())+
    theme(axis.text = element_text(size = 15,color = "black"))+
    theme(axis.text.x = element_text(size = 15,color = "black", angle = 0))+
    theme(plot.title = element_text(size = 20, face = "bold"))+
    theme(plot.subtitle=element_text(size = 20, hjust = 0, color="black"))+
    theme(axis.title.x=element_text(size = 18, hjust = 0.5, color="black"))+
    theme(axis.title.y=element_text(size = 18, hjust = 0.5, color="black"))
  
  
  df_frz <- df %>% dplyr::distinct(Precursor, Annotation, FragmentCharge)
  df_frz <- data.frame(table(df_frz$FragmentCharge)); colnames(df_frz) <- c('fr_z', 'count')
  df_frz$freq <- df_frz$count / sum(df_frz$count)
  # Plot I: Fragment ion charge distribution
  p9 <- ggplot(df_frz, aes(fr_z, freq))+ theme_bw()+
    geom_col(color = allColor[9], fill = allFill[9], width = 0.5)+
    geom_text(aes(label = paste0(sprintf('%.2f', freq * 100), '%')),
              size = 6, hjust = 0.5, vjust = 0, position = "stack", color = 'black')+
    labs(title = 'Fragment ion charge', x = 'Fragment ion charge', y = 'Frequency in percent (%)')+
    scale_y_continuous(labels = scales::percent)+
    theme_bw()+
    theme(legend.text = element_text(size = 18, color = "black"),legend.position = 'top',
          legend.title = element_text(size = 18, color="black") ,
          panel.grid.major =element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"))+
    theme(panel.grid =element_blank())+
    theme(axis.text = element_text(size = 15,color = "black"))+
    theme(axis.text.x = element_text(size = 15,color = "black", angle = 0))+
    theme(plot.title = element_text(size = 20, face = "bold"))+
    theme(plot.subtitle=element_text(size = 20, hjust = 0, color="black"))+
    theme(axis.title.x=element_text(size = 18, hjust = 0.5, color="black"))+
    theme(axis.title.y=element_text(size = 18, hjust = 0.5, color="black"))
  
  pdf(file = paste0(rlt_dir, '/', pdfName), width = 16, height = 16)
  p <- ggarrange(p1, p2, p3, p4, p5, p6, p7, p8, p9 , labels = c("A", "B", "C", "D", "E", "F", "G", "H", "I"), font.label = list(size = 14, face = "bold"), ncol = 3, nrow = 3)
  print(p)
  if (length(dev.list()) > 0){
    dev.off()
  }
}

# --------- main -------------
# read library
libPath <- 'isoformsemi_library_QC.tsv'
df <- read.delim(libPath, stringsAsFactors = F, check.names = F)

# set color style
#allColors = c("#071f3b", '#145774', '#965e11', '#823616', '#3f302d')
#allFills = c('#03172f', '#06445e', '#834d03', '#682407', '#281a18')
#allFills = c('#0d3b6f', '#1b749b', '#b97416', '#af491e', '#624c47')
allColors = c('#639d98', '#1e5d64', '#203d26', '#b49259', '#7d2a16')
allFills = c('#75bcb6', '#2b828c', '#2f5b38', '#e0b66e', '#a6381d')
image(x = 1:5, y = 1, z = as.matrix(1:5), col = colorRampPalette(allFills)(5))

# plotting
plotting(df, pdfName = 'isoformsemi_library_QC.pdf',
         allColor = allColors[c(1, 2, 3, 4, 5, 2, 3, 1, 4)],
         allFill = allFills[c(1, 2, 3, 4, 5, 2, 3, 1, 4)])



reviewedfull_library_QC.tsv
reviewedsemi_library_QC.tsv
isoformfull_library_QC.tsv
isoformsemi_library_QC.tsv

min(df$PrecursorMz)
max(df$PrecursorMz)
a<-unique(df$ModifiedPeptideSequence)
a[1:10]
