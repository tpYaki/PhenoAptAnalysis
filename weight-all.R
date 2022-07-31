install.packages("ggplot2")
install.packages("grid")
install.packages("scales")
install.packages("ggpubr")
install.packages("ggplotify")
install.packages("vcd")
install.packages("lattice")
setwd("/Users/liyaqi/Rproject/PhenoAPT_trial_rank")


##-----全局-------##
dat_0 = read.csv('/Users/liyaqi/PycharmProjects/PhenoAptAnalysis/hpo_id_8_rank_with_vcf.tsv', row.names = 2, header = T, sep = '\t')
dat_1 = read.csv('/Users/liyaqi/PycharmProjects/PhenoAptAnalysis/hpo_id_8_rr_with_vcf.tsv', row.names = 2, header = T, sep = '\t')
dat_2 = read.csv('/Users/liyaqi/PycharmProjects/PhenoAptAnalysis/Weight_8_rank_with_vcf.tsv', row.names = 2, header = T, sep = '\t')
dat_3 = read.csv('/Users/liyaqi/PycharmProjects/PhenoAptAnalysis/Weight_8_rr_with_vcf.tsv', row.names = 2, header = T, sep = '\t')

# data_rank = dat_0[c(5,13,21,29,37,45,46,47)]
##8种直接rank
data_rank = dat_0
data_rank = data_rank[,5:ncol(data_rank)]
data_rr = dat_1


##-----绘图区-------##
map=data_rank
mapname = names(map)
##---统计各个类别的数目----#
rowname = c('TOP1','TOP5','TOP10','TOP20','TOP50','TOP100','>100','filtered_by_REVEL_0.25','filtered_by_REVEL_0.5','filtered_by_REVEL_0.75','filtered_by_CADD_10','filtered_by_CADD_15','filtered_by_CADD_20','not ranked','not_in_filtered_tsv')
map_dat = data.frame()
map_dat = data.frame(row.names = rowname)
for (j in rowname){
  for (i in mapname){
    print(i)
    print(map[,i])
    map_dat[j,i]<-sum(map[,i]==j)
  }
}
map_dat$rank<-rowname
map_dat$rank<-factor(map_dat$rank,levels = rev(map_dat$rank))
##----搁在map_dat里----##
library(reshape2)
data_set<-melt(map_dat,id.vars='rank')
data_set$strategy<-data_set$variable
##-----代码粘贴区------###
data_set$strategy<-gsub('phrank_rank$','Pheno_only_tools',data_set$strategy)
data_set$strategy<-gsub('phrank_intersect_rank','PUMCHpipeline',data_set$strategy)
data_set$strategy<-gsub('phrank_rank_REVEL_0.25','REVEL_0.25',data_set$strategy)
data_set$strategy<-gsub('phrank_rank_REVEL_0.5','REVEL_0.5',data_set$strategy)
data_set$strategy<-gsub('phrank_rank_REVEL_0.75','REVEL_0.75',data_set$strategy)
data_set$strategy<-gsub('phrank_rank_CADD_10','CADD_10',data_set$strategy)
data_set$strategy<-gsub('phrank_rank_CADD_15','CADD_15',data_set$strategy)
data_set$strategy<-gsub('phrank_rank_CADD_20','CADD_20',data_set$strategy)
data_set$strategy<-gsub('phenolyzer_rank$','Pheno_only_tools',data_set$strategy)
data_set$strategy<-gsub('phenolyzer_intersect_rank','PUMCHpipeline',data_set$strategy)
data_set$strategy<-gsub('phenolyzer_rank_REVEL_0.25','REVEL_0.25',data_set$strategy)
data_set$strategy<-gsub('phenolyzer_rank_REVEL_0.5','REVEL_0.5',data_set$strategy)
data_set$strategy<-gsub('phenolyzer_rank_REVEL_0.75','REVEL_0.75',data_set$strategy)
data_set$strategy<-gsub('phenolyzer_rank_CADD_10','CADD_10',data_set$strategy)
data_set$strategy<-gsub('phenolyzer_rank_CADD_15','CADD_15',data_set$strategy)
data_set$strategy<-gsub('phenolyzer_rank_CADD_20','CADD_20',data_set$strategy)
data_set$strategy<-gsub('GADO_rank$','Pheno_only_tools',data_set$strategy)
data_set$strategy<-gsub('GADO_intersect_rank','PUMCHpipeline',data_set$strategy)
data_set$strategy<-gsub('GADO_rank_REVEL_0.25','REVEL_0.25',data_set$strategy)
data_set$strategy<-gsub('GADO_rank_REVEL_0.5','REVEL_0.5',data_set$strategy)
data_set$strategy<-gsub('GADO_rank_REVEL_0.75','REVEL_0.75',data_set$strategy)
data_set$strategy<-gsub('GADO_rank_CADD_10','CADD_10',data_set$strategy)
data_set$strategy<-gsub('GADO_rank_CADD_15','CADD_15',data_set$strategy)
data_set$strategy<-gsub('GADO_rank_CADD_20','CADD_20',data_set$strategy)
data_set$strategy<-gsub('phen2gene_rank$','Pheno_only_tools',data_set$strategy)
data_set$strategy<-gsub('phen2gene_intersect_rank','PUMCHpipeline',data_set$strategy)
data_set$strategy<-gsub('phen2gene_rank_REVEL_0.25','REVEL_0.25',data_set$strategy)
data_set$strategy<-gsub('phen2gene_rank_REVEL_0.5','REVEL_0.5',data_set$strategy)
data_set$strategy<-gsub('phen2gene_rank_REVEL_0.75','REVEL_0.75',data_set$strategy)
data_set$strategy<-gsub('phen2gene_rank_CADD_10','CADD_10',data_set$strategy)
data_set$strategy<-gsub('phen2gene_rank_CADD_15','CADD_15',data_set$strategy)
data_set$strategy<-gsub('phen2gene_rank_CADD_20','CADD_20',data_set$strategy)
data_set$strategy<-gsub('phenoapt_rank$','Pheno_only_tools',data_set$strategy)
data_set$strategy<-gsub('phenoapt_intersect_rank','PUMCHpipeline',data_set$strategy)
data_set$strategy<-gsub('phenoapt_rank_REVEL_0.25','REVEL_0.25',data_set$strategy)
data_set$strategy<-gsub('phenoapt_rank_REVEL_0.5','REVEL_0.5',data_set$strategy)
data_set$strategy<-gsub('phenoapt_rank_REVEL_0.75','REVEL_0.75',data_set$strategy)
data_set$strategy<-gsub('phenoapt_rank_CADD_10','CADD_10',data_set$strategy)
data_set$strategy<-gsub('phenoapt_rank_CADD_15','CADD_15',data_set$strategy)
data_set$strategy<-gsub('phenoapt_rank_CADD_20','CADD_20',data_set$strategy)
data_set$strategy<-gsub('LIRICAL_rank$','Integrated_tools',data_set$strategy)
data_set$strategy<-gsub('Exomiser_rank$','Integrated_tools',data_set$strategy)
data_set$strategy<-gsub('Phen_gen_rank$','Integrated_tools',data_set$strategy)
data_set$strategy<-factor(data_set$strategy,levels=c('Pheno_only_tools', 'PUMCHpipeline', 'REVEL_0.25', 'REVEL_0.5', 'REVEL_0.75', 'CADD_10', 'CADD_15', 'CADD_20', 'Integrated_tools'))

data_set$variable<-gsub('phrank_rank$','phrank',data_set$variable)
data_set$variable<-gsub('phrank_intersect_rank','phrank',data_set$variable)
data_set$variable<-gsub('phrank_rank_REVEL_0.25','phrank',data_set$variable)
data_set$variable<-gsub('phrank_rank_REVEL_0.5','phrank',data_set$variable)
data_set$variable<-gsub('phrank_rank_REVEL_0.75','phrank',data_set$variable)
data_set$variable<-gsub('phrank_rank_CADD_10','phrank',data_set$variable)
data_set$variable<-gsub('phrank_rank_CADD_15','phrank',data_set$variable)
data_set$variable<-gsub('phrank_rank_CADD_20','phrank',data_set$variable)
data_set$variable<-gsub('phenolyzer_rank$','phenolyzer',data_set$variable)
data_set$variable<-gsub('phenolyzer_intersect_rank','phenolyzer',data_set$variable)
data_set$variable<-gsub('phenolyzer_rank_REVEL_0.25','phenolyzer',data_set$variable)
data_set$variable<-gsub('phenolyzer_rank_REVEL_0.5','phenolyzer',data_set$variable)
data_set$variable<-gsub('phenolyzer_rank_REVEL_0.75','phenolyzer',data_set$variable)
data_set$variable<-gsub('phenolyzer_rank_CADD_10','phenolyzer',data_set$variable)
data_set$variable<-gsub('phenolyzer_rank_CADD_15','phenolyzer',data_set$variable)
data_set$variable<-gsub('phenolyzer_rank_CADD_20','phenolyzer',data_set$variable)
data_set$variable<-gsub('GADO_rank$','GADO',data_set$variable)
data_set$variable<-gsub('GADO_intersect_rank','GADO',data_set$variable)
data_set$variable<-gsub('GADO_rank_REVEL_0.25','GADO',data_set$variable)
data_set$variable<-gsub('GADO_rank_REVEL_0.5','GADO',data_set$variable)
data_set$variable<-gsub('GADO_rank_REVEL_0.75','GADO',data_set$variable)
data_set$variable<-gsub('GADO_rank_CADD_10','GADO',data_set$variable)
data_set$variable<-gsub('GADO_rank_CADD_15','GADO',data_set$variable)
data_set$variable<-gsub('GADO_rank_CADD_20','GADO',data_set$variable)
data_set$variable<-gsub('phen2gene_rank$','phen2gene',data_set$variable)
data_set$variable<-gsub('phen2gene_intersect_rank','phen2gene',data_set$variable)
data_set$variable<-gsub('phen2gene_rank_REVEL_0.25','phen2gene',data_set$variable)
data_set$variable<-gsub('phen2gene_rank_REVEL_0.5','phen2gene',data_set$variable)
data_set$variable<-gsub('phen2gene_rank_REVEL_0.75','phen2gene',data_set$variable)
data_set$variable<-gsub('phen2gene_rank_CADD_10','phen2gene',data_set$variable)
data_set$variable<-gsub('phen2gene_rank_CADD_15','phen2gene',data_set$variable)
data_set$variable<-gsub('phen2gene_rank_CADD_20','phen2gene',data_set$variable)
data_set$variable<-gsub('phenoapt_rank$','phenoapt',data_set$variable)
data_set$variable<-gsub('phenoapt_intersect_rank','phenoapt',data_set$variable)
data_set$variable<-gsub('phenoapt_rank_REVEL_0.25','phenoapt',data_set$variable)
data_set$variable<-gsub('phenoapt_rank_REVEL_0.5','phenoapt',data_set$variable)
data_set$variable<-gsub('phenoapt_rank_REVEL_0.75','phenoapt',data_set$variable)
data_set$variable<-gsub('phenoapt_rank_CADD_10','phenoapt',data_set$variable)
data_set$variable<-gsub('phenoapt_rank_CADD_15','phenoapt',data_set$variable)
data_set$variable<-gsub('phenoapt_rank_CADD_20','phenoapt',data_set$variable)
data_set$variable<-gsub('LIRICAL_rank$','LIRICAL',data_set$variable)
data_set$variable<-gsub('Exomiser_rank$','Exomiser',data_set$variable)
data_set$variable<-gsub('Phen_gen_rank$','Phen_gen',data_set$variable)
##------代码粘贴区--------##
write.csv(data_set,"data_Set.csv")
##--compareREVEL\CADD\intersect---###
library(scales)
library(ggplot2)
library(RColorBrewer)
library(RColorBrewer)
# Define the number of colors you want
nb.cols <- length(rowname)
mycolors <- colorRampPalette(brewer.pal(10, "Set3"))(nb.cols)
p <- ggplot(data = data_set, aes(x = variable, y = value, fill = rank)) + 
  
  theme_bw()+
  
  geom_bar(stat='identity',position = 'fill',width = 0.5) +#堆叠图，position = fill 表示堆叠图
  
  facet_grid(~strategy,scales='free_x',space = 'free')+
  
  theme(strip.text = element_text(size = 20))+
  
  labs(x = 'Tools',y = 'Percentage',fill = NULL)+ #定义坐标轴以及图例标题
  
  scale_fill_manual(values = mycolors) + ##自定义颜色，可通过`library(RColorBrewer);display.brewer.all()`来展示所有备选项；geom_point(colour =)也可以设定点颜色
  scale_y_continuous(n.breaks =10,labels = scales::percent) +
  
  guides(fill = guide_legend(ncol = 8, bycol = T, override.aes = list(size = 5))) +#定义图例的布局，1列，排序，图例中色块的大小增大5倍
  
  theme(axis.title.y = element_text(face = 'bold',color = 'black',size = 30),
        axis.title.x = element_text(face = 'bold',color = 'black',size = 30,vjust = -1.2),
        axis.text.y = element_text(face = 'bold',color = 'black',size = 30),
        axis.text.x = element_text(face = 'bold',color = 'black',size = 20,angle = 45,vjust = 1,hjust = 1), #x轴标签偏转45°，并下降0.5
        panel.grid = element_blank(),
        legend.position = 'top',
        legend.key.height = unit(0.6,'cm'),#定义图例中色块的高度
        legend.text = element_text(face = 'bold',color = 'black',size = 25),
        legend.title = element_text(face = 'bold',color = 'black',size = 25))
plot(p)

ggsave(filename = 'benchmark.png',width = 30, height = 10)

###-----三种加权方式是否有显著差异----###
library(ggpubr)
library(RColorBrewer)
dat_3_id <- cbind(newColName3 = rownames(dat_3), dat_3)
rownames(dat_3_id) <- 1:nrow(dat_3)
dat_1_id <- cbind(newColName1 = rownames(dat_1), dat_1)
rownames(dat_1_id) <- 1:nrow(dat_1)
df0_ori<-merge(dat_3_id[,c(1,46,54,38)],dat_1_id[,c(1,38)],by.x = 'newColName3',by.y = 'newColName1',all = TRUE)
# df0_ori = df0_ori[,2:ncol(df0_ori)]
names(df0_ori) = c('patients','Clinical_indication','Clinical_indication+weight','All_hpo+weight','All_hpo')
# df0_ori = rownames_to_column(df0_ori,'patients')
df0=data.frame()
for (i in 1:length(df0_ori[,1])){
  df0[i,'Weight_strategy'] = names(df0_ori)[2]
  df0[i,'Reciprocal_rank'] = df0_ori[i,2]
  df0[i,'patients'] = df0_ori[i,1]
}
for (i in (length(df0_ori[,1])+1):(length(df0_ori[,2])*2)){
  df0[i,'Weight_strategy'] = names(df0_ori)[3]
  df0[i,'Reciprocal_rank'] = df0_ori[i-length(df0_ori[,1]),3]
  df0[i,'patients'] = df0_ori[i-length(df0_ori[,1]),1]
}
for (i in (length(df0_ori[,1])*2+1):(length(df0_ori[,1])*3)){
  df0[i,'Weight_strategy'] = names(df0_ori)[4]
  df0[i,'Reciprocal_rank'] = df0_ori[i-length(df0_ori[,1])*2,4]
  df0[i,'patients'] = df0_ori[i-length(df0_ori[,1])*2,1]
}
for (i in (length(df0_ori[,1])*3+1):(length(df0_ori[,1])*4)){
  df0[i,'Weight_strategy'] = names(df0_ori)[5]
  df0[i,'Reciprocal_rank'] = df0_ori[i-length(df0_ori[,1])*3,5]
  df0[i,'patients'] = df0_ori[i-length(df0_ori[,1])*3,1]
}

###---箱线图----###
# p6<-
#   ggplot(df0,aes(x=Weight_strategy, y=rr, color=Weight_strategy, fill=Weight_strategy))+ 
#   geom_boxplot(alpha=0.2)+ # 如果这里的alpha写道了上面aes()里，则下面的jitter也会半透明
#   geom_jitter()+            # 反正要整活儿，那就多整点   
#   theme_bw()+ 
#   scale_fill_manual(values = c("#FFCCBC", "#9dd3cc", "#e1eeff"))+ #颜色数量和变量数量要一致，或更多
#   scale_color_manual(values = c("#BF360C", "#00695C", "#235784"))+
#   stat_compare_means(comparisons = my_comparisons, label = "p.signif")+#label这里表示选择显著性标记（星号） 
#   stat_compare_means(label.y = 1.29,label.x = 0.7)+
#   #ggtitle('Three different weight strategies on PhenoApt')+
#   theme(axis.title.y = element_text(face = 'bold',color = 'black',size = 15),
#       axis.title.x = element_text(face = 'bold',color = 'black',size = 15),
#       axis.text.y = element_text(face = 'bold',color = 'black',size = 15),
#       axis.text.x = element_text(face = 'bold',color = 'black',size = 7), 
#       panel.grid = element_blank(),
#       legend.position = 'top',
#       legend.key.height = unit(0.6,'cm'),#定义图例中色块的高度
#       legend.text = element_text(face = 'bold',color = 'black',size = 10))
# p6

###---小提琴图----###
comp<-unique(df0$Weight_strategy)
my_comparisons <- list(c(comp[1],comp[2]),c(comp[1],comp[3]),c(comp[1],comp[4]),c(comp[2],comp[3]),c(comp[2],comp[4]),c(comp[3],comp[4]))
p7=ggviolin(df0,x='Weight_strategy', y='Reciprocal_rank', fill='Weight_strategy',palette = c("#FFCCBC", "#9dd3cc", "#e1eeff","#d9d2e9"),trim = TRUE)+
  geom_jitter()+
  theme_bw()+
  #stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "t.test",ref.group = "0.5")+#label这里表示选择显著性标记（星号） 
  #stat_compare_means(label.y = 1.29,label.x = 0.7)+
  #ggtitle('Three different weight strategies on PhenoApt')+
  theme(axis.title.y = element_text(face = 'bold',color = 'black',size = 30),
        axis.title.x = element_text(face = 'bold',color = 'black',size = 30),
        axis.text.y = element_text(face = 'bold',color = 'black',size = 30),
        axis.text.x = element_text(face = 'bold',color = 'black',size = 12), 
        panel.grid = element_blank(),
        legend.position = 'top',
        legend.key.height = unit(0.6,'cm'),#定义图例中色块的高度
        legend.text = element_text(face = 'bold',color = 'black',size = 12),
        legend.title = element_text(face = 'bold',color = 'black',size = 12))
p7

# ##--老式pair图---####
# All_hpo_weight = df0_ori[,'All_hpo+weight']
# Clinical_indication = df0_ori[,"Clinical_indication"]
# Clinical_indication_weight = df0_ori[,"Clinical_indication+weight"]
# p1 = ggpaired(df0_ori, cond2 = "All_hpo+weight", cond1 = "Clinical_indication+weight",
#          fill = "condition", palette = c("","#FFCCBC"))
# p1
# p2 = ggpaired(df0_ori, cond2 = "Clinical_indication+weight", cond1 = "Clinical_indication",
#               fill = "condition", palette = "Set3")
# p2
# p3 = ggpaired(df0_ori, cond1 = "Clinical_indication", cond2 = "All_hpo+weight",
#               fill = "condition")
#               
# p3
# ###----####
##----新pair图----##
MRR = aggregate(df0$Reciprocal_rank, by=list(Category=df0$Weight_strategy), mean)
write.csv(MRR,'3 phenoapt MRR.csv')
df0$Weight_strategy <- as.character(df0$Weight_strategy)
df0$Weight_strategy <- factor(df0$Weight_strategy, levels=unique(df0$Weight_strategy))
Weight_strategy = df0[,'Weight_strategy']
Reciprocal_rank = df0[,'Reciprocal_rank']
p2 = ggplot(df0, aes(x = Weight_strategy, y = Reciprocal_rank)) +
  geom_boxplot(aes(fill = Weight_strategy), alpha = 0.2, col = c("#BF360C", "#00695C", "#235784","#d9d2e9")) +
  geom_point(aes(col = Weight_strategy)) +
  geom_line(aes(group = df0[,'patients']))+
  stat_summary(fun=mean, colour="darkred", geom="point", hape=18, size=3,show.guide = FALSE)+
  stat_summary(fun=mean, colour="darkred", geom="text", show.guide = FALSE, size=5,
               vjust=-0.7, aes( label=round(..y.., digits=1)))+
  theme_bw()+
  theme(axis.title.y = element_text(face = 'bold',color = 'black',size = 30),
      axis.title.x = element_text(face = 'bold',color = 'black',size = 30),
      axis.text.y = element_text(face = 'bold',color = 'black',size = 25),
      axis.text.x = element_text(face = 'bold',color = 'black',size = 12), 
      panel.grid = element_blank(),
      legend.position = 'top',
      legend.key.height = unit(0.6,'cm'),#定义图例中色块的高度
      legend.text = element_text(face = 'bold',color = 'black',size = 12),
      legend.title = element_text(face = 'bold',color = 'black',size = 12))
p2
###-----------###
#install.packages('patchwork')
library(patchwork)
(p7+p2)
myfilename ="weight_compare.png"
ggsave(filename = myfilename,width = 20, height = 10)

##---最老式的pair图--###
# install.packages("PairedData")
# install.packages("MASS")
# install.packages("dplyr")
# library(PairedData)
# library(MASS)
# library(dplyr)
# All_hpo_weight = df0_ori[,'All_hpo+weight']
# Clinical_indication = df0_ori[,"Clinical_indication"]
# Clinical_indication_weight = df0_ori[,"Clinical_indication+weight"]
# my_data <- data.frame( 
#   group = rep(c('All_hpo+weight',"Clinical_indication","Clinical_indication+weight"), each = 32),
#   Reciprocal_rank = c(All_hpo_weight,  Clinical_indication, Clinical_indication_weight)
# )
# my_data
# All_hpo_weight <- subset(my_data,  group == "All_hpo+weight", Reciprocal_rank,
#                  drop = TRUE)
# 
# Clinical_indication_weight <- subset(my_data,  group == "Clinical_indication+weight", Reciprocal_rank,
#                 drop = TRUE)
# Clinical_indication <- subset(my_data,  group == "Clinical_indication", Reciprocal_rank,
#                                      drop = TRUE)
# pd <- paired(Clinical_indication_weight,All_hpo_weight)
# p1 = plot(pd, type = "profile", col="blue") + theme_bw()
# 
# pd <- paired(Clinical_indication, Clinical_indication_weight,palette ='Set3')
# plot(pd,type = "profile") + theme_bw()
# ##-----###

##----hpo、gene分组云图-----###
library(ggplot2)
library(ggfortify)
library(devtools)
dat_comp = dat_1
dat_comp$hpo_id_organ_system_number =as.character(as.numeric(dat_comp$hpo_id_organ_system_number))
data = scale(dat_comp[,5:(ncol(dat_comp)-5)], center=T, scale=T)
data <- cbind(newColName = rownames(data), data)
rownames(data) <- 1:nrow(data)
dat_comp <- cbind(newColName = rownames(dat_comp), dat_comp)
rownames(dat_comp) <- 1:nrow(dat_comp)
data = merge(data,dat_comp[,c(1,5,53)],by.x='newColName',by.y='newColName')
autoplot(prcomp(as.data.frame(lapply(data[,2:(ncol(data)-2)],as.numeric))), data=data,colour = 'hpo_id_organ_system_number',label=FALSE,label.size=0.01,frame=TRUE,frame.type="norm")+theme_bw()
##数据太少，看不出来什么
myfilename ="PCA_hpo_id_organ_system_number_2.png"
ggsave(filename = myfilename,width = 5, height = 2.5)

dat_comp = dat_1
dat_comp$Symbol =as.character(dat_comp$Symbol)
data = scale(dat_comp[,5:(ncol(dat_comp)-5)], center=T, scale=T)
data <- cbind(newColName = rownames(data), data)
rownames(data) <- 1:nrow(data)
dat_comp <- cbind(newColName = rownames(dat_comp), dat_comp)
rownames(dat_comp) <- 1:nrow(dat_comp)
data = merge(data,dat_comp[,c(1,5,53)],by.x='newColName',by.y='newColName')
autoplot(prcomp(as.data.frame(lapply(data[,2:(ncol(data)-2)],as.numeric))), data=dat_comp,colour = 'Symbol',label=FALSE,label.size=1,frame=TRUE,frame.type="norm")+theme_bw()+
  theme(axis.title.y = element_text(face = 'bold',color = 'black',size = 30),
        axis.title.x = element_text(face = 'bold',color = 'black',size = 30),
        axis.text.y = element_text(face = 'bold',color = 'black',size = 25),
        axis.text.x = element_text(face = 'bold',color = 'black',size = 20), 
        panel.grid = element_blank(),
        legend.position = 'right',
        legend.key.height = unit(0.6,'cm'),#定义图例中色块的高度
        legend.text = element_text(face = 'bold',color = 'black',size = 20),
        legend.title = element_text(face = 'bold',color = 'black',size = 20))
##数据太少，看不出来什么
myfilename ="PCA_Symbol_2.png"
ggsave(filename = myfilename,width = 10, height = 8)

##---heatmap of MRR of different Gene--##
df = dat_1
LIRICAL_rank_MRR = aggregate(df$LIRICAL_rank, by=list(type=df$Symbol),mean)
install.packages('pmsampsize')
require('pmsampsize')
SD_d = aggregate(df0$Reciprocal_rank, by=list(Category=df0$Weight_strategy), sd)
pmsampsize(type = 'c',rsquared = 0.8, parameters = 30, intercept = 0.5, sd = 0.41)

##sanky
dat_subgroup_rr = dat_1
dat_subgroup_rr$Frequency <- cut(dat_subgroup_rr$Frequency, breaks = c(-Inf, 0.000001, 0.00001, 0.0001,0.001,0.01, Inf), labels = c("<1E-06","[1E-06~1E-05)","[1E-05~0.01%)","[0.01%~0.1%)","[0.1%~1%)","≥1%"), right=FALSE, include.lowest=TRUE)
dat_subgroup_rr$Gene = dat_subgroup_rr$Symbol
ggplot(data = dat_subgroup_rr,aes(axis1=Gene,axis2 = Variant_class, axis3 = Frequency,axis4 = hpo_id_organ_system_number))  +
  geom_alluvium(aes(fill=Gene), width = 2/5, discern = TRUE) +
  geom_stratum(width = 2/5, discern = TRUE,alpha=0.6) +
  geom_text(stat = "stratum", discern = TRUE, aes(label = after_stat(stratum)
                                                     ),color = "black",size = 6) + 
  geom_label(stat = "stratum", aes(label = after_stat(stratum)), min.y = 200) + 
  scale_x_discrete(limits = c("Gene","Variant_class","MAF","Affected_organ_system_number"), expand = c(.1, .1)) +
  labs(x='Demographic Features', y = 'Number of Patients (N=32)',fill = 'Gene')+
  theme_minimal()+
  theme_bw()+
  theme(axis.title.y = element_text(face = 'bold',color = 'black',size = 30),
        axis.title.x = element_text(face = 'bold',color = 'black',size = 30),
        axis.text.y = element_text(face = 'bold',color = 'black',size = 30),
        axis.text.x = element_text(face = 'bold',color = 'black',size = 20), 
        panel.grid = element_blank(),
        legend.position = 'right',
        legend.key.height = unit(0.6,'cm'),#定义图例中色块的高度
        legend.text = element_text(face = 'bold',color = 'black',size = 20),
        legend.title = element_text(face = 'bold',color = 'black',size = 20))
myfilename ="demographic_sanky_discernT.png"
ggsave(filename = myfilename,width = 20, height = 10)


