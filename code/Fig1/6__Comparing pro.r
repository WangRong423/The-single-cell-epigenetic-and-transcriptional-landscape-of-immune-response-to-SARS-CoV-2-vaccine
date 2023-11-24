library(tidyverse)
library(RColorBrewer)
library("do")
library(reshape2)

signac_pbmc <- readRDS("/database/wangrong/Results/0712_ATAC+RNA/HIPPO/lightHIPPO/6.13.signac_peak_to_gene.rds")
pro <- table(signac_pbmc$Participants,signac_pbmc$Timepoints,signac_pbmc$celltype)
pro_1 <- ftable(pro)
pro_ra <- prop.table(pro_1,1)
pro_2 <- as.data.frame(pro_ra)
pro_2$Var2 <- factor(pro_2$Var2,levels=c("Day0","Day1","Day3","Day6","Day14","Day30","Day31","Day33","Day36","Day44","Day171","Day185"))
pro_2$Var3 <- factor(pro_2$Var3,levels=c("Naïve B"," Intermediate B","Plasmablasts/Memory B","Naïve CD4+ T cells","CD4+ Tcm","Treg","Naïve CD8+ T cells","CD8+ Tem","MAIT","NK/NKT" ,"CD14+ Mono","CD16+ Mono","Plasmacytoid DC"))
b6<-brewer.pal(12,"Paired")
b7<-brewer.pal(8,"Set2")
c<-c(b6[c(1,2,3,4,5,7,8,9,10,11,12)],b7)
ggplot(pro_2, aes(x=Var1, y=Freq, fill=Var3)) + 
  facet_grid(.~Var2)+
  geom_bar(stat = "identity", width=0.5, col='black') +
  scale_fill_manual(values = c)+
  theme(panel.grid.major=element_blank())+RotatedAxis()+ 
  theme_bw()
antibody <- read.table("/home/wangrong/script/timepoint.csv", sep = ",", header = TRUE)
antibody$样本编号<-sub('LCL','P1',antibody$样本编号)
antibody$样本编号<-sub('LFW','P3',antibody$样本编号)
antibody$样本编号<-sub('LXJ','P4',antibody$样本编号)
antibody$样本编号<-sub('ZLH','P2',antibody$样本编号)
#antibody <- antibody[-45,]
antibody$时间<-mapvalues(antibody$时间, c("0","1","3","6","14","30","31","33","37","45","171","185"),c("Day0","Day1","Day3","Day6","Day14","Day30","Day31","Day33","Day36","Day44","Day171","Day185"))
antibody<-antibody[antibody$时间!="60",]
antibody1<-antibody[,c("样本编号","时间","生物所中和抗体","IgG浓度","NAb浓度")]
names(antibody1) <- c("样本编号","时间","浓度NAb (true virus)","浓度IgG","浓度NAb (in vitro)")
antibody1<-melt(antibody1,id=c("样本编号","时间"))
names(antibody1)=c('Var1','Var2','Var3','Freq')
an_df<-rbind(pro_2,antibody1)
an_df$Var2 <- factor(an_df$Var2,levels=c("Day0","Day1","Day3","Day6","Day14","Day30","Day31","Day33","Day36","Day44","Day171","Day185"))
an_df <- an_df[-666,]
an_df <- an_df[-665,]
ggplot(an_df,aes(x=Var2,y=Freq,colour=Var1,group=Var1,fill=Var1,shape=Var2)) +
  facet_wrap(Var3~.,scale="free_y")+
  scale_shape_manual(values = c(20,20,20,20,20,20,20,20,20,20,20,20))+
  geom_line(size =0.8)+
  geom_point(size=3)+
  theme(axis.line = element_line(colour = "black"))+
  labs( x = 'Cell types', y = 'Cell proportion')+
  scale_color_manual(values = b7)+
  theme_bw()+
  theme(axis.title = element_text(size=20),panel.grid=element_blank())+
  theme(axis.text.x = element_text(angle = -90,size = rel(1.2)))+
  theme(axis.text.y= element_text(size = rel(1.2)))+
  theme(strip.background = element_rect(color="black", fill=c13, size=1.5))
```
```{r}
library(ggplot2)
library(ggthemes)
library(tidyverse)
library(ggalluvial)
library(ggsci)
library(cowplot)
library(tastypie)
pro_3 <- pro_2[pro_2$Var2!="Day185",]
pro_3 <- pro_3[pro_3$Var2!="Day171",]
ggplot(pro_3, aes(x=Var1, y=Freq, fill=Var3, 
  stratum = Var3, alluvium = Var3)) +facet_grid(.~Var2)+
  geom_col(width = 0.4,color=NA)+
  geom_flow(width = 0.4,alpha = 0.2,knot.pos = 0)+ 
  geom_alluvium( width = 0.4,alpha = 0.2,knot.pos = 0)+ 
  scale_fill_manual(values = c)+
  theme_map()+
  theme(axis.text.x=element_text(size=20,vjust = 5),legend.position = 'none')+
  theme(panel.grid.major=element_blank())+
  theme_bw()

ggplot(pro_3, aes(x=Var1, y=Freq, fill=Var3, 
                stratum = Var3, alluvium = Var3)) +facet_grid(.~Var1)+
  geom_col(width = 0.4,color=NA)+
  geom_flow(width = 0.4,alpha = 0.2,knot.pos = 0)+ 
  geom_alluvium( width = 0.4,alpha = 0.2,knot.pos = 0)+ 
  scale_fill_manual(values = c)+
  theme_map()+
  theme(axis.text.x=element_text(size=10,vjust = 5),legend.position = 'none')+
  theme(panel.grid.major=element_blank())
 
ggplot(pro_3, aes(x=Var2, y=Freq, fill=Var3, 
                stratum = Var3, alluvium = Var3)) +facet_grid(.~Var1)+
  geom_col(width = 0.4,color=NA)+
  geom_flow(width = 0.4,alpha = 0.2,knot.pos = 0)+ 
  geom_alluvium( width = 0.4,alpha = 0.2,knot.pos = 0)+ 
  scale_fill_manual(values = c)+
  theme_map()+
  theme(axis.text.x=element_text(angle = 45,size=10),legend.position = 'none')+
  theme(panel.grid.major=element_blank())

save(pro_ra,pro_2,an_df,pro_3,antibody,file = "/database/wangrong/Results/0712_ATAC+RNA/HIPPO/lightHIPPO/6.21.celltype.pro.Rdata")
```

```{r}
#2.15 柱形图统计检验
library(ggplot2)
library(dplyr)
library(ggpubr)

celltypes = c('Naïve B',' Intermediate B','Plasmablasts/Memory B','Naïve CD4+ T cells','CD4+ Tcm','Treg','Naïve CD8+ T cells','CD8+ Tem','MAIT','NK/NKT','CD14+ Mono','CD16+ Mono','Plasmacytoid DC')

pplist = list()
for(i in celltypes){
  cellper_  = pro_3[pro_3$Var3==i,]
  colnames(cellper_) = c('participants','timepoints','celltype','percent')#对选择数据列命名
  cellper_$percent = as.numeric(cellper_$percent)#数值型数据
  cellper_ <- cellper_ %>% group_by(timepoints) %>% mutate(upper =  quantile(percent, 0.75), 
                                                      lower = quantile(percent, 0.25),
                                                      mean = mean(percent),
                                                      median = median(percent))#上下分位数
  #print(i)
  #print(cellper_$median)
  
  pp1 = ggplot(cellper_,aes(x=timepoints,y=percent)) + #ggplot作图
    geom_point(shape = 21,aes(fill=participants),width = 0.25) + 
    stat_summary(fun=mean, geom="point", color="grey60") +
    theme_cowplot() +
    theme(axis.text = element_text(size = 10),axis.title = element_text(size = 10),legend.text = element_text(size = 10),
          legend.title = element_text(size = 10),plot.title = element_text(size = 10,face = 'plain'),
          axis.text.x = element_text(size = 8,angle = 30,vjust = 0.85,hjust = 0.75)) +
    labs(title = i,y='Percentage') +
    geom_errorbar(aes(ymin = lower, ymax = upper),col = "grey60",width =  1)
    
    
  
  ###组间t检验分析
  labely = max(cellper_$percent)
  compare_means(percent ~ timepoints,  data = cellper_)
  my_comparisons <- list( c("Day0", "Day1"),c("Day0", "Day3"),c("Day0", "Day6"),c("Day0", "Day14"),c("Day0", "Day30"),c("Day0", "Day31"),c("Day0", "Day33"),c("Day0", "Day44"))
  pp1 = pp1 + stat_compare_means(comparisons = my_comparisons,size = 3,method = "t.test",label = "p.signif")
    
  pplist[[i]] = pp1
}

library(cowplot)
plot_grid(pplist[[1]],
          pplist[[2]],
          pplist[[3]],
          pplist[[5]],
          pplist[[5]],
          pplist[[6]],
          pplist[[7]],
          pplist[[8]],
          pplist[[9]],
          pplist[[10]],
          pplist[[11]],
          pplist[[12]],
          pplist[[13]])
