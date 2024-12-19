library(readr)
library(lme4)
library(lmerTest)
library(emmeans)
library(ggplot2)
library(readxl)
library(dplyr)
library(DescTools)
library(geomorph)
library(gridExtra)
library(graphics)

symmetry_metadata <- read_csv("metadata.csv")
linear_measure_average <- read_csv("linear_measure_average.csv")
symmetry_metadata_unbent <- read_csv("metadata_unbent.csv")
assymetry_pairwise_comparisons_body <- read_excel("assymetry pairwise comparisons body.xlsx")
assymetry_pairwise_comparisons_head <- read_excel("assymetry pairwise comparisons head.xlsx")

#---- LINEAR MEASUREMENT PLOTS AND STATS ----

#---- PECTORAL FIN LENGTH ----

#DIFFERENCES IN PECTORAL BETWEEN LEFT AND RIGHT
linear_measure_average$rearing <- factor(linear_measure_average$rearing, 
                                         levels = c("standard", "enriched", "semi-natural", "wild reared"))

ggplot(linear_measure_average, aes(x=rearing, y=pec_diff, color=rearing)) +
  scale_colour_manual(values = c("#D55E00","#E69F00","#56B4E9","#0072B2" ))+
  geom_boxplot(outlier.shape=NA)+theme_bw()+
  geom_jitter(data = linear_measure_average, aes(x=rearing, y=pec_diff),size=2)+
  xlab(NULL) +
  ylab(NULL) +
  theme(
    axis.text.x = element_text(size = 16), 
    axis.text.y = element_text(size = 16))

lm_pec_diff<-  lm(pec_diff ~rearing+length_cm, data=linear_measure_average)
anova(lm_pec_diff)
emmeans(lm_pec_diff, pairwise ~ rearing)

#REARING DIFFERENCES IN PECTORAL
ggplot(linear_measure_average, aes(x=length_cm, y=pec_mean,colour=rearing)) + 
  geom_point(size=2)+ scale_colour_manual(values = c("#D55E00","#E69F00", "#56B4E9","#0072B2"))+
  geom_smooth(method=lm)+theme_bw()+
  xlab(NULL) +
  ylab(NULL) +
  theme(
    axis.text.x = element_text(size = 16), 
    axis.text.y = element_text(size = 16))

pec_allometry<-  lm(pec_mean ~length_cm*rearing, data=linear_measure_average)
summary(pec_allometry)
anova(pec_allometry)
pec_allometry$coefficients
m.lst <- lstrends(pec_allometry, "rearing", var="length_cm")
pairs(m.lst)

#---- MOUTH LENGTH ----

#DIFFERENCES IN MOUTH BETWEEN LEFT AND RIGHT SIDE
ggplot(linear_measure_average, aes(x=rearing, y=mouth_diff, color=rearing)) +
  scale_colour_manual(values = c("#D55E00","#E69F00","#56B4E9","#0072B2" ))+
  geom_boxplot(outlier.shape=NA)+theme_bw()+
  geom_jitter(data = linear_measure_average, aes(x=rearing, y=mouth_diff),size=2)+
  xlab(NULL) +
  ylab(NULL) +
  theme(
    axis.text.x = element_text(size = 16), 
    axis.text.y = element_text(size = 16))

lm_mouth_diff<-  lm(mouth_diff ~rearing+length_cm, data=linear_measure_average)
anova(lm_mouth_diff)
emmeans(lm_mouth_diff, pairwise ~ rearing)

#REARING DIFFERENCES IN MOUTH
ggplot(linear_measure_average, aes(x=length_cm, y=mouth_mean,colour=rearing)) + 
  geom_point()+
  geom_point(size=2)+ scale_colour_manual(values = c("#D55E00","#E69F00", "#56B4E9","#0072B2"))+
  geom_smooth(method=lm)+theme_bw()+
  xlab(NULL) +
  ylab(NULL) +
  theme(
    axis.text.x = element_text(size = 16), 
    axis.text.y = element_text(size = 16))

mouth_allometry<-  lm(mouth_mean ~length_cm*rearing, data=linear_measure_average)
summary(mouth_allometry)
anova(mouth_allometry)

mouth_allometry$coefficients
m.lst <- lstrends(mouth_allometry, "rearing", var="length_cm")
pairs(m.lst)

#---- EYE LENGTH ----

#DIFFERENCES IN EYE BETWEEN LEFT AND RIGHT SIDE
ggplot(linear_measure_average, aes(x=rearing, y=eye_width_diff, color=rearing)) +
  scale_colour_manual(values = c("#D55E00","#E69F00","#56B4E9","#0072B2" ))+
  geom_boxplot(outlier.shape=NA)+theme_bw()+
  geom_jitter(data = linear_measure_average, aes(x=rearing, y=eye_width_diff),size=2)+
  xlab(NULL) +
  ylab(NULL) +
  theme(
    axis.text.x = element_text(size = 16), 
    axis.text.y = element_text(size = 16))

lm_eye_width_diff<-  lm(eye_width_diff ~rearing+length_cm, data=linear_measure_average)
anova(lm_eye_width_diff)
emmeans(lm_eye_width_diff, pairwise ~ rearing)

# REARING DIFFERENCES IN EYE WIDTH
ggplot(linear_measure_average, aes(x=length_cm, y=eye_width_mean,colour=rearing)) + 
  geom_point()+
  geom_point(size=2)+ scale_colour_manual(values = c("#D55E00","#E69F00", "#56B4E9","#0072B2"))+
  geom_smooth(method=lm)+theme_bw()+
  xlab(NULL) +
  ylab(NULL) +
  theme(
    axis.text.x = element_text(size = 16), 
    axis.text.y = element_text(size = 16))

eye_width_allometry<-  lm(eye_width_mean ~length_cm*rearing, data=linear_measure_average)
summary(eye_width_allometry)
anova(eye_width_allometry)

eye_width_allometry$coefficients
m.lst <- lstrends(eye_width_allometry, "rearing", var="length_cm")
pairs(m.lst)


#---- GEOMETRIC MORPHOMETRICS - BODY ----
asymmetry_geometric <- readland.tps("tpsleft_right_lm_unbent_19.TPS",specID=c("ID"))

corrected_asymmetry_geometric<-gpagen(asymmetry_geometric)
coords<-corrected_asymmetry_geometric$coords

plotAllSpecimens(corrected_asymmetry_geometric$coords, mean=TRUE, label=TRUE)
gm.prcomp_asymmetry_geometric<-gm.prcomp(corrected_asymmetry_geometric$coords)
plot(gm.prcomp_asymmetry_geometric)

PCA_asymmetry_geometric<-gm.prcomp_asymmetry_geometric
summary(PCA_asymmetry_geometric)

PCA_scores<-PCA_asymmetry_geometric$x
write.csv(PCA_scores,"PCA_scores.csv")
PCA_scores<-data.frame(PCA_scores)
PCA_scores<-PCA_scores %>% select(Comp1,Comp2,Comp3)

M <- mshape(corrected_asymmetry_geometric$coords)
plot(M)

#PLOTTING SHAPE OF PC1
PC1 <- PCA_asymmetry_geometric$x[,1]
preds <- shape.predictor(corrected_asymmetry_geometric$coords, x= PC1, Intercept = FALSE, 
                         pred1 = min(PC1), pred2 = max(PC1)) 
plotRefToTarget(M, preds$pred1)
plotRefToTarget(M, preds$pred2)

#PLOTTING SHAPE OF PC2
PC2 <- PCA_asymmetry_geometric$x[,2]
preds <- shape.predictor(corrected_asymmetry_geometric$coords, x= PC2, Intercept = FALSE, 
                         pred1 = min(PC2), pred2 = max(PC2)) 
plotRefToTarget(M, preds$pred1)#min
plotRefToTarget(M, preds$pred2)#max

#PLOTTING PCA SCORES
PCA_scores_meta<-cbind(PCA_scores,symmetry_metadata_unbent)
#write.csv(PCA_scores_meta,"PCA_scores_meta.csv")

PCA_scores_meta$side<-factor(PCA_scores_meta$side)
ggplot(PCA_scores_meta, aes(Comp1 ,Comp2, color=side))+
  geom_point()+
  stat_ellipse()+
  facet_grid(rows = vars(rearing))+theme_bw()+
  theme(text = element_text(size=25))

gdf<- geomorph.data.frame(corrected_asymmetry_geometric,
                          Length=PCA_scores_meta$length,
                          side=PCA_scores_meta$side,
                          rearing=PCA_scores_meta$rearing,
                          weight=PCA_scores_meta$weight,
                          rearing_side=PCA_scores_meta$rearing_side) 


fit11 <- procD.lm(coords ~ side*rearing*Length, data=gdf, iter=999, print.progress = T)
anova(fit11)

pleth.posthoc<- pairwise(fit11,groups=gdf$rearing_side,covariate=gdf$Length)
summary(pleth.posthoc,test.type = 'var')

MD <- morphol.disparity(coords ~ side*rearing*Length,groups= gdf$rearing_side,
                        data = gdf, iter = 999)
summary_MD<-summary(MD)
PV.dist<-summary_MD[["PV.dist"]]
PV.dist.Pval<-summary_MD[["PV.dist.Pval"]]

ggplot(assymetry_pairwise_comparisons_body ,aes(contrast,Z,fill=significance))+
  geom_col()+
  scale_fill_manual(values = c("grey","red"))+  theme_bw()+
  theme(axis.text.x = element_text(angle = 60,hjust = 1))+
  theme(
    axis.text.x = element_text(angle = 60, hjust = 1, size = 12), 
    axis.text.y = element_text(size = 20), 
    axis.title = element_text(size = 20), 
    legend.text = element_text(size = 20),
    legend.title = element_text(size = 20))


# ---- GEOMETRIC MORPHOMETRICS - HEAD ----

asymmetry_geometric <- readland.tps("tpsleft_right_lm_11.TPS",specID=c("ID"))

corrected_asymmetry_geometric<-gpagen(asymmetry_geometric)
coords<-corrected_asymmetry_geometric$coords

plotAllSpecimens(corrected_asymmetry_geometric$coords, mean=TRUE, label=TRUE)
gm.prcomp_asymmetry_geometric<-gm.prcomp(corrected_asymmetry_geometric$coords)
plot(gm.prcomp_asymmetry_geometric)

PCA_asymmetry_geometric<-gm.prcomp_asymmetry_geometric
summary(PCA_asymmetry_geometric)

PCA_scores<-PCA_asymmetry_geometric$x
write.csv(PCA_scores,"PCA_scores_head.csv")
PCA_scores<-data.frame(PCA_scores)
PCA_scores<-PCA_scores %>% select(Comp1,Comp2,Comp3)

M <- mshape(corrected_asymmetry_geometric$coords)
plot(M)
PC <- PCA_asymmetry_geometric$x[,2]
preds <- shape.predictor(corrected_asymmetry_geometric$coords, x= PC, Intercept = FALSE, 
                         pred1 = min(PC), pred2 = max(PC)) 
plotRefToTarget(M, preds$pred1)
plotRefToTarget(M, preds$pred2)

PCA_scores_meta<-cbind(PCA_scores,symmetry_metadata)
#write.csv(PCA_scores_meta,"PCA_scores_meta_head.csv")

PCA_scores_meta$side<-factor(PCA_scores_meta$side)
ggplot(PCA_scores_meta, aes(Comp1 ,Comp2, color=side))+
  geom_point()+
  stat_ellipse()+
  facet_grid(rows = vars(rearing))+theme_bw()+
  theme(text = element_text(size=25))

ggplot(data=symmetry_metadata,aes(x=rearing,y=length))+
  geom_boxplot(aes(fill=rearing),outlier.shape = NA)+
  geom_jitter(data=symmetry_metadata,aes(x=rearing,y=length),size=0.5)+
  scale_fill_manual(values = c("#56B4E9",
                               "#0072B2",
                               "#E69F00",
                               "#F0E442"))+
  theme(panel.background = element_blank())+theme(axis.line = element_line(colour = "black"))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5))

gdf<- geomorph.data.frame(corrected_asymmetry_geometric,
                          Length=PCA_scores_meta$length,
                          side=PCA_scores_meta$side,
                          rearing=PCA_scores_meta$rearing,
                          weight=PCA_scores_meta$weight,
                          rearing_side=PCA_scores_meta$rearing_side) 

fit11 <- procD.lm(coords ~ side*rearing*Length, data=gdf, iter=999, print.progress = T)
anova(fit11)

pleth.posthoc<- pairwise(fit11,groups=gdf$rearing_side,covariate=gdf$Length)
summary(pleth.posthoc,test.type = 'var')

MD <- morphol.disparity(coords ~ side*rearing*Length,groups= gdf$rearing_side,
                        data = gdf, iter = 999)

ggplot(assymetry_pairwise_comparisons_head ,aes(contrast,Z,fill=significance))+
  geom_col()+
  scale_fill_manual(values = c("red","grey"))+  theme_bw()+
  theme(axis.text.x = element_text(angle = 60,hjust = 1))+
  theme(
    axis.text.x = element_text(angle = 60, hjust = 1, size = 12), 
    axis.text.y = element_text(size = 20), 
    axis.title = element_text(size = 20), 
    legend.text = element_text(size = 20),
    legend.title = element_text(size = 20))