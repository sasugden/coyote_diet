#### PREPARE WORKSPACE ####
library(SIBER)
library(simmr)
library(ggplot2)
library(randomForest)
library(gridExtra)
library(car)
library(magrittr)
library(reshape2)

# Upload data.
data_single_point <- read.csv("raw_data/data_single_point.csv",
                              fileEncoding="UTF-8-BOM")
data_all <- read.csv("raw_data/data_all_SIA_samples.csv",
                     fileEncoding="UTF-8-BOM")

# Set graphical parameters.
arguments <- list()

arguments$community.hulls.args <- list(col = 1, lty = 1, lwd = 1)
arguments$group.ellipses.args  <- list(n = 100, p.interval=0.40, lty = 1, lwd = 2)
arguments$group.hull.args      <- list(lty = 2)

# Set the parameters for running JAGS.
parms <- list()
parms$n.iter <- 2 * 10^4   # number of iterations to run the model for
parms$n.burnin <- 1 * 10^3 # discard the first set of values
parms$n.thin <- 10     # thin the posterior by this many
parms$n.chains <- 2        # run this many chains

# Define the priors.
priors <- list()
priors$R <- 1 * diag(2)
priors$k <- 2
priors$tau.mu <- 1.0E-3

# Check the number of samples available in different groupings.
data_single_point %>% dplyr::group_by(HABITAT, SEX) %>% dplyr::count()
data_all %>% dplyr::group_by(HABITAT, PAIRED) %>% dplyr::count()
data_all %>% dplyr::group_by(HABITAT) %>% dplyr::count()

# Sex-based differences in isotopes, controlled for location.
shapiro.test(data_single_point$d13C) # p < 0.001
leveneTest(data_single_point$d13C ~ data_single_point$SEX)
ggplot(data_single_point, aes(x=d13C)) + geom_histogram()

ggplot(subset(data_single_point, is.na(SEX)=="FALSE"),
       aes(x=SEX, y=d13C)) + geom_boxplot() + 
  facet_wrap(~LOCATION) + theme_bw()
anova(lm(d13C ~ HABITAT + SEX, data_single_point))

shapiro.test(data_single_point$d15N) # p < 0.001
leveneTest(data_single_point$d15N ~ data_single_point$SEX)
ggplot(data_single_point, aes(x=d15N)) + geom_histogram()

ggplot(subset(data_single_point, is.na(SEX)=="FALSE"), 
       aes(x=SEX, y=d15N)) + geom_boxplot() +
  facet_wrap(~LOCATION) + theme_bw()
anova(lm(d15N ~ HABITAT + SEX, data_single_point))

#### HAIR-CLAW ISOTOPIC COMPARISON ####
data_hair_claw <- read.csv("raw_data/data_hair_claw.csv")

# Calculate means for each sample type.
with(data_hair_claw, tapply(d13C, Type, mean))
with(data_hair_claw, tapply(d15N, Type, mean))

# Hotelling's T2 test.
with(data_hair_claw, ICSNP::HotellingsT2(cbind(d13C, d15N) ~ Type))

############################################### PART ONE: LAND-USE BASED EFFECTS ###########################################
#### SIBER STANDARD ELLIPSES ####
# (a) Ellipses comparing coyotes by land use. ####
SIBER_landuse_data <- data_single_point[,c("d13C","d15N","LOCATION","HABITAT")]
colnames(SIBER_landuse_data) <- c("iso1","iso2","group","community")

SIBER_landuse_data$group <- as.integer(SIBER_landuse_data$community)
SIBER_landuse_data$community <- rep(1)
SIBER_landuse_data <- createSiberObject(SIBER_landuse_data)

groupMetricsML(SIBER_landuse_data) #1.5 = Urban, 1.2=Matrix, 1.1=Greenspace, 1.4=Suburban, 1.3=Rural
plotSiberObject(SIBER_landuse_data,
                ax.pad = 2, 
                hulls = F, arguments$community.hulls.args, 
                ellipses = T, arguments$group.ellipses.args,
                group.hulls = F, arguments$group.hull.args,
                bty = "L",
                iso.order = c(1,2),
                xlab=expression({delta}^13*C~'\u2030'),
                ylab=expression({delta}^15*N~'\u2030'),
                cex = 0.5)
plotGroupEllipses(SIBER_landuse_data, n = 100, p.interval = 0.4,
                  ci.mean = T, lty = 2, lwd = 2) 

SIBER_landuse_posterior <- siberMVN(SIBER_landuse_data, parms, priors)
SIBER_landuse_ellipses <- siberEllipses(SIBER_landuse_posterior)

siberDensityPlot(SIBER_landuse_ellipses,
                 xticklabels = c("Urban (no GPS)", "Matrix", "Greenspace", "Suburban", "Rural"), 
                 xlab = c("Location"),
                 ylab = expression("Standard Ellipse Area " ('\u2030' ^2) ),
                 bty = "L",
                 las = 1,
                 main = "SIBER ellipses - with land use data", ylim=c(0,15))

# Create a function to extract probability densities.
extract <- function(column){
  data <- hdrcde::hdr(column, prob = c(50,75,95))
  vector <- c(data$hdr[1,1], data$hdr[2,1], data$hdr[3,1], data$hdr[3,2], data$hdr[2,2], data$hdr[1,2], data$mode)
  names(vector) <- c("Min_95", "Min_75", "Min_50", "Max_50", "Max_75", "Max_95", "Mode")
  vector
}

# Save probability densities to a data frame for plotting.
SIBER_landuse_densities <- data.frame( #1.5 = Urban, 1.2=Matrix, 1.1=Greenspace, 1.4=Suburban, 1.3=Rural
  urban = extract(SIBER_landuse_ellipses[,1]),
  matrix = extract(SIBER_landuse_ellipses[,2]),
  greenspace = extract(SIBER_landuse_ellipses[,3]),
  suburban = extract(SIBER_landuse_ellipses[,4]),
  rural = extract(SIBER_landuse_ellipses[,5])
)

SIBER_landuse_densities <- as.data.frame(t(SIBER_landuse_densities))
SIBER_landuse_densities$Landuse <- rownames(SIBER_landuse_densities)
SIBER_landuse_densities$Landuse <- factor(SIBER_landuse_densities$Landuse, 
                                  levels=c("rural", "suburban", "greenspace", "matrix", "urban"))

# (b) Overlap among land use-based ellipses. ####
ellipse1 <- "1.1" # Greenspace
ellipse2 <- "1.2" # Matrix
ellipse3 <- "1.3" # Rural
ellipse4 <- "1.4" # Suburban

# 40% confidence interval ellipses.
overlap_mtrx_green <- bayesianOverlap(ellipse2, ellipse1, SIBER_landuse_posterior,
                                   draws = 100, p.interval = 0.40, n = 100)
overlap_mtrx_rur <- bayesianOverlap(ellipse2, ellipse3, SIBER_landuse_posterior,
                                   draws = 100, p.interval = 0.40, n = 100)
overlap_mtrx_sub <- bayesianOverlap(ellipse2, ellipse4, SIBER_landuse_posterior,
                                   draws = 100, p.interval = 0.40, n = 100)
overlap_green_rur <- bayesianOverlap(ellipse1, ellipse3, SIBER_landuse_posterior,
                                    draws = 100, p.interval = 0.40, n = 100)
overlap_green_sub <- bayesianOverlap(ellipse1, ellipse4, SIBER_landuse_posterior,
                                    draws = 100, p.interval = 0.40, n = 100)
overlap_rur_sub <- bayesianOverlap(ellipse3, ellipse4, SIBER_landuse_posterior,
                                     draws = 100, p.interval = 0.40, n = 100)

overlap_mtrx_green_prop <- (overlap_mtrx_green[,3] / (overlap_mtrx_green[,2] + 
                                                   overlap_mtrx_green[,1] -
                                                   overlap_mtrx_green[,3]))
overlap_mtrx_rur_prop <- (overlap_mtrx_rur[,3] / (overlap_mtrx_rur[,2] + 
                                                   overlap_mtrx_rur[,1] -
                                                   overlap_mtrx_rur[,3]))
overlap_mtrx_sub_prop <- (overlap_mtrx_sub[,3] / (overlap_mtrx_sub[,2] + 
                                                   overlap_mtrx_sub[,1] -
                                                   overlap_mtrx_sub[,3]))
overlap_green_rur_prop <- (overlap_green_rur[,3] / (overlap_green_rur[,2] + 
                                                    overlap_green_rur[,1] -
                                                    overlap_green_rur[,3]))
overlap_green_sub_prop <- (overlap_green_sub[,3] / (overlap_green_sub[,2] + 
                                                    overlap_green_sub[,1] -
                                                    overlap_green_sub[,3]))
overlap_rur_sub_prop <- (overlap_rur_sub[,3] / (overlap_rur_sub[,2] + 
                                                  overlap_rur_sub[,1] -
                                                  overlap_rur_sub[,3]))

vector <- c(mean(overlap_mtrx_green_prop), # 23.60%
            mean(overlap_mtrx_sub_prop), # 19.52%
            mean(overlap_mtrx_rur_prop), # 8.84%
            mean(overlap_green_sub_prop), # 44.96%
            mean(overlap_green_rur_prop), # 28.83%
            mean(overlap_rur_sub_prop)) # 41.48%

# Assemble data.
SIBER_landuse_overlap <- as.data.frame(100*vector)
colnames(SIBER_landuse_overlap) <- c("overlap_40")
rownames(SIBER_landuse_overlap) <- c("Mtrx_Green", "Mtrx_Sub", "Mtrx_Rur", "Green_Sub",
                                     "Green_Rur", "Rur_Sub")

rm(ellipse1, ellipse2, ellipse3, ellipse4, vector,
   overlap_mtrx_green, overlap_mtrx_rur, overlap_mtrx_sub, overlap_green_rur, overlap_green_sub, overlap_rur_sub,
   overlap_mtrx_green_prop, overlap_mtrx_sub_prop, overlap_mtrx_rur_prop, overlap_green_sub_prop,
   overlap_green_rur_prop, overlap_rur_sub_prop)

# (c) Layman metrics by land use. ####
SIBER_landuse_data <- data_single_point[,c("d13C","d15N","LOCATION","HABITAT")]
colnames(SIBER_landuse_data) <- c("iso1","iso2","group","community")
SIBER_landuse_data$community <- factor(SIBER_landuse_data$community,
                                       levels= c("rural", "suburban", "greenspace", "matrix", "urban"))
SIBER_landuse_data$community <- as.integer(SIBER_landuse_data$community)

layman_landuse <- data.frame(matrix(ncol=5, nrow=6))

for(i in c(1:5)){
  temp <- subset(SIBER_landuse_data, community==i) # Greenspace
  temp$group <- c(1:nrow(temp))
  temp <- createSiberObject(temp)
  temp <- communityMetricsML(temp)
  layman_landuse[,i] <- temp
  rm(temp)
}

colnames(layman_landuse) <- c("rural", "suburban", "greenspace", "matrix", "urban")
rownames(layman_landuse) <- c("NR", "NC", "TA", "CD", "NND", "SDNND")
layman_landuse$Measure <- rownames(layman_landuse)
layman_landuse$Measure <- forcats::fct_inorder(layman_landuse$Measure)

# (d) Mixing model: by land use (grouped, concentration-dependent (unadjusted), discrimination 1.0/3.5) ####
mix <- as.matrix(data_single_point[,c("d13C", "d15N")])
colnames(mix) <- c("d13C", "d15N")

sources <- read.csv("mixing_model_inputs/sources.csv")
s_names <- as.vector(sources$Source)
s_means <- matrix(c(sources$Meand13C, sources$Meand15N), nrow=7, ncol=2)
s_sds <- matrix(c(sources$SDd13C, sources$SDd15N), ncol=2, nrow=7)

discrimination <- read.csv("mixing_model_inputs/discrimination_2.0_3.5.csv")
c_means <- matrix(c(discrimination$Meand13C, discrimination$Meand15N), ncol=2, nrow=7)
c_sds <- matrix(c(discrimination$SDd13C, discrimination$SDd15N), ncol=2, nrow=7)

conc <- read.csv("mixing_model_inputs/concentrations.csv")
conc <- matrix(c(conc$d13C, conc$d15N), nrow=7, ncol=2)

grp <- data_single_point$HABITAT

simmr_model_load = simmr_load(mixtures=mix,
                          source_names=s_names,
                          source_means=s_means,
                          source_sds=s_sds,
                          correction_means=c_means,
                          correction_sds=c_sds,
                          group=grp)

plot(simmr_model_load, group=1:5, xlab=expression(paste(delta^13, "C (\u2030)",sep="")), 
     ylab=expression(paste(delta^15, "N (\u2030)",sep="")))

simmr_model_output = simmr_mcmc(simmr_model_load)

summary(simmr_model_output, type=c('quantiles', 'statistics'), group=c(1:5))

mm_landuse_main <- summary(simmr_model_output, type=c('statistics'), group=c(1:5))$statistics
for(i in 1:length(mm_landuse_main)){
  mm_landuse_main[[i]] <- mm_landuse_main[[i]][c(2:8),c(1:2)]
  mm_landuse_main[[i]] <- as.data.frame(mm_landuse_main[[i]])
  if(i==1){mm_landuse_main[[i]]$Site <- rep("greenspace")}
  if(i==2){mm_landuse_main[[i]]$Site <- rep("matrix")}
  if(i==3){mm_landuse_main[[i]]$Site <- rep("rural")}
  if(i==4){mm_landuse_main[[i]]$Site <- rep("suburban")}
  if(i==5){mm_landuse_main[[i]]$Site <- rep("urban")}
  mm_landuse_main[[i]]$Source <- rownames(mm_landuse_main[[i]])
}
mm_landuse_main <- rbind(mm_landuse_main[[1]],
                          mm_landuse_main[[2]],
                          mm_landuse_main[[3]],
                         mm_landuse_main[[4]],
                         mm_landuse_main[[5]])
mm_landuse_main[,c(1:2)] <- 100*mm_landuse_main[,c(1:2)]
mm_landuse_main <- mm_landuse_main[order(mm_landuse_main$Source), ]
rownames(mm_landuse_main) <- c(1:nrow(mm_landuse_main))

compare_groups(simmr_model_output, source_name = "Anthropogenic food", group=c(1:5))

rm(simmr_model_output, simmr_model_load)

mm_landuse_main$Source <- factor(mm_landuse_main$Source,
                                 levels=c("Anthropogenic food", "Domestic pets",
                                          "Insects", "Cricetid rodents ", "Small herbivores",
                                          "Ungulates", "Berries"))
mm_landuse_main$Site <- factor(mm_landuse_main$Site, levels=c("rural", "suburban", "greenspace",
                                                              "matrix", "urban"))

ggplot(mm_landuse_main, aes(x=Source, y=mean, fill=Site)) +
  geom_bar(stat="identity", position="dodge") +
  geom_errorbar(aes(x=Source, ymin=mean-sd, ymax=mean+sd), position=position_dodge(width=0.9), color="black", width=0.3) +
  theme_bw() +
  scale_fill_manual(values=c("springgreen3", "tomato1", "goldenrod", "deepskyblue2",  "grey")) +
  theme(axis.text.y=element_text(size=11, color="black"),
        axis.text.x=element_text(size=10, color="black", angle=45, hjust=1),
        axis.title=element_text(size=12, color="black"),
        legend.text=element_text(size=11, color="black"),
        legend.position="right") +
  labs(x="Source", y="Proportion")

#### SIGNIFICANT DIFFERENCES AMONG GROUPS ####
# Examine distributions.
ggplot(data_single_point, aes(d13C)) + geom_histogram()
ggplot(data_single_point, aes(d15N)) + geom_histogram()

# Check for normality.
shapiro.test(data_single_point$d13C) # p < 0.001
shapiro.test(data_single_point$d15N) # p < 0.001

# Check for homogeneity of dispersion
leveneTest(d13C ~ HABITAT, data_single_point) # p < 0.001 **
leveneTest(d15N ~ HABITAT, data_single_point) # p < 0.001 **

# Boxplots to visualize the differences.
ggplot(data_single_point, aes(x=HABITAT, y=d13C)) + geom_boxplot()
ggplot(data_single_point, aes(x=HABITAT, y=d15N)) + geom_boxplot()

# Statistical tests (ANOVA + posthoc test) for land use-based data.
stats_landuse <- list()
stats_landuse$anova <- data.frame()
stats_landuse$anova[1,1] <- oneway.test(d13C ~ HABITAT, data_single_point, var.equal=FALSE)$statistic
stats_landuse$anova[1,2] <- oneway.test(d13C ~ HABITAT, data_single_point, var.equal=FALSE)[["parameter"]][["num df"]]
stats_landuse$anova[1,3] <- oneway.test(d13C ~ HABITAT, data_single_point, var.equal=FALSE)[["parameter"]][["denom df"]]
stats_landuse$anova[1,4] <- oneway.test(d13C ~ HABITAT, data_single_point, var.equal=FALSE)$p.value

stats_landuse$anova[2,1] <- oneway.test(d15N ~ HABITAT, data_single_point, var.equal=FALSE)$statistic
stats_landuse$anova[2,2] <- oneway.test(d15N ~ HABITAT, data_single_point, var.equal=FALSE)[["parameter"]][["num df"]]
stats_landuse$anova[2,3] <- oneway.test(d15N ~ HABITAT, data_single_point, var.equal=FALSE)[["parameter"]][["denom df"]]
stats_landuse$anova[2,4] <- oneway.test(d15N ~ HABITAT, data_single_point, var.equal=FALSE)$p.value

colnames(stats_landuse$anova) <- c("F", "num df", "denom df", "p")
rownames(stats_landuse$anova) <- c("d13C", "d15N")
stats_landuse$anova <- as.data.frame(t(stats_landuse$anova))

pairwise <- userfriendlyscience::posthocTGH(data_single_point$d13C, data_single_point$HABITAT)
stats_landuse$d13C <- as.data.frame(pairwise[["output"]][["tukey"]])

pairwise <- userfriendlyscience::posthocTGH(data_single_point$d15N, data_single_point$HABITAT)
stats_landuse$d15N <- as.data.frame(pairwise[["output"]][["tukey"]])

rm(pairwise)

#### RE-TEST WITHOUT MANGE COYOTES ####
# (1a) Ellipses comparing all five land use types, without mange coyotes. ####
SIBER_no.mange_data <- subset(data_single_point, MANGE==0 | is.na(MANGE)=="TRUE")[,c("d13C","d15N","LOCATION","HABITAT")]
colnames(SIBER_no.mange_data) <- c("iso1","iso2","group","community")

SIBER_no.mange_data$group <- as.integer(SIBER_no.mange_data$community)
SIBER_no.mange_data$community <- rep(1)
SIBER_no.mange_data <- createSiberObject(SIBER_no.mange_data)

groupMetricsML(SIBER_no.mange_data) #1.2 = Matrix, 1.5 = Urban, 1.1=Greenspace, 1.4=Suburban, 1.3=Rural
plotSiberObject(SIBER_no.mange_data,
                ax.pad = 2, 
                hulls = F, community.hulls.args, 
                ellipses = T, group.ellipses.args,
                group.hulls = F, group.hull.args,
                bty = "L",
                iso.order = c(1,2),
                xlab=expression({delta}^13*C~'\u2030'),
                ylab=expression({delta}^15*N~'\u2030'),
                cex = 0.5)
plotGroupEllipses(SIBER_no.mange_data, n = 100, p.interval = 0.4,
                  ci.mean = T, lty = 2, lwd = 2) 

SIBER_no.mange_posterior <- siberMVN(SIBER_no.mange_data, parms, priors)
SIBER_no.mange_ellipses <- siberEllipses(SIBER_no.mange_posterior)

siberDensityPlot(SIBER_no.mange_ellipses,
                 xticklabels = c("Matrix", "Urban (no GPS)", "Greenspace", "Suburban", "Rural"), 
                 xlab = c("Location"),
                 ylab = expression("Standard Ellipse Area " ('\u2030' ^2) ),
                 bty = "L",
                 las = 1,
                 main = "SIBER ellipses - with land use data", ylim=c(0,35))

SIBER_no.mange_densities <- data.frame(  #1.2 = Matrix, 1.5 = Urban, 1.1=Greenspace, 1.4=Suburban, 1.3=Rural
  matrix = extract(SIBER_no.mange_ellipses[,1]),
  urban = extract(SIBER_no.mange_ellipses[,2]),
  greenspace = extract(SIBER_no.mange_ellipses[,3]),
  suburban = extract(SIBER_no.mange_ellipses[,4]),
  rural = extract(SIBER_no.mange_ellipses[,5])
)

SIBER_no.mange_densities <- as.data.frame(t(SIBER_no.mange_densities))
SIBER_no.mange_densities$Landuse <- rownames(SIBER_no.mange_densities)
SIBER_no.mange_densities$Landuse <- factor(SIBER_no.mange_densities$Landuse, 
                                  levels=c("rural", "suburban", "greenspace", "matrix", "urban"))

# (1b) Niche overlap - coyotes without mange, by location ####
ellipse1 <- "1.1" # Greenspace
ellipse2 <- "1.2" # Matrix
ellipse3 <- "1.3" # Rural
ellipse4 <- "1.4" # Suburban

# 40% confidence interval ellipses.
overlap_mtrx_green <- bayesianOverlap(ellipse2, ellipse1, SIBER_landuse_posterior,
                                      draws = 100, p.interval = 0.40, n = 100)
overlap_mtrx_rur <- bayesianOverlap(ellipse2, ellipse3, SIBER_landuse_posterior,
                                    draws = 100, p.interval = 0.40, n = 100)
overlap_mtrx_sub <- bayesianOverlap(ellipse2, ellipse4, SIBER_landuse_posterior,
                                    draws = 100, p.interval = 0.40, n = 100)
overlap_green_rur <- bayesianOverlap(ellipse1, ellipse3, SIBER_landuse_posterior,
                                     draws = 100, p.interval = 0.40, n = 100)
overlap_green_sub <- bayesianOverlap(ellipse1, ellipse4, SIBER_landuse_posterior,
                                     draws = 100, p.interval = 0.40, n = 100)
overlap_rur_sub <- bayesianOverlap(ellipse3, ellipse4, SIBER_landuse_posterior,
                                   draws = 100, p.interval = 0.40, n = 100)

overlap_mtrx_green_prop <- (overlap_mtrx_green[,3] / (overlap_mtrx_green[,2] + 
                                                        overlap_mtrx_green[,1] -
                                                        overlap_mtrx_green[,3]))
overlap_mtrx_rur_prop <- (overlap_mtrx_rur[,3] / (overlap_mtrx_rur[,2] + 
                                                    overlap_mtrx_rur[,1] -
                                                    overlap_mtrx_rur[,3]))
overlap_mtrx_sub_prop <- (overlap_mtrx_sub[,3] / (overlap_mtrx_sub[,2] + 
                                                    overlap_mtrx_sub[,1] -
                                                    overlap_mtrx_sub[,3]))
overlap_green_rur_prop <- (overlap_green_rur[,3] / (overlap_green_rur[,2] + 
                                                      overlap_green_rur[,1] -
                                                      overlap_green_rur[,3]))
overlap_green_sub_prop <- (overlap_green_sub[,3] / (overlap_green_sub[,2] + 
                                                      overlap_green_sub[,1] -
                                                      overlap_green_sub[,3]))
overlap_rur_sub_prop <- (overlap_rur_sub[,3] / (overlap_rur_sub[,2] + 
                                                  overlap_rur_sub[,1] -
                                                  overlap_rur_sub[,3]))

vector1 <- c(mean(overlap_mtrx_green_prop), # 23.60%
             mean(overlap_mtrx_sub_prop), # 19.52%
             mean(overlap_mtrx_rur_prop), # 8.84%
             mean(overlap_green_sub_prop), # 44.96%
             mean(overlap_green_rur_prop), # 28.83%
             mean(overlap_rur_sub_prop)) # 41.48%

# Assemble data.
SIBER_no.mange_overlap <- as.data.frame(100*vector1)
colnames(SIBER_no.mange_overlap) <- c("overlap_40")
rownames(SIBER_no.mange_overlap) <- c("Mtrx_Green", "Mtrx_Sub", "Mtrx_Rur", "Green_Sub",
                                     "Green_Rur", "Rur_Sub")

rm(ellipse1, ellipse2, ellipse3, ellipse4, vector1,
   overlap_mtrx_green, overlap_mtrx_rur, overlap_mtrx_sub, overlap_green_rur, overlap_green_sub, overlap_rur_sub,
   overlap_mtrx_green_prop, overlap_mtrx_sub_prop, overlap_mtrx_rur_prop, overlap_green_sub_prop,
   overlap_green_rur_prop, overlap_rur_sub_prop)

# (1c) Layman metrics comparing all five land use types, without mange coyotes. ####
SIBER_no.mange_data <- subset(data_single_point, MANGE==0 | is.na(MANGE)=="TRUE")[,c("d13C","d15N","LOCATION","HABITAT")]
colnames(SIBER_no.mange_data) <- c("iso1","iso2","group","community")
SIBER_no.mange_data$community <- factor(SIBER_no.mange_data$community,
                                       levels= c("rural", "suburban", "greenspace", "matrix", "urban"))
SIBER_no.mange_data$group <- as.integer(SIBER_no.mange_data$community)
SIBER_no.mange_data$community <- rep(1)

layman_no.mange <- data.frame(matrix(ncol=5, nrow=6))

for(i in c(1:5)){
  temp <- subset(SIBER_no.mange_data, group==i) # Greenspace
  temp$group <- c(1:nrow(temp))
  temp <- createSiberObject(temp)
  temp <- communityMetricsML(temp)
  layman_no.mange[,i] <- temp
  rm(temp)
}

colnames(layman_no.mange) <- c("rural", "suburban", "greenspace", "matrix", "urban")
rownames(layman_no.mange) <- c("NR", "NC", "TA", "CD", "NND", "SDNND")
layman_no.mange$Measure <- rownames(layman_no.mange)
layman_no.mange$Measure <- forcats::fct_inorder(layman_no.mange$Measure)

# (1d) Mixing models without manage coyotes ####
mix <- as.matrix(subset(data_single_point, MANGE==0 | is.na(MANGE)=="TRUE")[,c("d13C", "d15N")])
colnames(mix) <- c("d13C", "d15N")
grp <- subset(data_single_point, MANGE==0 | is.na(MANGE)=="TRUE")$HABITAT

simmr_model_load = simmr_load(mixtures=mix,
                                  source_names=s_names,
                                  source_means=s_means,
                                  source_sds=s_sds,
                                  correction_means=c_means,
                                  correction_sds=c_sds,
                                  group=grp)

plot(simmr_model_load, group=1:5, xlab=expression(paste(delta^13, "C (\u2030)",sep="")), 
     ylab=expression(paste(delta^15, "N (\u2030)",sep="")))

simmr_model_output = simmr_mcmc(simmr_model_load)

summary(simmr_model_output, type=c('quantiles', 'statistics'), group=c(1:5))

mm_no.mange_main <- summary(simmr_model_output, type=c('statistics'), group=c(1:5))$statistics
for(i in 1:length(mm_no.mange_main)){
  mm_no.mange_main[[i]] <- mm_no.mange_main[[i]][c(2:8),c(1:2)]
  mm_no.mange_main[[i]] <- as.data.frame(mm_no.mange_main[[i]])
  if(i==1){mm_no.mange_main[[i]]$Site <- rep("greenspace")}
  if(i==2){mm_no.mange_main[[i]]$Site <- rep("matrix")}
  if(i==3){mm_no.mange_main[[i]]$Site <- rep("rural")}
  if(i==4){mm_no.mange_main[[i]]$Site <- rep("suburban")}
  if(i==5){mm_no.mange_main[[i]]$Site <- rep("urban")}
  mm_no.mange_main[[i]]$Source <- rownames(mm_no.mange_main[[i]])
}
mm_no.mange_main <- rbind(mm_no.mange_main[[1]],
                         mm_no.mange_main[[2]],
                         mm_no.mange_main[[3]],
                         mm_no.mange_main[[4]],
                         mm_no.mange_main[[5]])
mm_no.mange_main[,c(1:2)] <- 100*mm_no.mange_main[,c(1:2)]
mm_no.mange_main <- mm_no.mange_main[order(mm_no.mange_main$Source), ]
rownames(mm_no.mange_main) <- c(1:nrow(mm_no.mange_main))

rm(simmr_model_output, simmr_model_load)

mm_no.mange_main$Source <- factor(mm_no.mange_main$Source,
                                 levels=c("Anthropogenic food", "Domestic pets",
                                          "Insects", "Cricetid rodents ", "Small herbivores", "Ungulates", "Berries"))
mm_no.mange_main$Site <- factor(mm_no.mange_main$Site, levels=c("rural", "suburban", "greenspace", "matrix", "urban"))

ggplot(mm_no.mange_main, aes(x=Source, y=mean, fill=Site)) +
  geom_bar(stat="identity", position="dodge") +
  geom_errorbar(aes(x=Source, ymin=mean-sd, ymax=mean+sd), position=position_dodge(width=0.9), color="black", width=0.3) +
  theme_bw() +
  scale_fill_manual(values=c("springgreen3", "tomato1", "goldenrod", "deepskyblue2",  "grey")) +
  theme(axis.text.y=element_text(size=11, color="black"),
        axis.text.x=element_text(size=10, color="black", angle=45, hjust=1),
        axis.title=element_text(size=12, color="black"),
        legend.text=element_text(size=11, color="black"),
        legend.position="right") +
  labs(x="Source", y="Proportion")

# (2a) Ellipses comparing all urban coyotes with and without mange. ####
SIBER_mange_data <- subset(data_single_point, HABITAT=="matrix")[,c("d13C","d15N","HABITAT","MANGE")]
colnames(SIBER_mange_data) <- c("iso1","iso2","group","community")

SIBER_mange_data$group <- as.integer(SIBER_mange_data$community)
SIBER_mange_data$community <- rep(1)
SIBER_mange_data <- createSiberObject(SIBER_mange_data)

groupMetricsML(SIBER_mange_data) #1.1 = Mange, 1.0 = No mange
plotSiberObject(SIBER_mange_data,
                ax.pad = 2, 
                hulls = F, community.hulls.args, 
                ellipses = T, group.ellipses.args,
                group.hulls = F, group.hull.args,
                bty = "L",
                iso.order = c(1,2),
                xlab=expression({delta}^13*C~'\u2030'),
                ylab=expression({delta}^15*N~'\u2030'),
                cex = 0.5)
plotGroupEllipses(SIBER_mange_data, n = 100, p.interval = 0.4,
                  ci.mean = T, lty = 2, lwd = 2) 

SIBER_mange_posterior <- siberMVN(SIBER_mange_data, parms, priors)
SIBER_mange_ellipses <- siberEllipses(SIBER_mange_posterior)

siberDensityPlot(SIBER_mange_ellipses,
                 xticklabels = c("Mange", "No mange"), 
                 xlab = c("Location"),
                 ylab = expression("Standard Ellipse Area " ('\u2030' ^2) ),
                 bty = "L",
                 las = 1,
                 main = "SIBER ellipses - with land use data", ylim=c(0,35))

SIBER_mange_densities <- data.frame( #1.1 = Mange, 1.0 = No mange
  Mange = extract(SIBER_mange_ellipses[,1]),
  No.Mange = extract(SIBER_mange_ellipses[,2])
)

SIBER_mange_densities <- as.data.frame(t(SIBER_mange_densities))
SIBER_mange_densities$Status <- rownames(SIBER_mange_densities)
SIBER_mange_densities$Status <- factor(SIBER_mange_densities$Status, 
                                           levels=c("Mange", "No.Mange"))



# (2b) Layman metrics, mange vs. no mange ####
SIBER_mange_data <- subset(data_single_point, HABITAT=="matrix")[,c("d13C","d15N","HABITAT","MANGE")]
colnames(SIBER_mange_data) <- c("iso1","iso2","group","community")
SIBER_mange_data$group <- as.integer(SIBER_mange_data$community)
SIBER_mange_data$community <- rep(1)

layman_mange <- data.frame(matrix(ncol=2, nrow=6))

for(i in c(0:1)){
  temp <- subset(SIBER_mange_data, group==i) # Greenspace
  temp$group <- c(1:nrow(temp))
  temp <- createSiberObject(temp)
  temp <- communityMetricsML(temp)
  layman_mange[,i+1] <- temp
  rm(temp)
}

colnames(layman_mange) <- c("no.mange", "mange")
rownames(layman_mange) <- c("NR", "NC", "TA", "CD", "NND", "SDNND")
layman_mange$Measure <- rownames(layman_mange)
layman_mange$Measure <- forcats::fct_inorder(layman_mange$Measure)

# (2c) Mixing models, mange vs. no mange ####
mix <- as.matrix(subset(data_single_point, LOCATION=="urban")[,c("d13C", "d15N")])
colnames(mix) <- c("d13C", "d15N")
grp <- subset(data_single_point,  LOCATION=="urban")$MANGE

simmr_model_load = simmr_load(mixtures=mix,
                                   source_names=s_names,
                                   source_means=s_means,
                                   source_sds=s_sds,
                                   correction_means=c_means,
                                   correction_sds=c_sds,
                                   group=grp)

plot(simmr_model_load, group=1:2, xlab=expression(paste(delta^13, "C (\u2030)",sep="")), 
     ylab=expression(paste(delta^15, "N (\u2030)",sep="")))

simmr_model_output = simmr_mcmc(simmr_model_load)

summary(simmr_model_output, type=c('quantiles', 'statistics'), group=c(1:2))

mm_mange_main <- summary(simmr_model_output, type=c('statistics'), group=c(1:2))$statistics
for(i in 1:length(mm_mange_main)){
  mm_mange_main[[i]] <- mm_mange_main[[i]][c(2:8),c(1:2)]
  mm_mange_main[[i]] <- as.data.frame(mm_mange_main[[i]])
  if(i==1){mm_mange_main[[i]]$Site <- rep("no.mange")}
  if(i==2){mm_mange_main[[i]]$Site <- rep("mange")}
  mm_mange_main[[i]]$Source <- rownames(mm_mange_main[[i]])
}
mm_mange_main <- rbind(mm_mange_main[[1]],
                          mm_mange_main[[2]])
mm_mange_main[,c(1:2)] <- 100*mm_mange_main[,c(1:2)]
mm_mange_main <- mm_mange_main[order(mm_mange_main$Source), ]
rownames(mm_mange_main) <- c(1:nrow(mm_mange_main))

rm(simmr_model_output, simmr_model_load)

mm_mange_main$Source <- factor(mm_mange_main$Source,
                                  levels=c("Anthropogenic food", "Domestic pets",
                                           "Insects", "Cricetid rodents ", "Small herbivores", "Ungulates", "Berries"))

ggplot(mm_mange_main, aes(x=Source, y=mean, fill=Site)) +
  geom_bar(stat="identity", position="dodge") +
  geom_errorbar(aes(x=Source, ymin=mean-sd, ymax=mean+sd), position=position_dodge(width=0.9), color="black", width=0.3) +
  theme_bw() +
  theme(axis.text.y=element_text(size=11, color="black"),
        axis.text.x=element_text(size=10, color="black", angle=45, hjust=1),
        axis.title=element_text(size=12, color="black"),
        legend.text=element_text(size=11, color="black"),
        legend.position="right") +
  labs(x="Source", y="Proportion")

############################################### PART TWO: SEASONAL COMPARISONS ############################################
# Data preparation ####
# Add a variable that combines location and season together.
data_all$Location_Season <- paste0(data_all$LOCATION, "_", data_all$WINTER)
data_all$Location_Season <- as.factor(data_all$Location_Season)
data_all %>% dplyr::group_by(Location_Season) %>% dplyr::count() # How many oberservations in each group?

# Select samples that have season information available. (Some samples had no collection date marked).
data_season <- subset(data_all, WINTER !="NA")
data_season$Location_Season <- factor(data_season$Location_Season,
                                      levels=c("rural_0", "rural_1",
                                               "suburban_0", "suburban_1",
                                               "urban_0", "urban_1"))

#### NICHE CHARACTERISTICS ####
# Niche width #####
SIBER_season_data <- data_season[,c("d13C","d15N","WINTER","LOCATION")]
SIBER_season_plot <- data_season[,c("d13C","d15N","WINTER","LOCATION","Location_Season")]
colnames(SIBER_season_data) <- c("iso1","iso2","group","community")
SIBER_season_data$community <- as.integer(SIBER_season_data$community)
SIBER_season_data$group <- as.integer(as.factor(SIBER_season_data$group))
SIBER_season_data <- createSiberObject(SIBER_season_data)

groupMetricsML(SIBER_season_data) # 1.1=rural/summer, 1.2=rural/winter,
                                  # 2.1=suburban/summer, 2.2=suburban/winter,
                                  # 3.1=urban/summer, 3.2=urban/winter
plotSiberObject(SIBER_season_data,
                ax.pad = 2, 
                hulls = F, arguments$community.hulls.args, 
                ellipses = T, arguments$group.ellipses.args,
                group.hulls = F, arguments$group.hull.args,
                bty = "L",
                iso.order = c(1,2),
                xlab=expression({delta}^13*C~'\u2030'),
                ylab=expression({delta}^15*N~'\u2030'),
                cex = 0.5)
plotGroupEllipses(SIBER_season_data, n = 100, p.interval = 0.99,
                  ci.mean = T, lty = 2, lwd = 2) 

SIBER_season_posterior <- siberMVN(SIBER_season_data, parms, priors)
SIBER_season_ellipses <- siberEllipses(SIBER_season_posterior)

siberDensityPlot(SIBER_season_ellipses,
                 xticklabels = c("Rur/Sum", "Rur/Wint", "Sub/Sum", "Sub/Wint", "Urb/Sum", "Urb/Wint"), 
                 xlab = c("Location/Season"),
                 ylab = expression("Standard Ellipse Area " ('\u2030' ^2) ),
                 bty = "L",
                 las = 1,
                 main = "SIBER ellipses - location + season (all samples)", ylim=c(0,15))

SIBER_season_densities <- data.frame(
  rur.sum = extract(SIBER_season_ellipses[,1]),
  rur.wint = extract(SIBER_season_ellipses[,2]),
  sub.sum = extract(SIBER_season_ellipses[,3]),
  sub.wint = extract(SIBER_season_ellipses[,4]),
  urb.sum = extract(SIBER_season_ellipses[,5]),
  urb.wint = extract(SIBER_season_ellipses[,6])
)

SIBER_season_densities <- as.data.frame(t(SIBER_season_densities))
SIBER_season_densities$Group <- rownames(SIBER_season_densities)
SIBER_season_densities$Group <- factor(SIBER_season_densities$Group, 
                                          levels=c("rur.sum", "rur.wint", "sub.sum", "sub.wint", 
                                                   "urb.sum", "urb.wint"))


# Niche overlap ####
# Pick one of the above options to run (a, b, or c) and then run niche overlap from there.
ellipse1 <- "1.1" # rural/summer
ellipse2 <- "1.2" # rural/winter
ellipse3 <- "2.1" # suburban/summer
ellipse4 <- "2.2" # suburban/winter
ellipse5 <- "3.1" # urban/summer
ellipse6 <- "3.2" # urban/winter

# 40% confidence interval ellipses.
overlap_rursum_rurwint <- bayesianOverlap(ellipse2, ellipse1, SIBER_season_posterior,
                                      draws = 100, p.interval = 0.40, n = 100)
overlap_subsum_subwint <- bayesianOverlap(ellipse4, ellipse3, SIBER_season_posterior,
                                    draws = 100, p.interval = 0.40, n = 100)
overlap_urbsum_urbwint <- bayesianOverlap(ellipse6, ellipse5, SIBER_season_posterior,
                                    draws = 100, p.interval = 0.40, n = 100)

overlap_rursum_subsum <- bayesianOverlap(ellipse1, ellipse3, SIBER_season_posterior,
                                     draws = 100, p.interval = 0.40, n = 100)
overlap_rursum_urbsum <- bayesianOverlap(ellipse1, ellipse5, SIBER_season_posterior,
                                     draws = 100, p.interval = 0.40, n = 100)
overlap_subsum_urbsum <- bayesianOverlap(ellipse3, ellipse5, SIBER_season_posterior,
                                   draws = 100, p.interval = 0.40, n = 100)

overlap_rurwint_subwint <- bayesianOverlap(ellipse2, ellipse4, SIBER_season_posterior,
                                         draws = 100, p.interval = 0.40, n = 100)
overlap_rurwint_urbwint <- bayesianOverlap(ellipse2, ellipse6, SIBER_season_posterior,
                                         draws = 100, p.interval = 0.40, n = 100)
overlap_subwint_urbwint <- bayesianOverlap(ellipse4, ellipse6, SIBER_season_posterior,
                                         draws = 100, p.interval = 0.40, n = 100)

overlap_rursum_rurwint_prop <- (overlap_rursum_rurwint[,3] / (overlap_rursum_rurwint[,2] + 
                                                                overlap_rursum_rurwint[,1] -
                                                                overlap_rursum_rurwint[,3]))
overlap_subsum_subwint_prop <- (overlap_subsum_subwint[,3] / (overlap_subsum_subwint[,2] + 
                                                                overlap_subsum_subwint[,1] -
                                                                overlap_subsum_subwint[,3]))
overlap_urbsum_urbwint_prop <- (overlap_urbsum_urbwint[,3] / (overlap_urbsum_urbwint[,2] + 
                                                                overlap_urbsum_urbwint[,1] -
                                                                overlap_urbsum_urbwint[,3]))

overlap_rursum_subsum_prop <- (overlap_rursum_subsum[,3] / (overlap_rursum_subsum[,2] + 
                                                         overlap_rursum_subsum[,1] -
                                                         overlap_rursum_subsum[,3]))
overlap_rursum_urbsum_prop <- (overlap_rursum_urbsum[,3] / (overlap_rursum_urbsum[,2] + 
                                                         overlap_rursum_urbsum[,1] -
                                                         overlap_rursum_urbsum[,3]))
overlap_subsum_urbsum_prop <- (overlap_subsum_urbsum[,3] / (overlap_subsum_urbsum[,2] + 
                                                         overlap_subsum_urbsum[,1] -
                                                         overlap_subsum_urbsum[,3]))

overlap_rurwint_subwint_prop <- (overlap_rurwint_subwint[,3] / (overlap_rurwint_subwint[,2] + 
                                                             overlap_rurwint_subwint[,1] -
                                                           overlap_rurwint_subwint[,3]))
overlap_rurwint_urbwint_prop <- (overlap_rurwint_urbwint[,3] / (overlap_rurwint_urbwint[,2] + 
                                                             overlap_rurwint_urbwint[,1] -
                                                             overlap_rurwint_urbwint[,3]))
overlap_subwint_urbwint_prop <- (overlap_subwint_urbwint[,3] / (overlap_subwint_urbwint[,2] + 
                                                             overlap_subwint_urbwint[,1] -
                                                             overlap_subwint_urbwint[,3]))

vector1 <- c(mean(overlap_rursum_rurwint_prop),
             mean(overlap_subsum_subwint_prop),
             mean(overlap_urbsum_urbwint_prop),
             mean(overlap_rursum_subsum_prop),
             mean(overlap_rursum_urbsum_prop),
             mean(overlap_subsum_urbsum_prop),
             mean(overlap_rurwint_subwint_prop),
             mean(overlap_rurwint_urbwint_prop),
             mean(overlap_subwint_urbwint_prop))

# Assemble data.
SIBER_season_overlap <- 100*as.data.frame(vector1)
colnames(SIBER_season_overlap) <- c("overlap_40")
rownames(SIBER_season_overlap) <- c("Rur_Sum.v.Wint", "Sub_Sum.v.Wint", "Urb_Sum.v.Wint",
                                     "Sum_Rur.v.Sub", "Sum_Rur.v.Urb", "Sum_Sub.v.Urb",
                                     "Wint_Rur.v.Sub", "Wint_Rur.v.Urb", "Wint_Sub.v.Urb")

rm(ellipse1, ellipse2, ellipse3, ellipse4, ellipse5, ellipse6, vector1,
   overlap_rursum_rurwint_prop, overlap_subsum_subwint_prop, overlap_urbsum_urbwint_prop,
   overlap_rursum_subsum_prop, overlap_rursum_urbsum_prop, overlap_subsum_urbsum_prop,
   overlap_rurwint_subwint_prop, overlap_rurwint_urbwint_prop, overlap_subwint_urbwint_prop,
   overlap_rursum_rurwint, overlap_subsum_subwint, overlap_urbsum_urbwint,
   overlap_rursum_subsum, overlap_rursum_urbsum, overlap_subsum_urbsum,
   overlap_rurwint_subwint, overlap_rurwint_urbwint, overlap_subwint_urbwint)

# Layman metrics ####
SIBER_season_data <- data_season[,c("d13C","d15N","WINTER","Location_Season")] # Choose data source here.
colnames(SIBER_season_data) <- c("iso1","iso2","group","community")
SIBER_season_data$community <- as.integer(SIBER_season_data$community)
SIBER_season_data$group <- rep(1)

layman_season <- data.frame(matrix(ncol=6, nrow=6))

for(i in c(1:6)){
  temp <- subset(SIBER_season_data, community==i) # Greenspace
  temp$group <- c(1:nrow(temp))
  temp <- createSiberObject(temp)
  temp <- communityMetricsML(temp)
  layman_season[,i] <- temp
  rm(temp)
}

colnames(layman_season) <- c("rural,summer", "rural,winter",
                             "suburban,summer", "suburban,winter", 
                             "urban,summer", "urban,winter")
rownames(layman_season) <- c("NR", "NC", "TA", "CD", "NND", "SDNND")
layman_season$Measure <- rownames(layman_season)
layman_season$Measure <- forcats::fct_inorder(layman_season$Measure)

#### MIXING MODELS (GROUPED, CONCENTRATION-DEPENDENT (UNADJUSTED), DISCRIMINATION 1.0/3.5) ####
# -Summer model (all sources available) ####
mix <- as.matrix(subset(data_season, WINTER==0)[,c("d13C", "d15N")])
colnames(mix) <- c("d13C", "d15N")
grp <- factor(subset(data_season, WINTER==0)$Location_Season)

sources <- read.csv("mixing_model_inputs/sources.csv")
s_names <- as.vector(sources$Source)
s_means <- matrix(c(sources$Meand13C, sources$Meand15N), nrow=7, ncol=2)
s_sds <- matrix(c(sources$SDd13C, sources$SDd15N), ncol=2, nrow=7)

discrimination <- read.csv("mixing_model_inputs/discrimination_2.0_3.5.csv")
c_means <- matrix(c(discrimination$Meand13C, discrimination$Meand15N), ncol=2, nrow=7)
c_sds <- matrix(c(discrimination$SDd13C, discrimination$SDd15N), ncol=2, nrow=7)

conc <- read.csv("mixing_model_inputs/concentrations.csv")
conc <- matrix(c(conc$d13C, conc$d15N), nrow=7, ncol=2)

simmr_model_load = simmr_load(mixtures=mix,
                          source_names=s_names,
                          source_means=s_means,
                          source_sds=s_sds,
                          correction_means=c_means,
                          correction_sds=c_sds,
                          concentration_means=conc,
                          group=grp)

plot(simmr_model_load, group=1:3, xlab=expression(paste(delta^13, "C (\u2030)",sep="")), 
     ylab=expression(paste(delta^15, "N (\u2030)",sep="")))

simmr_model_output = simmr_mcmc(simmr_model_load)

mm_season_summer <- summary(simmr_model_output, type=c('statistics'), group=c(1:3))$statistics
for(i in 1:length(mm_season_summer)){
  mm_season_summer[[i]] <- mm_season_summer[[i]][c(2:8),c(1:2)]
  mm_season_summer[[i]] <- as.data.frame(mm_season_summer[[i]])
  if(i==1){mm_season_summer[[i]]$Group <- rep("rural,summer")}
  if(i==2){mm_season_summer[[i]]$Group <- rep("suburban,summer")}
  if(i==3){mm_season_summer[[i]]$Group <- rep("urban,summer")}
  mm_season_summer[[i]]$Source <- rownames(mm_season_summer[[i]])
}
mm_season_summer <- rbind(mm_season_summer[[1]],
                          mm_season_summer[[2]],
                          mm_season_summer[[3]])
mm_season_summer[,c(1:2)] <- 100*mm_season_summer[,c(1:2)]
mm_season_summer <- mm_season_summer[order(mm_season_summer$Source), ]
rownames(mm_season_summer) <- c(1:nrow(mm_season_summer))

mm_season_summer$Source <- factor(mm_season_summer$Source,
                                levels=c("Anthropogenic food", "Domestic pets", "Insects", "Cricetid rodents ", "Small herbivores", "Ungulates", "Berries"))

mm_season_summer <- tidyr::separate(mm_season_summer, col=Group, into=c("Site","Season"), sep=",", remove=FALSE)

mm_season_summer$Site <- factor(mm_season_summer$Site, levels=c("rural", "suburban", "urban"))

# -Winter model (remove insects) ####
mix <- as.matrix(subset(data_season, WINTER==1)[,c("d13C", "d15N")])
colnames(mix) <- c("d13C", "d15N")
grp <- factor(subset(data_season, WINTER==1)$Location_Season)

sources <- read.csv("mixing_model_inputs/sources.csv")
sources <- subset(sources, Source !="Insects")
s_names <- as.vector(factor(sources$Source))
s_means <- matrix(c(sources$Meand13C, sources$Meand15N), nrow=6, ncol=2)
s_sds <- matrix(c(sources$SDd13C, sources$SDd15N), ncol=2, nrow=6)

discrimination <- read.csv("mixing_model_inputs/discrimination_2.0_3.5.csv")
discrimination <- subset(discrimination, X !="Insects")
c_means <- matrix(c(discrimination$Meand13C, discrimination$Meand15N), ncol=2, nrow=6)
c_sds <- matrix(c(discrimination$SDd13C, discrimination$SDd15N), ncol=2, nrow=6)

conc <- read.csv("mixing_model_inputs/concentrations.csv")
conc <- subset(conc, Source !="Insects")
conc <- matrix(c(conc$d13C, conc$d15N), nrow=6, ncol=2)

simmr_model_load = simmr_load(mixtures=mix,
                          source_names=s_names,
                          source_means=s_means,
                          source_sds=s_sds,
                          correction_means=c_means,
                          correction_sds=c_sds,
                          concentration_means=conc,
                          group=grp)

plot(simmr_model_load, group=1:3, xlab=expression(paste(delta^13, "C (\u2030)",sep="")), 
     ylab=expression(paste(delta^15, "N (\u2030)",sep="")))

simmr_model_output = simmr_mcmc(simmr_model_load)

mm_season_winter <- summary(simmr_model_output, type=c('statistics'), group=c(1:3))$statistics
for(i in 1:length(mm_season_winter)){
  mm_season_winter[[i]] <- mm_season_winter[[i]][c(2:7),c(1:2)]
  mm_season_winter[[i]] <- as.data.frame(mm_season_winter[[i]])
  if(i==1){mm_season_winter[[i]]$Group <- rep("rural,winter")}
  if(i==2){mm_season_winter[[i]]$Group <- rep("suburban,winter")}
  if(i==3){mm_season_winter[[i]]$Group <- rep("urban,winter")}
  mm_season_winter[[i]]$Source <- rownames(mm_season_winter[[i]])
}
mm_season_winter <- rbind(mm_season_winter[[1]],
                          mm_season_winter[[2]],
                          mm_season_winter[[3]])
mm_season_winter[,c(1:2)] <- 100*mm_season_winter[,c(1:2)]
mm_season_winter <- mm_season_winter[order(mm_season_winter$Source), ]
rownames(mm_season_winter) <- c(1:nrow(mm_season_winter))

mm_season_winter$Source <- factor(mm_season_winter$Source,
                                  levels=c("Anthropogenic food", "Domestic pets", "Cricetid rodents ", "Small herbivores", "Ungulates", "Berries"))

mm_season_winter <- tidyr::separate(mm_season_winter, col=Group, into=c("Site","Season"), sep=",", remove=FALSE)

mm_season_winter$Site <- factor(mm_season_winter$Site, levels=c("rural", "suburban", "urban"))

# -Concatenate and graph ####
mm_season_separated <- rbind(mm_season_summer, mm_season_winter)
mm_season_separated$Season <- factor(mm_season_separated$Season, levels=c("summer","winter"))

df <- data.frame(
  mean=c(0,0,0),
  sd=c(0,0,0),
  Group=c("rural,winter", "suburban,winter", "urban,winter"),
  Site=c("rural","suburban","urban"),
  Season=rep("winter"),
  Source=rep("Insects")
)

mm_season_separated <- rbind(mm_season_separated, df)
rm(df, simmr_model_load, simmr_model_output, mm_season_summer, mm_season_winter)

############################################### PART THREE: MODELING PROCEDURES ############################################
#### I. RANDOM FOREST MODELS: RECLASSIFY 'URBAN/UNKNOWN' COYOTES ####
# (a) Run the random forest model and make new predictions. ####
# Define the model predictors as the isotope values of matrix and greenspace coyotes.
predictors <- subset(data_single_point, HABITAT %in% c("matrix", "greenspace"))[,c("d13C", "d15N")]

# Define the model response as either matrix or greenspace.
response <- as.character(subset(data_single_point, HABITAT %in% c("matrix", "greenspace"))$HABITAT)

# Combine the predictors and response into a single data frame.
rf.data <- data.frame(response, predictors)

# Perform random forest model.
set.seed(2)
rf.raw.results <- randomForest(response~., data = rf.data, ntree = 1000)
print(rf.raw.results)

# Predict matrix/greenspace assignments for unclassified urban coyotes.
newdata <- subset(data_single_point, HABITAT=="urban")[,c("COYOTEID","d13C","d15N")]
rownames(newdata) <- newdata$COYOTEID
newdata$COYOTEID <- NULL
rf.predictions <- as.data.frame(predict(rf.raw.results, newdata, type="response"))
colnames(rf.predictions) <- "HABITAT"
rf.predictions$COYOTEID <- rownames(rf.predictions)
data_predicted <- data_single_point
data_predicted$HABITAT[match(rf.predictions$COYOTEID, data_predicted$COYOTEID)] <- rf.predictions$HABITAT

rm(predictors, response, rf.data, rf.predictions, newdata)

# (b) Calculate and plot SIBER ellipses for predicted data. ####
SIBER_predicted_data <- data_predicted[,c("d13C","d15N","LOCATION","HABITAT")]
SIBER_predicted_plot <- SIBER_predicted_data
SIBER_predicted_plot$HABITAT <- factor(SIBER_predicted_plot$HABITAT)
colnames(SIBER_predicted_data) <- c("iso1","iso2","group","community")
SIBER_predicted_data$group <- as.integer(SIBER_predicted_data$community)
SIBER_predicted_data$community <- rep(1)
SIBER_predicted_data <- createSiberObject(SIBER_predicted_data)

groupMetricsML(SIBER_predicted_data) # 1.4=Matrix, 1.3=Greenspace, 1.2=Suburban, 1.1=Rural

SIBER_predicted_posterior <- siberMVN(SIBER_predicted_data, parms, priors)
SIBER_predicted_ellipses <- siberEllipses(SIBER_predicted_posterior)

siberDensityPlot(SIBER_predicted_ellipses,
                 xticklabels = c("Matrix", "Greenspace", "Suburban", "Rural"), 
                 xlab = c("Location"),
                 ylab = expression("Standard Ellipse Area " ('\u2030' ^2) ),
                 bty = "L",
                 las = 1,
                 main = "SIBER ellipses - with random-forest-modeled data", ylim=c(0,20))

SIBER_predicted_densities <- data.frame( #1.4=Matrix, 1.3=Greenspace, 1.2=Suburban, 1.1=Rural
  matrix = extract(SIBER_predicted_ellipses[,1]),
  greenspace = extract(SIBER_predicted_ellipses[,2]),
  suburban = extract(SIBER_predicted_ellipses[,3]),
  rural = extract(SIBER_predicted_ellipses[,4])
)

SIBER_predicted_densities <- as.data.frame(t(SIBER_predicted_densities))
SIBER_predicted_densities$Landuse <- rownames(SIBER_predicted_densities)
SIBER_predicted_densities$Landuse <- factor(SIBER_predicted_densities$Landuse, 
                                  levels=c("rural", "suburban", "greenspace", "matrix"))

# (c) Calculate niche overlap for the random forest modeled data. ####
ellipse1 <- "1.1" # Rural
ellipse2 <- "1.2" # Suburban
ellipse3 <- "1.3" # Greenspace
ellipse4 <- "1.4" # Matrix

# 40% confidence ellipses
overlap_mtrx_green <- bayesianOverlap(ellipse4, ellipse3, SIBER_predicted_posterior,
                                      draws = 100, p.interval = 0.40, n = 100)
overlap_mtrx_rur <- bayesianOverlap(ellipse4, ellipse1, SIBER_predicted_posterior,
                                    draws = 100, p.interval = 0.40, n = 100)
overlap_mtrx_sub <- bayesianOverlap(ellipse4, ellipse2, SIBER_predicted_posterior,
                                    draws = 100, p.interval = 0.40, n = 100)
overlap_green_rur <- bayesianOverlap(ellipse3, ellipse1, SIBER_predicted_posterior,
                                     draws = 100, p.interval = 0.40, n = 100)
overlap_green_sub <- bayesianOverlap(ellipse3, ellipse2, SIBER_predicted_posterior,
                                     draws = 100, p.interval = 0.40, n = 100)
overlap_rur_sub <- bayesianOverlap(ellipse1, ellipse2, SIBER_predicted_posterior,
                                   draws = 100, p.interval = 0.40, n = 100)

overlap_mtrx_green_prop <- (overlap_mtrx_green[,3] / (overlap_mtrx_green[,2] + 
                                                        overlap_mtrx_green[,1] -
                                                        overlap_mtrx_green[,3]))

# Area of overlap / (area of ellipse 1 + area of ellipse 2 - overlap area)
overlap_mtrx_rur_prop <- (overlap_mtrx_rur[,3] / (overlap_mtrx_rur[,2] + 
                                                    overlap_mtrx_rur[,1] -
                                                    overlap_mtrx_rur[,3]))
overlap_mtrx_sub_prop <- (overlap_mtrx_sub[,3] / (overlap_mtrx_sub[,2] + 
                                                    overlap_mtrx_sub[,1] -
                                                    overlap_mtrx_sub[,3]))
overlap_green_rur_prop <- (overlap_green_rur[,3] / (overlap_green_rur[,2] + 
                                                      overlap_green_rur[,1] -
                                                      overlap_green_rur[,3]))
overlap_green_sub_prop <- (overlap_green_sub[,3] / (overlap_green_sub[,2] + 
                                                      overlap_green_sub[,1] -
                                                      overlap_green_sub[,3]))
overlap_rur_sub_prop <- (overlap_rur_sub[,3] / (overlap_rur_sub[,2] + 
                                                  overlap_rur_sub[,1] -
                                                  overlap_rur_sub[,3]))

vector1 <- c(mean(overlap_mtrx_green_prop), # 3.76%
             mean(overlap_mtrx_sub_prop), # 0.913%
             mean(overlap_mtrx_rur_prop), # 1.26%
             mean(overlap_green_sub_prop), # 27.91%
             mean(overlap_green_rur_prop), # 5.78%
             mean(overlap_rur_sub_prop)) # 17.70%

SIBER_predicted_overlap <- 100*as.data.frame(vector1)
colnames(SIBER_predicted_overlap) <- c("overlap_40")
rownames(SIBER_predicted_overlap) <- c("Mtrx_Green", "Mtrx_Sub", "Mtrx_Rur", "Green_Sub",
                                     "Green_Rur", "Rur_Sub")

rm(ellipse1, ellipse2, ellipse3, ellipse4, vector1,
   overlap_mtrx_green, overlap_mtrx_rur, overlap_mtrx_sub, overlap_green_rur, overlap_green_sub, overlap_rur_sub,
   overlap_mtrx_green_prop, overlap_mtrx_sub_prop, overlap_mtrx_rur_prop, overlap_green_sub_prop,
   overlap_green_rur_prop, overlap_rur_sub_prop)

# (d) Calculate layman metrics for random forest modeled data. ####
SIBER_predicted_data <- data_predicted[,c("d13C","d15N","LOCATION","HABITAT")]
colnames(SIBER_predicted_data) <- c("iso1","iso2","group","community")
SIBER_predicted_data$community <- factor(SIBER_predicted_data$community,
                                       levels= c("rural", "suburban", "greenspace", "matrix"))
SIBER_predicted_data$community <- as.integer(SIBER_predicted_data$community)

layman_predicted <- data.frame(matrix(ncol=4, nrow=6))

for(i in c(1:4)){
  temp <- subset(SIBER_predicted_data, community==i) # Greenspace
  temp$group <- c(1:nrow(temp))
  temp <- createSiberObject(temp)
  temp <- communityMetricsML(temp)
  layman_predicted[,i] <- temp
  rm(temp)
}

colnames(layman_predicted) <- c("rural", "suburban", "greenspace", "matrix")
rownames(layman_predicted) <- c("NR", "NC", "TA", "CD", "NND", "SDNND")
layman_predicted$Measure <- rownames(layman_predicted)
layman_predicted$Measure <- forcats::fct_inorder(layman_predicted$Measure)


# (e) Mixing model: predicted data (grouped, concentration-dependent (unadjusted), discrimination 1.0/3.5 ####
data_predicted$HABITAT <- factor(data_predicted$HABITAT,
                                  levels=c("rural", "suburban", "greenspace", "matrix"))
mix <- as.matrix(data_predicted[,c("d13C", "d15N")])
colnames(mix) <- c("d13C", "d15N")

sources <- read.csv("~/Downloads/SIA_project/sources.csv")
s_names <- as.vector(sources$Source)
s_means <- matrix(c(sources$Meand13C, sources$Meand15N), nrow=7, ncol=2)
s_sds <- matrix(c(sources$SDd13C, sources$SDd15N), ncol=2, nrow=7)

discrimination <- read.csv("~/Downloads/SIA_project/discrimination_1.0_3.5.csv")
c_means <- matrix(c(discrimination$Meand13C, discrimination$Meand15N), ncol=2, nrow=7)
c_sds <- matrix(c(discrimination$SDd13C, discrimination$SDd15N), ncol=2, nrow=7)

conc <- read.csv("~/Downloads/SIA_project/concentrations.csv")
conc <- matrix(c(conc$d13C, conc$d15N), nrow=7, ncol=2)

grp <- data_predicted$HABITAT

simmr_model_load = simmr_load(mixtures=mix,
                                  source_names=s_names,
                                  source_means=s_means,
                                  source_sds=s_sds,
                                  correction_means=c_means,
                                  correction_sds=c_sds,
                                  group=grp)

plot(simmr_model_load, group=1:4, xlab=expression(paste(delta^13, "C (\u2030)",sep="")), 
     ylab=expression(paste(delta^15, "N (\u2030)",sep="")))

simmr_model_output = simmr_mcmc(simmr_model_load)

summary(simmr_model_output, type=c('quantiles', 'statistics'), group=c(1:4))

mm_predicted_main <- summary(simmr_model_output, type=c('statistics'), group=c(1:4))$statistics
for(i in 1:length(mm_predicted_main)){
  mm_predicted_main[[i]] <- mm_predicted_main[[i]][c(2:8),c(1:2)]
  mm_predicted_main[[i]] <- as.data.frame(mm_predicted_main[[i]])
  if(i==1){mm_predicted_main[[i]]$Site <- rep("rural")}
  if(i==2){mm_predicted_main[[i]]$Site <- rep("suburban")}
  if(i==3){mm_predicted_main[[i]]$Site <- rep("greenspace")}
  if(i==4){mm_predicted_main[[i]]$Site <- rep("matrix")}
  mm_predicted_main[[i]]$Source <- rownames(mm_predicted_main[[i]])
}
mm_predicted_main <- rbind(mm_predicted_main[[1]],
                         mm_predicted_main[[2]],
                         mm_predicted_main[[3]],
                         mm_predicted_main[[4]])
mm_predicted_main[,c(1:2)] <- 100*mm_predicted_main[,c(1:2)]
mm_predicted_main <- mm_predicted_main[order(mm_predicted_main$Source), ]
rownames(mm_predicted_main) <- c(1:nrow(mm_predicted_main))

rm(simmr_model_output, simmr_model_load)

mm_predicted_main$Source <- factor(mm_predicted_main$Source,
                            levels=c("Anthropogenic food", "Domestic pets", "Insects", "Cricetid rodents ", "Small herbivores", "Ungulates", "Berries"))
mm_predicted_main$Site <- factor(mm_predicted_main$Site,
                                 levels=c("rural", "suburban", "greenspace", "matrix"))

#### II. GLMs WITH LAND USE, FROM RANDOM FOREST PREDICTIONS ####
model_data <- subset(data_predicted,
                     is.na(WINTER)=="FALSE" & is.na(SEX)=="FALSE") # With random forest-modeled data

model_data$MANGE[is.na(model_data$MANGE)=="TRUE"] <- 0
model_data$ADULT[is.na(model_data$ADULT)=="TRUE"] <- 0

ggplot(model_data, aes(x=as.factor(MANGE), y=d13C)) + geom_boxplot() + facet_wrap(~HABITAT)
ggplot(model_data, aes(x=as.factor(MANGE), y=d15N)) + geom_boxplot() + facet_wrap(~HABITAT)

model_data$MANGE[is.na(model_data$MANGE)=="TRUE"] <- 0
model_data$ADULT[is.na(model_data$ADULT)=="TRUE"] <- 0

# Model 1: d13C IN ALL COYOTES, LAND USE, WITHOUT COVARIATES ####
model <- lm(d13C ~ SEX + HABITAT + WINTER + ADULT, 
            model_data, na.action="na.fail")
model_d13C_predict <- dredge(model) # Location, mange, and winter are all top predictors
temp.dredge <- as.data.frame(model_d13C_predict)

model_d13C_predict <- model.avg(model_d13C_predict, subset=delta < 2)
model_d13C_predict <- cbind(confint(model_d13C_predict, level=0.95, full=TRUE)[,1],
              confint(model_d13C_predict, level=0.5, full=TRUE)[,1],
              model_d13C_predict[["coefficients"]][1,],
              confint(model_d13C_predict, level=0.5, full=TRUE)[,2],
              confint(model_d13C_predict, level=0.95, full=TRUE)[,2])
model_d13C_predict <- as.data.frame(model_d13C_predict)
colnames(model_d13C_predict) <- c("psd_2.5", "psd_25", "psd_Coef", "psd_75", "psd_97.5")
model_d13C_predict$merge <- rownames(model_d13C_predict)
model_d13C_predict <- subset(model_d13C_predict, !(merge %in% c("(Intercept)", "RESEARCHERSS")))

weights <- data.frame(matrix(nrow=0, ncol=1))
for(j in c(2:6)){
  temp <- temp.dredge
  temp <- subset(temp, is.na(temp[,j])=="FALSE")
  weights[j-1,1] <- sum(temp[,ncol(temp)])
  rm(temp)
}
rownames(weights) <- colnames(temp.dredge[,c(2:6)])
colnames(weights) <- "cum_weight"
weights[("HABITATsuburban"),] <- weights[("HABITAT"),]
weights[("HABITATgreenspace"),] <- weights[("HABITAT"),]
weights[("HABITATmatrix"),] <- weights[("HABITAT"),]
weights$merge <- rownames(weights)

model_d13C_predict <- merge(model_d13C_predict, weights, by="merge", all=FALSE)

ggplot(model_d13C_predict) +
  ggstance::geom_linerangeh(aes(y=merge, xmin=psd_25, xmax=psd_75),
                            size=2, color="black") + 
  ggstance::geom_linerangeh(aes(y=merge, xmin=psd_2.5, xmax=psd_97.5),
                            color="black", size=0.5) +
  geom_point(aes(y=merge, x=psd_Coef, size=cum_weight),
             fill="white", shape=21, stroke=1, show.legend=TRUE) +
  coord_flip() + theme_bw() +
  labs(x="Coefficient", y="Predictor", title="d13C: all coyotes, land use, random forest predictions") +
  geom_vline(xintercept=0, linetype="dashed", color="black") +
  theme(axis.text.x=element_text(size=9, color="black", angle=45, vjust=1, hjust=1),
        axis.text.y=element_text(size=9, color="black"),
        axis.title=element_text(size=10, color="black"),
        plot.title=element_text(size=10, color="black"))

# Model 2: d15N IN ALL COYOTES, LAND USE, WITHOUT COVARIATES ####
model <- lm(d15N ~ SEX + HABITAT + WINTER + ADULT, 
            model_data, na.action="na.fail")
model_d15N_predict <- dredge(model) # Location, mange, and winter are all top predictors
temp.dredge <- as.data.frame(model_d15N_predict)

model_d15N_predict <- model.avg(model_d15N_predict, subset=delta < 2)
model_d15N_predict <- cbind(confint(model_d15N_predict, level=0.95, full=TRUE)[,1],
              confint(model_d15N_predict, level=0.5, full=TRUE)[,1],
              model_d15N_predict[["coefficients"]][1,],
              confint(model_d15N_predict, level=0.5, full=TRUE)[,2],
              confint(model_d15N_predict, level=0.95, full=TRUE)[,2])
model_d15N_predict <- as.data.frame(model_d15N_predict)
colnames(model_d15N_predict) <- c("psd_2.5", "psd_25", "psd_Coef", "psd_75", "psd_97.5")
model_d15N_predict$merge <- rownames(model_d15N_predict)
model_d15N_predict <- subset(model_d15N_predict, !(merge %in% c("(Intercept)", "RESEARCHERSS")))

weights <- data.frame(matrix(nrow=0, ncol=1))
for(j in c(2:6)){
  temp <- temp.dredge
  temp <- subset(temp, is.na(temp[,j])=="FALSE")
  weights[j-1,1] <- sum(temp[,ncol(temp)])
  rm(temp)
}
rownames(weights) <- colnames(temp.dredge[,c(2:6)])
colnames(weights) <- "cum_weight"
weights[("HABITATsuburban"),] <- weights[("HABITAT"),]
weights[("HABITATgreenspace"),] <- weights[("HABITAT"),]
weights[("HABITATmatrix"),] <- weights[("HABITAT"),]
weights$merge <- rownames(weights)

model_d15N_predict <- merge(model_d15N_predict, weights, by="merge", all=FALSE)

ggplot(model_d15N_predict) +
  ggstance::geom_linerangeh(aes(y=merge, xmin=psd_25, xmax=psd_75),
                            size=2, color="black") + 
  ggstance::geom_linerangeh(aes(y=merge, xmin=psd_2.5, xmax=psd_97.5),
                            color="black", size=0.5) +
  geom_point(aes(y=merge, x=psd_Coef),
             fill="white", shape=21, stroke=1, show.legend=TRUE) +
  coord_flip() + theme_bw() +
  labs(x="Coefficient", y="Predictor", title="d15N: all coyotes, land use, random forest predictions") +
  geom_vline(xintercept=0, linetype="dashed", color="black") +
  theme(axis.text.x=element_text(size=9, color="black", angle=45, vjust=1, hjust=1),
        axis.text.y=element_text(size=9, color="black"),
        axis.title=element_text(size=10, color="black"),
        plot.title=element_text(size=10, color="black"))

## MODELS WITH LAND USE, IGNORING UNCLASSIFIED COYOTES ####
model_data <- subset(data_single_point, is.na(SEX)=="FALSE" & HABITAT != "urban" & is.na(WINTER)=="FALSE") 
model_data$HABITAT <- factor(model_data$HABITAT, levels=c("rural", "suburban", "greenspace", "matrix"))

model_data$MANGE[is.na(model_data$MANGE)=="TRUE"] <- 0
model_data$ADULT[is.na(model_data$ADULT)=="TRUE"] <- 0

# Model 1: d13C IN ALL COYOTES, LAND USE, WITHOUT COVARIATES ####
model <- lm(d13C ~ SEX + HABITAT + WINTER + ADULT, 
            model_data, na.action="na.fail")
model_d13C_raw <- dredge(model) 
temp.dredge <- model_d13C_raw

model_d13C_raw <- model.avg(model_d13C_raw, subset=delta < 2)
model_d13C_raw <- cbind(confint(model_d13C_raw, level=0.95, full=TRUE)[,1],
              confint(model_d13C_raw, level=0.5, full=TRUE)[,1],
              model_d13C_raw[["coefficients"]][1,],
              confint(model_d13C_raw, level=0.5, full=TRUE)[,2],
              confint(model_d13C_raw, level=0.95, full=TRUE)[,2])
model_d13C_raw <- as.data.frame(model_d13C_raw)
colnames(model_d13C_raw) <- c("psd_2.5", "psd_25", "psd_Coef", "psd_75", "psd_97.5")
model_d13C_raw$merge <- rownames(model_d13C_raw)
model_d13C_raw <- subset(model_d13C_raw, !(merge %in% c("(Intercept)", "RESEARCHERSS")))

weights <- data.frame(matrix(nrow=0, ncol=1))
for(j in c(2:6)){
  temp <- temp.dredge
  temp <- subset(temp, is.na(temp[,j])=="FALSE")
  weights[j-1,1] <- sum(temp[,ncol(temp)])
  rm(temp)
}
rownames(weights) <- colnames(temp.dredge[,c(2:6)])
colnames(weights) <- "cum_weight"
weights[("HABITATsuburban"),] <- weights[("HABITAT"),]
weights[("HABITATgreenspace"),] <- weights[("HABITAT"),]
weights[("HABITATmatrix"),] <- weights[("HABITAT"),]
weights$merge <- rownames(weights)

model_d13C_raw <- merge(model_d13C_raw, weights, by="merge", all=FALSE)

ggplot(model_d13C_raw) +
  ggstance::geom_linerangeh(aes(y=merge, xmin=psd_25, xmax=psd_75),
                            size=2, color="black") + 
  ggstance::geom_linerangeh(aes(y=merge, xmin=psd_2.5, xmax=psd_97.5),
                            color="black", size=0.5) +
  geom_point(aes(y=merge, x=psd_Coef),
             fill="white", shape=21, stroke=1, show.legend=TRUE) +
  coord_flip() + theme_bw() +
  labs(x="Coefficient", y="Predictor", title="d13C: by land use, all coyotes, excluding urban-unclassified") +
  geom_vline(xintercept=0, linetype="dashed", color="black") +
  theme(axis.text.x=element_text(size=9, color="black", angle=45, vjust=1, hjust=1),
        axis.text.y=element_text(size=9, color="black"),
        axis.title=element_text(size=10, color="black"),
        plot.title=element_text(size=10, color="black"))

# Model 2: d15N IN ALL COYOTES, LAND USE, WITHOUT COVARIATES ####
model <- lm(d15N ~ SEX + HABITAT + WINTER + ADULT, 
            model_data, na.action="na.fail")
model_d15N_raw <- dredge(model) 

model_d15N_raw <- model.avg(model_d15N_raw, subset=delta < 2)
model_d15N_raw <- cbind(confint(model_d15N_raw, level=0.95, full=TRUE)[,1],
              confint(model_d15N_raw, level=0.5, full=TRUE)[,1],
              model_d15N_raw[["coefficients"]][1,],
              confint(model_d15N_raw, level=0.5, full=TRUE)[,2],
              confint(model_d15N_raw, level=0.95, full=TRUE)[,2])
model_d15N_raw <- as.data.frame(model_d15N_raw)
colnames(model_d15N_raw) <- c("psd_2.5", "psd_25", "psd_Coef", "psd_75", "psd_97.5")
model_d15N_raw$merge <- rownames(model_d15N_raw)
model_d15N_raw <- subset(model_d15N_raw, !(merge %in% c("(Intercept)", "RESEARCHERSS")))

weights <- data.frame(matrix(nrow=0, ncol=1))
for(j in c(2:6)){
  temp <- temp.dredge
  temp <- subset(temp, is.na(temp[,j])=="FALSE")
  weights[j-1,1] <- sum(temp[,ncol(temp)])
  rm(temp)
}
rownames(weights) <- colnames(temp.dredge[,c(2:6)])
colnames(weights) <- "cum_weight"
weights[("HABITATsuburban"),] <- weights[("HABITAT"),]
weights[("HABITATgreenspace"),] <- weights[("HABITAT"),]
weights[("HABITATmatrix"),] <- weights[("HABITAT"),]
weights$merge <- rownames(weights)

model_d15N_raw <- merge(model_d15N_raw, weights, by="merge", all=FALSE)

ggplot(model_d15N_raw) +
  ggstance::geom_linerangeh(aes(y=merge, xmin=psd_25, xmax=psd_75),
                            size=2, color="black") + 
  ggstance::geom_linerangeh(aes(y=merge, xmin=psd_2.5, xmax=psd_97.5),
                            color="black", size=0.5) +
  geom_point(aes(y=merge, x=psd_Coef),
             fill="white", shape=21, stroke=1, show.legend=TRUE) +
  coord_flip() + theme_bw() +
  labs(x="Coefficient", y="Predictor", title="d15N: by land use, all coyotes, excluding urban-unclassified") +
  geom_vline(xintercept=0, linetype="dashed", color="black") +
  theme(axis.text.x=element_text(size=9, color="black", angle=45, vjust=1, hjust=1),
        axis.text.y=element_text(size=9, color="black"),
        axis.title=element_text(size=10, color="black"),
        plot.title=element_text(size=10, color="black"))

############################################### PART FOUR: MIXING MODEL SENSITIVITY #########################################
sensitivity_results <- list()

mix <- as.matrix(data_single_point[,c("d13C", "d15N")])
colnames(mix) <- c("d13C", "d15N")
sources <- read.csv("mixing_model_inputs/sources.csv")
s_names <- as.vector(sources$Source)
s_means <- matrix(c(sources$Meand13C, sources$Meand15N), nrow=7, ncol=2)
s_sds <- matrix(c(sources$SDd13C, sources$SDd15N), ncol=2, nrow=7)

discrimination <- read.csv("mixing_model_inputs/discrimination_2.0_3.0.csv")
c_means <- matrix(c(discrimination$Meand13C, discrimination$Meand15N), ncol=2, nrow=7)
c_sds <- matrix(c(discrimination$SDd13C, discrimination$SDd15N), ncol=2, nrow=7)

conc <- read.csv("mixing_model_inputs/concentrations.csv")
conc <- matrix(c(conc$d13C, conc$d15N), nrow=7, ncol=2)

conc_digest <- read.csv("mixing_model_inputs/concentrations_digest.csv")
conc_digest <- matrix(c(conc_digest$d13C, conc_digest$d15N), nrow=7, ncol=2)

grp <- data_single_point$HABITAT

#### (1-3) 2.0 / 3.0, grouped, independent ####
for(p in c(1:3)){
  if(p==1){simmr_groups = simmr_load(mixtures=mix,
                                     source_names=s_names,
                                     source_means=s_means,
                                     source_sds=s_sds,
                                     correction_means=c_means,
                                     correction_sds=c_sds,
                                     #concentration_means=conc,
                                     group=grp)}
  if(p==2){simmr_groups = simmr_load(mixtures=mix,
                                     source_names=s_names,
                                     source_means=s_means,
                                     source_sds=s_sds,
                                     correction_means=c_means,
                                     correction_sds=c_sds,
                                     concentration_means=conc,
                                     group=grp)}
  if(p==3){simmr_groups = simmr_load(mixtures=mix,
                                     source_names=s_names,
                                     source_means=s_means,
                                     source_sds=s_sds,
                                     correction_means=c_means,
                                     correction_sds=c_sds,
                                     concentration_means=conc_digest,
                                     group=grp)}
  
  #plot(simmr_groups, group=1:3, xlab=expression(paste(delta^13, "C (\u2030)",sep="")), 
  #     ylab=expression(paste(delta^15, "N (\u2030)",sep="")))
  
  simmr_groups_out = simmr_mcmc(simmr_groups)
  
  summary(simmr_groups_out, type=c('statistics'), group=c(1:5))
  
  mm_plot_data <- summary(simmr_groups_out, type=c('statistics'), group=c(1:5))$statistics
  for(i in 1:length(mm_plot_data)){
    mm_plot_data[[i]] <- mm_plot_data[[i]][c(2:8),c(1:2)]
    mm_plot_data[[i]] <- as.data.frame(mm_plot_data[[i]])
    if(i==1){mm_plot_data[[i]]$Location <- rep("greenspace")}
    if(i==2){mm_plot_data[[i]]$Location <- rep("matrix")}
    if(i==3){mm_plot_data[[i]]$Location <- rep("rural")}
    if(i==4){mm_plot_data[[i]]$Location <- rep("suburban")}
    if(i==5){mm_plot_data[[i]]$Location <- rep("urban")}
    mm_plot_data[[i]]$Source <- rownames(mm_plot_data[[i]])
  }
  
  mm_plot_data <- rbind(mm_plot_data[[1]],
                        mm_plot_data[[2]],
                        mm_plot_data[[3]],
                        mm_plot_data[[4]],
                        mm_plot_data[[5]])
  mm_plot_data[,c(1:2)] <- 100*mm_plot_data[,c(1:2)]
  mm_plot_data <- mm_plot_data[order(mm_plot_data$Source), ]
  
  sensitivity_results[[p]] <- as.data.frame(mm_plot_data)
}


#### (4-6): 2.0 / 3.5, grouped, independent ####
mix <- as.matrix(data_single_point[,c("d13C", "d15N")])
colnames(mix) <- c("d13C", "d15N")
discrimination <- read.csv("mixing_model_inputs/discrimination_2.0_3.5.csv") 
c_means <- matrix(c(discrimination$Meand13C, discrimination$Meand15N), ncol=2, nrow=7)
c_sds <- matrix(c(discrimination$SDd13C, discrimination$SDd15N), ncol=2, nrow=7)
grp <- data_single_point$HABITAT

for(p in c(1:3)){
  if(p==1){simmr_groups = simmr_load(mixtures=mix,
                                     source_names=s_names,
                                     source_means=s_means,
                                     source_sds=s_sds,
                                     correction_means=c_means,
                                     correction_sds=c_sds,
                                     #concentration_means=conc,
                                     group=grp)}
  if(p==2){simmr_groups = simmr_load(mixtures=mix,
                                     source_names=s_names,
                                     source_means=s_means,
                                     source_sds=s_sds,
                                     correction_means=c_means,
                                     correction_sds=c_sds,
                                     concentration_means=conc,
                                     group=grp)}
  if(p==3){simmr_groups = simmr_load(mixtures=mix,
                                     source_names=s_names,
                                     source_means=s_means,
                                     source_sds=s_sds,
                                     correction_means=c_means,
                                     correction_sds=c_sds,
                                     concentration_means=conc_digest,
                                     group=grp)}
  
  #plot(simmr_groups, group=1:3, xlab=expression(paste(delta^13, "C (\u2030)",sep="")), 
  #     ylab=expression(paste(delta^15, "N (\u2030)",sep="")))
  
  simmr_groups_out = simmr_mcmc(simmr_groups)
  
  summary(simmr_groups_out, type=c('statistics'), group=c(1:5))
  
  mm_plot_data <- summary(simmr_groups_out, type=c('statistics'), group=c(1:5))$statistics
  for(i in 1:length(mm_plot_data)){
    mm_plot_data[[i]] <- mm_plot_data[[i]][c(2:8),c(1:2)]
    mm_plot_data[[i]] <- as.data.frame(mm_plot_data[[i]])
    if(i==1){mm_plot_data[[i]]$Location <- rep("greenspace")}
    if(i==2){mm_plot_data[[i]]$Location <- rep("matrix")}
    if(i==3){mm_plot_data[[i]]$Location <- rep("rural")}
    if(i==4){mm_plot_data[[i]]$Location <- rep("suburban")}
    if(i==5){mm_plot_data[[i]]$Location <- rep("urban")}
    mm_plot_data[[i]]$Source <- rownames(mm_plot_data[[i]])
  }
  
  mm_plot_data <- rbind(mm_plot_data[[1]],
                        mm_plot_data[[2]],
                        mm_plot_data[[3]],
                        mm_plot_data[[4]],
                        mm_plot_data[[5]])
  mm_plot_data[,c(1:2)] <- 100*mm_plot_data[,c(1:2)]
  mm_plot_data <- mm_plot_data[order(mm_plot_data$Source), ]
  
  sensitivity_results[[p+3]] <- as.data.frame(mm_plot_data)
}

#### (7-9): 2.5 / 3.5, grouped, independent ####
mix <- as.matrix(data_single_point[,c("d13C", "d15N")])
colnames(mix) <- c("d13C", "d15N")
discrimination <- read.csv("mixing_model_inputs/discrimination_2.5_3.5.csv") 
c_means <- matrix(c(discrimination$Meand13C, discrimination$Meand15N), ncol=2, nrow=7)
c_sds <- matrix(c(discrimination$SDd13C, discrimination$SDd15N), ncol=2, nrow=7)
grp <- data_single_point$HABITAT

for(p in c(1:3)){
  if(p==1){simmr_groups = simmr_load(mixtures=mix,
                                     source_names=s_names,
                                     source_means=s_means,
                                     source_sds=s_sds,
                                     correction_means=c_means,
                                     correction_sds=c_sds,
                                     #concentration_means=conc,
                                     group=grp)}
  if(p==2){simmr_groups = simmr_load(mixtures=mix,
                                     source_names=s_names,
                                     source_means=s_means,
                                     source_sds=s_sds,
                                     correction_means=c_means,
                                     correction_sds=c_sds,
                                     concentration_means=conc,
                                     group=grp)}
  if(p==3){simmr_groups = simmr_load(mixtures=mix,
                                     source_names=s_names,
                                     source_means=s_means,
                                     source_sds=s_sds,
                                     correction_means=c_means,
                                     correction_sds=c_sds,
                                     concentration_means=conc_digest,
                                     group=grp)}
  
  #plot(simmr_groups, group=1:3, xlab=expression(paste(delta^13, "C (\u2030)",sep="")), 
  #     ylab=expression(paste(delta^15, "N (\u2030)",sep="")))
  
  simmr_groups_out = simmr_mcmc(simmr_groups)
  
  summary(simmr_groups_out, type=c('statistics'), group=c(1:5))
  
  mm_plot_data <- summary(simmr_groups_out, type=c('statistics'), group=c(1:5))$statistics
  for(i in 1:length(mm_plot_data)){
    mm_plot_data[[i]] <- mm_plot_data[[i]][c(2:8),c(1:2)]
    mm_plot_data[[i]] <- as.data.frame(mm_plot_data[[i]])
    if(i==1){mm_plot_data[[i]]$Location <- rep("greenspace")}
    if(i==2){mm_plot_data[[i]]$Location <- rep("matrix")}
    if(i==3){mm_plot_data[[i]]$Location <- rep("rural")}
    if(i==4){mm_plot_data[[i]]$Location <- rep("suburban")}
    if(i==5){mm_plot_data[[i]]$Location <- rep("urban")}
    mm_plot_data[[i]]$Source <- rownames(mm_plot_data[[i]])
  }
  
  mm_plot_data <- rbind(mm_plot_data[[1]],
                        mm_plot_data[[2]],
                        mm_plot_data[[3]],
                        mm_plot_data[[4]],
                        mm_plot_data[[5]])
  mm_plot_data[,c(1:2)] <- 100*mm_plot_data[,c(1:2)]
  mm_plot_data <- mm_plot_data[order(mm_plot_data$Source), ]
  
  sensitivity_results[[p+6]] <- as.data.frame(mm_plot_data)
}

#### (10-12): 0.5 / 3.5, grouped, independent ####
mix <- as.matrix(data_single_point[,c("d13C", "d15N")])
colnames(mix) <- c("d13C", "d15N")
discrimination <- read.csv("mixing_model_inputs/discrimination_0.5_3.5.csv") # Derbridge et al.
c_means <- matrix(c(discrimination$Meand13C, discrimination$Meand15N), ncol=2, nrow=7)
c_sds <- matrix(c(discrimination$SDd13C, discrimination$SDd15N), ncol=2, nrow=7)
grp <- data_single_point$HABITAT

for(p in c(1:3)){
  if(p==1){simmr_groups = simmr_load(mixtures=mix,
                                     source_names=s_names,
                                     source_means=s_means,
                                     source_sds=s_sds,
                                     correction_means=c_means,
                                     correction_sds=c_sds,
                                     #concentration_means=conc,
                                     group=grp)}
  if(p==2){simmr_groups = simmr_load(mixtures=mix,
                                     source_names=s_names,
                                     source_means=s_means,
                                     source_sds=s_sds,
                                     correction_means=c_means,
                                     correction_sds=c_sds,
                                     concentration_means=conc,
                                     group=grp)}
  if(p==3){simmr_groups = simmr_load(mixtures=mix,
                                     source_names=s_names,
                                     source_means=s_means,
                                     source_sds=s_sds,
                                     correction_means=c_means,
                                     correction_sds=c_sds,
                                     concentration_means=conc_digest,
                                     group=grp)}
  
  #plot(simmr_groups, group=1:3, xlab=expression(paste(delta^13, "C (\u2030)",sep="")), 
  #     ylab=expression(paste(delta^15, "N (\u2030)",sep="")))
  
  simmr_groups_out = simmr_mcmc(simmr_groups)
  
  summary(simmr_groups_out, type=c('statistics'), group=c(1:5))
  
  mm_plot_data <- summary(simmr_groups_out, type=c('statistics'), group=c(1:5))$statistics
  for(i in 1:length(mm_plot_data)){
    mm_plot_data[[i]] <- mm_plot_data[[i]][c(2:8),c(1:2)]
    mm_plot_data[[i]] <- as.data.frame(mm_plot_data[[i]])
    if(i==1){mm_plot_data[[i]]$Location <- rep("greenspace")}
    if(i==2){mm_plot_data[[i]]$Location <- rep("matrix")}
    if(i==3){mm_plot_data[[i]]$Location <- rep("rural")}
    if(i==4){mm_plot_data[[i]]$Location <- rep("suburban")}
    if(i==5){mm_plot_data[[i]]$Location <- rep("urban")}
    mm_plot_data[[i]]$Source <- rownames(mm_plot_data[[i]])
  }
  
  mm_plot_data <- rbind(mm_plot_data[[1]],
                        mm_plot_data[[2]],
                        mm_plot_data[[3]],
                        mm_plot_data[[4]],
                        mm_plot_data[[5]])
  mm_plot_data[,c(1:2)] <- 100*mm_plot_data[,c(1:2)]
  mm_plot_data <- mm_plot_data[order(mm_plot_data$Source), ]
  
  sensitivity_results[[p+9]] <- as.data.frame(mm_plot_data)
}


#### (13-15): 4.25 / 3.05, grouped, independent ####
mix <- as.matrix(data_single_point[,c("d13C", "d15N")])
colnames(mix) <- c("d13C", "d15N")
discrimination <- read.csv("mixing_model_inputs/discrimination_4.25_3.05.csv") # Derbridge et al.
c_means <- matrix(c(discrimination$Meand13C, discrimination$Meand15N), ncol=2, nrow=7)
c_sds <- matrix(c(discrimination$SDd13C, discrimination$SDd15N), ncol=2, nrow=7)
grp <- data_single_point$HABITAT

for(p in c(1:3)){
  if(p==1){simmr_groups = simmr_load(mixtures=mix,
                                     source_names=s_names,
                                     source_means=s_means,
                                     source_sds=s_sds,
                                     correction_means=c_means,
                                     correction_sds=c_sds,
                                     #concentration_means=conc,
                                     group=grp)}
  if(p==2){simmr_groups = simmr_load(mixtures=mix,
                                     source_names=s_names,
                                     source_means=s_means,
                                     source_sds=s_sds,
                                     correction_means=c_means,
                                     correction_sds=c_sds,
                                     concentration_means=conc,
                                     group=grp)}
  if(p==3){simmr_groups = simmr_load(mixtures=mix,
                                     source_names=s_names,
                                     source_means=s_means,
                                     source_sds=s_sds,
                                     correction_means=c_means,
                                     correction_sds=c_sds,
                                     concentration_means=conc_digest,
                                     group=grp)}
  
  #plot(simmr_groups, group=1:3, xlab=expression(paste(delta^13, "C (\u2030)",sep="")), 
  #     ylab=expression(paste(delta^15, "N (\u2030)",sep="")))
  
  simmr_groups_out = simmr_mcmc(simmr_groups)
  
  summary(simmr_groups_out, type=c('statistics'), group=c(1:5))
  
  mm_plot_data <- summary(simmr_groups_out, type=c('statistics'), group=c(1:5))$statistics
  for(i in 1:length(mm_plot_data)){
    mm_plot_data[[i]] <- mm_plot_data[[i]][c(2:8),c(1:2)]
    mm_plot_data[[i]] <- as.data.frame(mm_plot_data[[i]])
    if(i==1){mm_plot_data[[i]]$Location <- rep("greenspace")}
    if(i==2){mm_plot_data[[i]]$Location <- rep("matrix")}
    if(i==3){mm_plot_data[[i]]$Location <- rep("rural")}
    if(i==4){mm_plot_data[[i]]$Location <- rep("suburban")}
    if(i==5){mm_plot_data[[i]]$Location <- rep("urban")}
    mm_plot_data[[i]]$Source <- rownames(mm_plot_data[[i]])
  }
  
  mm_plot_data <- rbind(mm_plot_data[[1]],
                        mm_plot_data[[2]],
                        mm_plot_data[[3]],
                        mm_plot_data[[4]],
                        mm_plot_data[[5]])
  mm_plot_data[,c(1:2)] <- 100*mm_plot_data[,c(1:2)]
  mm_plot_data <- mm_plot_data[order(mm_plot_data$Source), ]
  
  sensitivity_results[[p+12]] <- as.data.frame(mm_plot_data)
}

#### (16-18): 3.4 / 3.0, grouped, independent ####
mix <- as.matrix(data_single_point[,c("d13C", "d15N")])
colnames(mix) <- c("d13C", "d15N")
discrimination <- read.csv("mixing_model_inputs/discrimination_3.4_3.0.csv") # Derbridge et al.
c_means <- matrix(c(discrimination$Meand13C, discrimination$Meand15N), ncol=2, nrow=7)
c_sds <- matrix(c(discrimination$SDd13C, discrimination$SDd15N), ncol=2, nrow=7)
grp <- data_single_point$HABITAT

for(p in c(1:3)){
  if(p==1){simmr_groups = simmr_load(mixtures=mix,
                                     source_names=s_names,
                                     source_means=s_means,
                                     source_sds=s_sds,
                                     correction_means=c_means,
                                     correction_sds=c_sds,
                                     #concentration_means=conc,
                                     group=grp)}
  if(p==2){simmr_groups = simmr_load(mixtures=mix,
                                     source_names=s_names,
                                     source_means=s_means,
                                     source_sds=s_sds,
                                     correction_means=c_means,
                                     correction_sds=c_sds,
                                     concentration_means=conc,
                                     group=grp)}
  if(p==3){simmr_groups = simmr_load(mixtures=mix,
                                     source_names=s_names,
                                     source_means=s_means,
                                     source_sds=s_sds,
                                     correction_means=c_means,
                                     correction_sds=c_sds,
                                     concentration_means=conc_digest,
                                     group=grp)}
  
  #plot(simmr_groups, group=1:3, xlab=expression(paste(delta^13, "C (\u2030)",sep="")), 
  #     ylab=expression(paste(delta^15, "N (\u2030)",sep="")))
  
  simmr_groups_out = simmr_mcmc(simmr_groups)
  
  summary(simmr_groups_out, type=c('statistics'), group=c(1:5))
  
  mm_plot_data <- summary(simmr_groups_out, type=c('statistics'), group=c(1:5))$statistics
  for(i in 1:length(mm_plot_data)){
    mm_plot_data[[i]] <- mm_plot_data[[i]][c(2:8),c(1:2)]
    mm_plot_data[[i]] <- as.data.frame(mm_plot_data[[i]])
    if(i==1){mm_plot_data[[i]]$Location <- rep("greenspace")}
    if(i==2){mm_plot_data[[i]]$Location <- rep("matrix")}
    if(i==3){mm_plot_data[[i]]$Location <- rep("rural")}
    if(i==4){mm_plot_data[[i]]$Location <- rep("suburban")}
    if(i==5){mm_plot_data[[i]]$Location <- rep("urban")}
    mm_plot_data[[i]]$Source <- rownames(mm_plot_data[[i]])
  }
  
  mm_plot_data <- rbind(mm_plot_data[[1]],
                        mm_plot_data[[2]],
                        mm_plot_data[[3]],
                        mm_plot_data[[4]],
                        mm_plot_data[[5]])
  mm_plot_data[,c(1:2)] <- 100*mm_plot_data[,c(1:2)]
  mm_plot_data <- mm_plot_data[order(mm_plot_data$Source), ]
  
  sensitivity_results[[p+15]] <- as.data.frame(mm_plot_data)
}


#### CONCATENATE RESULTS ####
openxlsx::write.xlsx(sensitivity_results, "~/sensitivity_results.xlsx", row.names=TRUE)

sensitivity_df <- data.frame(
  Source = sensitivity_results[[1]]$Source,
  Location = sensitivity_results[[2]]$Location)

for(i in 1:length(sensitivity_results)){
  sensitivity_df <- cbind(sensitivity_df, sensitivity_results[[i]]$mean)
}

colnames(sensitivity_df) <- c("Source", "Location",
                              "grp_ind_2.0_3.0", "grp_dep_2.0_3.0", "grp_dig_2.0_3.0",
                              "grp_ind_2.0_3.5", "grp_dep_2.0_3.5", "grp_dig_2.0_3.5",
                              "grp_ind_2.5_3.5", "grp_dep_2.5_3.5", "grp_dig_2.5_3.5",
                              "grp_ind_0.5_3.5", "grp_dep_0.5_3.5", "grp_dig_0.5_3.5",
                              "grp_ind_4.25_3.05", "grp_dep_4.25_3.05", "grp_dig_4.25_3.05",
                              "grp_ind_3.4_3.0", "grp_dep_3.4_3.0", "grp_dig_3.4_3.0")

sensitivity_df_melt <- melt(sensitivity_df)

sensitivity_corr <- Hmisc::rcorr(as.matrix(sensitivity_df[,c(3:ncol(sensitivity_df))]), type="spearman")
sensitivity_corr_res <- as.data.frame(sensitivity_corr[["r"]])
sensitivity_corr_p <- as.data.frame(sensitivity_corr[["P"]])

sensitivity_corr_res$var1 <- rownames(sensitivity_corr_res)
sensitivity_corr_res <- melt(sensitivity_corr_res)

sensitivity_corr_p$var1 <- rownames(sensitivity_corr_p)
sensitivity_corr_p <- melt(sensitivity_corr_p)

colnames(sensitivity_corr_res) <- c("Model1", "Model2", "R")
colnames(sensitivity_corr_p) <- c("Model1", "Model2", "p")

sensitivity_corr_res <- merge(sensitivity_corr_res, sensitivity_corr_p, by=c("Model1", "Model2"), all=TRUE)
sensitivity_corr <- subset(sensitivity_corr_res, is.na(R)=="FALSE")
rm(sensitivity_corr_res, sensitivity_corr_p)

write.csv(sensitivity_df, "~/sensitivity.df.csv")

sensitivity_corr <- subset(sensitivity_corr, is.na(R)=="FALSE")
sensitivity_corr <- subset(sensitivity_corr, R < 1)
