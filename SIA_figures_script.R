# Define plot theme ####
plot_theme <- theme(plot.tag = element_text(color="black", face="bold"),
                    strip.text = element_text(size=9),
                    panel.grid = element_blank(),
                    strip.background=element_rect(color="black", size=1),
                    panel.border = element_rect(color="black", fill=NA, size=1),
                    axis.text = element_text(color= "black", size=8),
                    axis.title = element_text(color="black", size=9),
                    legend.title = element_text(color="black", size=9),
                    plot.title=element_text(color="black", size=10),
                    legend.text = element_text(color="black", size=8))

#### FIG. 1: MAP ####
# The map was created using qGIS.
#### FIG. 2: LAND-USE ISOTOPIC RESULTS ####
### (a) Isoscape showing ellipses
plot_data <- data_single_point
plot_data$HABITAT <- factor(plot_data$HABITAT,
                             levels=c("rural", "suburban", "greenspace", "matrix", "urban"))

isoscape <- ggplot() +
  # Add the coyotes
  geom_point(data=plot_data, aes(x=d13C, y=d15N, color = HABITAT, shape=HABITAT), size = 2) +
  stat_ellipse(data=plot_data,
               aes(x=d13C, y=d15N, group = HABITAT, 
                   color = HABITAT), size=1, 
               level = 0.4, alpha=0,
               type = "norm", geom = "polygon") + 
  scale_color_manual(values=c("springgreen3", "goldenrod", "tomato1", "deepskyblue2", "orchid1")) +

  theme_bw() + plot_theme +
  ylab(expression(paste(delta^{15}, "N (\u2030)")))+
  xlab(expression(paste(delta^{13}, "C (\u2030)"))) + 
  labs(tag="a") +
  theme(legend.position="none")

### (b) Niche width
niche.widths <- ggplot() +
  geom_boxplot(data=SIBER_landuse_densities, 
               aes(x=Landuse, ymin=Min_95, ymax=Max_95, 
                   lower=Min_95, middle=Mode, upper=Max_95),
               stat="identity", width=0.4, fatten=NULL, fill="lightgrey") +
  geom_boxplot(data=SIBER_landuse_densities,
               aes(x=Landuse, ymin=Min_75, ymax=Max_75, 
                   lower=Min_75, middle=Mode, upper=Max_75),
               stat="identity", width=0.6, fatten=NULL, fill="grey") +
  geom_boxplot(data=SIBER_landuse_densities, 
               aes(x=Landuse, ymin=Min_50, ymax=Max_50, 
                   lower=Min_50, middle=Mode, upper=Max_50),
               stat="identity", width=0.8, fatten=NULL, fill="darkgrey") +
  geom_point(data=SIBER_landuse_densities, aes(x=Landuse, y=Mode)) +
  theme_bw() + plot_theme +
  labs(x="\n Habitat use", y="Standard ellipse area", tag="b") +
  theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1))

### (c) Niche overlap
overlap.data <- read.csv("overlap2.csv") # Convert overlaps into bar graph format (done manually)
overlap.data$axis <- factor(overlap.data$axis, levels=c("rural", "suburban", "greenspace", "matrix"))
overlap.data$fill <- factor(overlap.data$fill, levels=c("rural", "suburban", "greenspace", "matrix"))

overlap <- ggplot(overlap.data, aes(x=axis, y=overlap_40, fill=fill)) +
  geom_bar(stat="identity") +
  scale_fill_manual(values=c("springgreen3", "goldenrod", "tomato1", "deepskyblue2")) +
  theme_bw() + plot_theme +
  theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1),
        legend.position="none") +
  labs(x="\n Habitat use", y="Ellipse overlap", tag="c") +
  scale_y_continuous(expand = expansion(mult = c(0, .1)))

### (d) Layman metrics
layman_landuse <- reshape2::melt(layman_landuse, id.vars=c("Measure"))
colnames(layman_landuse) <- c("Measure", "Land.use", "value")
layman_landuse$Land.use <- factor(layman_landuse$Land.use,
                                  levels=c("rural", "suburban", "greenspace", "matrix", "urban"))
layman_landuse$Measure

layman <- ggplot(layman_landuse, aes(x=Land.use, y=value, fill=Land.use)) + geom_bar(stat="identity", position="dodge") +
  facet_wrap(~Measure, scales="free_y", nrow=3, ncol=2) +
  theme_bw() + plot_theme +
  scale_fill_manual(values=c("springgreen3", "goldenrod", "tomato1", "deepskyblue2",  "orchid1")) +
  scale_y_continuous(expand = expansion(mult = c(0, .1))) +
  theme(legend.position="none",
        axis.text.x=element_text(angle=45, hjust=1, vjust=1),
        strip.text = element_blank()) +
  labs(x="Habitat use", y="Value", tag="d")

### (e) Mixing model
mm_landuse_main$Source <- factor(mm_landuse_main$Source)
levels(mm_landuse_main$Source)[levels(mm_landuse_main$Source)=="Anthropogenic food"] <- "Anthropogenic\nfood"
levels(mm_landuse_main$Source)[levels(mm_landuse_main$Source)=="Domestic pets"] <- "Domestic\npets"
levels(mm_landuse_main$Source)[levels(mm_landuse_main$Source)=="Cricetid rodents "] <- "Cricetid\nrodents"
levels(mm_landuse_main$Source)[levels(mm_landuse_main$Source)=="Small herbivores"] <- "Small\nherbivores"
mm_landuse_main$Site <- factor(mm_landuse_main$Site,
                               levels=c("rural", "suburban", "greenspace", "matrix", "urban"))
mm_landuse_main$Source <- factor(mm_landuse_main$Source,
                                 levels=c("Anthropogenic\nfood", "Domestic\npets", "Insects", "Cricetid\nrodents",
                                          "Small\nherbivores", "Ungulates", "Berries"))

model <- ggplot(mm_landuse_main, aes(x=Source, y=mean, fill=Site)) +
  geom_bar(stat="identity", position="dodge") +
  geom_errorbar(aes(x=Source, ymin=mean-sd, ymax=mean+sd), position=position_dodge(width=0.9), color="black", width=0.3) +
  theme_bw() + plot_theme +
  scale_fill_manual(values=c("springgreen3", "goldenrod", "tomato1", "deepskyblue2",  "orchid1")) +
  theme(legend.position="none",
        axis.text.x=element_text(angle=45, hjust=1, vjust=1)) +
  labs(x="Source", y="Proportion", tag="e") +
  scale_y_continuous(expand = expansion(mult = c(0, .1)))

### Assemble figure.
grid.arrange(isoscape, niche.widths, overlap, ncol=3, widths=c(0.5,0.25,0.2))
top <- arrangeGrob(isoscape, niche.widths, overlap, ncol=3, widths=c(0.4,0.25,0.2))

layman <- ggplotGrob(layman)
model <- ggplotGrob(model)
bottom <- cowplot::plot_grid(layman, model,
                             align = 'h', axis = 'b', rel_widths=c(0.31,0.6))

grid.arrange(top, bottom, nrow=2)

levels(plot_data$HABITAT)[levels(plot_data$HABITAT)=="urban"] <- "urban (no GPS)"

legend <- cowplot::get_legend(ggplot() +
                                geom_point(data=plot_data, aes(x=d13C, y=d15N, color = HABITAT, 
                                                                        shape=HABITAT), size = 2) +
                                stat_ellipse(data=plot_data,
                                             aes(x=d13C, y=d15N, group = HABITAT, 
                                                 color = HABITAT), size=1, 
                                             level = 0.4,
                                             alpha=0,
                                             type = "norm",
                                             geom = "polygon") + 
                                scale_color_manual(values=c("springgreen3", "goldenrod", 
                                                            "tomato1", "deepskyblue2",  "orchid1")) +
                                guides(color=guide_legend(title="Land Use", ncol=2), 
                                       linetype=guide_legend(title="Land Use", ncol=2), 
                                       shape=guide_legend(title="Land Use", ncol=2)) +
                                theme_bw() + plot_theme +
                                ylab(expression(paste(delta^{15}, "N (\u2030)")))+
                                xlab(expression(paste(delta^{13}, "C (\u2030)"))) + 
                                labs(tag="a") +
                                theme(legend.position="bottom", legend.title=element_blank()))

grid.arrange(top, bottom, legend, nrow=3, heights=c(0.4, 0.4, 0.1))

Fig2 <- arrangeGrob(top, bottom, legend, nrow=3, heights=c(0.4, 0.4, 0.1))
ggsave("publication_figures/Fig2.jpg", Fig2, width=6.5, height=6.5, units="in", dpi=300)

rm(top, bottom, legend, Fig2, overlap.data, overlap, niche.widths, layman, model, isoscape)

#### FIG. 3: SEASONAL ISOTOPIC RESULTS ####
### (a) Isoscape showing ellipses
plot_data <- data_season

levels(plot_data$Location_Season)[levels(plot_data$Location_Season)=="rural_0"] <- c("rural, summer")
levels(plot_data$Location_Season)[levels(plot_data$Location_Season)=="rural_1"] <- c("rural, winter")
levels(plot_data$Location_Season)[levels(plot_data$Location_Season)=="suburban_0"] <- c("suburban, summer")
levels(plot_data$Location_Season)[levels(plot_data$Location_Season)=="suburban_1"] <- c("suburban, winter")
levels(plot_data$Location_Season)[levels(plot_data$Location_Season)=="urban_0"] <- c("urban, summer")
levels(plot_data$Location_Season)[levels(plot_data$Location_Season)=="urban_1"] <- c("urban, winter")

isoscape <- ggplot() +
  geom_point(data=plot_data, aes(x=d13C, y=d15N, color = Location_Season), size = 1) +
  stat_ellipse(data=subset(plot_data, WINTER==0),
               aes(x=d13C, y=d15N, group = Location_Season, 
                   color = Location_Season), size=1, 
               level = 0.4, alpha=0, type = "norm", geom = "polygon", show.legend=TRUE) + 
  stat_ellipse(data=subset(plot_data, WINTER==1),
               aes(x=d13C, y=d15N, group = Location_Season, 
                   color = Location_Season), size=1, linetype="dashed",
               level = 0.4, alpha=0, type = "norm", geom = "polygon", show.legend=TRUE) + 
  theme_bw() + plot_theme +
  scale_color_manual(values=c("springgreen3", "darkgreen",
                              "goldenrod1", "darkgoldenrod3", 
                              "orchid1", "darkorchid3")) +
  ylab(expression(paste(delta^{15}, "N (\u2030)")))+
  xlab(expression(paste(delta^{13}, "C (\u2030)"))) + 
  labs(tag="a") +
  theme(legend.position="none") +
  guides(color=guide_legend(title="Group"), linetype=guide_legend(title="Group"), shape=guide_legend(title="Group"))

### (b) Niche width
levels(SIBER_season_densities$Group)[levels(SIBER_season_densities$Group)=="rur.sum"] <- "rur,sum"
levels(SIBER_season_densities$Group)[levels(SIBER_season_densities$Group)=="rur.wint"] <- "rur,wint"
levels(SIBER_season_densities$Group)[levels(SIBER_season_densities$Group)=="sub.sum"] <- "sub,sum"
levels(SIBER_season_densities$Group)[levels(SIBER_season_densities$Group)=="sub.wint"] <- "sub,wint"
levels(SIBER_season_densities$Group)[levels(SIBER_season_densities$Group)=="urb.sum"] <- "urb,sumr"
levels(SIBER_season_densities$Group)[levels(SIBER_season_densities$Group)=="urb.wint"] <- "urb,wint"

niche.widths <- ggplot() +
  geom_boxplot(data=SIBER_season_densities, aes(x=Group, ymin=Min_95, ymax=Max_95, 
                                                lower=Min_95, middle=Mode, upper=Max_95), stat="identity", width=0.4, fatten=NULL, fill="lightgrey") +
  geom_boxplot(data=SIBER_season_densities, aes(x=Group, ymin=Min_75, ymax=Max_75, 
                                                lower=Min_75, middle=Mode, upper=Max_75), stat="identity", width=0.6, fatten=NULL, fill="grey") +
  geom_boxplot(data=SIBER_season_densities, aes(x=Group, ymin=Min_50, ymax=Max_50, 
                                                lower=Min_50, middle=Mode, upper=Max_50), stat="identity", width=0.8, fatten=NULL, fill="darkgrey") +
  geom_point(data=SIBER_season_densities, aes(x=Group, y=Mode)) +
  theme_bw() + plot_theme +
  labs(x="\n Location, season", y="Standard ellipse area", tag="b") +
  theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1))

### (c) Layman metrics
layman_season <- reshape2::melt(layman_season)
colnames(layman_season) <- c("Measure", "Group", "Value")

layman_season <- tidyr::separate(layman_season, col=Group, into=c("Site","Season"), sep=",", remove=FALSE)

layman_season$Site <- factor(layman_season$Site,
                             levels=c("rural", "suburban", "urban"))
layman_season$Season <- factor(layman_season$Season,
                               levels=c("summer", "winter"))

layman_plot <- ggplot(layman_season, aes(x=Group, y=Value, fill=Group)) + geom_bar(stat="identity", position="dodge") +
  facet_wrap(~Measure, scales="free_y", nrow=3, ncol=2) +
  scale_y_continuous(expand = expansion(mult = c(0, .1))) +
  scale_fill_manual(values=c("springgreen3", "goldenrod1", "orchid1",
                             "darkgreen", "darkgoldenrod3", "darkorchid3")) +
  theme_bw() + plot_theme +
  theme(legend.position="none",
        axis.text.x=element_blank(),
        strip.text = element_blank()) +
  labs(x="\nLocation, season", y="Value", tag="c")

### (d) Mixing models
levels(mm_season_separated$Source)[levels(mm_season_separated$Source)=="Anthropogenic food"] <- "Anthropogenic\nfood"
levels(mm_season_separated$Source)[levels(mm_season_separated$Source)=="Domestic pets"] <- "Domestic\npets"
levels(mm_season_separated$Source)[levels(mm_season_separated$Source)=="Cricetid rodents "] <- "Cricetid\nrodents"
levels(mm_season_separated$Source)[levels(mm_season_separated$Source)=="Small herbivores"] <- "Small\nherbivores"

mm_season_separated$Group <- factor(mm_season_separated$Group)

model <- ggplot(mm_season_separated, aes(x=Source, y=mean, fill=Group)) +
  geom_bar(stat="identity", position="dodge") +
  geom_errorbar(aes(x=Source, ymin=mean-sd, ymax=mean+sd), position=position_dodge(width=0.9), color="black", width=0.3) +
  theme_bw() + plot_theme +
  scale_fill_manual(values=c("springgreen3", "darkgreen",
                              "goldenrod1", "darkgoldenrod3", 
                              "orchid1", "darkorchid3")) +
  facet_grid(~Season) +
  theme(legend.position="none",
        axis.text.x=element_text(angle=45, hjust=1, vjust=1)) +
  labs(x="Source", y="Proportion", tag="d") +
  scale_y_continuous(expand = expansion(mult = c(0, .1)))

### Assembly
grid.arrange(isoscape, niche.widths, layman_plot,
             model,
             layout_matrix=rbind(c(1,2,3),c(4,4,4)))

legend <- cowplot::get_legend(ggplot() +
                                geom_point(data=plot_data, aes(x=d13C, y=d15N, color = Location_Season), size = 2) +
                                stat_ellipse(data=subset(plot_data, WINTER==0),
                                             aes(x=d13C, y=d15N, group = Location_Season, 
                                                 color = Location_Season), size=1, 
                                             level = 0.4, alpha=0, type = "norm", geom = "polygon", show.legend=TRUE) + 
                                stat_ellipse(data=subset(plot_data, WINTER==1),
                                             aes(x=d13C, y=d15N, group = Location_Season, 
                                                 color = Location_Season), size=1, linetype="dashed",
                                             level = 0.4, alpha=0, type = "norm", geom = "polygon", show.legend=TRUE) + 
                                theme_bw() + plot_theme +
                                scale_color_manual(values=c("springgreen3", "darkgreen",
                                                            "goldenrod1", "darkgoldenrod3", 
                                                            "orchid1", "darkorchid3")) +
                                theme(legend.position="bottom", legend.title=element_blank()))

grid.arrange(isoscape, niche.widths, layman_plot,
             model,
             legend,
             layout_matrix=rbind(c(1,2,3),c(4,4,4),c(5,5,5)),
             heights=c(0.4, 0.4, 0.1))

top <- arrangeGrob(isoscape, niche.widths, layman_plot,
                   widths=c(0.4,0.3,0.4), ncol=3)

grid.arrange(top, model, legend, nrow=3, heights=c(0.4,0.4,0.1))

Fig3 <- arrangeGrob(top, model, legend, nrow=3, heights=c(0.4,0.4,0.1))

ggsave("publication_figures/Fig3.jpg", Fig3, width=6.5, height=6.5, units="in", dpi=300)

rm(Fig3, top, isoscape, niche.widths, layman_plot, model, legend, plot_data)

#### FIG. 4: MODELS PREDICTING ISOTOPE VALUES ####
model_d13C_predict$merge[model_d13C_predict$merge=="HABITATsuburban"] <- "Habitat use\nsuburban"
model_d13C_predict$merge[model_d13C_predict$merge=="HABITATgreenspace"] <- "Habitat use\ngreenspace"
model_d13C_predict$merge[model_d13C_predict$merge=="HABITATmatrix"] <- "Habitat use\nmatrix"
model_d13C_predict$merge[model_d13C_predict$merge=="ADULT"] <- "Adult (Y)"
model_d13C_predict$merge[model_d13C_predict$merge=="WINTER"] <- "Winter (Y)"

model_d13C_predict$merge <- factor(model_d13C_predict$merge,
                                   levels=c("Habitat use\nsuburban", "Habitat use\ngreenspace", "Habitat use\nmatrix",
                                            "Winter (Y)"))

plot1 <- ggplot(model_d13C_predict) +
  ggstance::geom_linerangeh(aes(y=merge, xmin=psd_25, xmax=psd_75),
                            size=2, color="black") + 
  ggstance::geom_linerangeh(aes(y=merge, xmin=psd_2.5, xmax=psd_97.5),
                            color="black", size=0.5) +
  geom_point(aes(y=merge, x=psd_Coef),
             fill="white", shape=21, stroke=1, show.legend=TRUE) +
  coord_flip() + theme_bw() + plot_theme +
  labs(x="Coefficient", y="Predictor", tag="a", title=expression({delta}^13*C)) +
  geom_vline(xintercept=0, linetype="dashed", color="black") +
  theme(axis.text.x=element_text(angle=45, vjust=1, hjust=1),
        panel.grid.major.y=element_line(color="grey", size=0.25))


model_d15N_predict$merge[model_d15N_predict$merge=="HABITATsuburban"] <- "Habitat use\nsuburban"
model_d15N_predict$merge[model_d15N_predict$merge=="HABITATgreenspace"] <- "Habitat use\ngreenspace"
model_d15N_predict$merge[model_d15N_predict$merge=="HABITATmatrix"] <- "Habitat use\nmatrix"
model_d15N_predict$merge[model_d15N_predict$merge=="ADULT"] <- "Adult (Y)"
model_d15N_predict$merge[model_d15N_predict$merge=="WINTER"] <- "Winter (Y)"

model_d15N_predict$merge <- factor(model_d15N_predict$merge,
                                   levels=c("Habitat use\nsuburban", "Habitat use\ngreenspace", "Habitat use\nmatrix",
                                            "Adult (Y)", "Winter (Y)"))

plot2 <- ggplot(model_d15N_predict) +
  ggstance::geom_linerangeh(aes(y=merge, xmin=psd_25, xmax=psd_75),
                            size=2, color="black") + 
  ggstance::geom_linerangeh(aes(y=merge, xmin=psd_2.5, xmax=psd_97.5),
                            color="black", size=0.5) +
  geom_point(aes(y=merge, x=psd_Coef),
             fill="white", shape=21, stroke=1, show.legend=TRUE) +
  coord_flip() + theme_bw() + plot_theme + 
  labs(x="Coefficient", y="Predictor", tag="b", title=expression({delta}^15*N)) +
  geom_vline(xintercept=0, linetype="dashed", color="black") +
  theme(axis.text.x=element_text(angle=45, vjust=1, hjust=1),
        panel.grid.major.y=element_line(color="grey", size=0.25))

grid.arrange(plot1, plot2, ncol=2)

models.plot <- arrangeGrob(plot1, plot2, ncol=2)
ggsave("models.plot.jpg", models.plot, width=6, height=4, units="in", dpi=300)

#### FIG. 5: DIFFERENT MIXING MODEL FORMULATIONS ####
sensitivity_df_melt <- subset(sensitivity_df_melt, 
                              !variable %in% c("grp_ind_0.5_3.5", "grp_dep_0.5_3.5", "grp_dig_0.5_3.5"))

sensitivity_df_melt$Location <- factor(sensitivity_df_melt$Location,
                                       levels=c("rural", "suburban", "greenspace", "matrix", "urban"))

levels(sensitivity_df_melt$Location)[levels(sensitivity_df_melt$Location)=="urban"] <- "urban\n(no GPS)"

levels(sensitivity_df_melt$Source)[levels(sensitivity_df_melt$Source)=="Anthropogenic food"] <- "Anthropogenic\nfood"
levels(sensitivity_df_melt$Source)[levels(sensitivity_df_melt$Source)=="Domestic pets"] <- "Domestic\npets"
levels(sensitivity_df_melt$Source)[levels(sensitivity_df_melt$Source)=="Cricetid rodents "] <- "Cricetid\nrodents"
levels(sensitivity_df_melt$Source)[levels(sensitivity_df_melt$Source)=="Small herbivores"] <- "Small\nherbivores"

sensitivity_df_melt$Source <- factor(sensitivity_df_melt$Source,
                                     levels=c("Anthropogenic\nfood", "Domestic\npets", "Insects", "Cricetid\nrodents",
                                              "Small\nherbivores", "Ungulates", "Berries"))

sensitivity_df_melt <- tidyr::separate(sensitivity_df_melt, col=variable,
                                    into=c("grp","form","c1", "c2"), sep="_", remove=FALSE)
sensitivity_df_melt$c1[sensitivity_df_melt$c1==4.25] <- 4.3
sensitivity_df_melt$c2[sensitivity_df_melt$c2==3.05] <- 3.1
sensitivity_df_melt$discr <- paste0(sensitivity_df_melt$c1, " / ", sensitivity_df_melt$c2)
sensitivity_df_melt$form <- factor(sensitivity_df_melt$form, levels=c("ind", "dep", "dig"))
levels(sensitivity_df_melt$form)[levels(sensitivity_df_melt$form)=="ind"] <- "concentration-independent"
levels(sensitivity_df_melt$form)[levels(sensitivity_df_melt$form)=="dep"] <- "concentration-dependent"
levels(sensitivity_df_melt$form)[levels(sensitivity_df_melt$form)=="dig"] <- "concentration-dependent,\nadjusted for digestibility"



plot.b <- ggplot(sensitivity_df_melt, aes(x=Location, y=value, color=discr, shape=form)) + 
  geom_point(size=2) + 
  #coord_flip() + 
  facet_wrap(~Source, scales="free_y") +
  #scale_x_discrete(limits=rev(levels(sensitivity_df_melt$Location))) +
  #scale_color_manual(values=c("springgreen3", "tomato1", "goldenrod", "deepskyblue2", "orchid1")) +
  scale_shape_manual(values=c(16, 6, 17)) +
  guides(color=guide_legend(title.position="top", title.hjust=0.5, title="Discrimination\nfactors", order=1),
         shape=guide_legend(title.position="top", title.hjust=0.5, title="Model type", order=2)) +
  theme_bw() + plot_theme +
  theme(panel.grid.major.y=element_line(color="lightgrey", size=0.25),
        axis.text.x=element_text(angle=45, vjust=1, hjust=1),
        legend.position=c(1,0),
        legend.box="horizontal",
        legend.justification=c(1.2,0.2)) +
  labs(y="Dietary proportion (%)", x="Habitat use")

ggsave("publication_figures/Fig5_detailed.jpg", plot.b, width=6.5, height=7, units="in", dpi=300)

ggsave("publication_figures/spearman_matrix.jpg", plot.a, width=4, height=4, units="in", dpi=300)

grid.arrange(plot.a + labs(tag="a"), plot.b, ncol=2)
Fig4 <- arrangeGrob(plot.a, plot.b, ncol=2)
ggsave("publication_figures/Fig4_DRAFT.jpg", Fig4, width=8, height=5, units="in", dpi=300)

sensitivity_df_melt$discr[sensitivity_df_melt$discr=="2.0_3.0"] <- "2.0, 3.0"
sensitivity_df_melt$discr[sensitivity_df_melt$discr=="2.0_3.5"] <- "2.0, 3.5"
sensitivity_df_melt$discr[sensitivity_df_melt$discr=="2.5_3.5"] <- "2.5, 3.5"
sensitivity_df_melt$discr[sensitivity_df_melt$discr=="3.4_3.0"] <- "3.4, 3.0"
sensitivity_df_melt$discr[sensitivity_df_melt$discr=="4.25_3.05"] <- "4.25, 3.05"

legend <- ggplot(sensitivity_df_melt, aes(x=Location, y=value)) + 
  geom_point(aes(shape=discr), size=2) + 
  coord_flip() + 
  facet_wrap(~Source) +
  scale_x_discrete(limits=rev(levels(sensitivity_df_melt$Location))) +
  scale_color_manual(values=c("springgreen3", "tomato1", "goldenrod", "deepskyblue2", "orchid1")) +
  scale_shape_manual(values=c(0, 16, 2, 18, 13),
                     guide=guide_legend(title="Discrimination\nfactors", title.hjust=0.5)) +
  theme_bw() + plot_theme +
  theme(panel.grid.major.x=element_line(color="lightgrey", size=0.25),
        legend.position="right") +
  labs(y="Dietary proportion (%)", x="Source", tag="b")

ggsave("publication_figures/Fig4_LEGEND.jpg", legend, width=5, height=5, units="in", dpi=300)

temp <- sensitivity_corr
temp <- subset(temp, R < 1)
temp <- tidyr::separate(temp, col=Model1,
                        into=c("grp","m1_form","m1_d13c", "m1_d15N"), sep="_", remove=FALSE)
temp$grp <- NULL
temp <- tidyr::separate(temp, col=Model2,
                        into=c("grp","m2_form","m2_d13c", "m2_d15N"), sep="_", remove=FALSE)
temp$grp <- NULL
temp$m1_discr <- paste0(temp$m1_d13c, temp$m1_d15N)
temp$m2_discr <- paste0(temp$m2_d13c, temp$m2_d15N)

temp$same_discr <- temp$m1_discr == temp$m2_discr
temp$same_conc <- temp$m1_form == temp$m2_form

temp1 <- subset(temp, temp$same_discr=="TRUE")
temp2 <- subset(temp, temp$same_conc=="TRUE")

mean(temp1$R)
sd(temp1$R)

mean(temp2$R)
sd(temp2$R)

#### FIG. S1: HAIR/CLAW ISOTOPE COMPARISONS ####
data_hair_claw$CoyoteID <- as.factor(data_hair_claw$CoyoteID)

temp1 <- subset(data_hair_claw, Type=="claw")
colnames(temp1) <- c("CoyoteID", "Type", "SampleID", "d13C_Claw", "d15N_Claw")
temp1$SampleID <- NULL
temp1$Type <- NULL

temp2 <- subset(data_hair_claw, Type=="hair")
colnames(temp2) <- c("CoyoteID", "Type", "SampleID", "d13C_Hair", "d15N_Hair")
temp2$SampleID <- NULL
temp2$Type <- NULL

data_hair_claw_plot <- merge(temp1, temp2, by="CoyoteID", all=TRUE)
data_hair_claw_plot$d13C <- rep("d13C")
data_hair_claw_plot$d15N <- rep("d15N")

plot1 <- ggplot(data_hair_claw_plot, aes(x=d13C_Claw, y=d13C_Hair)) +
  geom_point() +
  geom_smooth(method="lm") +
  theme_bw() + plot_theme +
  facet_wrap(~d13C, labeller = label_bquote(cols = {delta}^13*C)) +
  theme(axis.text=element_text(color="black", size=10),
        axis.title=element_text(color="black", size=11),
        plot.title=element_text(color="black", size=13)) +
  labs(x=expression({delta}^13*C~'\u2030'~" (claw)"),
       y=expression({delta}^13*C~'\u2030'~" (hair)"))

plot2 <- ggplot(data_hair_claw_plot, aes(x=d15N_Claw, y=d15N_Hair)) +
  geom_point() +
  geom_smooth(method="lm") +
  theme_bw() + plot_theme +
  facet_wrap(~d15N, labeller = label_bquote(cols = {delta}^15*N)) +
  scale_y_continuous(breaks=c(8,8.5,9,9.5,10,10.5,11)) +
  scale_x_continuous(breaks=c(8.5, 9.0, 9.5, 10.0)) +
  theme(axis.text=element_text(color="black", size=10),
        axis.title=element_text(color="black", size=11),
        plot.title=element_text(color="black", size=13)) +
  labs(x=expression({delta}^15*N~'\u2030'~" (claw)"),
       y=expression({delta}^15*N~'\u2030'~" (hair)"))

grid.arrange(plot1, plot2, nrow=1)
hair_claw <- arrangeGrob(plot1, plot2, nrow=1)

ggsave("publication_figures/hair_claw_regressions.jpg", hair_claw, width=6, height=3, units="in", dpi=300)

rm(plot1, plot2, hair_claw, data_hair_claw_plot, temp1, temp2)

#### FIG. S2: LAND-USE ISOTOPIC RESULTS WITHOUT MANGE ####
### (a) Isoscape showing ellipses
plot_data <- subset(data_single_point, MANGE==0 | is.na(MANGE)=="TRUE")[,c("d13C","d15N","LOCATION","HABITAT")]
plot_data$HABITAT <- factor(plot_data$HABITAT,
                             levels=c("rural", "suburban", "greenspace", "matrix", "urban"))

isoscape <- ggplot() +
  # Add the coyotes
  geom_point(data=plot_data, aes(x=d13C, y=d15N, color = HABITAT, shape=HABITAT), size = 2) +
  stat_ellipse(data=plot_data,
               aes(x=d13C, y=d15N, group = HABITAT, 
                   color = HABITAT), size=1, 
               level = 0.4, alpha=0,
               type = "norm", geom = "polygon") + 
  scale_color_manual(values=c("springgreen3", "goldenrod", "tomato1", "deepskyblue2", "orchid1")) +
  
  theme_bw() + plot_theme +
  ylab(expression(paste(delta^{15}, "N (\u2030)")))+
  xlab(expression(paste(delta^{13}, "C (\u2030)"))) + 
  labs(tag="a") +
  theme(legend.position="none")

### (b) Niche width
niche.widths <- ggplot() +
  geom_boxplot(data=SIBER_no.mange_densities, 
               aes(x=Landuse, ymin=Min_95, ymax=Max_95, 
                   lower=Min_95, middle=Mode, upper=Max_95),
               stat="identity", width=0.4, fatten=NULL, fill="lightgrey") +
  geom_boxplot(data=SIBER_no.mange_densities,
               aes(x=Landuse, ymin=Min_75, ymax=Max_75, 
                   lower=Min_75, middle=Mode, upper=Max_75),
               stat="identity", width=0.6, fatten=NULL, fill="grey") +
  geom_boxplot(data=SIBER_no.mange_densities, 
               aes(x=Landuse, ymin=Min_50, ymax=Max_50, 
                   lower=Min_50, middle=Mode, upper=Max_50),
               stat="identity", width=0.8, fatten=NULL, fill="darkgrey") +
  geom_point(data=SIBER_no.mange_densities, aes(x=Landuse, y=Mode)) +
  theme_bw() + plot_theme +
  labs(x="\n Habitat use", y="Standard ellipse area", tag="b") +
  theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1))

### (c) Niche overlap
write.csv(SIBER_no.mange_overlap, "overlap.csv")
overlap.data <- read.csv("overlap_no_mange.csv")
overlap.data$axis <- factor(overlap.data$axis, levels=c("rural", "suburban", "greenspace", "matrix"))
overlap.data$fill <- factor(overlap.data$fill, levels=c("rural", "suburban", "greenspace", "matrix"))

overlap <- ggplot(overlap.data, aes(x=axis, y=overlap_40, fill=fill)) +
  geom_bar(stat="identity") +
  scale_fill_manual(values=c("springgreen3", "goldenrod", "tomato1", "deepskyblue2")) +
  theme_bw() + plot_theme +
  theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1),
        legend.position="none") +
  labs(x="\n Habitat use", y="Ellipse overlap", tag="c") +
  scale_y_continuous(expand = expansion(mult = c(0, .1)))

### (d) Layman metrics
layman_no.mange <- reshape2::melt(layman_no.mange, id.vars="Measure")
colnames(layman_no.mange) <- c("Measure", "Land.use", "Value")
layman_no.mange$Land.use <- factor(layman_no.mange$Land.use,
                                  levels=c("rural", "suburban", "greenspace", "matrix", "urban"))

layman <- ggplot(layman_no.mange, aes(x=Land.use, y=Value, fill=Land.use)) + geom_bar(stat="identity", position="dodge") +
  facet_wrap(~Measure, scales="free_y", nrow=3, ncol=2) +
  theme_bw() + plot_theme +
  scale_fill_manual(values=c("springgreen3", "goldenrod", "tomato1", "deepskyblue2",  "orchid1")) +
  scale_y_continuous(expand = expansion(mult = c(0, .1))) +
  theme(legend.position="none",
        axis.text.x=element_text(angle=45, hjust=1, vjust=1),
        strip.text = element_blank()) +
  labs(x="Habitat use", y="Value", tag="d")

### (e) Mixing model
levels(mm_no.mange_main$Source)[levels(mm_no.mange_main$Source)=="Anthropogenic food"] <- "Anthropogenic\nfood"
levels(mm_no.mange_main$Source)[levels(mm_no.mange_main$Source)=="Domestic pets"] <- "Domestic\npets"
levels(mm_no.mange_main$Source)[levels(mm_no.mange_main$Source)=="Cricetid rodents "] <- "Cricetid\nrodents"
levels(mm_no.mange_main$Source)[levels(mm_no.mange_main$Source)=="Small herbivores"] <- "Small\nherbivores"

model <- ggplot(mm_no.mange_main, aes(x=Source, y=mean, fill=Site)) +
  geom_bar(stat="identity", position="dodge") +
  geom_errorbar(aes(x=Source, ymin=mean-sd, ymax=mean+sd), position=position_dodge(width=0.9), color="black", width=0.3) +
  theme_bw() + plot_theme +
  scale_fill_manual(values=c("springgreen3", "goldenrod", "tomato1", "deepskyblue2",  "orchid1")) +
  theme(legend.position="none",
        axis.text.x=element_text(angle=45, hjust=1, vjust=1)) +
  labs(x="Source", y="Proportion", tag="e") +
  scale_y_continuous(expand = expansion(mult = c(0, .1)))

### Assembly
grid.arrange(isoscape, niche.widths, overlap, ncol=3, widths=c(0.5,0.25,0.2))
top <- arrangeGrob(isoscape, niche.widths, overlap, ncol=3, widths=c(0.4,0.25,0.2))

layman <- ggplotGrob(layman)
model <- ggplotGrob(model)
bottom <- cowplot::plot_grid(layman, model,
                             align = 'h', axis = 'b', rel_widths=c(0.31,0.6))

grid.arrange(top, bottom, nrow=2)

levels(plot_data$HABITAT)[levels(plot_data$HABITAT)=="urban"] <- "urban (no GPS)"

legend <- cowplot::get_legend(ggplot() +
                                geom_point(data=plot_data, aes(x=d13C, y=d15N, color = HABITAT, 
                                                               shape=HABITAT), size = 2) +
                                stat_ellipse(data=plot_data,
                                             aes(x=d13C, y=d15N, group = HABITAT, 
                                                 color = HABITAT), size=1, 
                                             level = 0.4,
                                             alpha=0,
                                             type = "norm",
                                             geom = "polygon") + 
                                scale_color_manual(values=c("springgreen3", "goldenrod", 
                                                            "tomato1", "deepskyblue2",  "orchid1")) +
                                guides(color=guide_legend(title="Land Use", ncol=2), 
                                       linetype=guide_legend(title="Land Use", ncol=2), 
                                       shape=guide_legend(title="Land Use", ncol=2)) +
                                theme_bw() + plot_theme +
                                ylab(expression(paste(delta^{15}, "N (\u2030)")))+
                                xlab(expression(paste(delta^{13}, "C (\u2030)"))) + 
                                labs(tag="a") +
                                theme(legend.position="bottom", legend.title=element_blank()))

grid.arrange(top, bottom, legend, nrow=3, heights=c(0.4, 0.4, 0.1))

FigS2 <- arrangeGrob(top, bottom, legend, nrow=3, heights=c(0.4, 0.4, 0.1))
ggsave("publication_figures/FigS2_no_mange.jpg", FigS2, width=6.5, height=6.5, units="in", dpi=300)

rm(top, bottom, legend, FigS2, overlap.data, overlap, niche.widths, layman, model, isoscape)

#### FIG. S3: NICHE OVERLAP AMONG SEASONS ####
overlap.season <- SIBER_season_overlap
overlap.season$Group <- rownames(overlap.season)
overlap.season <- tidyr::separate(overlap.season, col=Group, into=c("constant","comparisons"), sep="_", remove=TRUE)
overlap.season <- tidyr::separate(overlap.season, col=comparisons,
                                  into=c("axis","fill"), sep=".v.", remove=TRUE)

# Suburban/rural comparison needs to get reproduced.
df <- subset(overlap.season, axis=="Rur" & fill=="Sub")
df$axis <- rep("Sub")
df$fill <- rep("Rur")
rownames(df) <- c(1, 2)
overlap.season <- rbind(overlap.season, df)

df <- subset(overlap.season, fill=="Urb")
df$fill <- df$axis
df$axis <- rep("Urb")
rownames(df) <- c(1:4)
overlap.season <- rbind(overlap.season, df)

overlap.season$facet <- rep("habitat")
for(i in 1:nrow(overlap.season)){
  if(overlap.season[i,"constant"] == "Sum"){
    overlap.season[i,"facet"] <- "summer"
  }
  if(overlap.season[i,"constant"] == "Wint"){
    overlap.season[i,"facet"] <- "winter"
  }
}

df <- subset(overlap.season, facet=="habitat")
df$axis <- c("rural", "suburban", "urban")
df$fill <- c("neutral", "neutral", "neutral")

overlap.season <- subset(overlap.season, facet !="habitat")
overlap.season <- rbind(overlap.season, df)
overlap.season$facet <- factor(overlap.season$facet)
overlap.season$axis[overlap.season$axis=="Rur"] <- "rural"
overlap.season$axis[overlap.season$axis=="Sub"] <- "suburban"
overlap.season$axis[overlap.season$axis=="Urb"] <- "urban"
overlap.season$fill[overlap.season$fill=="Rur"] <- "rural"
overlap.season$fill[overlap.season$fill=="Sub"] <- "suburban"
overlap.season$fill[overlap.season$fill=="Urb"] <- "urban"


facet1 <- subset(overlap.season, facet=="habitat")
levels(facet1$facet)[levels(facet1$facet)=="habitat"] <- "habitats:\nbetween seasons"

facet2 <- subset(overlap.season, facet=="summer")
levels(facet2$facet)[levels(facet2$facet)=="summer"] <- "summer:\nbetween habitats"

facet3 <- subset(overlap.season, facet=="winter")
levels(facet3$facet)[levels(facet3$facet)=="winter"] <- "winter:\nbetween habitats"

plot1 <- ggplot(facet1, aes(x=axis, y=overlap_40, fill=fill)) +
  geom_bar(stat="identity") +
  theme_bw() + plot_theme +
  facet_wrap(~facet) +
  scale_y_continuous(expand = expansion(mult = c(0, .1)), limits=c(0,40)) +
  theme(legend.position="none",
        axis.text.x=element_text(angle=45, hjust=1, vjust=1)) +
  scale_fill_manual(values=c("grey")) +
  labs(x="\nHabitat use", y="Percent overlap")

plot2 <- ggplot(facet2, aes(x=axis, y=overlap_40, fill=fill)) +
  geom_bar(stat="identity") +
  theme_bw() + plot_theme +
  facet_wrap(~facet) +
  scale_y_continuous(expand = expansion(mult = c(0, .1)), limits=c(0,40)) +
  theme(legend.position="none",
        axis.text.x=element_text(angle=45, hjust=1, vjust=1)) +
  scale_fill_manual(values=c("springgreen3", "goldenrod1", "orchid1"))+
  labs(x="\nHabitat use", y="Percent overlap")

plot3 <- ggplot(facet3, aes(x=axis, y=overlap_40, fill=fill)) +
  geom_bar(stat="identity") +
  theme_bw() + plot_theme +
  facet_wrap(~facet) +
  scale_y_continuous(expand = expansion(mult = c(0, .1)), limits=c(0,40))  +
  theme(legend.position="none",
        axis.text.x=element_text(angle=45, hjust=1, vjust=1)) +
  scale_fill_manual(values=c("darkgreen", "darkgoldenrod3", "darkorchid3")) +
  labs(x="\nHabitat use", y="Percent overlap")

grid.arrange(plot1, plot2, plot3, ncol=3)

df <- data.frame(
  group = c("rural,summer", "rural,winter", "suburban,summer","suburban,winter","urban,summer","urban,winter"),
  values = c(1,2,3,4,5,6)
)

legend <- cowplot::get_legend(ggplot(df, aes(x=group, y=values, fill=group)) +
                                geom_bar(stat="identity") +
                                guides(guide_legend(ncol=3)) +
                                theme_bw() + plot_theme +
                                theme(legend.position="bottom",
                                      legend.title=element_blank()) +
                                scale_fill_manual(values=c("springgreen3", "darkgreen", 
                                                           "goldenrod1", "darkgoldenrod3", 
                                                           "orchid1", "darkorchid3"),
                                                  guide=guide_legend(ncol=3)) +
                                labs(x="Habitat use", y="Percent overlap"))

grid.arrange(plot1, plot2, plot3, legend, layout_matrix=rbind(c(1,2,3),c(4,4,4)), heights=c(0.8,0.1))

season.overlap <- arrangeGrob(plot1, plot2, plot3, legend, 
                              layout_matrix=rbind(c(1,2,3),c(4,4,4)), heights=c(0.8,0.1))

ggsave("publication_figures/supp.season.overlap.jpg", season.overlap, width=5.5, height=5, units="in", dpi=300)

#### FIG. S4: RANDOM FOREST MODEL PREDICTED DATA #### 
SIBER_predicted_plot$HABITAT <- factor(SIBER_predicted_plot$HABITAT,
                                        levels=c("rural", "suburban", "greenspace", "matrix")) 

### (a) Isoscape with ellipses
isoscape <- ggplot() +
  geom_point(data=SIBER_predicted_plot, aes(x=d13C, y=d15N, color = HABITAT, shape=HABITAT), size = 2) +
  stat_ellipse(data=SIBER_predicted_plot,
               aes(x=d13C, y=d15N, group = HABITAT, 
                   color = HABITAT), size=1, 
               level = 0.4,
               alpha=0,
               type = "norm",
               geom = "polygon") + 
  scale_color_manual(values=c("springgreen3", "goldenrod", "tomato1", "deepskyblue2")) +
  
  theme_bw() + plot_theme +
  ylab(expression(paste(delta^{15}, "N (\u2030)")))+
  xlab(expression(paste(delta^{13}, "C (\u2030)"))) + 
  labs(tag="a") +
  theme(legend.position="none")

### (b) Niche widths
niche.widths <- ggplot() +
  geom_boxplot(data=SIBER_predicted_densities, aes(x=Landuse, ymin=Min_95, ymax=Max_95, 
                                         lower=Min_95, middle=Mode, upper=Max_95), stat="identity", width=0.4, fatten=NULL, fill="lightgrey") +
  geom_boxplot(data=SIBER_predicted_densities, aes(x=Landuse, ymin=Min_75, ymax=Max_75, 
                                         lower=Min_75, middle=Mode, upper=Max_75), stat="identity", width=0.6, fatten=NULL, fill="grey") +
  geom_boxplot(data=SIBER_predicted_densities, aes(x=Landuse, ymin=Min_50, ymax=Max_50, 
                                         lower=Min_50, middle=Mode, upper=Max_50), stat="identity", width=0.8, fatten=NULL, fill="darkgrey") +
  geom_point(data=SIBER_predicted_densities, aes(x=Landuse, y=Mode)) +
  theme_bw() + plot_theme +
  labs(x="\n Habitat use", y="Standard ellipse area", tag="b") +
  theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1))

### (c) Niche overlap
write.csv(SIBER_predicted_overlap, "predict.overlap.csv")
SIBER_predicted_overlap <- read.csv("predict.overlap2.csv")

SIBER_predicted_overlap$axis <- factor(SIBER_predicted_overlap$axis, levels=c("rural", "suburban", "greenspace", "matrix"))
SIBER_predicted_overlap$fill <- factor(SIBER_predicted_overlap$fill, levels=c("rural", "suburban", "greenspace", "matrix"))

overlap <- ggplot(SIBER_predicted_overlap, aes(x=axis, y=overlap_40, fill=fill)) +
  geom_bar(stat="identity") +
  scale_fill_manual(values=c("springgreen3", "goldenrod", "tomato1", "deepskyblue2")) +
  theme_bw() + plot_theme +
  theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1),
        legend.position="none") +
  labs(x="\n Habitat use", y="Ellipse overlap", tag="c") +
  
  scale_y_continuous(expand = expansion(mult = c(0, .1)))

### (d) Layman metrics
layman_predicted <- reshape2::melt(layman_predicted, id.vars="Measure")
colnames(layman_predicted) <- c("Measure", "Land.use", "Value")
layman_predicted$Land.use <- factor(layman_predicted$Land.use,
                                  levels=c("rural", "suburban", "greenspace", "matrix"))

layman <- ggplot(layman_predicted, aes(x=Land.use, y=Value, fill=Land.use)) +
  geom_bar(stat="identity", position="dodge") +
  facet_wrap(~Measure, scales="free_y", nrow=3, ncol=2) +
  theme_bw() + plot_theme +
  scale_fill_manual(values=c("springgreen3", "goldenrod", "tomato1", "deepskyblue2")) +
  scale_y_continuous(expand = expansion(mult = c(0, .1))) +
  theme(legend.position="none",
        axis.text.x=element_text(angle=45, hjust=1, vjust=1),
        strip.text = element_blank()) +
  labs(x="Habitat use", y="Value", tag="d")

#### (e) Mising models
levels(mm_predicted_main$Source)[levels(mm_predicted_main$Source)=="Anthropogenic food"] <- "Anthropogenic\nfood"
levels(mm_predicted_main$Source)[levels(mm_predicted_main$Source)=="Domestic pets"] <- "Domestic\npets"
levels(mm_predicted_main$Source)[levels(mm_predicted_main$Source)=="Cricetid rodents "] <- "Cricetid\nrodents"
levels(mm_predicted_main$Source)[levels(mm_predicted_main$Source)=="Small herbivores"] <- "Small\nherbivores"

model <- ggplot(mm_predicted_main, aes(x=Source, y=mean, fill=Site)) +
  geom_bar(stat="identity", position="dodge") +
  geom_errorbar(aes(x=Source, ymin=mean-sd, ymax=mean+sd), position=position_dodge(width=0.9), color="black", width=0.3) +
  theme_bw() + plot_theme +
  scale_fill_manual(values=c("springgreen3", "goldenrod", "tomato1", "deepskyblue2")) +
  theme(legend.position="none",
        axis.text.x=element_text(angle=45, hjust=1, vjust=1)) +
  labs(x="Source", y="Proportion", tag="e") +
  scale_y_continuous(expand = expansion(mult = c(0, .1)))

### Assembly 
grid.arrange(isoscape, niche.widths, overlap, ncol=3, widths=c(0.5,0.25,0.2))
top <- arrangeGrob(isoscape, niche.widths, overlap, ncol=3, widths=c(0.4,0.25,0.2))

layman <- ggplotGrob(layman)
model <- ggplotGrob(model)
bottom <- cowplot::plot_grid(layman, model,
                             align = 'h', axis = 'b', rel_widths=c(0.31,0.6))

grid.arrange(top, bottom, nrow=2)

legend <- cowplot::get_legend(ggplot() +
                                geom_point(data=data_predicted, aes(x=d13C, y=d15N, color = HABITAT, 
                                                                       shape=HABITAT), size = 2) +
                                stat_ellipse(data=data_predicted,
                                             aes(x=d13C, y=d15N, group = HABITAT, 
                                                 color = HABITAT), size=1, 
                                             level = 0.4,
                                             alpha=0,
                                             type = "norm",
                                             geom = "polygon") + 
                                scale_color_manual(values=c("springgreen3", "goldenrod", 
                                                            "tomato1", "deepskyblue2")) +
                                guides(color=guide_legend(title="Land Use", ncol=2), 
                                       linetype=guide_legend(title="Land Use", ncol=2), 
                                       shape=guide_legend(title="Land Use", ncol=2)) +
                                theme_bw() + plot_theme +
                                ylab(expression(paste(delta^{15}, "N (\u2030)")))+
                                xlab(expression(paste(delta^{13}, "C (\u2030)"))) + 
                                labs(tag="a") +
                                theme(legend.position="bottom", legend.title=element_blank()))

grid.arrange(top, bottom, legend, nrow=3, heights=c(0.4, 0.4, 0.1))

FigS <- arrangeGrob(top, bottom, legend, nrow=3, heights=c(0.4, 0.4, 0.1))
ggsave("publication_figures/Fig_random_forest.jpg", FigS, width=6.5, height=6.5, units="in", dpi=300)

rm(top, bottom, legend, FigS, overlap, niche.widths, layman, model, isoscape)

#### FIG. S5: MODELS PREDICTING ISOTOPE VALUES (coyotes w/o GPS data excluded) ####
model_d13C_raw$merge[model_d13C_raw$merge=="HABITATsuburban"] <- "Habitat use\nsuburban"
model_d13C_raw$merge[model_d13C_raw$merge=="HABITATgreenspace"] <- "Habitat use\ngreenspace"
model_d13C_raw$merge[model_d13C_raw$merge=="HABITATmatrix"] <- "Habitat use\nmatrix"
model_d13C_raw$merge[model_d13C_raw$merge=="ADULT"] <- "Adult (Y)"
model_d13C_raw$merge[model_d13C_raw$merge=="WINTER"] <- "Winter (Y)"

model_d13C_raw$merge <- factor(model_d13C_raw$merge,
                               levels=c("Habitat use\nsuburban", "Habitat use\ngreenspace", "Habitat use\nmatrix",
                                        "Adult (Y)", "Winter (Y)"))

plot1 <- ggplot(model_d13C_raw) +
  ggstance::geom_linerangeh(aes(y=merge, xmin=psd_25, xmax=psd_75),
                            size=2, color="black") + 
  ggstance::geom_linerangeh(aes(y=merge, xmin=psd_2.5, xmax=psd_97.5),
                            color="black", size=0.5) +
  geom_point(aes(y=merge, x=psd_Coef),
             fill="white", shape=21, stroke=1, show.legend=TRUE) +
  coord_flip() + theme_bw() + plot_theme +
  labs(x="Coefficient", y="Predictor", tag="a", title=expression({delta}^13*C)) +
  geom_vline(xintercept=0, linetype="dashed", color="black") +
  theme(axis.text.x=element_text(angle=45, vjust=1, hjust=1),
        panel.grid.major.y=element_line(color="grey", size=0.25))


model_d15N_raw$merge[model_d15N_raw$merge=="HABITATsuburban"] <- "Habitat use\nsuburban"
model_d15N_raw$merge[model_d15N_raw$merge=="HABITATgreenspace"] <- "Habitat use\ngreenspace"
model_d15N_raw$merge[model_d15N_raw$merge=="HABITATmatrix"] <- "Habitat use\nmatrix"
model_d15N_raw$merge[model_d15N_raw$merge=="WINTER"] <- "Winter (Y)"

model_d15N_raw$merge <- factor(model_d15N_raw$merge,
                               levels=c("Habitat use\nsuburban", "Habitat use\ngreenspace", "Habitat use\nmatrix",
                                        "Winter (Y)"))

plot2 <- ggplot(model_d15N_raw) +
  ggstance::geom_linerangeh(aes(y=merge, xmin=psd_25, xmax=psd_75),
                            size=2, color="black") + 
  ggstance::geom_linerangeh(aes(y=merge, xmin=psd_2.5, xmax=psd_97.5),
                            color="black", size=0.5) +
  geom_point(aes(y=merge, x=psd_Coef),
             fill="white", shape=21, stroke=1, show.legend=TRUE) +
  coord_flip() + theme_bw() + plot_theme + 
  labs(x="Coefficient", y="Predictor", tag="b", title=expression({delta}^15*N)) +
  geom_vline(xintercept=0, linetype="dashed", color="black") +
  theme(axis.text.x=element_text(angle=45, vjust=1, hjust=1),
        panel.grid.major.y=element_line(color="grey", size=0.25))

grid.arrange(plot1, plot2, ncol=2)

models.plot.raw <- arrangeGrob(plot1, plot2, ncol=2)
ggsave("models.plot.raw.jpg", models.plot.raw, width=6, height=4, units="in", dpi=300)

#### FIG. S6: SPEARMAN CORRELATION MATRIX ####
sensitivity_corr <- subset(sensitivity_corr, !(Model1 %in% c("grp_ind_0.5_3.5", "grp_dep_0.5_3.5", "grp_dig_0.5_3.5")))
sensitivity_corr <- subset(sensitivity_corr, !(Model2 %in% c("grp_ind_0.5_3.5", "grp_dep_0.5_3.5", "grp_dig_0.5_3.5")))

sensitivity_corr$Model2 <- factor(sensitivity_corr$Model2)
sensitivity_corr$Model1 <- factor(sensitivity_corr$Model1, levels=levels(sensitivity_corr$Model2))

sensitivity_corr$Model1 <- factor(sensitivity_corr$Model1, 
                                  levels=c(
                                    "grp_ind_2.0_3.0", "grp_dep_2.0_3.0", "grp_dig_2.0_3.0",
                                    "grp_ind_2.0_3.5", "grp_dep_2.0_3.5", "grp_dig_2.0_3.5",
                                    "grp_ind_2.5_3.5", "grp_dep_2.5_3.5", "grp_dig_2.5_3.5",
                                    "grp_ind_3.4_3.0", "grp_dep_3.4_3.0", "grp_dig_3.4_3.0",
                                    "grp_ind_4.25_3.05", "grp_dep_4.25_3.05", "grp_dig_4.25_3.05"))
sensitivity_corr$Model2 <- factor(sensitivity_corr$Model2, levels=levels(sensitivity_corr$Model1))

# Remove duplicates, so that I only have one half of the correlation matrix.
for(i in 1:nrow(sensitivity_corr)){
  if(as.numeric(sensitivity_corr[i,"Model2"]) < as.numeric(sensitivity_corr[i,"Model1"])){
    sensitivity_corr[i,"R"] <- NA
  }
}
sensitivity_corr <- subset(sensitivity_corr, is.na(R)=="FALSE")

sensitivity_corr$Model1num <- as.numeric(sensitivity_corr$Model1)
sensitivity_corr$Model2 <- factor(sensitivity_corr$Model2, levels=rev(levels(sensitivity_corr$Model2)))
sensitivity_corr$Model2num <- as.numeric(sensitivity_corr$Model2)

rectangles <- read.csv("rectangles.csv")

plot.a <- ggplot() + 
  geom_tile(data = sensitivity_corr, aes(x=Model1num+1.5, y=Model2num+1.5, fill=R)) +
  scale_fill_gradient2(low="blue", high="red", mid="white", midpoint=0, na.value="white",
                       guide=guide_colorbar(title="Spearman's R\n", title.hjust=0.5)) +
  theme_bw() + plot_theme +
  theme(axis.text=element_blank(),
        axis.ticks=element_blank(),
        axis.title=element_blank(),
        legend.position=c(0.8,0.8),
        panel.border=element_blank()) +
  geom_rect(data=subset(rectangles, location %in% c("edge", "bottom")),
            mapping=aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2),
            color="black", fill="white", alpha=0.5) +
  geom_segment(data=subset(rectangles, location=="segments"),
               aes(x=x1, xend=x2, y=y1, yend=y2), color="black") +
  geom_text(data=subset(rectangles, location=="bottom"),
            aes(x=x1+(x2-x1)/2, y=y1+(y2-y1)/2, label=label), size=2.75) +
  geom_text(data=subset(rectangles, location=="edge" & type=="discr"),
            aes(x=x1+(x2-x1)/2, y=y1+(y2-y1)/2, label=label), size=2.75, angle=90) +
  geom_text(data=subset(rectangles, location=="edge" & type=="type"),
            aes(x=x1+(x2-x1)/2, y=y1+(y2-y1)/2, label=label), size=2.75)

#### FIG. S7: ISOSCAPE WITH SOURCES LABELED ####
# 2.0 / 3.0 ####
sources_data <- read.csv("mixing_model_inputs/sources.csv")
discrimination_data <- read.csv("mixing_model_inputs/discrimination_2.0_3.0.csv")

source_plot <- data.frame(matrix(ncol=0, nrow=7))
source_plot$d13C <- (sources_data$Meand13C + discrimination_data$Meand13C)
source_plot$d15N <- (sources_data$Meand15N + discrimination_data$Meand15N)
source_plot$SD13C <- sources_data$SDd13C
source_plot$SD15N <- sources_data$SDd15N
source_plot$Source <- sources_data$Source

data_plot <- data_single_point
data_plot$HABITAT <- factor(data_plot$HABITAT,
                             levels=c("rural", "suburban", "greenspace", "matrix", "urban"))

source_plot_2.0_3.0 <- source_plot

isoscape2.0_3.0<- ggplot() +
  
  # Add the coyotes
  geom_point(data=data_plot, aes(x=d13C, y=d15N, color = HABITAT, shape=HABITAT), size = 2) +
  geom_segment(aes(x=subset(source_plot_2.0_3.0, Source=="Berries")$d13C,
                   xend=subset(source_plot_2.0_3.0, Source=="Anthropogenic food")$d13C,
                   y=subset(source_plot_2.0_3.0, Source=="Berries")$d15N,
                   yend=subset(source_plot_2.0_3.0, Source=="Anthropogenic food")$d15N), linetype="dashed") +
  geom_segment(aes(x=subset(source_plot_2.0_3.0, Source=="Anthropogenic food")$d13C,
                   xend=subset(source_plot_2.0_3.0, Source=="Domestic pets")$d13C,
                   y=subset(source_plot_2.0_3.0, Source=="Anthropogenic food")$d15N,
                   yend=subset(source_plot_2.0_3.0, Source=="Domestic pets")$d15N), linetype="dashed") +
  geom_segment(aes(x=subset(source_plot_2.0_3.0, Source=="Insects")$d13C,
                   xend=subset(source_plot_2.0_3.0, Source=="Domestic pets")$d13C,
                   y=subset(source_plot_2.0_3.0, Source=="Insects")$d15N,
                   yend=subset(source_plot_2.0_3.0, Source=="Domestic pets")$d15N), linetype="dashed") +
  geom_segment(aes(x=subset(source_plot_2.0_3.0, Source=="Berries")$d13C,
                   xend=subset(source_plot_2.0_3.0, Source=="Insects")$d13C,
                   y=subset(source_plot_2.0_3.0, Source=="Berries")$d15N,
                   yend=subset(source_plot_2.0_3.0, Source=="Insects")$d15N), linetype="dashed") +
  
  scale_color_manual(values=c("springgreen3", "goldenrod", "tomato1", "deepskyblue2", "orchid1")) +
  guides(color=guide_legend(title="Habitat use"), 
         linetype=guide_legend(title="Habitat use"), 
         shape=guide_legend(title="Habitat use")) +
  
  # Add the sources
  geom_errorbar(data=source_plot_2.0_3.0, aes(x=d13C, ymin=d15N-SD15N, ymax=d15N+SD15N), width=.2,
                position=position_dodge(.9)) +
  geom_errorbarh(data=source_plot_2.0_3.0, aes(y=d15N, xmin=d13C-SD13C, xmax=d13C+SD13C),
                 position=position_dodge(.9)) +
  geom_point(data=source_plot_2.0_3.0, aes(x=d13C, y=d15N), color="black", size=2) +
  theme_bw() + plot_theme +
  ylab(expression(paste(delta^{15}, "N (\u2030)")))+
  xlab(expression(paste(delta^{13}, "C (\u2030)"))) + 
  labs(title="1.0 / 3.0") +
  theme(legend.position="none")

# 2.0 / 3.5 ####
sources_data <- read.csv("mixing_model_inputs/sources.csv")
discrimination_data <- read.csv("mixing_model_inputs/discrimination_2.0_3.5.csv")

source_plot <- data.frame(matrix(ncol=0, nrow=7))
source_plot$d13C <- (sources_data$Meand13C + discrimination_data$Meand13C)
source_plot$d15N <- (sources_data$Meand15N + discrimination_data$Meand15N)
source_plot$SD13C <- sources_data$SDd13C
source_plot$SD15N <- sources_data$SDd15N
source_plot$Source <- sources_data$Source

source_plot_2.0_3.5 <- source_plot

isoscape2.0_3.5<- ggplot() +
  # Add the coyotes
  geom_point(data=data_plot, aes(x=d13C, y=d15N, color = HABITAT, shape=HABITAT), size = 2) +
  geom_segment(aes(x=subset(source_plot_2.0_3.5, Source=="Berries")$d13C,
                   xend=subset(source_plot_2.0_3.5, Source=="Anthropogenic food")$d13C,
                   y=subset(source_plot_2.0_3.5, Source=="Berries")$d15N,
                   yend=subset(source_plot_2.0_3.5, Source=="Anthropogenic food")$d15N), linetype="dashed") +
  geom_segment(aes(x=subset(source_plot_2.0_3.5, Source=="Anthropogenic food")$d13C,
                   xend=subset(source_plot_2.0_3.5, Source=="Domestic pets")$d13C,
                   y=subset(source_plot_2.0_3.5, Source=="Anthropogenic food")$d15N,
                   yend=subset(source_plot_2.0_3.5, Source=="Domestic pets")$d15N), linetype="dashed") +
  geom_segment(aes(x=subset(source_plot_2.0_3.5, Source=="Insects")$d13C,
                   xend=subset(source_plot_2.0_3.5, Source=="Domestic pets")$d13C,
                   y=subset(source_plot_2.0_3.5, Source=="Insects")$d15N,
                   yend=subset(source_plot_2.0_3.5, Source=="Domestic pets")$d15N), linetype="dashed") +
  geom_segment(aes(x=subset(source_plot_2.0_3.5, Source=="Berries")$d13C,
                   xend=subset(source_plot_2.0_3.5, Source=="Insects")$d13C,
                   y=subset(source_plot_2.0_3.5, Source=="Berries")$d15N,
                   yend=subset(source_plot_2.0_3.5, Source=="Insects")$d15N), linetype="dashed") +
  
  scale_color_manual(values=c("springgreen3", "goldenrod", "tomato1", "deepskyblue2", "orchid1")) +
  guides(color=guide_legend(title="Habitat use"), 
         linetype=guide_legend(title="Habitat use"), 
         shape=guide_legend(title="Habitat use")) +
  
  # Add the sources
  geom_errorbar(data=source_plot_2.0_3.5, aes(x=d13C, ymin=d15N-SD15N, ymax=d15N+SD15N), width=.2,
                position=position_dodge(.9)) +
  geom_errorbarh(data=source_plot_2.0_3.5, aes(y=d15N, xmin=d13C-SD13C, xmax=d13C+SD13C),
                 position=position_dodge(.9)) +
  geom_point(data=source_plot_2.0_3.5, aes(x=d13C, y=d15N), color="black", size=2) +
  theme_bw() + plot_theme +
  ylab(expression(paste(delta^{15}, "N (\u2030)")))+
  xlab(expression(paste(delta^{13}, "C (\u2030)"))) + 
  labs(title="1.0 / 3.0") +
  theme(legend.position="none")

# 2.5 / 3.5 ####
sources_data <- read.csv("mixing_model_inputs/sources.csv")
discrimination_data <- read.csv("mixing_model_inputs/discrimination_2.5_3.5.csv")

source_plot <- data.frame(matrix(ncol=0, nrow=7))
source_plot$d13C <- (sources_data$Meand13C + discrimination_data$Meand13C)
source_plot$d15N <- (sources_data$Meand15N + discrimination_data$Meand15N)
source_plot$SD13C <- sources_data$SDd13C
source_plot$SD15N <- sources_data$SDd15N
source_plot$Source <- sources_data$Source

source_plot_2.5_3.5 <- source_plot

isoscape2.5_3.5<- ggplot() +
  # Add the coyotes
  geom_point(data=data_plot, aes(x=d13C, y=d15N, color = HABITAT, shape=HABITAT), size = 2) +
  geom_segment(aes(x=subset(source_plot_2.5_3.5, Source=="Berries")$d13C,
                   xend=subset(source_plot_2.5_3.5, Source=="Anthropogenic food")$d13C,
                   y=subset(source_plot_2.5_3.5, Source=="Berries")$d15N,
                   yend=subset(source_plot_2.5_3.5, Source=="Anthropogenic food")$d15N), linetype="dashed") +
  geom_segment(aes(x=subset(source_plot_2.5_3.5, Source=="Anthropogenic food")$d13C,
                   xend=subset(source_plot_2.5_3.5, Source=="Domestic pets")$d13C,
                   y=subset(source_plot_2.5_3.5, Source=="Anthropogenic food")$d15N,
                   yend=subset(source_plot_2.5_3.5, Source=="Domestic pets")$d15N), linetype="dashed") +
  geom_segment(aes(x=subset(source_plot_2.5_3.5, Source=="Insects")$d13C,
                   xend=subset(source_plot_2.5_3.5, Source=="Domestic pets")$d13C,
                   y=subset(source_plot_2.5_3.5, Source=="Insects")$d15N,
                   yend=subset(source_plot_2.5_3.5, Source=="Domestic pets")$d15N), linetype="dashed") +
  geom_segment(aes(x=subset(source_plot_2.5_3.5, Source=="Berries")$d13C,
                   xend=subset(source_plot_2.5_3.5, Source=="Insects")$d13C,
                   y=subset(source_plot_2.5_3.5, Source=="Berries")$d15N,
                   yend=subset(source_plot_2.5_3.5, Source=="Insects")$d15N), linetype="dashed") +
  
  scale_color_manual(values=c("springgreen3", "goldenrod", "tomato1", "deepskyblue2", "orchid1")) +
  guides(color=guide_legend(title="Habitat use"), 
         linetype=guide_legend(title="Habitat use"), 
         shape=guide_legend(title="Habitat use")) +
  
  # Add the sources
  geom_errorbar(data=source_plot_2.5_3.5, aes(x=d13C, ymin=d15N-SD15N, ymax=d15N+SD15N), width=.2,
                position=position_dodge(.9)) +
  geom_errorbarh(data=source_plot_2.5_3.5, aes(y=d15N, xmin=d13C-SD13C, xmax=d13C+SD13C),
                 position=position_dodge(.9)) +
  geom_point(data=source_plot_2.5_3.5, aes(x=d13C, y=d15N), color="black", size=2) +
  theme_bw() + plot_theme +
  ylab(expression(paste(delta^{15}, "N (\u2030)")))+
  xlab(expression(paste(delta^{13}, "C (\u2030)"))) + 
  labs(title="1.0 / 3.0") +
  theme(legend.position="none")

# 3.4 / 3.0 ####
sources_data <- read.csv("mixing_model_inputs/sources.csv")
discrimination_data <- read.csv("mixing_model_inputs/discrimination_3.4_3.0.csv")

source_plot <- data.frame(matrix(ncol=0, nrow=7))
source_plot$d13C <- (sources_data$Meand13C + discrimination_data$Meand13C)
source_plot$d15N <- (sources_data$Meand15N + discrimination_data$Meand15N)
source_plot$SD13C <- sources_data$SDd13C
source_plot$SD15N <- sources_data$SDd15N
source_plot$Source <- sources_data$Source

source_plot_3.4_3.0 <- source_plot

isoscape3.4_3.0<- ggplot() +
  # Add the coyotes
  geom_point(data=data_plot, aes(x=d13C, y=d15N, color = HABITAT, shape=HABITAT), size = 2) +
  geom_segment(aes(x=subset(source_plot_3.4_3.0, Source=="Berries")$d13C,
                   xend=subset(source_plot_3.4_3.0, Source=="Anthropogenic food")$d13C,
                   y=subset(source_plot_3.4_3.0, Source=="Berries")$d15N,
                   yend=subset(source_plot_3.4_3.0, Source=="Anthropogenic food")$d15N), linetype="dashed") +
  geom_segment(aes(x=subset(source_plot_3.4_3.0, Source=="Anthropogenic food")$d13C,
                   xend=subset(source_plot_3.4_3.0, Source=="Domestic pets")$d13C,
                   y=subset(source_plot_3.4_3.0, Source=="Anthropogenic food")$d15N,
                   yend=subset(source_plot_3.4_3.0, Source=="Domestic pets")$d15N), linetype="dashed") +
  geom_segment(aes(x=subset(source_plot_3.4_3.0, Source=="Insects")$d13C,
                   xend=subset(source_plot_3.4_3.0, Source=="Domestic pets")$d13C,
                   y=subset(source_plot_3.4_3.0, Source=="Insects")$d15N,
                   yend=subset(source_plot_3.4_3.0, Source=="Domestic pets")$d15N), linetype="dashed") +
  geom_segment(aes(x=subset(source_plot_3.4_3.0, Source=="Berries")$d13C,
                   xend=subset(source_plot_3.4_3.0, Source=="Insects")$d13C,
                   y=subset(source_plot_3.4_3.0, Source=="Berries")$d15N,
                   yend=subset(source_plot_3.4_3.0, Source=="Insects")$d15N), linetype="dashed") +
  
  scale_color_manual(values=c("springgreen3", "goldenrod", "tomato1", "deepskyblue2", "orchid1")) +
  guides(color=guide_legend(title="Habitat use"), 
         linetype=guide_legend(title="Habitat use"), 
         shape=guide_legend(title="Habitat use")) +
  
  # Add the sources
  geom_errorbar(data=source_plot_3.4_3.0, aes(x=d13C, ymin=d15N-SD15N, ymax=d15N+SD15N), width=.2,
                position=position_dodge(.9)) +
  geom_errorbarh(data=source_plot_3.4_3.0, aes(y=d15N, xmin=d13C-SD13C, xmax=d13C+SD13C),
                 position=position_dodge(.9)) +
  geom_point(data=source_plot_3.4_3.0, aes(x=d13C, y=d15N), color="black", size=2) +
  theme_bw() + plot_theme +
  ylab(expression(paste(delta^{15}, "N (\u2030)")))+
  xlab(expression(paste(delta^{13}, "C (\u2030)"))) + 
  labs(title="1.0 / 3.0") +
  theme(legend.position="none")

# 4.25 / 3.05 ####
sources_data <- read.csv("mixing_model_inputs/sources.csv")
discrimination_data <- read.csv("mixing_model_inputs/discrimination_4.25_3.05.csv")

source_plot <- data.frame(matrix(ncol=0, nrow=7))
source_plot$d13C <- (sources_data$Meand13C + discrimination_data$Meand13C)
source_plot$d15N <- (sources_data$Meand15N + discrimination_data$Meand15N)
source_plot$SD13C <- sources_data$SDd13C
source_plot$SD15N <- sources_data$SDd15N
source_plot$Source <- sources_data$Source

source_plot_4.25_3.05 <- source_plot

isoscape4.25_3.05<- ggplot() +
  # Add the coyotes
  geom_point(data=data_plot, aes(x=d13C, y=d15N, color = HABITAT, shape=HABITAT), size = 2) +
  geom_segment(aes(x=subset(source_plot_4.25_3.05, Source=="Berries")$d13C,
                   xend=subset(source_plot_4.25_3.05, Source=="Anthropogenic food")$d13C,
                   y=subset(source_plot_4.25_3.05, Source=="Berries")$d15N,
                   yend=subset(source_plot_4.25_3.05, Source=="Anthropogenic food")$d15N), linetype="dashed") +
  geom_segment(aes(x=subset(source_plot_4.25_3.05, Source=="Anthropogenic food")$d13C,
                   xend=subset(source_plot_4.25_3.05, Source=="Domestic pets")$d13C,
                   y=subset(source_plot_4.25_3.05, Source=="Anthropogenic food")$d15N,
                   yend=subset(source_plot_4.25_3.05, Source=="Domestic pets")$d15N), linetype="dashed") +
  geom_segment(aes(x=subset(source_plot_4.25_3.05, Source=="Insects")$d13C,
                   xend=subset(source_plot_4.25_3.05, Source=="Domestic pets")$d13C,
                   y=subset(source_plot_4.25_3.05, Source=="Insects")$d15N,
                   yend=subset(source_plot_4.25_3.05, Source=="Domestic pets")$d15N), linetype="dashed") +
  geom_segment(aes(x=subset(source_plot_4.25_3.05, Source=="Berries")$d13C,
                   xend=subset(source_plot_4.25_3.05, Source=="Insects")$d13C,
                   y=subset(source_plot_4.25_3.05, Source=="Berries")$d15N,
                   yend=subset(source_plot_4.25_3.05, Source=="Insects")$d15N), linetype="dashed") +
  
  scale_color_manual(values=c("springgreen3", "goldenrod", "tomato1", "deepskyblue2", "orchid1")) +
  guides(color=guide_legend(title="Habitat use"), 
         linetype=guide_legend(title="Habitat use"), 
         shape=guide_legend(title="Habitat use")) +
  
  # Add the sources
  geom_errorbar(data=source_plot_4.25_3.05, aes(x=d13C, ymin=d15N-SD15N, ymax=d15N+SD15N), width=.2,
                position=position_dodge(.9)) +
  geom_errorbarh(data=source_plot_4.25_3.05, aes(y=d15N, xmin=d13C-SD13C, xmax=d13C+SD13C),
                 position=position_dodge(.9)) +
  geom_point(data=source_plot_4.25_3.05, aes(x=d13C, y=d15N), color="black", size=2) +
  theme_bw() + plot_theme +
  ylab(expression(paste(delta^{15}, "N (\u2030)")))+
  xlab(expression(paste(delta^{13}, "C (\u2030)"))) + 
  labs(title="1.0 / 3.0") +
  theme(legend.position="none")


# Assemble figure ####
levels(data_plot$HABITAT)[levels(data_plot$HABITAT)=="urban"] <- "urban (no GPS)"

legend <- cowplot::get_legend(ggplot() +
                                # Add the coyotes
                                geom_point(data=data_plot, aes(x=d13C, y=d15N, color = HABITAT, shape=HABITAT), size = 2) +
                                
                                scale_color_manual(values=c("springgreen3", "goldenrod", "tomato1", "deepskyblue2", "orchid1")) +
                                guides(color=guide_legend(title="Habitat use"), 
                                       linetype=guide_legend(title="Habitat use"), 
                                       shape=guide_legend(title="Habitat use")) +
                                
                                theme_bw() + plot_theme +
                                ylab(expression(paste(delta^{15}, "N (\u2030)")))+
                                xlab(expression(paste(delta^{13}, "C (\u2030)"))) + 
                                labs(title="1.0 / 3.0") +
                                theme(legend.position="right",
                                      legend.title=element_text(hjust=0.5)))

grid.arrange(isoscape2.0_3.0 + labs(title="2.0 / 3.0") + theme(axis.title=element_blank(), plot.title=element_text(size=9)),
             isoscape2.0_3.5 + labs(title="2.0 / 3.5") + theme(axis.title=element_blank(), plot.title=element_text(size=9)),
             isoscape2.5_3.5 + labs(title="2.5 / 3.5") + theme(axis.title=element_blank(), plot.title=element_text(size=9)),
             isoscape3.4_3.0 + labs(title="3.4 / 3.0") + theme(axis.title=element_blank(), plot.title=element_text(size=9)),
             isoscape4.25_3.05 + labs(title="4.3 / 3.1") + theme(axis.title=element_blank(), plot.title=element_text(size=9)),
             legend,
             bottom=textGrob(expression(paste(delta^{13}, "C (\u2030)")), hjust=0.5, gp=gpar(cex=0.75)),
             left=textGrob(expression(paste(delta^{15}, "N (\u2030)")), hjust=0.5, rot=90, gp=gpar(cex=0.75)),
             layout_matrix=rbind(c(1,1,2,2,3,3),
                                 c(4,4,5,5,6,6)))

isoscape.panel <- arrangeGrob(isoscape2.0_3.0 + labs(title="2.0 / 3.0") + theme(axis.title=element_blank(), plot.title=element_text(size=9)),
                              isoscape2.0_3.5 + labs(title="2.0 / 3.5") + theme(axis.title=element_blank(), plot.title=element_text(size=9)),
                              isoscape2.5_3.5 + labs(title="2.5 / 3.5") + theme(axis.title=element_blank(), plot.title=element_text(size=9)),
                              isoscape3.4_3.0 + labs(title="3.4 / 3.0") + theme(axis.title=element_blank(), plot.title=element_text(size=9)),
                              isoscape4.25_3.05 + labs(title="4.3 / 3.1") + theme(axis.title=element_blank(), plot.title=element_text(size=9)),
                              legend,
                              bottom=textGrob(expression(paste(delta^{13}, "C (\u2030)")), hjust=0.5, gp=gpar(cex=0.75)),
                              left=textGrob(expression(paste(delta^{15}, "N (\u2030)")), hjust=0.5, rot=90, gp=gpar(cex=0.75)),
                              layout_matrix=rbind(c(1,1,2,2,3,3),
                                                  c(4,4,5,5,6,6)))

ggsave("publication_figures/isoscape_panel.jpg", isoscape.panel, width=7, height=5.5, units="in", dpi=300)
rm(sources_data, discrimination_data, source_plot, data_plot,
   isoscape2.0_3.0, isoscape2.0_3.5, isoscape2.5_3.5, isoscape0.5_3.5, isoscape4.25_3.05,
   source_plot_2.0_3.0, source_plot_2.0_3.5, source_plot_2.5_3.5, source_plot_0.5_3.5, source_plot_4.25_3.05)
