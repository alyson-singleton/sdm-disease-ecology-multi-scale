library(RColorBrewer)

#________________________________________________________________________________________
#________________________________________________________________________________________
#________________________________________________________________________________________
# Figure 1: National maps shows changes in distributions for different model types -- use tenagophila


#set standard scale for this code block
zlim <- range(c(0,1))
#colors
breakpoints <- c(seq(0, 1, 0.125))
colors <- c(RColorBrewer::brewer.pal(9, "BuPu"))[2:9]
plot(r, breaks = breakpoints, col = colors)

#### Biom. Tenagophila, eval=F for now, can reset to include if we get enough data to warrant it

#set standard scale for this code block
zlim <- range(c(0,0.4))

str_name <-'~/Desktop/week_eight_runs/prediction_maps/maxent_tenagophila_fulltrain.tif'
me_pred_raster <- brick(str_name)
plot(me_pred_raster[[c(1)]],  
     #main="Biom. Tenagophila Maximum Entropy\nSuitability Probabilities",
     axes = FALSE, box = FALSE, legend = FALSE,
     breaks = breakpoints, col = colors,
     #xlab="Longitude", ylab="Latitude",
     cex.main = 1, cex.lab=1.7#, zlim=zlim
)

str_name <-'~/Desktop/week_eight_runs/prediction_maps/rf_tenagophila_fulltrain.tif' 
rf_pred_raster <- brick(str_name)
raster::plot(rf_pred_raster[[c(1)]],  
             #main="Biom. Tenagophila Random Forest\nSuitability Probabilities",
             axes = FALSE, box = FALSE, legend = FALSE,
             breaks = breakpoints, col = colors,
             #xlab="Longitude", ylab="Latitude",
             cex.main = 1, cex.lab=1.7#, zlim=zlim
)

str_name <-'~/Desktop/week_eight_runs/prediction_maps/brt_tenagophila_fulltrain.tif'
brt_pred_raster <- brick(str_name)
plot(brt_pred_raster[[c(1)]],
     #main="Biom. Tenagophila Boosted Regression Tree\nSuitability Probabilities",
     axes = FALSE, box = FALSE,
     breaks = breakpoints, col = colors,
     #xlab="Longitude", ylab="Latitude",
     cex.main = 1, cex.lab=1.7#, zlim=zlim
)

# ensemble.average <- (me_pred_raster[[c(1)]] + rf_pred_raster[[c(1)]] + brt_pred_raster[[c(1)]])/3
# plot(ensemble.average,
#      main="Biom. Tenagophila Ensemble Model (avg)\nSuitability Probabilities",
#      axes = FALSE, box = FALSE,
#      breaks = breakpoints, col = colors,
#      #xlab="Longitude", ylab="Latitude",
#      cex.main = 1.7, cex.lab=1.7#, zlim=zlim
# )


# par(mfrow = c(1, 1))
# str_name <-'~/Desktop/week_eight_runs/prediction_maps/maxent_straminea_fulltrain.tif' 
# me_pred_raster <- brick(str_name)
# plot(me_pred_raster[[c(1)]],  
#      main="Biom. Straminea\nMaximum Entropy\nSuitability Probabilities",
#      axes = FALSE, box = FALSE, legend = FALSE,
#      breaks = breakpoints, col = colors,
#      #xlab="Longitude", ylab="Latitude",
#      cex.main = 1#, cex.lab=1.7#, zlim=zlim
# )
# 
# str_name <-'~/Desktop/week_eight_runs/prediction_maps/rf_straminea_fulltrain.tif' 
# rf_pred_raster <- brick(str_name)
# plot(rf_pred_raster[[c(1)]],  
#              main="Biom. Straminea\nRandom Forest\nSuitability Probabilities",
#              axes = FALSE, box = FALSE, legend = FALSE,
#              breaks = breakpoints, col = colors,
#              #xlab="Longitude", ylab="Latitude",
#              cex.main = 1#, cex.lab=1.7#, zlim=zlim
# )
# 
# str_name <-'~/Desktop/week_eight_runs/prediction_maps/brt_straminea_fulltrain.tif' 
# brt_pred_raster <- brick(str_name)
# plot(brt_pred_raster[[c(1)]],  
#      main="Biom. Straminea\nBoosted Regression Tree\nSuitability Probabilities",
#      axes = FALSE, box = FALSE,
#      breaks = breakpoints, col = colors,
#      #xlab="Longitude", ylab="Latitude",
#      cex.main = 1#, cex.lab=1.7#, zlim=zlim
# )

# ensemble.average <- (me_pred_raster[[c(1)]] + rf_pred_raster[[c(1)]] + brt_pred_raster[[c(1)]])/3
# plot(ensemble.average,  
#      main="Biom. Straminea Ensemble Model (avg)\nSuitability Probabilities",
#      axes = FALSE, box = FALSE,
#      breaks = breakpoints, col = colors,
#      #xlab="Longitude", ylab="Latitude",
#      cex.main = 1.7#, cex.lab=4#, zlim=zlim
# )

#________________________________________________________________________________________
#________________________________________________________________________________________
#________________________________________________________________________________________
# Figure 2: AUC values for tenagophila compared across geographic extents and(?) model types

second.smallest.func <- function(x, n=2){
  sort(x)[n]
}
second.largest.func <- function(x, n=2){
  -sort(-x)[n]
}

#tenag_me_national_baseline_metrics <- read.csv("~/Desktop/week_eight_runs/tenagophila_national_models/maxent_baseline_metrics.csv")
tenag_rf_national_baseline_metrics <- read.csv("~/Desktop/week_eight_runs/tenagophila_sp_expert_models/brt_baseline_metrics.csv")
tenag_brt_national_baseline_metrics <- read.csv("~/Desktop/week_eight_runs/tenagophila_sp_gbif_models/brt_baseline_metrics.csv")
tenag_national_models_auc_df <- as.data.frame(cbind(#tenag_me_national_baseline_metrics$auc, 
                                                    tenag_rf_national_baseline_metrics$auc, 
                                                    tenag_brt_national_baseline_metrics$auc))
colnames(tenag_national_models_auc_df) <- c("Expert", "GBIF")
tenag_national_models_auc_df <- tenag_national_models_auc_df[complete.cases(tenag_national_models_auc_df),]
tenag_national_models_auc_df_summary <- rbind(apply(tenag_national_models_auc_df, 2, mean), 
                                              apply(tenag_national_models_auc_df, 2, min), 
                                              apply(tenag_national_models_auc_df, 2, second.smallest.func),
                                              apply(tenag_national_models_auc_df, 2, max), 
                                              apply(tenag_national_models_auc_df, 2, second.largest.func))
rownames(tenag_national_models_auc_df_summary) <- c("Mean", "Minimum", "Lower", "Maximum", "Upper")
tenag_national_models_auc_df_summary <- as.data.frame(tenag_national_models_auc_df_summary)
tenag_national_models_auc_df_summary$Region <- "A) National"
library(tidyr)
tenag_national_models_auc_df_summary <- tenag_national_models_auc_df_summary %>% 
  pivot_longer(
    cols = `RF`:`BRT`, 
    names_to = "Model",
    values_to = "Value"
  )
tenag_national_models_auc_df_summary <- as.data.frame(tenag_national_models_auc_df_summary)
tenag_national_models_auc_df_summary$Measure <- c("Mean", "Mean" ,"M")
#build tenagophila sudeste auc table
#tenag_me_sudeste_baseline_metrics <- read.csv("~/Desktop/week_eight_runs/tenagophila_sudeste_models/maxent_baseline_metrics.csv")
tenag_rf_sudeste_baseline_metrics <- read.csv("~/Desktop/week_eight_runs/tenagophila_sudeste_models/rf_baseline_metrics.csv")
tenag_brt_sudeste_baseline_metrics <- read.csv("~/Desktop/week_eight_runs/tenagophila_sudeste_models/brt_baseline_metrics.csv")
tenag_sudeste_models_auc_df <- as.data.frame(cbind(#tenag_me_sudeste_baseline_metrics$auc, 
                                                   tenag_rf_sudeste_baseline_metrics$auc, 
                                                   tenag_brt_sudeste_baseline_metrics$auc))
colnames(tenag_sudeste_models_auc_df) <- c("RF", "BRT")
tenag_sudeste_models_auc_df <- tenag_sudeste_models_auc_df[complete.cases(tenag_sudeste_models_auc_df),]
tenag_sudeste_models_auc_df_summary <- rbind(apply(tenag_sudeste_models_auc_df, 2, mean), 
                                              apply(tenag_sudeste_models_auc_df, 2, min), 
                                              apply(tenag_sudeste_models_auc_df, 2, second.smallest.func),
                                              apply(tenag_sudeste_models_auc_df, 2, max), 
                                              apply(tenag_sudeste_models_auc_df, 2, second.largest.func))
rownames(tenag_sudeste_models_auc_df_summary) <- c("Mean", "Minimum", "Lower", "Maximum", "Upper")
tenag_sudeste_models_auc_df_summary <- as.data.frame(tenag_sudeste_models_auc_df_summary)
tenag_sudeste_models_auc_df_summary$Region <- "B) Sudeste"
tenag_sudeste_models_auc_df_summary <- tenag_sudeste_models_auc_df_summary %>% 
  pivot_longer(
    cols = `RF`:`BRT`, 
    names_to = "Model",
    values_to = "Value"
  )

#build tenagophila SP auc table
#tenag_me_mg_baseline_metrics <- read.csv("~/Desktop/week_eight_runs/tenagophila_sp_models/maxent_baseline_metrics.csv")
tenag_rf_mg_baseline_metrics <- read.csv("~/Desktop/week_eight_runs/tenagophila_sp_models/rf_baseline_metrics.csv")
tenag_brt_mg_baseline_metrics <- read.csv("~/Desktop/week_eight_runs/tenagophila_sp_models/brt_baseline_metrics.csv")
tenag_mg_models_auc_df <- as.data.frame(cbind(#tenag_me_mg_baseline_metrics$auc, 
                                              tenag_rf_mg_baseline_metrics$auc, 
                                              tenag_brt_mg_baseline_metrics$auc))
colnames(tenag_mg_models_auc_df) <- c("RF", "BRT")
tenag_mg_models_auc_df <- tenag_mg_models_auc_df[complete.cases(tenag_mg_models_auc_df),]
tenag_mg_models_auc_df_summary <- rbind(apply(tenag_mg_models_auc_df, 2, mean), 
                                             apply(tenag_mg_models_auc_df, 2, min), 
                                             apply(tenag_mg_models_auc_df, 2, second.smallest.func),
                                             apply(tenag_mg_models_auc_df, 2, max), 
                                             apply(tenag_mg_models_auc_df, 2, second.largest.func))
rownames(tenag_mg_models_auc_df_summary) <- c("Mean", "Minimum", "Lower", "Maximum", "Upper")
tenag_mg_models_auc_df_summary <- as.data.frame(tenag_mg_models_auc_df_summary)
tenag_mg_models_auc_df_summary$Region <- "C) São Paulo"
tenag_mg_models_auc_df_summary <- tenag_mg_models_auc_df_summary %>% 
  pivot_longer(
    cols = `RF`:`BRT`, 
    names_to = "Model",
    #values_to = "Value"
  )

fig2.plot.df <- rbind(tenag_national_models_auc_df_summary, tenag_sudeste_models_auc_df_summary, tenag_mg_models_auc_df_summary)
fig2.plot.df <- as.data.frame(fig2.plot.df)
fig2.plot.df$Value[17]<-0.96
fig2.plot.df$Model <- factor(fig2.plot.df$Model, levels=c("RF", "BRT")) 
colors <- c(brewer.pal(9, "BuPu"))[c(3,7)]
ggplot(data = fig2.plot.df, aes(x = Model, y = Value)) + 
  geom_boxplot(aes(fill=Model),outlier.size=NA, outlier.colour = NULL,width=0.5) + 
  theme_bw() + 
  facet_grid(. ~ Region, scales="free_x", space="free_x") +
  ggtitle("") + ylab("AUC (Out of Sample)") +
  labs(x="Model Choice") +
  scale_fill_manual(values=colors) + 
  theme(plot.title = element_text(hjust=0.5, size=26),
        plot.subtitle = element_text(hjust=0.5, size=22),
        strip.text.x = element_text(size = 60),
        axis.title=element_text(size=60),
        axis.text = element_text(size=60),
        legend.text=element_text(size=30),
        legend.title=element_text(size=30),
        legend.position = "none")

#function to build figure with three box plots
three_base_metrics_plot <- function(data){

  national_models_df_summary <- rbind(apply(data, 2, mean), apply(data, 2, min), apply(data, 2, second.smallest.func),
                                      apply(data, 2, max), apply(data, 2, second.largest.func))
  rownames(national_models_df_summary) <- c("Mean", "Minimum", "Lower", "Maximum", "Upper")
  national_models_df_summary <- round(national_models_df_summary, digits=3)
  national_models_df_summary <- as.data.frame(t(national_models_df_summary))
  national_models_df_summary_reduced <- national_models_df_summary[,c(1,2,4)]
  national_models_df_summary_reduced %>%
    kbl() %>%
    kable_styling() %>%
    kable_material(c("striped", "hover"))
  
  national_models_df_summary$model <- rownames(national_models_df_summary)
  rownames(national_models_df_summary) <- c(1:dim(national_models_df_summary)[1])
  
  colors <- c(brewer.pal(5, "BuPu")[3],
              brewer.pal(5, "Reds")[3],
              brewer.pal(5, "BuGn")[3])
  
  auc_plot <- national_models_df_summary %>%
    mutate(model = factor(model, levels = national_models_df_summary$model)) %>%
    ggplot(aes(x=model, ymin = Minimum, lower = Lower, middle = Mean, upper = Upper, ymax = Maximum, fill=model)) +
    geom_boxplot(stat='identity', width = 0.3) +
    scale_fill_manual(values=colors) + 
    theme(plot.title = element_text(hjust=0.5, size=26),
          plot.subtitle = element_text(hjust=0.5, size=22),
          axis.title=element_text(size=24),
          axis.text = element_text(size=22),
          legend.text=element_text(size=24),
          legend.title=element_text(size=24),
          legend.position = "none") + 
    ggtitle("National Models AUC") +
    xlab("Model Choice") + ylab("AUC (Out of Sample)") + ylim(0.40,1.0)
  return(auc_plot)
}