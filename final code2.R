
packages <- c("data.table", "tidyverse", "sf", "sp", "spdep",'ggplot2',
              "geofacet", "ggpubr", "ggthemes", 'readxl', 'lubridate','dplyr','readr',"extRemes","ggcorrplot")
lapply(packages, library, character.only = TRUE)

#############################Intercity connectivity enhances disaster resilience in the absence of mobility constraints####################

overall$tergdp<-overall$tergdp/100
overall$total_prec<-overall$total_prec/100
overall$mean_elevation<-overall$mean_elevation/1000

current<-read_csv('network23.csv')
iso<-read_csv('networkiso.csv')
weak<-read_csv('networkweak.csv')

##compared with connectivity baselines in isolated network
overall1 <- overall %>%
  left_join(
    iso %>% dplyr::select(citycode, interaction_unevenness(mad), economic-tier_mobility_concentration, rich_neighbor_preference), 
    by = "citycode"
  ) %>%
  mutate(
    interaction_unevenness(mad) = coalesce(interaction_unevenness(mad).y, interaction_unevenness(mad).x),
    economic-tier_mobility_concentration = coalesce(economic-tier_mobility_concentration.y, economic-tier_mobility_concentration.x),
    rich_neighbor_preference= coalesce(rich_neighbor_preference.y, rich_neighbor_preference.x)
  ) %>%
  dplyr::select(-interaction_unevenness(mad).x, -interaction_unevenness(mad).y, 
         -economic-tier_mobility_concentration.x, -economic-tier_mobility_concentration.y, 
         -rich_neighbor_preference.x, -rich_neighbor_preference.y)

hlm_model_test1 <- lm(
  total_intrac ~ total_prec + max_wind_diff + mean_elevation +proportion_10+proportion_60+
    proportion_90 + damdes + max_Rpolicy + max_Cpolicy  +mobility_per_capita+interaction_unevenness(mad)+economic-tier_mobility_concentration+rich_neighbor_preference,
  data = overall)

summary(hlm_model_test1)

overall1$predicted_total_intrac <- predict(hlm_model_test1, newdata = overall1)

overall1$netimproveratio<-(overall1$total_intrac-overall1$predicted_total_intrac)/abs(overall1$predicted_total_intrac)

###extended data fig 2c, mapped in arcgis
overall1_fige2<-overall1%>% 
  group_by(citycode)%>%
  summarise(n=n(),
    netimproveratio0=100*median(netimproveratio))


##compared with connectivity baselines in weak-control network
overall2 <- overall %>%
  left_join(
    weak %>% dplyr::select(citycode, interaction_unevenness(mad), economic-tier_mobility_concentration, rich_neighbor_preference), 
    by = "citycode"
  ) %>%
  mutate(
    interaction_unevenness(mad) = coalesce(interaction_unevenness(mad).y, interaction_unevenness(mad).x),
    economic-tier_mobility_concentration = coalesce(economic-tier_mobility_concentration.y, economic-tier_mobility_concentration.x),
    rich_neighbor_preference = coalesce(rich_neighbor_preference.y, rich_neighbor_preference.x)
  ) %>%
  dplyr::select(-interaction_unevenness(mad).x, -interaction_unevenness(mad).y, 
         -economic-tier_mobility_concentration.x, -economic-tier_mobility_concentration.y, 
         -rich_neighbor_preference.x, -rich_neighbor_preference.y)

hlm_model_test1 <- lm(
  total_intrac ~ total_prec + max_wind_diff + mean_elevation +proportion_10+proportion_60+
    proportion_90 + damdes + max_Rpolicy + max_Cpolicy  +mobility_per_capita+interaction_unevenness(mad)+economic-tier_mobility_concentration+rich_neighbor_preference,
  data = overall)

summary(hlm_model_test1)

overall2$predicted_total_intrac <- predict(hlm_model_test1, newdata = overall2)
overall2$netimprove<-(overall2$total_intrac-overall2$predicted_total_intrac)/abs(overall2$predicted_total_intrac)


###extended data fig 2b, mapped in arcgis
overall2_fige2<-overall2%>%
  group_by(citycode)%>%
  summarise(n=n(),
            netimprove=100*median(netimprove))


###each city in 100-year-return period rainfall
##compared with connectivity baselines in isolated network
maximpact<-read_csv("maximpact.csv")
maximpact<-maximpact%>%
  filter(max_return==100)

hlm_model_test1 <- lm(
  total_intrac ~ total_prec + max_wind_diff + mean_elevation +proportion_10+proportion_60+
    proportion_90 + damdes + max_Rpolicy + max_Cpolicy  +mobility_per_capita+interaction_unevenness(mad)+economic-tier_mobility_concentration+rich_neighbor_preference,
  data = overall)

summary(hlm_model_test1)

maximpact$predicted_total_intrac <- predict(hlm_model_test1, newdata = maximpact)

maximpact_updated <- maximpact %>%
  left_join(
    iso %>% dplyr::select(citycode, interaction_unevenness(mad), economic-tier_mobility_concentration, rich_neighbor_preference), 
    by = "citycode"
  ) %>%
  mutate(
    interaction_unevenness(mad) = coalesce(interaction_unevenness(mad).y, interaction_unevenness(mad).x),
    economic-tier_mobility_concentration = coalesce(economic-tier_mobility_concentration.y, economic-tier_mobility_concentration.x),
    rich_neighbor_preference = coalesce(rich_neighbor_preference.y, rich_neighbor_preference.x)
  ) %>%
  dplyr::select(-interaction_unevenness(mad).x, -interaction_unevenness(mad).y, 
         -economic-tier_mobility_concentration.x, -economic-tier_mobility_concentration.y, 
         -rich_neighbor_preference.x, -rich_neighbor_preference.y)


maximpact_updated$predicted_total_intrac0 <- predict(hlm_model_test1, newdata = maximpact_updated)
maximpact_updated<-maximpact_updated%>%
  dplyr::select(citycode,max_return,predicted_total_intrac0,interaction_unevenness(mad),economic-tier_mobility_concentration,rich_neighbor_preference)

maximpact<-merge(maximpact,maximpact_updated,by=c('citycode',"max_return"))
####fig 2
maximpact$interaction_unevenness(mad)<-(maximpact$interaction_unevenness(mad).x-maximpact$interaction_unevenness(mad).y)/maximpact$interaction_unevenness(mad).x
maximpact$economic-tier_mobility_concentration<-(maximpact$economic-tier_mobility_concentration.x-maximpact$economic-tier_mobility_concentration.y)/maximpact$economic-tier_mobility_concentration.x
maximpact$rich_neighbor_preference<-(maximpact$rich_neighbor_preference.x-maximpact$rich_neighbor_preference.y)/maximpact$rich_neighbor_preference.x

maximpact$netpredit<-maximpact$predicted_total_intrac-maximpact$predicted_total_intrac0
maximpact$netpreditratio<-maximpact$netpredit/abs(maximpact$predicted_total_intrac0)

maximpact<-maximpact%>%
  filter(!is.na(netpreditratio))

summary(maximpact$netpreditratio)

summary(maximpact$interaction_unevenness(mad))
summary(maximpact$economic-tier_mobility_concentration)
summary(maximpact$rich_neighbor_preference)

####fig s3
a <- sum(-172.6227408 * (maximpact$interaction_unevenness(mad).x - maximpact$interaction_unevenness(mad).y) / abs(maximpact$predicted_total_intrac), na.rm = TRUE)
b <- sum(-0.1470388 * (maximpact$economic-tier_mobility_concentration.x - maximpact$economic-tier_mobility_concentration.y) / abs(maximpact$predicted_total_intrac), na.rm = TRUE)
c <- sum(0.0845488 * (maximpact$rich_neighbor_preference.x - maximpact$rich_neighbor_preference.y) / abs(maximpact$predicted_total_intrac), na.rm = TRUE)


total_contribution <- a + b + c
a_ratio <- a / total_contribution
b_ratio <- b / total_contribution
c_ratio <- c / total_contribution

contribution_data <- data.frame(
  Variable = c("Mean Absolute Deviation", "economic-tier_mobility_concentration", "Dis_ex"),
  Contribution_Ratio = c(a_ratio, b_ratio, c_ratio)
)


p <- ggplot(contribution_data, aes(x = Variable, y = Contribution_Ratio, fill = Variable)) +
  geom_bar(stat = "identity", width = 0.7, color = "black") +
  labs(
    x = "Variable",
    y = "Contribution Ratio",
  ) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +  # 转换为百分比格式
  theme_bw() +
  theme(
    text = element_text(size = 7),
    plot.title = element_text(hjust = 0.5),
    legend.position = "none"
  )





##compared with connectivity baselines in weak-control network

maximpact1<-maximpact1%>%
  filter(max_return==100)
hlm_model_test1 <- lm(
  total_intrac ~ total_prec + max_wind_diff + mean_elevation +proportion_10+proportion_60+
    proportion_90 + damdes + max_Rpolicy + max_Cpolicy  +mobility_per_capita+interaction_unevenness(mad)+economic-tier_mobility_concentration+rich_neighbor_preference,
  data = overall)



maximpact1$predicted_total_intrac <- predict(hlm_model_test1, newdata = maximpact1)

maximpact_updated1 <- maximpact1 %>%
  left_join(
    weak %>% dplyr::select(citycode, interaction_unevenness(mad), economic-tier_mobility_concentration, rich_neighbor_preference), 
    by = "citycode"
  ) %>%
  mutate(
    interaction_unevenness(mad) = coalesce(interaction_unevenness(mad).y, interaction_unevenness(mad).x),
    economic-tier_mobility_concentration = coalesce(economic-tier_mobility_concentration.y, economic-tier_mobility_concentration.x),
    rich_neighbor_preference = coalesce(rich_neighbor_preference.y, rich_neighbor_preference.x)
  ) %>%
  dplyr::select(-interaction_unevenness(mad).x, -interaction_unevenness(mad).y, 
         -economic-tier_mobility_concentration.x, -economic-tier_mobility_concentration.y, 
         -rich_neighbor_preference.x, -rich_neighbor_preference.y)


maximpact_updated1$predicted_total_intrac0 <- predict(hlm_model_test1, newdata = maximpact_updated1)
maximpact_updated1<-maximpact_updated1%>%
  dplyr::select(citycode,max_return,predicted_total_intrac0,interaction_unevenness(mad),economic-tier_mobility_concentration,rich_neighbor_preference)

maximpact1<-merge(maximpact1,maximpact_updated1,by=c('citycode',"max_return"))

maximpact1$interaction_unevenness(mad)<-(maximpact1$interaction_unevenness(mad).x-maximpact1$interaction_unevenness(mad).y)/maximpact1$interaction_unevenness(mad).x
maximpact1$economic-tier_mobility_concentration<-(maximpact1$economic-tier_mobility_concentration.x-maximpact1$economic-tier_mobility_concentration.y)/maximpact1$economic-tier_mobility_concentration.x
maximpact1$rich_neighbor_preference<-(maximpact1$rich_neighbor_preference.x-maximpact1$rich_neighbor_preference.y)/maximpact1$rich_neighbor_preference.x

maximpact1$netpredit<-maximpact1$predicted_total_intrac-maximpact1$predicted_total_intrac0
maximpact1$netpreditratio<-maximpact1$netpredit/abs(maximpact1$predicted_total_intrac)


maximpact1<-maximpact1%>%
  filter(!is.na(netpreditratio))
summary(maximpact1$interaction_unevenness(mad))
summary(maximpact1$economic-tier_mobility_concentration)
summary(maximpact1$rich_neighbor_preference)

summary(maximpact1$netpreditratio)


a <- sum(-172.6227408 * (maximpact1$interaction_unevenness(mad).x - maximpact1$interaction_unevenness(mad).y) / abs(maximpact1$predicted_total_intrac), na.rm = TRUE)
b <- sum(-0.1470388 * (maximpact1$economic-tier_mobility_concentration.x - maximpact1$economic-tier_mobility_concentration.y) / abs(maximpact1$predicted_total_intrac), na.rm = TRUE)
c <- sum(0.0845488 * (maximpact1$rich_neighbor_preference.x - maximpact1$rich_neighbor_preference.y) / abs(maximpact1$predicted_total_intrac), na.rm = TRUE)



total_contribution <- a + b + c
a_ratio <- a / total_contribution
b_ratio <- b / total_contribution
c_ratio <- c / total_contribution


contribution_data <- data.frame(
  Variable = c("Mean Absolute Deviation", "economic-tier_mobility_concentration", "Dis_ex"),
  Contribution_Ratio = c(a_ratio, b_ratio, c_ratio)
)


p <- ggplot(contribution_data, aes(x = Variable, y = Contribution_Ratio, fill = Variable)) +
  geom_bar(stat = "identity", width = 0.7, color = "black") +
  labs(
    x = "Variable",
    y = "Contribution Ratio",
    title = "Contribution Ratios of Variables"
  ) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +  # 转换为百分比格式
  theme_minimal() +
  theme(
    text = element_text(size = 12),
    plot.title = element_text(hjust = 0.5),
    legend.position = "none"
  )
