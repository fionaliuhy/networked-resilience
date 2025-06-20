packages <- c("data.table", "tidyverse", "sf", "sp", "spdep",
              "dlnm", "tsModel", "hydroGOF","RColorBrewer", 
              "geofacet", "ggpubr", "ggthemes","extRemes")

# load packages
lapply(packages, library, character.only = TRUE)

options(digits = 8)  
options(scipen = 999)  

#############################Human mobility resilience under extreme rainfall####################

###Provide return-level estimates (e.g., 100-year extreme rainfall) for each city, used in Fig. 1a and Fig. 3
prec<-read_csv("daily_average_precipitation.csv")
prec$mean<-prec1$mean*1000
prec$time <- as.integer(format(prec1$date, "%Y%m%d"))
colnames(prec)<-c('citycode','date','prec','time')

prec <- prec %>%
  distinct(citycode, date, .keep_all = TRUE)

precbase<-prec%>%
  filter(year(date)>=1981&year(date)<=2020)

return_periods <- c(2, 3, 4, 5, 6, 7, 8, 9, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100)

return_levels <- precbase %>%
  group_by(citycode) %>%
  do({
    city_data <- na.omit(.$prec)  
    
    # use 95th
    threshold <- quantile(city_data, 0.95, na.rm = TRUE)
    if (sum(city_data > threshold) < 20) {
      return(data.frame(citycode = unique(.$citycode), return_period = return_periods, return_level = NA))
    }
    fit <- tryCatch({
      fevd(city_data, threshold = threshold, type = "GP")
    }, error = function(e) {
      NULL
    })
    if (is.null(fit)) {
      return(data.frame(citycode = unique(.$citycode), return_period = return_periods, return_level = NA))
    }
    
    return_levels <- sapply(return_periods, function(period) {
      tryCatch(as.numeric(return.level(fit, return.period = period)), error = function(e) NA)
    })
    
    data.frame(citycode = unique(.$citycode), return_period = return_periods, return_level = return_levels)
  }) %>%
  ungroup()


###identify extreme rainfall day
data<- read_csv("data.csv")

events0 <- data %>%
  filter((prec>=50|returnperiod>=3)) %>%
  filter(year==2023)

##identify extreme rainfall event
events1 <- events0 %>%
  arrange(citycode, date) %>%  
  group_by(citycode) %>%
  mutate(
    date_diff = date - lag(date, default = first(date) - 1), 
    group_id = cumsum(date_diff > 1)
  ) %>%
  group_by(citycode, group_id) %>%
  slice_max(prec, n = 1, with_ties = FALSE) %>% 
  ungroup() %>%
  dplyr::select(-date_diff, -group_id) 

adjacency_distance <- read_csv("adjacency_distance.csv") ##cities' distant to neighbours
adjacency_distance <- adjacency_distance %>%
  distinct(ocitycode, adjacent_codes) 
adjacency_map <- adjacency_distance %>%
  separate_rows(adjacent_codes, sep = ",") %>%
  mutate(adjacent_codes = as.numeric(adjacent_codes)) %>%
  rename(citycode = ocitycode, adjacent_citycode = adjacent_codes)

dis <- read_csv("distance.csv")##cities' distant to all other cities

events0 <- events0 %>%
  arrange(citycode, date) %>%
  group_by(citycode) %>%
  mutate(
    event1 = cumsum(c(1, diff(date)) > 1)  
  ) %>%
  ungroup()

# Process adjacency relationship data
adjacency_map <- adjacency_distance %>%
  separate_rows(adjacent_codes, sep = ",") %>%
  mutate(adjacent_citycode = as.numeric(adjacent_codes)) %>%
  rename(citycode = ocitycode)


events0$event2 <- NA_integer_
event_counter <- 1

# Based on temporal and spatial continuity to identiy the same event
for (i in seq_len(nrow(events0))) {
  if (is.na(events0$event2[i])) {
    current_event1 <- events0$event1[i]
    current_city <- events0$citycode[i]
    
    same_event1_records <- events0 %>%
      filter(
        citycode == current_city,
        event1 == current_event1
      )
    
    queue <- same_event1_records %>%
      dplyr::select(date, citycode) %>%
      distinct()
    
    checked <- data.frame(date = as.Date(character()), citycode = numeric())
    
    while (nrow(queue) > 0) {
      current <- queue[1, ]
      queue <- queue[-1, ]
      
      if (nrow(checked %>% filter(date == current$date & citycode == current$citycode)) > 0) {
        next
      }
      
      checked <- rbind(checked, current)
      
      current_city_event1 <- events0 %>%
        filter(citycode == current$citycode) %>%
        filter(event1 == events0$event1[events0$citycode == current$citycode & events0$date == current$date])
      
      events0$event2[events0$citycode %in% current_city_event1$citycode & 
                       events0$event1 %in% current_city_event1$event1] <- event_counter
      
      # Find neighboring cities in the same event
      adjacent_cities <- adjacency_map %>%
        filter(citycode == current$citycode) %>%
        pull(adjacent_citycode)
      
      # Find cities within 200 km in the same event
      close_cities <- dis %>%
        filter(ocitycode == current$citycode & distance <= 200) %>%
        pull(dcitycode)
      
      related_cities <- unique(c(adjacent_cities, close_cities))
      
      for (related_city in related_cities) {
        related_records <- events0 %>%
          filter(
            citycode == related_city,
            date == current$date,
            is.na(event2)
          )
        
        if (nrow(related_records) > 0) {
          related_event1_records <- events0 %>%
            filter(
              citycode == related_city,
              event1 == related_records$event1[1],
              is.na(event2)
            ) %>%
            dplyr::select(date, citycode) %>%
            distinct()
          
          queue <- rbind(queue, related_event1_records)
          queue <- unique(queue)
        }
      }
    }
    
    event_counter <- event_counter + 1
  }
}


###fig 1a, mapped in arcgis
times <- events0 %>%
  group_by(citycode) %>%
  summarize(event_count = n(), .groups = "drop")

### number of extreme-rainfall day in the same city-event
event_day_count <- events0 %>%
  group_by(event2,citycode) %>%
  summarize(event_day_count = n(), .groups = "drop")

### number of cities affected in the same city-event
event_counts <- events0 %>%
  group_by(event2) %>%
  summarise(unique_city_count = n_distinct(citycode), .groups = "drop")

com<-merge(events0,event_counts,by='event2')

##fig 1b
event2_counts <- com %>%
  group_by(unique_city_count) %>%
  summarize(
    event2_count = n_distinct(event2),   
    row_count = n(),     
    .groups = "drop"
  )


library(extrafont)
font_import(pattern = "Arial")
loadfonts(device = "pdf")

p1 <- ggplot(event2_counts, aes(x = unique_city_count, y = event2_count)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  labs(
    x = "Number of affected cities in the same event",
    y = "Number of event"
  ) +
  theme_classic() +
  theme(
    text = element_text(family = "Arial", size = 7),
    axis.title = element_text(size = 6.5),
    axis.text = element_text(size = 6.5),  
    plot.title = element_text(hjust = 0.5)
  )+coord_flip()

p2 <- ggplot(event2_counts, aes(x = unique_city_count, y = row_count)) +
  geom_bar(stat = "identity", fill = "darkorange") +
  labs(
    x = "Number of affected cities in the same event",
    y = "Number of city-day"
  ) +
  theme_classic() +
  theme(
    text = element_text(family = "Arial", size = 7),
    axis.title = element_text(size = 6.5),
    axis.text = element_text(size = 6.5),  
    plot.title = element_text(hjust = 0.5)
  )+coord_flip()


###fig. 1c
events0<-merge(events0,return_levels3,by='citycode',all.x=T)

comx<-com[,c('citycode','date','event2','unique_city_count')]

event_days <- events0 %>%
  group_by(citycode, event2) %>%
  summarize(extreme_day_count = n(), .groups = "drop")

ggplot(event_days, aes(x = extreme_day_count)) +
  geom_histogram(aes(y = after_stat(density)), binwidth = 1, color = "black", fill = "lightblue") +
  labs(x = "Number of extreme-rainfall days", y = "Number of city-event")+
  theme_classic() +
  theme(
    text = element_text(family = "Arial", size = 7),
    axis.title = element_text(size = 6.5),
    axis.text = element_text(size = 6.5),  
    plot.title = element_text(hjust = 0.5))+
  coord_flip()
    
freq_table <- event_days %>%
  count(extreme_day_count) %>%
  mutate(freq = n / sum(n))

###fig 1d --get the baselines and calculate the relative changes of intracity mobility intensity
events <- events0 %>%
  filter((prec>=50|returnperiod>=3)&holiday==0) %>%
  filter(year==2023) %>%
  arrange(citycode, date) %>% 
  group_by(citycode) %>%
  mutate(
    event_group = cumsum(c(1, diff(as.Date(date)) > 1))
  ) %>%
  group_by(citycode, event_group) %>%
  summarise(
    date = min(date), 
    returnperiod = first(returnperiod), 
    prec=max(prec),
    event2=event2,
    .groups = "drop"
  ) %>%
  mutate(eventid = row_number()) 


data1<-merge(data,events, by=c('citycode',"date","returnperiod"),all.x=T)


data0 <- data1 %>%
  mutate(date = as.Date(date)) %>%
  arrange(citycode, date) %>% 
  filter(!is.na(eventid)) %>%
  group_by(citycode, eventid,event2) %>%
  summarise(
    start_date = min(date),  
    .groups = "drop"
  ) %>%
  arrange(citycode, start_date) %>%  
  group_by(citycode) %>%
  mutate(
    next_start_date = lead(start_date), 
    gap = as.integer(next_start_date - start_date)  
  ) %>%
  rowwise() %>%
  mutate(
    range = ifelse(
      is.na(next_start_date) | gap > 28, 
      list(seq(start_date - 14, start_date + 14, by = "day")),
      list(seq(start_date - 14, start_date + floor(gap / 2) - 1, by = "day"))
    )
  ) %>%
  unnest(range) %>%
  rename(date = range) %>%
  ungroup() %>%
  mutate(date = as.Date(date)) %>%
  left_join(
    data %>% mutate(date = as.Date(date)),  
    by = c("citycode", "date")
  ) %>%
  group_by(citycode, eventid) %>%
  mutate(
    exposure = as.integer(date - start_date)  
  ) %>%
  ungroup()


base <- data0 %>%
  filter(prec <= 2 & exposure < 0 & holiday == 0 ) %>%
  group_by(citycode, eventid,event2) %>%
  summarise(
    windbase=mean(wind_speed),
    tempbase=mean(max_temp),
    base_day_6_and_0 = ifelse(
      !is.na(mean(intra[day %in% c(6, 0)], na.rm = TRUE)),
      mean(intra[day %in% c(6, 0)], na.rm = TRUE),
      mean(intra[day %in% 1:5], na.rm = TRUE)
    ),
    base_day_1_to_5 = ifelse(
      !is.na(mean(intra[day %in% 1:5], na.rm = TRUE)),
      mean(intra[day %in% 1:5], na.rm = TRUE),
      mean(intra[day %in% c(6, 0)], na.rm = TRUE)
    )
  )


data0 <- data0 %>%
  mutate(
    selected_base = ifelse(day %in% 1:5, "base_day_1_to_5",
                           ifelse(day %in% c(6, 0), "base_day_6_and_0", NA))
  )

data0 <- data0 %>%
  left_join(base %>% pivot_longer(cols = starts_with("base_day"), 
                                  names_to = "selected_base", 
                                  values_to = "base_value"), 
            by = c("citycode", "eventid", "selected_base"))


data0$intrac<-data0$intra/data0$base_value-1

data0<-merge(data0,event_day_count,by.x=c('citycode','event2.x'),by.y=c('citycode','event2'))
data0<-merge(data0,event_counts,by.x=c('event2.x'),by.y=c('event2'),all.x=T)

avg_all <- data0 %>%
  group_by(exposure) %>%
  summarise(
    mean_intrac = mean(intrac, na.rm = TRUE),
    q25 = quantile(intrac, 0.25, na.rm = TRUE),
    q75 = quantile(intrac, 0.75, na.rm = TRUE)
  )


avg_d1 <- data0 %>%
  filter(extreme_day_count == 1) %>%
  group_by(exposure) %>%
  summarise(
    mean_intrac = mean(intrac, na.rm = TRUE),
    q25 = quantile(intrac, 0.25, na.rm = TRUE),
    q75 = quantile(intrac, 0.75, na.rm = TRUE)
  )


avg_d2 <- data0 %>%
  filter(extreme_day_count == 2) %>%
  group_by(exposure) %>%
  summarise(
    mean_intrac = mean(intrac, na.rm = TRUE),
    q25 = quantile(intrac, 0.25, na.rm = TRUE),
    q75 = quantile(intrac, 0.75, na.rm = TRUE)
  )

p3 <- ggplot(data0, aes(x = exposure, y = intrac, group = interaction(citycode, eventid))) +
  geom_line(color = "blue", alpha = 0.05) 
  
  # 所有数据平均线 + 阴影
  geom_ribbon(data = avg_all, aes(x = exposure, ymin = q25, ymax = q75), 
              fill = "gray70", alpha = 0.3, inherit.aes = FALSE) +
  geom_line(data = avg_all, aes(x = exposure, y = mean_intrac), 
            color = "black", size = 0.6, inherit.aes = FALSE) +
  
  # extreme_day_count == 1
  geom_ribbon(data = avg_d1, aes(x = exposure, ymin = q25, ymax = q75), 
              fill = "green", alpha = 0.3, inherit.aes = FALSE) +
  geom_line(data = avg_d1, aes(x = exposure, y = mean_intrac), 
            color = "green4", size = 0.6, inherit.aes = FALSE) +
  
  # extreme_day_count == 2
  geom_ribbon(data = avg_d2, aes(x = exposure, ymin = q25, ymax = q75), 
              fill = "orange", alpha = 0.3, inherit.aes = FALSE) +
  geom_line(data = avg_d2, aes(x = exposure, y = mean_intrac), 
            color = "orange3", size = 0.6, inherit.aes = FALSE) +
  
  labs(
    x = "Exposure (day)",
    y = "Baseline (%)"
  ) +
  xlim(-7, 7) +
  ylim(-0.52, 0.27) +
  theme_bw() +
  theme(
    text = element_text(family = "Arial", size = 7),
    axis.title = element_text(size = 6.5),
    axis.text = element_text(size = 6.5),
    plot.title = element_text(hjust = 0.5)
  )

print(p3)

####fig 1e quantify the human mobility resilience 
data2<-data0[data0$exposure>=-3,]

result <- data2 %>%
  group_by(citycode, event2.x) %>%  
  summarise(
    base = mean(base_value, na.rm = TRUE),
    
    date_exposure_zero = min(date[exposure == 0], na.rm = TRUE),
    
    date_before = min(date[exposure <= 0 & exposure >= -3 & intrac <= 0], na.rm = TRUE),

    date_after = min(date[exposure >= 0 & exposure <= 14 & intrac < 0], na.rm = TRUE),
    
    period = ifelse(!is.na(date_before) & !is.na(date_after),
                    as.numeric(date_after - date_before), NA),
    
    total_intrac = ifelse(!is.na(date_before) & !is.na(date_after),
                          sum(intrac[date >= date_before & date <= date_after], na.rm = TRUE),
                          NA),

    total_prec = ifelse(!is.na(date_before) & !is.na(date_after),
                        sum(prec[date >= date_before & date <= date_after], na.rm = TRUE),
                        NA),
    
    max_prec = ifelse(!is.na(date_before) & !is.na(date_after),
                      max(prec[date >= date_before & date <= date_after], na.rm = TRUE),
                      NA),
    
    max_Rpolicy = ifelse(!is.na(date_before) & !is.na(date_after),
                         max(Rpolicy[date >= date_before & date <= date_after], na.rm = TRUE),
                         NA),
    
    max_Cpolicy = ifelse(!is.na(date_before) & !is.na(date_after),
                         max(Cpolicy[date >= date_before & date <= date_after], na.rm = TRUE),
                         NA),

    max_return = ifelse(!is.na(date_before) & !is.na(date_after),
                        max(returnperiod[date >= date_before & date <= date_after], na.rm = TRUE),
                        NA),

    max_temp_diff = ifelse(!is.na(date_before) & !is.na(date_after),
                           max(max_temp[date >= date_before & date <= date_after] - 
                                 tempbase[date >= date_before & date <= date_after], na.rm = TRUE),
                           NA),
    
    max_wind_diff = ifelse(!is.na(date_before) & !is.na(date_after),
                           max(wind_speed[date >= date_before & date <= date_after] - 
                                 windbase[date >= date_before & date <= date_after], na.rm = TRUE),
                           NA)
  ) %>%
  ungroup()


overall<-merge(overall,event_day_count,by=c('citycode',"event2.x"),by.y=c('citycode','event2'),all.x=TRUE)
overall<-merge(overall,event_counts,by=c('event2.x'),by.y=c('event2'),all.x=TRUE)


overall_long <- bind_rows(
  overall %>% mutate(group = "All"),
  overall %>% filter(event_day_count == 1) %>% mutate(group = "1 day"),
  overall %>% filter(event_day_count > 1) %>% mutate(group = ">1 day"),
  overall %>% filter(event_day_count == 1 & unique_city_count > 1) %>% mutate(group = "1 day_>1 city"),
  overall %>% filter(event_day_count == 1 & unique_city_count <= 1) %>% mutate(group = "1 day_1 city"),
  overall %>% filter(event_day_count == 1 & unique_city_count > 2) %>% mutate(group = "1 day_>2 city"),
  overall %>% filter( event_day_count== 1 & unique_city_count <= 2) %>% mutate(group = "1 day_1-2 city")
)

p4<-ggplot(overall_long, aes(x = group, y = total_intrac, fill = group)) +
  geom_half_violin(side = "r", alpha = 0.5, color = NA,width = 1.0) +
  geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.8) +
  stat_summary(fun = mean, geom = "point", shape = 21, size = 1.5, fill = "red", color = "black") +
  labs(x = "", y = "Human mobility resilience") +
  ylim(-0.55, 0) +
  theme_bw() +
  theme(
    text = element_text(family = "Arial", size = 6.5),
    axis.text.x = element_text(size = 6.5),
    legend.position = "none"
  )

summary(overall$total_intrac)
group1 <- overall$total_intrac[overall$unique_city_count>1]
group2 <- overall$total_intrac[overall$unique_city_count<=1]

group1 <- overall$total_intrac[overall$event_day_count==1]
group2 <- overall$total_intrac[overall$event_day_count>1]


group1 <- overall$total_intrac[overall$event_day_count==1&overall$unique_city_count <=2]
group2 <- overall$total_intrac[overall$event_day_count==1&overall$unique_city_count >2]

group1 <- overall$max_prec[overall$event_day_count==1]
group2 <- overall$max_prec[overall$event_day_count>1]
####p=0.0158
group1 <- overall$max_prec[overall$event_day_count==1&overall$unique_city_count >1]
group2 <- overall$max_prec[overall$event_day_count==1&overall$unique_city_count <=1]
####p= 0.000034328
group1 <- overall$max_prec[overall$event_day_count==1&overall$unique_city_count >2]
group2 <- overall$max_prec[overall$event_day_count==1&overall$unique_city_count <=2]
#p-value = 0.0000033043

wilcox_test_result <- wilcox.test(group1, group2)
print(wilcox_test_result)

##fig s1
p5<-ggplot(overall_long, aes(x = group, y =max_prec , fill = group)) +
  geom_half_violin(side = "r", alpha = 0.5, color = NA,width = 1.0) +
  geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.8) +
  stat_summary(fun = mean, geom = "point", shape = 21, size = 1.5, fill = "red", color = "black") +
  labs(x = "", y = "Total Intracity Mobility") +
  theme_bw() +
  theme(
    text = element_text(family = "Arial", size = 6.5),
    axis.text.x = element_text(size = 6.5),
    legend.position = "none"
  )


#############################Intercity connectivity affects resilience####################


data <- read_csv("2023.csv")  
data0<-data

###zhengzhou single-city event   
od <- data0 %>%
  filter(ocitycode %in% c(410100) | dcitycode %in% c(410100))

od <- od %>%
  mutate(mobility = replace_na(mobility, 0))

base_410100 <- od %>%
  filter(dcitycode == 410100, time >= "20230722" & time <="20230804") %>%
  group_by(ocitycode)%>%
  summarize(base = mean(mobility),
            normalratio=base/704.53805)

#sum(base_410100$base) =857.31936
in_410100s <- od %>%
  filter(dcitycode == 410100, time >= "20230805" & time <= "20230807") %>%
  group_by(ocitycode) %>%
  summarize(avg_mobility = mean(mobility), .groups = "drop") %>%
  left_join(base_410100, by = "ocitycode") %>%
  mutate(
    mobility_change = avg_mobility - base,
    changeratio=mobility_change/857.31936
  )

in_410100s$normalratio<-in_410100s$normalratio*100   ###fig 2a mapped in arcgis
in_410100s$changeratio<-in_410100s$changeratio*100   ###fig 2b mapped in arcgis

cor.test(in_410100s$changeratio,in_410100s$normalratio)


###zhengzhou multi-city event   
base_410100 <- od %>%
  filter(dcitycode == 410100, time >= "20230812" & time <="20230825") %>%
  group_by(ocitycode)%>%
  summarize(base = mean(mobility),
            normalratio=base/704.53805)

in_410100m <- od %>%
  filter(dcitycode == 410100, time >= "20230826" & time <= "20230828") %>%
  group_by(ocitycode) %>%
  summarize(avg_mobility = mean(mobility), .groups = "drop") %>%
  left_join(base_410100, by = "ocitycode") %>%
  mutate(
    mobility_change = avg_mobility - base,
    changeratio=mobility_change/704.53805
  )

in_410100m$normalratio<-in_410100m$normalratio*100  
in_410100m$changeratio<-in_410100m$changeratio*100 ###fig 2c mapped in arcgis


event17_cities <- events %>%
  filter(event2 == 17) %>%
  pull(citycode) %>%
  unique()

# see wether connected cities are affected in the same event
in_410100m <- in_410100m %>%
  mutate(group = ifelse(ocitycode %in% event17_cities, "In Event", "Not in Event"))

cor.test(in_410100m$changeratio[in_410100m$group=="In Event"],in_410100m$normalratio[in_410100m$group=="In Event"])
cor.test(in_410100m$changeratio[in_410100m$group=="Not in Event"],in_410100m$normalratio[in_410100m$group=="Not in Event"])

combined_data <- bind_rows(
  in_410100m %>% mutate(source = "in_410100m"),
  in_410100s %>% mutate(source = "in_410100s")
)

# fig 2g
ggplot(combined_data, aes(x = normalratio, y = changeratio)) +
  geom_point(data = subset(combined_data, source == "in_410100s"), 
             color = "#1b7837", size = 2, alpha = 0.6) + 
  geom_point(data = subset(combined_data, source == "in_410100m" & group == "In Event 17 Cities"), 
             color = "#92c5de", size = 2, alpha = 0.6) + 
  geom_point(data = subset(combined_data, source == "in_410100m" & group == "Not in Event 17 Cities"), 
             color = "#0571b0", size = 2, alpha = 0.6) + 
  geom_smooth(data = subset(combined_data, source == "in_410100s"), 
              method = "lm", se = TRUE, color = "#1b7837") + 
  geom_smooth(data = subset(combined_data, source == "in_410100m" & group == "In Event 17 Cities"), 
              method = "lm", se = TRUE, color = "#92c5de") + 
  geom_smooth(data = subset(combined_data, source == "in_410100m" & group == "Not in Event 17 Cities"), 
              method = "lm", se = TRUE, color = "#0571b0") +  
  labs(
    x = "normalratio", 
    y = "changeratio"
  ) +
  theme_bw()

t.test(in_410100s$changeratio, in_410100m$changeratio, paired = FALSE)
wilcox.test(in_410100s$changeratio, in_410100m$changeratio, paired = FALSE)

df <- data.frame(
  x = in_410100s$normalratio,
  y = in_410100m$normalratio
)

##fig s2
figs21<-ggplot(df, aes(x = x, y = y)) +
  geom_point() +  # 散点图
  geom_smooth(method = "lm", se = FALSE, color = "blue") +  # 拟合线
  labs(x = "Baseline inflow in single-city event",
       y = "Baseline inflow in multi-city event",
       title = "zhengzhou") +
  theme_bw()+
  theme(
    text = element_text(size = 7)  
  )

cor.test(in_410100s$normalratio,in_410100m$normalratio)
###r=0.99924051 ,p-value < 0.000000000000000222

###others cities are the same with above --fig. extended data fig.1

#######regression analysis --table 1, extended data table 2 & table s3

overall<-overall%>%
  filter(year==2023)

overall$schoolnum[is.na(overall$schoolnum)]<-0
overall$hosnum[is.na(overall$hosnum)]<-0
overall$gdppercapita<-as.numeric(overall$gdppercapita)


library(glmmTMB)
library(glmnet)
library(glmnetUtils)
library(olsrr)
library(corrplot)

####only city-level factors
X0 <- as.matrix(overall[, c("total_prec", "max_temp_diff", "max_wind_diff", "mean_elevation","mean_slope","slope_dif",
                            "urbanization", "old65",  "primarygdp","tergdp", "pop", "gdppercapita",
                           "proportion_10","proportion_30","proportion_50","proportion_60",
                           "proportion_80","proportion_90","ndvi","roaddes","roadleg","hosnum","schoolnum",
                           "damdes", "max_Rpolicy", "max_Cpolicy")])
y <- overall$total_intrac
library(MASS)

data <- cbind(y = y, X0)  

full_model <- lm(y ~ ., data = as.data.frame(data))
null_model <- lm(y ~ 1, data = as.data.frame(data))

#stepwise
stepwise_model <- stepAIC(full_model, direction = "both", trace = FALSE)

selected_vars <- names(coef(stepwise_model))[-1]
print(selected_vars)

##ols
best_ols_model <- lm(y ~ ., data = as.data.frame(data[, c("y", selected_vars)]))
summary(best_ols_model) 

best_data <- data[, c("y", selected_vars2)]  
correlation_matrix_best <- cor(best_data)  

# correlation
corrplot(correlation_matrix_best, method = "circle", type = "upper", tl.col = "black", tl.srt = 45)


vars <- overall %>%
  dplyr::select(total_prec, max_wind_diff, gdppercapita, mean_elevation,
                proportion_10, proportion_60, proportion_90, damdes,
                max_Rpolicy, max_Cpolicy)
cor_matrix <- cor(vars, use = "complete.obs")

par(family = "Arial") 

corrplot(
  cor_matrix,
  method = "color",
  type = "upper",
  tl.col = "black",
  tl.cex = 6.5 / 12,
  addCoef.col = "black"
)

# fig s9
model <- lm(100 * total_prec ~ max_prec, data = overall)
summary(model)
coef <- coefficients(model)
intercept <- round(coef[1], 2)
slope <- round(coef[2], 2)

r2 <- round(summary(model)$r.squared, 3)

figs9 <- ggplot(overall, aes(x = max_prec, y = 100 * total_prec)) +
  geom_point(alpha = 0.6, shape = 16, size = 1) +  
  geom_smooth(method = "lm", color = "blue", se = TRUE) +
  labs(x = "Max Daily Rainfall (mm)", y = "Cumulative Rainfall (%)") +
  theme_bw() +
  annotate("text", 
           x = Inf, y = -Inf, hjust = 1.1, vjust = -1,
           label = paste0("y = ", slope, "x + ", intercept, 
                          "\nR² = ", r2),
           size = 6.5)


library(lme4)
library(lmerTest)
library(QuantPsyc)
##overall is the data for each city-event in 2023
overall$tergdp<-overall$tergdp/100
overall$total_prec<-overall$total_prec/100
overall$mean_elevation<-overall$mean_elevation/1000

###table 1
hlm_model_test1 <- lm(
  total_intrac ~ total_prec + max_wind_diff + gdppercapita+mean_elevation +proportion_10+proportion_60+
    proportion_90 + damdes + max_Rpolicy + max_Cpolicy ,
  data = overall)

summary(hlm_model_test1)
hlm_summary <- summary(hlm_model_test1)
fixed_effects <- as.data.frame(hlm_summary$coefficients)


hlm_model_test2 <- lmer(
  total_intrac ~ total_prec + max_wind_diff +  mean_elevation + tergdp+proportion_10+proportion_60+
    proportion_90 + damdes + max_Rpolicy + max_Cpolicy  + (1|prov_code),
  data = overall
)

fixed_effects <- as.data.frame(hlm_summary$coefficients)
random_effects <- as.data.frame(hlm_summary$varcor)


hlm_model_test3 <- lmer(
  total_intrac ~ total_prec + max_wind_diff +  mean_elevation + tergdp+proportion_10+proportion_60+
    proportion_90 + damdes + max_Rpolicy + max_Cpolicy  + (1|unique_city_count),
  data = overall
)

fixed_effects <- as.data.frame(hlm_summary$coefficients)
random_effects <- as.data.frame(hlm_summary$varcor)


library(MuMIn)
r2_values <- r.squaredGLMM(hlm_model_test1)
print(r2_values)
r2_values <- r.squaredGLMM(hlm_model_test2)
print(r2_values)
r2_values <- r.squaredGLMM(hlm_model_test3)
print(r2_values)

predictions1 <- predict(hlm_model_test1, newdata = overall)
actuals1 <- overall$total_intrac
rmse1 <- sqrt(mean((predictions1 - actuals1)^2))
print(paste("RMSE for Model 1 (Mixed Model):", rmse1))
predictions2 <- predict(hlm_model_test2, newdata = overall, re.form = NA) 
rmse2 <- sqrt(mean((predictions2 - actuals1)^2))
print(paste("RMSE for Model 2 (Mixed Model):", rmse2))
predictions3 <- predict(hlm_model_test3, newdata = overall, re.form = NA) 
rmse3 <- sqrt(mean((predictions3 - actuals1)^2))
print(paste("RMSE for Model 3 (Mixed Model):", rmse3))


####add intercity connectivity
cor.test(overall$total_intrac,overall$logweightedeco)

X <- as.matrix(overall[, c("total_prec", "max_wind_diff", "mean_elevation", "tergdp"  ,"proportion_10","proportion_60","proportion_90",
                           "damdes", "max_Rpolicy", "max_Cpolicy", 
                           "Per_capita_interaction",
                           "betweenness", "eigenvector", "pagerank", "closeness_distance",
                           "Mobility-weighted_distance","Near-neighbor preference","withinprovince",'withneighbor',
                           "Mobility-weighted_GDP_percapita","rich_neighbor_preference",'Economic_status_relative_to_neighbors',
                           'economic-tier_mobility_concentration',"max_ratio","interaction_unevenness(mad)","normalized_entropy","gini")])


y <- overall$total_intrac
library(MASS)

data <- cbind(y = y, X) 
full_model <- lm(y ~ ., data = as.data.frame(data))
null_model <- lm(y ~ 1, data = as.data.frame(data))

# stepwise
stepwise_model <- stepAIC(full_model, direction = "both", trace = FALSE)

selected_vars <- names(coef(stepwise_model))[-1]  # 除去截距
print(selected_vars)

##ols
best_ols_model <- lm(y ~ ., data = as.data.frame(data[, c("y", selected_vars)]))
summary(best_ols_model) 

best_data <- data[, c("y", selected_vars2)]  ##factors in the ols model

###fig s8
cor_matrix <- cor(data[, selected_vars2], use = "pairwise.complete.obs")

print(cor_matrix)
library(corrplot)
corrplot(cor_matrix, method = "color", type = "upper", order = "hclust",
         tl.col = "black", tl.cex = 0.8, addCoef.col = "black")

corrplot(
  cor_matrix,
  method = "color",
  type = "upper",
  tl.col = "black",
  tl.cex = 6.5 / 12,
  addCoef.col = "black"
)


hlm_model_test1 <- lm(
  total_intrac ~ total_prec + max_wind_diff  + mean_elevation + proportion_10+proportion_60+
    proportion_90 + damdes + max_Rpolicy + Per_capita_interaction+max_Cpolicy+ interaction_unevenness(mad)+economic-tier_mobility_concentration+rich_neighbor_preference,
  data = overall)
summary(hlm_model_test1)

hlm_summary <- summary(hlm_model_test1)
fixed_effects <- as.data.frame(hlm_summary$coefficients)
random_effects <- as.data.frame(hlm_summary$varcor)


hlm_model_test2 <- lmer(
  total_intrac ~ total_prec + max_wind_diff  + mean_elevation + proportion_10+proportion_60+
    proportion_90 + damdes + max_Rpolicy + max_Cpolicy +Per_capita_interaction+ interaction_unevenness(mad)+economic-tier_mobility_concentration+rich_neighbor_preference +(1 | prov_code),
  data = overall
)
summary(hlm_model_test2)
hlm_summary <- summary(hlm_model_test2)
fixed_effects <- as.data.frame(hlm_summary$coefficients)
random_effects <- as.data.frame(hlm_summary$varcor)


hlm_model_test3 <- lmer(
  total_intrac ~ total_prec + max_wind_diff  + mean_elevation + proportion_10+proportion_60+
    proportion_90 + damdes + max_Rpolicy + max_Cpolicy +Per_capita_interaction+ interaction_unevenness(mad)+economic-tier_mobility_concentration+rich_neighbor_preference +(1|unique_city_count),
  data = overall
)

summary(hlm_model_test3)
hlm_summary <- summary(hlm_model_test3)
fixed_effects <- as.data.frame(hlm_summary$coefficients)
random_effects <- as.data.frame(hlm_summary$varcor)


r2_values <- r.squaredGLMM(hlm_model_test1)
print(r2_values)
r2_values <- r.squaredGLMM(hlm_model_test2)
print(r2_values)
r2_values <- r.squaredGLMM(hlm_model_test3)
print(r2_values)


predictions1 <- predict(hlm_model_test1, newdata = overall)
actuals1 <- overall$total_intrac
rmse1 <- sqrt(mean((predictions1 - actuals1)^2))
print(paste("RMSE for Model 1 (Mixed Model):", rmse1))
predictions2 <- predict(hlm_model_test2, newdata = overall, re.form = NA) 
rmse2 <- sqrt(mean((predictions2 - actuals1)^2))
print(paste("RMSE for Model 2 (Mixed Model):", rmse2))
predictions3 <- predict(hlm_model_test3, newdata = overall, re.form = NA) 
rmse3 <- sqrt(mean((predictions3 - actuals1)^2))
print(paste("RMSE for Model 3 (Mixed Model):", rmse3))

###extended data table 2

overall1<-overall[overall$unique_city_count<=2,]
hlm_model_test1 <- lm(
  total_intrac ~ total_prec + max_wind_diff  + mean_elevation + tergdp+proportion_10+proportion_60+
    proportion_90 + damdes + max_Rpolicy + max_Cpolicy ,
  data = overall1)
summary(hlm_model_test1)
summary(hlm_model_test1)

hlm_summary <- summary(hlm_model_test1)
fixed_effects <- as.data.frame(hlm_summary$coefficients)
random_effects <- as.data.frame(hlm_summary$varcor)


hlm_model_test2 <- lm(
  total_intrac ~ total_prec + max_wind_diff  + mean_elevation + proportion_10+proportion_60+
    proportion_90 + damdes + max_Rpolicy + max_Cpolicy +Per_capita_interaction+ interaction_unevenness(mad)+economic-tier_mobility_concentration+rich_neighbor_preference,
  data = overall1)
summary(hlm_model_test2)
hlm_summary <- summary(hlm_model_test2)
fixed_effects <- as.data.frame(hlm_summary$coefficients)
random_effects <- as.data.frame(hlm_summary$varcor)


predictions1 <- predict(hlm_model_test1, newdata = overall1)
actuals1 <- overall1$total_intrac
rmse1 <- sqrt(mean((predictions1 - actuals1)^2))
print(paste("RMSE for Model 1 (Mixed Model):", rmse1))
predictions2 <- predict(hlm_model_test2, newdata = overall1, re.form = NA) 
rmse2 <- sqrt(mean((predictions2 - actuals1)^2))
print(paste("RMSE for Model 2 (Mixed Model):", rmse2))

overall1<-overall[overall$unique_city_count>2,]
hlm_model_test1 <- lm(
  total_intrac ~ total_prec + max_wind_diff  + mean_elevation + tergdp+proportion_10+proportion_60+
    proportion_90 + damdes + max_Rpolicy + max_Cpolicy ,
  data = overall1)
summary(hlm_model_test1)
summary(hlm_model_test1)

hlm_summary <- summary(hlm_model_test1)
fixed_effects <- as.data.frame(hlm_summary$coefficients)
random_effects <- as.data.frame(hlm_summary$varcor)


hlm_model_test2 <- lm(
  total_intrac ~ total_prec + max_wind_diff  + mean_elevation + proportion_10+proportion_60+
    proportion_90 + damdes + max_Rpolicy + max_Cpolicy +Per_capita_interaction+ interaction_unevenness(mad)+economic-tier_mobility_concentration+rich_neighbor_preference,
  data = overall1)
summary(hlm_model_test2)
hlm_summary <- summary(hlm_model_test2)
fixed_effects <- as.data.frame(hlm_summary$coefficients)
random_effects <- as.data.frame(hlm_summary$varcor)


predictions1 <- predict(hlm_model_test1, newdata = overall1)
actuals1 <- overall1$total_intrac
rmse1 <- sqrt(mean((predictions1 - actuals1)^2))
print(paste("RMSE for Model 1 (Mixed Model):", rmse1))
predictions2 <- predict(hlm_model_test2, newdata = overall1, re.form = NA) 
rmse2 <- sqrt(mean((predictions2 - actuals1)^2))
print(paste("RMSE for Model 2 (Mixed Model):", rmse2))


