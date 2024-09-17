library(tidyverse)


metadata<-read.csv("~/Desktop/NoyesLab/OREI_Calfs/simple_cow_metadata.csv")

metadata$birth<-ifelse(metadata$DIM>=0,"Post birth","Pre birth")
View(metadata)
dim(metadata)
cow_group<-metadata %>% group_by(CowId,birth) %>%
  mutate(group = cur_group_id() ) %>% ungroup()

cow_group_filter<-cow_group %>% group_by(group) %>%
  summarize(count=n()) %>% filter(count>1) %>% 
  left_join(cow_group,by="group") %>% select(-count)


write.csv(cow_group_filter,"~/Desktop/NoyesLab/OREI_Calfs/cow_coassembly_groups.csv")

View(cow_group)
View(cow_group_filter)

