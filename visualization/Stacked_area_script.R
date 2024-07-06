library(readxl)
library(tidyverse)
library(writexl)

Chr3_file <-"Chr3_data_June2022.xlsx"

data <- read_excel(Chr3_file,
                   sheet = "1µg_mL",
                   col_types = c("text", "text", "numeric", 
                                 "numeric", "numeric", "numeric", 
                                 "numeric", "numeric", "numeric", 
                                 "numeric", "numeric", "numeric", 
                                 "numeric", "numeric"))

glimpse(data)

data_long <- data %>% 
  pivot_longer(cols = -c(Passage, Population), 
               names_to="well_ID", 
               values_to="percentage")

data_long <- data_long %>% 
  mutate(Population = factor(Population),
         Passage = as.integer(Passage),
         well_ID=factor(well_ID, levels = c("Lineage 1","Lineage 2","Lineage 3",
                                            "Lineage 4","Lineage 5","Lineage 6",
                                            "Lineage 7","Lineage 8","Lineage 9",
                                            "Lineage 10","Lineage 11","Lineage 12")))

data_long <- data_long %>% 
  group_by(Passage, well_ID) %>% 
  mutate(total_pop = sum(percentage),
         norm_percent = 100 * (percentage/total_pop))

write_xlsx(data_long, paste0(data.dir,"/Chr3_YPAD.xlsx"))

plot <- data_long %>% 
  ggplot(aes(x=Passage, y=norm_percent, fill = Population)) +
  geom_area(size = .1, color = "white") + 
  xlab("Passage Number") +
  ylab("Population (%)") +
  scale_x_continuous(breaks =c(0, 1, 5, 10, 12, 15), 
                     expand=c(0,0)) +
  ggtitle("1 ug/ml") +
  geom_vline(xintercept = 10, 
             linetype = "dashed",
             color = "black",
             size = 1) +
  scale_fill_manual("Population", values = c("1:1 ratio" = "grey",
                                             "CNV BFP" = "darkcyan",
                                             "CNV GFP" = "yellow",
                                             "GFP Only" = "#55BA2B",
                                             "BFP Only" = "deepskyblue")) +
  theme_classic() +
  theme(axis.text.x=element_text(size=rel(1.5)),
        axis.text.y=element_text(size=rel(1.5)),
        axis.text = element_text(face="bold",color = "black"),
        axis.title.x=element_text(size=rel(2)),
        axis.title.y=element_text(size=rel(2)),
        plot.title=element_text(size=rel(2.5), hjust=.7),
        strip.text.x= element_text(size=12),
        panel.spacing = unit(0.2, "lines")) +
  facet_wrap(~well_ID)