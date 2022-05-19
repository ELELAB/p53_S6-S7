library(ggplot2)
library(extrafont)
df <- read.csv("std_for_selected_21.csv", header = TRUE, row.names = 1)
#df <- df[ ,-c(3:5) ]
#df

plot1 <- ggplot(data = df, mapping = aes(x=factor(mutation, levels = mutation), y = ddG)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  theme(aspect.ratio=9/16) +
  theme(text=element_text(family="sans")) +
  geom_bar( aes(x=factor(mutation, levels = mutation), y=ddG), stat="identity", fill="#238A8DFF", alpha=0.7) +
  geom_errorbar( aes(x=factor(mutation, levels = mutation), ymin=ddG-std, ymax=ddG+std), width=0.4, colour="black", alpha=0.9) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(x = "Mutations")

ggsave(plot1, filename = "bar_DNA_mut.pdf", height=6, width=10)
