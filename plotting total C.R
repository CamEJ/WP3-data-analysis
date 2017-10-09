
# plotting total C
pd = read.table("totCN.txt", sep = "\t", header=T)
head(pd)



AvsDiv <- pd %>% 
  group_by(timepoint, Treatment) %>%
  summarise(Nav = mean(Nitrogen..), 
            Cav = mean(Carbon..),
            Csdv = sd(Carbon..),
            total = n())
AvsDiv[["timepoint"]] <- setFactorOrder(AvsDiv[["timepoint"]], c("T0", "T2", "T3",  "T5","T7", "T8", "T11", "T13"))

# order of treatments. 
AvsDiv[["Treatment"]] <- setFactorOrder(AvsDiv[["Treatment"]], c("Control", "Slurry", "Flood", "Flood+Slurry"))



#### and now plot :
p <- ggplot(AvsDiv, aes(x=timepoint, y=Cav, ymin=Cav-Csdv, ymax=Cav+Csdv))+
  geom_pointrange(size=0.75)+
  #geom_hline(yintercept = 204, linetype=2)+ # y intercept= Av of all control resp
  coord_flip() + #ylim = c(0, 1000)
  xlab('Time Point')+ 
  geom_point(colour="black", shape=21, size = 10) +
  aes(fill = factor(Treatment)) + 
  scale_fill_manual(values=c("black", "chocolate4", "slateblue", "olivedrab"))
p

p1 = p + ylab(expression(paste('% carbon' )))


# legend on side with heading
p2 <- p1 + theme_bw() + 
  theme(axis.text.y=element_text(size=16, colour="black"),
        axis.text.x=element_text(size=16, colour="black"),
        axis.title.x=element_text(size=16, colour="black"),
        legend.text = element_text(size=16),
        legend.title = element_text(size=17),
        axis.title.y=element_blank())+
  labs(fill="    Treatment") 
p2


# put legend inside plot 
PDRplot = p2 + theme(legend.justification = c(1, 1), legend.position = c(1, 1),
                     legend.background = element_rect(colour = "black", fill = "white"))

PDRplot
