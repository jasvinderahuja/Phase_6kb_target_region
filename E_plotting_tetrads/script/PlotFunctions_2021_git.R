require(readr)
clover <- "#008F01"

plotpertetradshiny <- function(tetidd, from=-4300, to=4600, tetradData1=tetradData, DSBval=TRUE, label_pre=NA) {
  tetrads<-as.list(sort(unique(tetradData1$tetrad)))
  snp_stops = sort(unique(tetradData1$snp_start))
  snp_stops_label = snp_stops[!(snp_stops %in% c(-3966, 3534, 4034, 4534, 4734))]
  #needs==Phased.strands.spore2.NH<-readRDS("Phased.strands.spore2.NH.RDS")
  xsc<-scale_x_continuous(breaks = unique(as.numeric(tetradData1$snp_start)),
                          name = NULL,
                          limits=c(from,to),
                          )
  ysc <- scale_y_continuous(limits=c(0,4.2), breaks = c(0.9, 1.9, 2.9, 3.9), labels = NULL, name=NULL)
  colours = c("H" = "green", "S" = "red", "W" = "blue", "UseLabel"="black")
  thm<-theme_bw()+theme(panel.grid.major = element_blank(),
                        panel.grid.minor.y = element_blank(),
                        plot.title = element_text(size=10, face="bold"),
                        axis.text=element_text(size=3),
                        axis.title=element_text(size=10),
                        axis.text.x = element_text(angle = 90, hjust = 1),
                        panel.grid.minor.x = element_line(colour="gray80", size=0.1, linetype = "dashed"),
                        legend.position="none",
                        text = element_text(size=5))
  ddcateg <- subset(tetradData1, tetrad == tetidd)
   labels1 <- ddcateg %>%
               select(spore, y_axis, nReads) %>%
               distinct() %>%
               mutate(P="UseLabel") %>%
               group_by(spore) %>%
               dplyr::mutate(y_axis=mean(y_axis), label1 = paste0(spore, unique(nReads), collapse=","))
  labeled<-geom_text(data=labels1, aes(y=y_axis, x=-3900, label=label1), size=2)
  plots <- ggplot(data=ddcateg)
    if(DSBval==TRUE){
      plots <- plots+
        ggplot2::annotate("rect", xmin = -35, ymin = -Inf, xmax = +22, ymax = +Inf, alpha=1, fill = clover)+
        ggplot2::annotate("rect", xmin = -55, ymin = -Inf, xmax = +78, ymax = +Inf, alpha=0.2, fill = clover)
      #the dark clover colored line covers the DSB cluster (~80% DSBs)
    }
  plots <- plots+
    aes(x = as.numeric(snp_start), y = y_axis, group=P, color = P) +
    geom_point(size=1) + labs(title = paste0(label_pre, unique(ddcateg$tetrad), collapse=""))
      #                           subtitle = unique(ddcateg$MLtake))
  plot=plots+
    labeled+
    scale_color_manual(values = colours)+xsc+ysc+thm
  a<-c(tetidd, plot)
  return(plot)
}





