library(ggplot2)
library("scales")
##Time course data
head(dd)
dd = read.csv("data/immigration_death_k1_1_k2_1.csv")
full = read.csv("data/immigration_death_full_k1_1_k2_1.csv")
g1 = ggplot(full, aes(x=times, y=n)) + 
    geom_step() + 
    xlab("Time") + ylab("Population")

pdf("figure1.pdf", width=4, height=2.5)
print(g1)
sink=dev.off()
system("pdfcrop figure1.pdf")

pdf("figure1b.pdf", width=4, height=2.5)
print(g1 +geom_point(data=dd, col=2, size=3))
sink=dev.off()
system("pdfcrop figure1b.pdf")

###likelihood for ID process
g = ggplot(df, aes(fill=exp(z), x=x, y=y)) + 
    geom_raster() + 
    guides(fill=FALSE) + 
    xlab(expression(k[1])) + ylab(expression(k[2])) +
    scale_fill_continuous(low="white", high="red") + 
    xlim(c(0,10)) + ylim(c(0,10))
g = g + scale_y_continuous(expand=c(0, 0)) + 
    scale_x_continuous(expand=c(0, 0))
g = g + stat_contour(aes(z=exp(z)), size=0.2)
pdf("figure2.pdf", width=4, height=2.5)
print(g)
sink=dev.off()
system("pdfcrop figure2.pdf")



