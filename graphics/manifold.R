require(ggplot2)
require(reshape2)
require(grid)
library(gridExtra)
theme_set(theme_bw(base_size = 11))
vplayout = function(x, y) 
    viewport(layout.pos.row = x, layout.pos.col = y)

std = function(value) {
    value/max(value)
}

dd_e = read.csv("data/exact_fisher.csv", header=TRUE)
dd_e$type = "Exact"
colnames(dd_e)[1:3] = c("d.dk1k1", "d.dk1k2", "d.dk2k2")

dd_mc = read.csv("data/mc_fisher.csv", header=TRUE)
dd_mc$type = "Moment~Closure"

# 
# dd_mc[3] = std(dd_mc[3])
# dd_mc[4] = std(dd_mc[4])
# dd_mc[5] = std(dd_mc[5])
# dd_e[1] = std(dd_e[1])
# dd_e[2] = std(dd_e[2])
# dd_e[3] = std(dd_e[3])
# apply(dd_e[1:3], 2, range)
# apply(dd_mc[3:5], 2, range)

dd = rbind(dd_mc, dd_e)
dd=dd[dd$k1 > 1.2 & dd$k2 > 1, ]
dd$close11 = "d^2/dk[1] * dk[1]"
dd$close12 = "d^2/dk[1] * dk[2]"
dd$close22 = "d^2/dk[2] * dk[2]"


dd = expand.grid(x=1:10, y=1:10)
dd = data.frame(dd, type=rep(LETTERS[1:2], each=100), 
                var =rep(c("C", "D"), each=200) )
dd$z = rnorm(400, rep(c(0, 100), each=200))


p1 <- ggplot(subset(dd, var=="C"), aes(x,y))+
    geom_raster(aes(fill=z)) + facet_grid(type ~ var) + 
    theme(legend.position="bottom", plot.margin = unit(c(1,-1,1,0.2), "line"))
p2 <- ggplot(subset(dd, var=="D"), aes(x,y))+
    geom_raster(aes(fill=z)) + facet_grid(type ~ var) + 
    theme(legend.position="bottom", 
          plot.margin = unit(c(1,1,1,-0.8), "line"),
          axis.text.y = element_blank(), axis.ticks.y = element_blank()) + ylab(NULL)
grid.arrange(arrangeGrob(p1, p2, nrow = 1))

# 
# dd = melt(dd, c("k1", "k2", "type"))
head(dd)

g = ggplot(dd,aes(k1, k2)) + 
    scale_x_continuous(expand=c(0, 0.01)) + 
    scale_y_continuous(expand=c(0, 0)) + 
    xlab(expression(k[1])) + 
    theme(legend.position="bottom") + 
    scale_fill_continuous(name = "") 
     

g1 =g + geom_raster(aes(x=k1, y=k2, fill=d.dk1k1))  + 
    stat_contour(aes(z=d.dk1k1)) +
    facet_grid(type~close11, labeller= label_parsed) 

g1b = g1 + theme(plot.margin = unit(c(1,-0.8,1,0.2), "line"))
    

g2 = g + geom_raster(aes(fill=d.dk1k2))  + 
    stat_contour(aes(z=d.dk1k2)) +
    facet_grid(type~close11, labeller= label_parsed) + 
    ylab(NULL)

g2b = g2 +  theme(legend.position="bottom", 
          plot.margin = unit(c(1,-0.8,1, 0), "line"),
          axis.text.y = element_blank(), axis.ticks.y = element_blank())
    

g3 = g + geom_raster(aes(fill=d.dk2k2))  + 
    stat_contour(aes(z=d.dk2k2)) +
    facet_grid(type~close11, labeller= label_parsed) + 
    ylab(NULL) 

g3b = g3 + theme(legend.position="bottom", 
                   plot.margin = unit(c(1,1,1,-0.25), "line"),
                   axis.text.y = element_blank(), axis.ticks.y = element_blank())



# 
# grid.newpage()
# pushViewport(viewport(layout = grid.layout(1, 3)))
# print(g1, vp = vplayout(1, 1))
# print(g2, vp = vplayout(1, 2))
# print(g3, vp = vplayout(1, 3))
# 



pdf("graphics/figure2.pdf", width=10, height=7)
grid.arrange(arrangeGrob(g1b, g2b, g3b, nrow = 1))
dev.off()

