library(ggplot2)
require(reshape2)
require(grid)
theme_set(theme_bw(base_size = 10))

dd_dis = read.csv("data/immigration_death.csv", header=TRUE)
dd_full = read.csv("data/full_immigration_death.csv", header=TRUE)

g = ggplot(dd_full, aes(times, n)) + 
    geom_step() + 
    geom_point(data=dd_dis, aes(times, n), size=3, colour="red") + 
    xlab("Time") + ylab("Population")

pdf("graphics/figure1.pdf", width=3, height=3)
print(g)
dev.off()