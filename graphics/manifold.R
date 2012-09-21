dd = NULL
k1 = seq(0.5, 1.5, 0.01/3)
files = list.files(pattern="^[0-9]{1,3}.csv$")

for(i in seq_along(files)) {
  dd_tmp = read.delim(files[i], header=FALSE)[1:616,]

  dd = rbind(dd, dd_tmp)
}
dd$k2 = seq(0.07, 1.3, 0.006/3)
i = as.numeric(sub(".csv", "", files))

dd$k1 = rep(k1[i], each=616)
write.csv(dd, file="exact_fisher.csv", row.names=FALSE)

#1.5-1 by 0.01/3
#0.07-1.3 by 0.006/3

?heatmap
require(ggplot2)
head(dd)

ggplot(dd) + geom_raster(aes(x=k1, y=k2, fill=V3))
