library(Mcomp)
library(forecast)
library(dplyr)
library(ggplot2)
library(plyr)

setwd("~/thesis/thesis-master/R")

csv <- list.files()
csv <- csv[grepl(pattern = ".csv", x = csv, fixed = T)]
series <- csv[!grepl(pattern = "N", x = csv, fixed = T)]
all_files <- Reduce(rbind.fill, lapply(series, read.csv))

N0049 <- read.csv("N1345.csv")
truth <- cbind("TRUTH", t(data.frame(as.numeric(M3$N1345$xx))))
truth <- as.data.frame(truth)
names(truth) <- c("method", paste0("V", 1:(ncol(truth)-1)))

N0049 %>% 
  select(-1,-2,-3) %>% 
  rbind(truth) %>% 
  melt(id.vars = 'method') -> N0049_melt

N0049_melt$variable <- gsub(pattern = "V", fixed = T, x = N0049_melt$variable, replacement = "")
names(N0049_melt)[1] <- "Method"
N0049_melt$value <- as.numeric(N0049_melt$value)
N0049_melt$variable <- as.numeric(N0049_melt$variable)
ggplot(N0049_melt, aes(x = variable, y = value)) + geom_line(aes(color = Method, group = Method, linetype = Method)) + 
  scale_linetype_manual(values = c(2:5,1:5,1:5,1)) +
  scale_color_manual(values=c(c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7"),c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2"),1)) +
  xlab("Horizon") +
  ylab("N1345") +
  ggtitle("Series N1345 (Demographic, Quarterly) Forecast Horizon") -> p

autoplot(M3$N1345$x) +
  ylab('N1345') +
  ggtitle("N1345 (Demographic, Quarterly)") -> q

#ggsave(p, width = 7, units = "in")
#ggsave(q, width = 7, units = "in")
