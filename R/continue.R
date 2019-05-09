setwd("~/thesis/thesis-master/R")
library(tidyverse)
library(plyr)
library(data.table)
library(forecast)
library(Mcomp)
library(reshape2)

csv <- list.files()
acc <- csv[grepl(pattern = "N", x = csv, fixed = T)]

all_files <- Reduce(rbind.fill, lapply(acc, read.csv))
fwrite(all_files, "allForecasts.csv", row.names = F)

################################################################################################
### MAPE
################################################################################################

output <- list()
for(i in 1:nrow(all_files)){
  cat(i, "\n")
  Series <- as.character(all_files$Series[i])
  Period <- as.character(all_files$Period[i])
  Type <- as.character(all_files$Type[i])
  Method <- as.character(all_files$method[i])
  
  test <-  as.numeric(M3[[all_files$Series[i]]]$xx)
  len <- length(test)
  forecast <- as.numeric(all_files[i, 5:(5 + length(M3[[all_files$Series[i]]]$xx)-1)])
  
  mapeVec <- vector()
  for(j in 1:len){
    mape <- accuracy(forecast[j], test[j])
    mapeVec <- c(mapeVec, mape[length(mape)])
  }
  
  while(length(mapeVec) < 18) mapeVec <- c(mapeVec, NA)
  
  out <- data.frame(Series = as.character(Series),
              Period = as.character(Period),
              Type = as.character(Type),
              Method = as.character(Method),
              MAPE_1 = mapeVec[1],
              MAPE_2 = mapeVec[2],
              MAPE_3 = mapeVec[3],
              MAPE_4 = mapeVec[4],
              MAPE_5 = mapeVec[5],
              MAPE_6 = mapeVec[6],
              MAPE_7 = mapeVec[7],
              MAPE_8 = mapeVec[8],
              MAPE_9 = mapeVec[9],
              MAPE_10 = mapeVec[10],
              MAPE_11 = mapeVec[11],
              MAPE_12 = mapeVec[12],
              MAPE_13 = mapeVec[13],
              MAPE_14 = mapeVec[14],
              MAPE_15 = mapeVec[15],
              MAPE_16 = mapeVec[16],
              MAPE_17 = mapeVec[17],
              MAPE_18 = mapeVec[18])
  
  output[[i]] <- out
}

mapeDF <- do.call(rbind,output)
mapeDF %>% 
  group_by(Method) %>% 
  select(-1,-2,-3,-4) %>% 
  dplyr::summarise(mean(MAPE_1, na.rm = T),
                   mean(MAPE_2, na.rm = T),
                   mean(MAPE_3, na.rm = T),
                   mean(MAPE_4, na.rm = T),
                   mean(MAPE_5, na.rm = T),
                   mean(MAPE_6, na.rm = T),
                   mean(MAPE_7, na.rm = T),
                   mean(MAPE_8, na.rm = T),
                   mean(MAPE_9, na.rm = T),
                   mean(MAPE_10, na.rm = T),
                   mean(MAPE_11, na.rm = T),
                   mean(MAPE_12, na.rm = T),
                   mean(MAPE_13, na.rm = T),
                   mean(MAPE_14, na.rm = T),
                   mean(MAPE_15, na.rm = T),
                   mean(MAPE_16, na.rm = T),
                   mean(MAPE_17, na.rm = T),
                   mean(MAPE_18, na.rm = T)) -> mapeSummary

names(mapeSummary) <- c("Method", paste0("h", 1:18))
df_melted = melt(mapeSummary, id.vars = 'Method')
ggplot(df_melted, aes(x = variable, y = value)) + geom_line(aes(color = Method, group = Method, linetype = Method)) + ylim(9,111)

#Yearly
mapeDF %>% 
  filter(Period == "YEARLY") %>% 
  group_by(Method) %>% 
  select(-1,-2,-3,-4, 5:10) %>% 
  dplyr::summarise(mean(MAPE_1, na.rm = T),
                   mean(MAPE_2, na.rm = T),
                   mean(MAPE_3, na.rm = T),
                   mean(MAPE_4, na.rm = T),
                   mean(MAPE_5, na.rm = T),
                   mean(MAPE_6, na.rm = T)
                   ) -> mapeSummary

names(mapeSummary) <- c("Method", paste0("h", 1:6))
df_melted = melt(mapeSummary, id.vars = 'Method')
ggplot(df_melted, aes(x = variable, y = value)) + geom_line(aes(color = Method, group = Method, linetype = Method)) +
  scale_linetype_manual(values = c(1:5,1:5,1:4))

# Quarterly
mapeDF %>% 
  filter(Period == "QUARTERLY") %>% 
  group_by(Method) %>% 
  select(-1,-2,-3,-4) %>% 
  dplyr::summarise(mean(MAPE_1, na.rm = T),
                   mean(MAPE_2, na.rm = T),
                   mean(MAPE_3, na.rm = T),
                   mean(MAPE_4, na.rm = T),
                   mean(MAPE_5, na.rm = T),
                   mean(MAPE_6, na.rm = T),
                   mean(MAPE_7, na.rm = T),
                   mean(MAPE_8, na.rm = T)
                  ) -> mapeSummary

names(mapeSummary) <- c("Method", paste0("h", 1:8))
df_melted = melt(mapeSummary, id.vars = 'Method')
ggplot(df_melted, aes(x = variable, y = value)) + geom_line(aes(color = Method, group = Method, linetype = Method)) +
  scale_linetype_manual(values = c(1:5,1:5,1:4))

#Monthly
mapeDF %>% 
  filter(Period == "MONTHLY") %>% 
  group_by(Method) %>% 
  select(-1,-2,-3,-4) %>% 
  dplyr::summarise(mean(MAPE_1, na.rm = T),
                   mean(MAPE_2, na.rm = T),
                   mean(MAPE_3, na.rm = T),
                   mean(MAPE_4, na.rm = T),
                   mean(MAPE_5, na.rm = T),
                   mean(MAPE_6, na.rm = T),
                   mean(MAPE_7, na.rm = T),
                   mean(MAPE_8, na.rm = T),
                   mean(MAPE_9, na.rm = T),
                   mean(MAPE_10, na.rm = T),
                   mean(MAPE_11, na.rm = T),
                   mean(MAPE_12, na.rm = T),
                   mean(MAPE_13, na.rm = T),
                   mean(MAPE_14, na.rm = T),
                   mean(MAPE_15, na.rm = T),
                   mean(MAPE_16, na.rm = T),
                   mean(MAPE_17, na.rm = T),
                   mean(MAPE_18, na.rm = T)) -> mapeSummary

names(mapeSummary) <- c("Method", paste0("h", 1:18))
df_melted = melt(mapeSummary, id.vars = 'Method')
ggplot(df_melted, aes(x = variable, y = value)) + geom_line(aes(color = Method, group = Method, linetype = Method)) +
  scale_linetype_manual(values = c(1:5,1:5,1:4))

#OTHER
mapeDF %>% 
  filter(Period == "OTHER") %>% 
  group_by(Method) %>% 
  select(-1,-2,-3,-4) %>% 
  dplyr::summarise(mean(MAPE_1, na.rm = T),
                   mean(MAPE_2, na.rm = T),
                   mean(MAPE_3, na.rm = T),
                   mean(MAPE_4, na.rm = T),
                   mean(MAPE_5, na.rm = T),
                   mean(MAPE_6, na.rm = T),
                   mean(MAPE_7, na.rm = T),
                   mean(MAPE_8, na.rm = T)) -> mapeSummary

names(mapeSummary) <- c("Method", paste0("h", 1:8))
df_melted = melt(mapeSummary, id.vars = 'Method')
ggplot(df_melted, aes(x = variable, y = value)) + 
  geom_line(aes(color = Method, group = Method, linetype = Method)) +
  scale_linetype_manual(values = c(1:5,1:5,1:4))


################################################################################################
### nRMSE
################################################################################################

output <- list()
for(i in 1:nrow(all_files)){
  cat(i, "\n")
  Series <- as.character(all_files$Series[i])
  Period <- as.character(all_files$Period[i])
  Type <- as.character(all_files$Type[i])
  Method <- as.character(all_files$method[i])
  mean <- mean(M3[[all_files$Series[i]]]$x)
  test <-  as.numeric(M3[[all_files$Series[i]]]$xx)/mean
  len <- length(test)
  forecast <- as.numeric(all_files[i, 5:(5 + length(M3[[all_files$Series[i]]]$xx)-1)])/mean
  
  nrmseVec <- vector()
  for(j in 1:len){
    nrmse <- accuracy(forecast[j], test[j])
    nrmseVec <- c(nrmseVec, nrmse[2])
  }
  
  while(length(nrmseVec) < 18) nrmseVec <- c(nrmseVec, NA)
  
  out <- data.frame(Series = as.character(Series),
                    Period = as.character(Period),
                    Type = as.character(Type),
                    Method = as.character(Method),
                    rmse_1 = nrmseVec[1],
                    rmse_2 = nrmseVec[2],
                    rmse_3 = nrmseVec[3],
                    rmse_4 = nrmseVec[4],
                    rmse_5 = nrmseVec[5],
                    rmse_6 = nrmseVec[6],
                    rmse_7 = nrmseVec[7],
                    rmse_8 = nrmseVec[8],
                    rmse_9 = nrmseVec[9],
                    rmse_10 = nrmseVec[10],
                    rmse_11 = nrmseVec[11],
                    rmse_12 = nrmseVec[12],
                    rmse_13 = nrmseVec[13],
                    rmse_14 = nrmseVec[14],
                    rmse_15 = nrmseVec[15],
                    rmse_16 = nrmseVec[16],
                    rmse_17 = nrmseVec[17],
                    rmse_18 = nrmseVec[18])
  
  output[[i]] <- out
}

rmseDF <- do.call(rbind,output)
rmseDF %>% 
  group_by(Method) %>% 
  select(-1,-2,-3,-4) %>% 
  dplyr::summarise(mean(rmse_1, na.rm = T),
                   mean(rmse_2, na.rm = T),
                   mean(rmse_3, na.rm = T),
                   mean(rmse_4, na.rm = T),
                   mean(rmse_5, na.rm = T),
                   mean(rmse_6, na.rm = T),
                   mean(rmse_7, na.rm = T),
                   mean(rmse_8, na.rm = T),
                   mean(rmse_9, na.rm = T),
                   mean(rmse_10, na.rm = T),
                   mean(rmse_11, na.rm = T),
                   mean(rmse_12, na.rm = T),
                   mean(rmse_13, na.rm = T),
                   mean(rmse_14, na.rm = T),
                   mean(rmse_15, na.rm = T),
                   mean(rmse_16, na.rm = T),
                   mean(rmse_17, na.rm = T),
                   mean(rmse_18, na.rm = T)) -> rmseSummary

names(rmseSummary) <- c("Method", paste0("h", 1:18))
df_melted = melt(rmseSummary, id.vars = 'Method')
ggplot(df_melted, aes(x = variable, y = value)) + geom_line(aes(color = Method, group = Method, linetype = Method))+
  scale_linetype_manual(values = c(1:5,1:5,1:4))

#Yearly
rmseDF %>% 
  filter(Period == "YEARLY") %>% 
  group_by(Method) %>% 
  select(-1,-2,-3,-4, 5:10) %>% 
  dplyr::summarise(mean(rmse_1, na.rm = T),
                   mean(rmse_2, na.rm = T),
                   mean(rmse_3, na.rm = T),
                   mean(rmse_4, na.rm = T),
                   mean(rmse_5, na.rm = T),
                   mean(rmse_6, na.rm = T)
  ) -> rmseSummary

names(rmseSummary) <- c("Method", paste0("h", 1:6))
df_melted = melt(rmseSummary, id.vars = 'Method')
ggplot(df_melted, aes(x = variable, y = value)) + geom_line(aes(color = Method, group = Method, linetype = Method)) +
  scale_linetype_manual(values = c(1:5,1:5,1:4))

# Quarterly
rmseDF %>% 
  filter(Period == "QUARTERLY") %>% 
  group_by(Method) %>% 
  select(-1,-2,-3,-4) %>% 
  dplyr::summarise(mean(rmse_1, na.rm = T),
                   mean(rmse_2, na.rm = T),
                   mean(rmse_3, na.rm = T),
                   mean(rmse_4, na.rm = T),
                   mean(rmse_5, na.rm = T),
                   mean(rmse_6, na.rm = T),
                   mean(rmse_7, na.rm = T),
                   mean(rmse_8, na.rm = T)
  ) -> rmseSummary

names(rmseSummary) <- c("Method", paste0("h", 1:8))
df_melted = melt(rmseSummary, id.vars = 'Method')
ggplot(df_melted, aes(x = variable, y = value)) + geom_line(aes(color = Method, group = Method, linetype = Method)) +
  scale_linetype_manual(values = c(1:5,1:5,1:4))

#Monthly
rmseDF %>% 
  filter(Period == "MONTHLY") %>% 
  group_by(Method) %>% 
  select(-1,-2,-3,-4) %>% 
  dplyr::summarise(mean(rmse_1, na.rm = T),
                   mean(rmse_2, na.rm = T),
                   mean(rmse_3, na.rm = T),
                   mean(rmse_4, na.rm = T),
                   mean(rmse_5, na.rm = T),
                   mean(rmse_6, na.rm = T),
                   mean(rmse_7, na.rm = T),
                   mean(rmse_8, na.rm = T),
                   mean(rmse_9, na.rm = T),
                   mean(rmse_10, na.rm = T),
                   mean(rmse_11, na.rm = T),
                   mean(rmse_12, na.rm = T),
                   mean(rmse_13, na.rm = T),
                   mean(rmse_14, na.rm = T),
                   mean(rmse_15, na.rm = T),
                   mean(rmse_16, na.rm = T),
                   mean(rmse_17, na.rm = T),
                   mean(rmse_18, na.rm = T)) -> rmseSummary

names(rmseSummary) <- c("Method", paste0("h", 1:18))
df_melted = melt(rmseSummary, id.vars = 'Method')
ggplot(df_melted, aes(x = variable, y = value)) + geom_line(aes(color = Method, group = Method, linetype = Method)) +
  scale_linetype_manual(values = c(1:5,1:5,1:4))

#OTHER
rmseDF %>% 
  filter(Period == "OTHER") %>% 
  group_by(Method) %>% 
  select(-1,-2,-3,-4) %>% 
  dplyr::summarise(mean(rmse_1, na.rm = T),
                   mean(rmse_2, na.rm = T),
                   mean(rmse_3, na.rm = T),
                   mean(rmse_4, na.rm = T),
                   mean(rmse_5, na.rm = T),
                   mean(rmse_6, na.rm = T),
                   mean(rmse_7, na.rm = T),
                   mean(rmse_8, na.rm = T)) -> rmseSummary

names(rmseSummary) <- c("Method", paste0("h", 1:8))
df_melted = melt(rmseSummary, id.vars = 'Method')
ggplot(df_melted, aes(x = variable, y = value)) + 
  geom_line(aes(color = Method, group = Method, linetype = Method)) +
  scale_linetype_manual(values = c(1:5,1:5,1:4))    

################################################################################################
### nMAE
################################################################################################

output <- list()
for(i in 1:nrow(all_files)){
  cat(i, "\n")
  Series <- as.character(all_files$Series[i])
  Period <- as.character(all_files$Period[i])
  Type <- as.character(all_files$Type[i])
  Method <- as.character(all_files$method[i])
  mean <- mean(M3[[all_files$Series[i]]]$x)
  test <-  as.numeric(M3[[all_files$Series[i]]]$xx)/mean
  len <- length(test)
  forecast <- as.numeric(all_files[i, 5:(5 + length(M3[[all_files$Series[i]]]$xx)-1)])/mean
  
  nmaeVec <- vector()
  for(j in 1:len){
    nmae <- accuracy(forecast[j], test[j])
    nmaeVec <- c(nmaeVec, nmae[3])
  }
  
  while(length(nmaeVec) < 18) nmaeVec <- c(nmaeVec, NA)
  
  out <- data.frame(Series = as.character(Series),
                    Period = as.character(Period),
                    Type = as.character(Type),
                    Method = as.character(Method),
                    mae_1 = nmaeVec[1],
                    mae_2 = nmaeVec[2],
                    mae_3 = nmaeVec[3],
                    mae_4 = nmaeVec[4],
                    mae_5 = nmaeVec[5],
                    mae_6 = nmaeVec[6],
                    mae_7 = nmaeVec[7],
                    mae_8 = nmaeVec[8],
                    mae_9 = nmaeVec[9],
                    mae_10 = nmaeVec[10],
                    mae_11 = nmaeVec[11],
                    mae_12 = nmaeVec[12],
                    mae_13 = nmaeVec[13],
                    mae_14 = nmaeVec[14],
                    mae_15 = nmaeVec[15],
                    mae_16 = nmaeVec[16],
                    mae_17 = nmaeVec[17],
                    mae_18 = nmaeVec[18])
  
  output[[i]] <- out
}

maeDF <- do.call(rbind,output)
maeDF %>% 
  group_by(Method) %>% 
  select(-1,-2,-3,-4) %>% 
  dplyr::summarise(mean(mae_1, na.rm = T),
                   mean(mae_2, na.rm = T),
                   mean(mae_3, na.rm = T),
                   mean(mae_4, na.rm = T),
                   mean(mae_5, na.rm = T),
                   mean(mae_6, na.rm = T),
                   mean(mae_7, na.rm = T),
                   mean(mae_8, na.rm = T),
                   mean(mae_9, na.rm = T),
                   mean(mae_10, na.rm = T),
                   mean(mae_11, na.rm = T),
                   mean(mae_12, na.rm = T),
                   mean(mae_13, na.rm = T),
                   mean(mae_14, na.rm = T),
                   mean(mae_15, na.rm = T),
                   mean(mae_16, na.rm = T),
                   mean(mae_17, na.rm = T),
                   mean(mae_18, na.rm = T)) -> maeSummary

names(maeSummary) <- c("Method", paste0("h", 1:18))
df_melted = melt(maeSummary, id.vars = 'Method')
ggplot(df_melted, aes(x = variable, y = value)) + geom_line(aes(color = Method, group = Method, linetype = Method))+
  scale_linetype_manual(values = c(1:5,1:5,1:4))

#Yearly
maeDF %>% 
  filter(Period == "YEARLY") %>% 
  group_by(Method) %>% 
  select(-1,-2,-3,-4, 5:10) %>% 
  dplyr::summarise(mean(mae_1, na.rm = T),
                   mean(mae_2, na.rm = T),
                   mean(mae_3, na.rm = T),
                   mean(mae_4, na.rm = T),
                   mean(mae_5, na.rm = T),
                   mean(mae_6, na.rm = T)
  ) -> maeSummary

names(maeSummary) <- c("Method", paste0("h", 1:6))
df_melted = melt(maeSummary, id.vars = 'Method')
ggplot(df_melted, aes(x = variable, y = value)) + geom_line(aes(color = Method, group = Method, linetype = Method)) +
  scale_linetype_manual(values = c(1:5,1:5,1:4))

# Quarterly
maeDF %>% 
  filter(Period == "QUARTERLY") %>% 
  group_by(Method) %>% 
  select(-1,-2,-3,-4) %>% 
  dplyr::summarise(mean(mae_1, na.rm = T),
                   mean(mae_2, na.rm = T),
                   mean(mae_3, na.rm = T),
                   mean(mae_4, na.rm = T),
                   mean(mae_5, na.rm = T),
                   mean(mae_6, na.rm = T),
                   mean(mae_7, na.rm = T),
                   mean(mae_8, na.rm = T)
  ) -> maeSummary

names(maeSummary) <- c("Method", paste0("h", 1:8))
df_melted = melt(maeSummary, id.vars = 'Method')
ggplot(df_melted, aes(x = variable, y = value)) + geom_line(aes(color = Method, group = Method, linetype = Method)) +
  scale_linetype_manual(values = c(1:5,1:5,1:4))

#Monthly
maeDF %>% 
  filter(Period == "MONTHLY") %>% 
  group_by(Method) %>% 
  select(-1,-2,-3,-4) %>% 
  dplyr::summarise(mean(mae_1, na.rm = T),
                   mean(mae_2, na.rm = T),
                   mean(mae_3, na.rm = T),
                   mean(mae_4, na.rm = T),
                   mean(mae_5, na.rm = T),
                   mean(mae_6, na.rm = T),
                   mean(mae_7, na.rm = T),
                   mean(mae_8, na.rm = T),
                   mean(mae_9, na.rm = T),
                   mean(mae_10, na.rm = T),
                   mean(mae_11, na.rm = T),
                   mean(mae_12, na.rm = T),
                   mean(mae_13, na.rm = T),
                   mean(mae_14, na.rm = T),
                   mean(mae_15, na.rm = T),
                   mean(mae_16, na.rm = T),
                   mean(mae_17, na.rm = T),
                   mean(mae_18, na.rm = T)) -> maeSummary

names(maeSummary) <- c("Method", paste0("h", 1:18))
df_melted = melt(maeSummary, id.vars = 'Method')
ggplot(df_melted, aes(x = variable, y = value)) + geom_line(aes(color = Method, group = Method, linetype = Method)) +
  scale_linetype_manual(values = c(1:5,1:5,1:4))

#OTHER
maeDF %>% 
  filter(Period == "OTHER") %>% 
  group_by(Method) %>% 
  select(-1,-2,-3,-4) %>% 
  dplyr::summarise(mean(mae_1, na.rm = T),
                   mean(mae_2, na.rm = T),
                   mean(mae_3, na.rm = T),
                   mean(mae_4, na.rm = T),
                   mean(mae_5, na.rm = T),
                   mean(mae_6, na.rm = T),
                   mean(mae_7, na.rm = T),
                   mean(mae_8, na.rm = T)) -> maeSummary

names(maeSummary) <- c("Method", paste0("h", 1:8))
df_melted = melt(maeSummary, id.vars = 'Method')
ggplot(df_melted, aes(x = variable, y = value)) + 
  geom_line(aes(color = Method, group = Method, linetype = Method)) +
  scale_linetype_manual(values = c(1:5,1:5,1:4)) 
