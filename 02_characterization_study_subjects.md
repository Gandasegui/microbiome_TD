# Now, let's analyse the data to create the Table1

```R
library(tidyverse)
library(data.table)

#Loeading metadata file
metadata <- read_tsv(file = "metadata_TD.tsv")
metadata <- metadata[-1, ]
metadata_tib <- as_tibble(metadata)
metadata_tib$bristol <- as.numeric(metadata_tib$bristol)
```

For the three-group-criteria - variable is codified as "clingroup"
For the two-group-criteria - variable is codified as "clingroup_2"

#### For the three-group-criteria

```R
x <- c('1A', '1B', '2')
#Number of participants
for (i in x){
  q <- metadata_tib %>% filter(clingroup == i)
  print(paste0('Group ', i, ', N = ', nrow(q)))
}
#Age
metadata_tib$age <- as.numeric(metadata_tib$age)
for (i in x){
  q <- metadata_tib %>% filter(clingroup == i)
  print(paste0('Group ', i, ', median_age = ', median(na.omit(q$age)),
               ', IQR = ', quantile(na.omit(q$age))[2], '-', quantile(na.omit(q$age))[4]
  ))
}
kruskal.test(age ~ clingroup, data = metadata_tib)#p-value = 0.064

#Sex
metadata_tib$sex <- as.factor(metadata_tib$sex)
for (i in x){
  q <- metadata_tib %>% filter(clingroup == i)
  p <- q %>% filter(sex == 'female')
  print(paste0('Group ', i, ', %_women = ', nrow(p)/nrow(q)*100))
}
kruskal.test(sex ~ clingroup, data = metadata_tib)#p-value = 0.229

#Bristol
metadata_tib$bristol <- as.numeric(metadata_tib$bristol)
for (i in x){
  q <- metadata_tib %>% filter(clingroup == i)
  print(paste0('Group ', i, ', median_bristol = ', median(na.omit(q$bristol)),
               ', IQR = ', quantile(q$bristol)[2], '-', quantile(q$bristol)[4]
               ))
}
kruskal.test(bristol ~ clingroup, data = metadata_tib)#p-value = 5.328e-08

#Nºstools
metadata_tib$stoolno <- as.numeric(metadata_tib$stoolno)
for (i in x){
  q <- metadata_tib %>% filter(clingroup == i)
  print(paste0('Group ', i, ', median_Nstool = ', median(na.omit(q$stoolno)),
               ', IQR = ', quantile(na.omit(q$stoolno))[2], '-', quantile(na.omit(q$stoolno))[4]
               ))
}
kruskal.test(stoolno ~ clingroup, data = metadata_tib)#p-value = 0.07526

#LCN2
metadata_tib$lcn2 <- as.numeric(metadata_tib$lcn2)
for (i in x){
  q <- metadata_tib %>% filter(clingroup == i)
  print(paste0('Group ', i, ', median_LCN2 = ', median(na.omit(q$lcn2)),
               ', IQR = ', quantile(q$lcn2)[2], '-', quantile(q$lcn2)[4]
  ))
}
kruskal.test(lcn2 ~ clingroup, data = metadata_tib)#p-value = 0.2551
```

#### For the two-group-criteria

```R
#Number of participants
for (i in (1:2)){
  q <- metadata_tib %>% filter(clingroup_2 == i)
  print(paste0('Group ', i, ', N = ', nrow(q)))
}
#Age
metadata_tib$age <- as.numeric(metadata_tib$age)
for (i in (1:2)){
  q <- metadata_tib %>% filter(clingroup_2 == i)
  print(paste0('Group ', i, ', median_age = ', median(na.omit(q$age)),
               ', IQR = ', quantile(na.omit(q$age))[2], '-', quantile(na.omit(q$age))[4]
  ))
}
kruskal.test(age ~ clingroup_2, data = metadata_tib)#p-value = 0.784

#Sex
metadata_tib$sex <- as.factor(metadata_tib$sex)
for (i in (1:2)){
  q <- metadata_tib %>% filter(clingroup_2 == i)
  p <- q %>% filter(sex == 'female')
  print(paste0('Group ', i, ', %_women = ', nrow(p)/nrow(q)*100))
}
kruskal.test(sex ~ clingroup_2, data = metadata_tib)#p-value = 0.7563

#Bristol
metadata_tib$bristol <- as.numeric(metadata_tib$bristol)
for (i in (1:2)){
  q <- metadata_tib %>% filter(clingroup_2 == i)
  print(paste0('Group ', i, ', median_bristol = ', median(na.omit(q$bristol)),
               ', IQR = ', quantile(q$bristol)[2], '-', quantile(q$bristol)[4]
  ))
}
kruskal.test(bristol ~ clingroup_2, data = metadata_tib)#p-value = 7.484e-09
#Nºstools
metadata_tib$stoolno <- as.numeric(metadata_tib$stoolno)
for (i in (1:2)){
  q <- metadata_tib %>% filter(clingroup_2 == i)
  print(paste0('Group ', i, ', median_Nstool = ', median(na.omit(q$stoolno)),
               ', IQR = ', quantile(na.omit(q$stoolno))[2], '-', quantile(na.omit(q$stoolno))[4]
  ))
}
kruskal.test(stoolno ~ clingroup_2, data = metadata_tib)#p-value = 0.03901
#LCN2
metadata_tib$lcn2 <- as.numeric(metadata_tib$lcn2)
for (i in (1:2)){
  q <- metadata_tib %>% filter(clingroup_2 == i)
  print(paste0('Group ', i, ', median_LCN2 = ', median(na.omit(q$lcn2)),
               ', IQR = ', quantile(q$lcn2)[2], '-', quantile(q$lcn2)[4]
  ))
}
kruskal.test(lcn2 ~ clingroup_2, data = metadata_tib)#p-value = 0.9086
```
