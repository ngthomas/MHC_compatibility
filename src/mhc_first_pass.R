# MHC's raw first pass through
# written by Thomas Ng

library(tidyverse)

parent.pair.tbl <- read.csv("mhcdatafilesforanalysis/parents_MHC.csv",
                            stringsAsFactors = F) %>% tbl_df()
parent.sex.tbl <- read.csv("mhcdatafilesforanalysis/Adult_sex_year_MHC.csv",
                           stringsAsFactors = F) %>% tbl_df()

entry.label <- read.table("mhcdnaandaafiles/adults_AA_MHC_936_steelhead_274bp.fas",
                          stringsAsFactors = F) %>%
  tbl_df() %>%
  filter(row_number() %%2==1) %>%
  rename(entry.name= V1)

indiv.aa <- read.table("mhcdnaandaafiles/adults_AA_MHC_936_steelhead_274bp.fas",
                       stringsAsFactors = F) %>%
  tbl_df() %>%
  filter(row_number() %%2==0) %>%
  rename(aa.seq= V1)

aa.tbl <- cbind.data.frame(entry.label, indiv.aa) %>%
  group_by(entry.name) %>%
  mutate(
    indiv = strsplit(as.character(entry.name),"_|>")[[1]][2],
    allele = strsplit(as.character(entry.name), "_|>")[[1]][3]) %>%
  ungroup() %>%
  mutate(aa.factor = factor(aa.seq))

# including sex and year information to the aa.seq
aa.join.tbl <- left_join(aa.tbl,parent.sex.tbl %>% rename(indiv=ID), by="indiv")

aa.seq.factor <- factor(levels(aa.tbl$aa.factor))

edit.distance <- combn(aa.seq.factor, 2) %>% 
  apply(.,2, function(x) c(x[1],
                           x[2],
                           adist(aa.seq.factor[x[1]],
                                 aa.seq.factor[x[2]]))) %>%
  t() %>%
  tbl_df() 

colnames(edit.distance) <- c("aa.1","aa.2","edit")
length.aa <- str_length(aa.seq.factor[1])

edit.distance <- edit.distance %>% 
  ungroup() %>%
  mutate(p.distance = edit/length.aa)


# sanity check: making sure that the pairing show up only once 
parent.pair.reduce <- parent.pair.tbl %>% mutate(pa.ma = paste0(Pa, Ma)) %>%
  group_by(pa.ma, SpawnYear) %>%
  summarise(pa = Pa[1],
            ma = Ma[1])

just.sex <- aa.join.tbl %>% group_by(indiv) %>% summarise(Sex = Sex[1])

# founding out that parent.pair.tbl are not exactly sorted properly by sex
reassign.sex <- left_join(parent.pair.reduce %>%
                            rename(p1 = pa,
                                   p2 = ma),
                          just.sex %>% rename(p1=indiv), by="p1" ) %>%
  rename(sex.p1 = Sex) %>%
  left_join(., 
            just.sex %>% rename(p2=indiv), by="p2" ) %>%
  rename(sex.p2 = Sex) 

#write.csv(reassign.sex %>% filter(sex.pa!="Male" | sex.ma!="Female") %>% ungroup %>% select(-pa.ma) , file ="weird_parentpairs.csv", quote = F,row.names = F)


parent.pair.cleanup <- reassign.sex %>% 
  filter((sex.p1=="Male" & sex.p2=="Female") || (sex.p2=="Male" & sex.p1=="Female")) %>%
  ungroup() %>%
  mutate(pa = ifelse(sex.p1=="Male", p1, p2),
         ma = ifelse(sex.p2=="Female", p2, p1)) %>%
  group_by(pa, ma, SpawnYear) %>%
  summarise(ct = n())

pa.count <- parent.pair.cleanup %>%
  group_by(pa) %>%
  summarise(n.partners = n()) %>%
  ungroup() %>%
  mutate(n.cat = n()) %>%
  group_by(n.partners) %>%
  summarise(fraction = n()/n.cat[1],
            perspective="male")

ma.count <- parent.pair.cleanup %>%
  group_by(ma) %>%
  summarise(n.partners = n()) %>%
  ungroup() %>%
  mutate(n.cat = n()) %>%
  group_by(n.partners) %>%
  summarise(fraction = n()/n.cat[1],
            perspective="female")

parent.count <- bind_rows(pa.count, ma.count)
parent.count$n.partners <- factor(parent.count$n.partners, levels=unique(parent.count$n.partners) %>% sort(decreasing = T))

pa.count <- parent.pair.cleanup %>%
  group_by(pa, SpawnYear) %>%
  summarise(n.partners = n()) %>%
  group_by(SpawnYear) %>%
  mutate(n.cat = n()) %>%
  group_by(n.partners,SpawnYear) %>%
  summarise(fraction = n()/n.cat[1],
            perspective="male",
            count = n())

ma.count <- parent.pair.cleanup %>%
  group_by(ma, SpawnYear) %>%
  summarise(n.partners = n()) %>%
  group_by(SpawnYear) %>%
  mutate(n.cat = n()) %>%
  group_by(n.partners, SpawnYear) %>%
  summarise(fraction = n()/n.cat[1],
            perspective="female",
            count=n())

parent.count <- bind_rows(pa.count, ma.count)
parent.count$n.partners <- factor(parent.count$n.partners, levels=unique(parent.count$n.partners) %>% sort(decreasing = T))


# precompute edit distance matrix 
n.aa <- length(aa.seq.factor)

aa.matrix <- matrix(0, nrow=n.aa, ncol=n.aa)
aa.matrix[cbind(edit.distance$aa.1, edit.distance$aa.2)] <- edit.distance$p.distance
aa.matrix[cbind(edit.distance$aa.2, edit.distance$aa.1)] <- edit.distance$p.distance


# Boostraping sampling based on all observed males and female (excluding unknown - it is easier to manipulate and manage; at the same time, there are plenty of collected unpaired known male and female);
# another approach is to only randomized the ordering of the observed pa and ma (it would be good to see the difference)
# These two approaches would make a difference whether we assume that MHC profile has an effect on whether an individual will be driven to look for a mate. (For the null model, we assume that 1)everyone has an equal opportunity to find a partner or all of them are matured to mate 2) and come from the same population group)

# the second layer is that whether the number of mates will be influenced by MHC (is it intrinsic or extrinsic properties?) (for now, let's assume that it's an extrinsic random property such that the number of mating partners assigned to individual is randomized)
# separate each year (filter), 
# randomize the ordering of males and give randomized number of tickets to male for number of partners that could be chose (sample & rep), then randomized the ordering of female (sample),
# sum all pairwise MHC hamming's distance (possibly do a min(four #) for the next try)

set.seed(939204760)

SimulateParentPairs <- function(){
  
  n.iter <- 10000
  year.ls <- parent.pair.cleanup$SpawnYear %>% unique
  full.pdist.ls <- lapply(year.ls, function(year.i) {
    
    parent.pair.yr <- parent.pair.cleanup %>% filter(SpawnYear == year.i)
    n.mating.pair <- parent.pair.yr %>% dim %>% .[1]
    n.male <- parent.pair.yr$pa %>% unique %>% length
    n.female <- parent.pair.yr$ma %>% unique %>% length
    n.partner.male <- parent.pair.yr$pa %>% table %>% as.numeric()
    n.partner.female <- parent.pair.yr$ma %>% table %>% as.numeric()
    
    indiv.aa <- aa.join.tbl %>%
      mutate(aa.numeric = as.numeric(aa.factor)) %>%
      filter(Sex != "Unk", Year == year.i) %>%
      select(indiv, Sex, aa.numeric, allele)
    
    n.obs.male <- indiv.aa %>% filter(Sex=="Male", allele ==1 )%>% dim %>% .[1]
    n.obs.female <- indiv.aa %>% filter(Sex=="Female", allele ==1 ) %>% dim %>% .[1]
    
    indiv.aa.male.1 <- indiv.aa %>% filter(Sex=="Male", allele == 1) %>% select(aa.numeric) %>% unlist %>% as.integer()
    indiv.aa.male.2 <- indiv.aa %>% filter(Sex=="Male", allele == 2) %>% select(aa.numeric) %>% unlist %>% as.integer()
    indiv.aa.female.1 <- indiv.aa %>% filter(Sex=="Female", allele == 1) %>% select(aa.numeric) %>% unlist %>% as.integer()
    indiv.aa.female.2 <- indiv.aa %>% filter(Sex=="Female", allele == 2) %>% select(aa.numeric) %>% unlist %>% as.integer()
    
    mean.sample.pdist <- sapply(1:n.iter, function(i){
      order.dad <- sample(1:n.obs.male, n.male)
      order.mates.male <- sample(n.partner.male)
      order.mates.female <- sample(n.partner.female)
      order.mom <- sample(1:n.obs.female, n.female)
      
      (aa.matrix[cbind(rep(indiv.aa.male.1[order.dad], order.mates.male),
                       rep(indiv.aa.female.1[order.mom], order.mates.female))] +
          aa.matrix[cbind(rep(indiv.aa.male.1[order.dad], order.mates.male),
                          rep(indiv.aa.female.2[order.mom], order.mates.female))] +
          aa.matrix[cbind(rep(indiv.aa.male.2[order.dad], order.mates.male),
                          rep(indiv.aa.female.1[order.mom], order.mates.female))] +
          aa.matrix[cbind(rep(indiv.aa.male.2[order.dad], order.mates.male),
                          rep(indiv.aa.female.2[order.mom], order.mates.female))]) %>%
        mean
    })
    
    
    # true hamming distance of the observed set
    pa.obs.indx <- left_join(parent.pair.yr,
                             indiv.aa %>%
                               rename(pa=indiv)) %>%
      ungroup()
    
    ma.obs.indx <- left_join(parent.pair.yr,
                             indiv.aa %>%
                               rename(ma=indiv)) %>%
      ungroup()
    
    pa.obs.aa.1 <- pa.obs.indx %>% filter(allele==1) %>% select(aa.numeric) %>% unlist %>% as.integer()
    pa.obs.aa.2 <- pa.obs.indx %>% filter(allele==2) %>% select(aa.numeric) %>% unlist %>% as.integer()
    ma.obs.aa.1 <- ma.obs.indx %>% filter(allele==1) %>% select(aa.numeric) %>% unlist %>% as.integer()
    ma.obs.aa.2 <- ma.obs.indx %>% filter(allele==2) %>% select(aa.numeric) %>% unlist %>% as.integer()
    
    obs.sum.pdist <- (aa.matrix[cbind(pa.obs.aa.1, ma.obs.aa.1)] +
                        aa.matrix[cbind(pa.obs.aa.1, ma.obs.aa.2)] +
                        aa.matrix[cbind(pa.obs.aa.2, ma.obs.aa.1)]+
                        aa.matrix[cbind(pa.obs.aa.2, ma.obs.aa.2)])
    
    mean.obs.pdist <- mean(obs.sum.pdist)
    sd.mean.obs.pdist <- sd(obs.sum.pdist)
    
    min.sample.pdist <- sapply(1:n.iter, function(i){
      order.dad <- sample(1:n.obs.male, n.male)
      order.mates.male <- sample(n.partner.male)
      order.mates.female <- sample(n.partner.female)
      order.mom <- sample(1:n.obs.female, n.female)
      
      pmin(aa.matrix[cbind(rep(indiv.aa.male.1[order.dad], order.mates.male),
                           rep(indiv.aa.female.1[order.mom], order.mates.female))] ,
           aa.matrix[cbind(rep(indiv.aa.male.1[order.dad], order.mates.male),
                           rep(indiv.aa.female.2[order.mom], order.mates.female))] ,
           aa.matrix[cbind(rep(indiv.aa.male.2[order.dad], order.mates.male),
                           rep(indiv.aa.female.1[order.mom], order.mates.female))] ,
           aa.matrix[cbind(rep(indiv.aa.male.2[order.dad], order.mates.male),
                           rep(indiv.aa.female.2[order.mom], order.mates.female))]) %>%
        mean
    })
    
    obs.min.pdist <- pmin(aa.matrix[cbind(pa.obs.aa.1, ma.obs.aa.1)],
                          aa.matrix[cbind(pa.obs.aa.1, ma.obs.aa.2)],
                          aa.matrix[cbind(pa.obs.aa.2, ma.obs.aa.1)],
                          aa.matrix[cbind(pa.obs.aa.2, ma.obs.aa.2)])
    
    mean.obs.min.pdist <- mean(obs.min.pdist)
    sd.obs.min.pdist <- sd(obs.min.pdist)
    
    max.sample.pdist <- sapply(1:n.iter, function(i){
      order.dad <- sample(1:n.obs.male, n.male)
      order.mates.male <- sample(n.partner.male)
      order.mates.female <- sample(n.partner.female)
      order.mom <- sample(1:n.obs.female, n.female)
      
      pmax(aa.matrix[cbind(rep(indiv.aa.male.1[order.dad], order.mates.male),
                           rep(indiv.aa.female.1[order.mom], order.mates.female))] ,
           aa.matrix[cbind(rep(indiv.aa.male.1[order.dad], order.mates.male),
                           rep(indiv.aa.female.2[order.mom], order.mates.female))] ,
           aa.matrix[cbind(rep(indiv.aa.male.2[order.dad], order.mates.male),
                           rep(indiv.aa.female.1[order.mom], order.mates.female))] ,
           aa.matrix[cbind(rep(indiv.aa.male.2[order.dad], order.mates.male),
                           rep(indiv.aa.female.2[order.mom], order.mates.female))]) %>%
        mean
    })
    
    obs.max.pdist <- pmax(aa.matrix[cbind(pa.obs.aa.1, ma.obs.aa.1)],
                          aa.matrix[cbind(pa.obs.aa.1, ma.obs.aa.2)],
                          aa.matrix[cbind(pa.obs.aa.2, ma.obs.aa.1)],
                          aa.matrix[cbind(pa.obs.aa.2, ma.obs.aa.2)])
    
    mean.obs.max.pdist <- mean(obs.max.pdist)
    sd.obs.max.pdist <- sd(obs.max.pdist)
    
    list(year = year.i,
         
         pa=parent.pair.yr$pa,
         ma=parent.pair.yr$ma,
         
         sum.mean.sample = mean.sample.pdist,
         sum.pdist.obs = obs.sum.pdist,
         mean.sum.obs = mean.obs.pdist,
         sd.sum.obs = sd.mean.obs.pdist,
         
         min.mean.sample = min.sample.pdist,
         min.pdist.obs = obs.min.pdist,
         mean.min.obs = mean.obs.min.pdist,
         sd.min.obs = sd.obs.min.pdist,
         
         max.mean.sample = max.sample.pdist,
         max.pdist.obs = obs.max.pdist,
         mean.max.obs = mean.obs.max.pdist,
         sd.max.obs = sd.obs.max.pdist
    )
    
    
  })
  return(full.pdist.ls)
}
full.pdist.ls <- SimulateParentPairs()

write.all.df <- function(full.pdist.ls, output.suffix ="allaa") {
  sample.df <- lapply(full.pdist.ls,function(i) data.frame(year=i$year,
                                                           sum.mean.sample=i$sum.mean.sample,
                                                           min.mean.sample=i$min.mean.sample,
                                                           max.mean.sample=i$max.mean.sample
  )) %>% bind_rows()
  
  obs.df <- lapply(full.pdist.ls, function(i) data.frame(year=i$year,
                                                         pa=i$pa,
                                                         ma=i$ma,  
                                                         sum.pdist.obs=i$sum.pdist.obs,
                                                         min.pdist.obs=i$min.pdist.obs,
                                                         max.pdist.obs=i$max.pdist.obs
  )) %>% bind_rows()
  
  mean.obs.df <- lapply(full.pdist.ls, function(i) data.frame(year=i$year,
                                                              mean.sum.obs=i$mean.sum.obs,
                                                              mean.min.obs=i$mean.min.obs,
                                                              mean.max.obs=i$mean.max.obs
  )) %>% bind_rows()
  
  sd.obs.df <- lapply(full.pdist.ls, function(i) data.frame(year=i$year,
                                                            sd.sum.obs=i$sd.sum.obs,
                                                            sd.min.obs=i$sd.min.obs,
                                                            sd.max.obs=i$sd.max.obs
  )) %>% bind_rows()
  
  write.csv(sample.df, file= paste0("out/sample_pdist_mean_",output.suffix,".csv"), row.names = F, quote=F)
  write.csv(obs.df, file= paste0("out/obs_",output.suffix,".csv"), row.names = F, quote = F)
  write.csv(mean.obs.df, file= paste0("out/obs_pdist_mean_",output.suffix,".csv"), row.names = F, quote = F)
  write.csv(sd.obs.df, paste0("out/obs_pdist_sd_",output.suffix,".csv"), row.names = F, quote = F)
  
}

write.all.df(full.pdist.ls)

# reformat the data list into ggplot-friendly data frame
PlottingMean <- function(full.pdist.ls){
  
  sample.mean.df <- lapply(full.pdist.ls, function(i) data.frame(year=i$year, sample.mean=i$sum.mean.sample)) %>%
    bind_rows()
  obs.mean.df <- lapply(full.pdist.ls, function(i) data.frame(year=i$year, obs.mean=i$mean.sum.obs)) %>%
    bind_rows()
  
  sample.mean.stat.df <- sample.mean.df %>%
    group_by(year) %>%
    summarise(
      q0.01 = quantile(sample.mean, 0.01),
      q0.05 = quantile(sample.mean, 0.05),
      q0.1 = quantile(sample.mean, 0.1),
      q0.2 = quantile(sample.mean, 0.2),
      q0.5 = quantile(sample.mean, 0.5)
    ) %>% gather(.,"quantile","value",2:6)
  
  sample.mean.stat.df$quantile <- factor(sample.mean.stat.df$quantile,
                                         levels=sort(levels(factor(sample.mean.stat.df$quantile )),decreasing = T))
  
  ggplot() +
    geom_violin(data=sample.mean.df, aes(x=1,y=sample.mean), alpha=0.1)+
    geom_point(data=sample.mean.stat.df, aes(x=1, y=value, color=quantile), shape=95, cex=8)+
    geom_point(data=obs.mean.df, aes(x=1, y=obs.mean), shape=42, cex=8)+
    facet_grid(~year)+
    theme_bw()+
    xlab("")+
    ylab("sum of pairwise p-distance")+
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank())
}
#PlottingMean(full.pdist.ls)
PlottingMin <- function(){
  
  sample.min.df <- lapply(full.pdist.ls, function(i) data.frame(year=i$year, sample.min=i$min.mean.sample)) %>%
    bind_rows()
  
  obs.min.df <- lapply(full.pdist.ls, function(i) data.frame(year=i$year, obs.min=i$mean.min.obs)) %>%
    bind_rows()
  
  sample.min.stat.df <- sample.min.df %>%
    group_by(year) %>%
    summarise(
      q0.01 = quantile(sample.min, 0.01),
      q0.05 = quantile(sample.min, 0.05),
      q0.1 = quantile(sample.min, 0.1),
      q0.2 = quantile(sample.min, 0.2),
      q0.5 = quantile(sample.min, 0.5)
    ) %>% gather(.,"quantile","value",2:6)
  
  sample.min.stat.df$quantile <- factor(sample.min.stat.df$quantile,
                                        levels=sort(levels(factor(sample.min.stat.df$quantile )),decreasing = T))
  
  ggplot() +
    geom_violin(data=sample.min.df, aes(x=1,y=sample.min), alpha=0.1)+
    geom_point(data=sample.min.stat.df, aes(x=1, y=value, color=quantile), shape=95, cex=8)+
    geom_point(data=obs.min.df, aes(x=1, y=obs.min), shape=42, cex=8)+
    facet_grid(~year)+
    theme_bw()+
    xlab("")+
    ylab("minimum of pairwise p-distance")+
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank())
}
#PlottingMin()


#### isolate AA to only select functional sites

#The information about which aa functional sites are important are highlighted from an email sent by Kerry. - see"AA's of PBS for Steelhead"

#The AA's thought to play a role in peptide binding are the following:

#9, 11, 13, 28, 30, 32, 37, 38, 47, 56, 60 61, 65, 68, 70, 71, 74, 78, 81, 82, 85, 86, 88, 89

important.aa.indx <- c(9, 11, 13, 28, 30, 32, 37, 38, 47, 56, 60, 61, 65, 68, 70, 71, 74, 78, 81, 82, 85, 86, 88, 89)

indiv.aa.tbl <- read.table("mhcdnaandaafiles/adults_AA_MHC_936_steelhead_274bp.fas",
                           stringsAsFactors = F) %>%
  tbl_df() %>%
  filter(row_number() %%2==0) %>%
  rename(aa.seq= V1) 


aa.tbl <- cbind.data.frame(entry.label, indiv.aa.tbl) %>%
  group_by(entry.name) %>%
  mutate(
    indiv = strsplit(as.character(entry.name),"_|>")[[1]][2],
    allele = strsplit(as.character(entry.name), "_|>")[[1]][3],
    aa.seq = paste0(strsplit(as.character(aa.seq),"")[[1]][important.aa.indx], collapse="")
  ) %>%
  ungroup() %>%
  mutate(aa.factor = factor(aa.seq))

# including sex and year information to the aa.seq
aa.join.tbl <- left_join(aa.tbl,parent.sex.tbl %>% rename(indiv=ID), by="indiv")

aa.seq.factor <- factor(levels(aa.tbl$aa.factor))

edit.distance <- combn(aa.seq.factor, 2) %>% 
  apply(.,2, function(x) c(x[1],
                           x[2],
                           adist(aa.seq.factor[x[1]],
                                 aa.seq.factor[x[2]]))) %>%
  t() %>%
  tbl_df() 

colnames(edit.distance) <- c("aa.1","aa.2","edit")
length.aa <- str_length(aa.seq.factor[1])

edit.distance <- edit.distance %>% 
  ungroup() %>%
  mutate(p.distance = edit/length.aa)

# precompute edit distance matrix 
n.aa <- length(aa.seq.factor)

aa.matrix <- matrix(0, nrow=n.aa, ncol=n.aa)
aa.matrix[cbind(edit.distance$aa.1, edit.distance$aa.2)] <- edit.distance$p.distance
aa.matrix[cbind(edit.distance$aa.2, edit.distance$aa.1)] <- edit.distance$p.distance

full.pdist.ls <- SimulateParentPairs()
PlottingMean(full.pdist.ls)
write.all.df(full.pdist.ls,output.suffix ="fnaa")
PlottingMin()
