library(DEsingle)

#------------------spin1 vs spin 2---------------------
#params
grp1=69
grp2=54

#files
f_in <- "./out/a06_deg_00_pp/RNAlater_spin-root_69_54.csv"
f_out <- "./out/a06_deg_01_de/RNAlater_spin-root.csv"

#1. load
counts<- read.table(f_in, sep=",", header=TRUE, row.names = 1)
group <- factor(c(rep(1,grp1), rep(2,grp2)))

#2. analysis
results <- DEsingle(counts = counts, group = group, parallel = TRUE)
results.classified <- DEtype(results = results, threshold = 0.05)

#3. write
write.csv(results.classified, f_out)

##------------------spin vs others---------------------
##params
#grp1=280
#grp2=5391

##files
#f_in <- "./out/a06_deg_00_pp/spin_others.csv"
#f_out <- "./out/a06_deg_01_de/spin_others.csv"

##1. load
#counts<- read.table(f_in, sep=",", header=TRUE, row.names = 1)
#group <- factor(c(rep(1,grp1), rep(2,grp2)))

##2. analysis
#results <- DEsingle(counts = counts, group = group, parallel = TRUE)
#results.classified <- DEtype(results = results, threshold = 0.05)

##3. write
#write.csv(results.classified, f_out)
