library(ggplot2)
library(jsonlite)
# list all sample names
snames <- c("MATS_0001_Ly_P_WG_1",
            "MATS_0001_Ly_P_WG_2",
            "MATS_0002_Ly_P_WG_1",
            "MATS_0002_Ly_P_WG_2",
            "MATS_0003_Ly_P_WG_1",
            "MATS_0003_Ly_P_WG_2",
            "MATS_0004_Ly_P_WG",
            "MATS_0004_Ly_P_WG_2",
            "MATS_0004_Ly_P_WG_MATS10_140133_Blasts",
            "MATS_0004_Ly_P_WG_MATS10_140133_HSC",
            "MATS_0004_Ly_P_WG_MATS10_140133_MPP1",
            "MATS_0004_Ly_P_WG_MATS10_140133_MPP2",
            "MATS_0009_Ly_P_WG_AML_01_S2",
            "MATS_0009_Ly_P_WG_MPN_01_S1",
            "MATS_0010_Ly_P_WG_AML_02_S5",
            "MATS_0010_Ly_P_WG_MPN_02_S4",
            "MATS_0011_Ly_P_WG_AML_03_S8",
            "MATS_0011_Ly_P_WG_MPN_03_S7",
            "MATS_0012_Ly_P_WG_AML_04_S3",
            "MATS_0012_Ly_P_WG_MPN_04_S2",
            "MATS_0013_Ly_P_WG_AML_06_S1",
            "MATS_0013_Ly_P_WG_MPN_06_S8",
            "MATS_0015_Ly_P_WG_AML_05_S6",
            "MATS_0015_Ly_P_WG_MATS05_140943_GMP",
            "MATS_0015_Ly_P_WG_MATS05_140943_HSC",
            "MATS_0015_Ly_P_WG_MATS05_140943_Myeloid",
            "MATS_0015_Ly_P_WG_MATS05_150472_HSC1",
            "MATS_0015_Ly_P_WG_MPN_05_S5")

for (sname in snames){
# read json file
json.file <- paste0("/Volumes/gsiprojects/dicklab/MATS/data/WG/sequenza/sols/",sname, "_alternative_solutions.json")
json.reads <- jsonlite::fromJSON(json.file)

gammas <- names(json.reads)

json.df<- c("gamma", "cellularity", "ploidy", "segs")
for (g in gammas){
  segs <- read.csv(paste0("/Volumes/gsiprojects/dicklab/MATS/data/WG/sequenza/segs/",sname, "/gammas/", g, "/", sname,"_Total_CN.seg"),
                   sep = "\t", as.is = T)
  s <- dim(segs)[1]
  c <- json.reads[[g]]$cellularity
  p <- json.reads[[g]]$ploidy
  # for (x in 1:length(c)){
  entr <- c(g, c[1], p[1], s)
    json.df <- rbind(json.df, entr)
  # }
}
colnames(json.df) <- json.df[1,]
json.df <- data.frame(json.df[-1,])

json.df <- json.df[order(as.numeric(json.df$gamma)),]
row.names(json.df) <- c()

json.df$gamma <- as.numeric(json.df$gamma)
json.df$cellularity <- as.numeric(json.df$cellularity)
json.df$ploidy <- as.numeric(json.df$ploidy)
json.df$segs <- as.numeric(json.df$segs)


# plot 
json.df.m <- melt(json.df, id = "gamma")
plt <- ggplot(json.df.m, aes(x = as.factor(gamma), y = value)) +
  geom_point() +
  geom_line(group = 1) +
  facet_grid(variable ~ ., scales = "free") + 
  theme_bw() + xlab ("gamma") #+ scale_x_continuous(gamma, labels = as.numeric(gamma), breaks = gamma)
plt

ggsave(paste0("~/Projects/GSI/Projects/MATS/WG/Sequenza/", sname, "_gamma_selector.png"), plt)
}



# settle the gamma for the steady value

# select solution

