## Script author: Gabriel E. B. de Barros
## Git: https://github.com/gabrielbarea

## Scripts in 'A review of the glacial environment arthropod trace fossils 
##  Umfolozia and Warvichnium with the description of new ichnotaxa'

## Scripts and data used for statistical analysis and graphical production of 
##  the paper.

#### Header ####
setwd() ## set directory
getwd() ## get the directory
Sys.setenv(LANG = "en") ## set the Rstudio to English

if(!require(devtools)) install.packages("devtools")
if(!require(pacman)) install.packages("pacman")
if(!require(magrittr)) install.packages("magrittr")
if(!require(MASS)) install.packages("MASS")
if(!require(factoextra)) install.packages("factoextra")
if(!require(broom)) install.packages("broom")
if(!require(dplyr)) install.packages("dplyr")
if(!require(plyr)) install.packages("plyr")
if(!require(readxl)) install.packages("readxl") 
if(!require(ggplot2)) install.packages("ggplot2")
if(!require(ggforce)) install.packages("ggforce")
if(!require(ggbiplot)) install.packages("ggbiplot")
if(!require(ggpubr)) install.packages("ggpubr")
if(!require(cowplot)) install.packages("cowplot")
library(devtools)
devtools::install_github('Mikata-Project/ggthemr', force = TRUE)
devtools::install_github("vqv/ggbiplot", force = TRUE)

pacman::p_load(magrittr, MASS, factoextra, broom, dplyr, readxl, ggplot2, 
               ggforce, ggbiplot, ggpubr, cowplot, ggthemr, plyr)

#### Reading the Data ####

data_warv <- read_excel("measurements_um_wa.xlsx", 
                        trim_ws = FALSE, sheet = "Warvichnium")

data_umfo <- read_excel("measurements_um_wa.xlsx", 
                        trim_ws = FALSE, sheet = "Umfolozia")

umfo_csv <- read.csv("umfolozia.csv")

warv_csv <- read.csv("warvichnium.csv")

ggthemr('flat')

#### Umfolozia CVA and PCA ####
#### CVA Umfolozia
lda.umfo <- lda(morph ~ TWEW + TRW + TAM + IW + ESL + ETRW + REIW + REELS, 
                data = umfo_csv)

lda.umfo$scaling

umfo.sub <- umfo_csv %>% dplyr::select(-morph) %>% # drop Species column 
    as.matrix # cast to matrix for calculations

# calculate CV scores
CVA.scores <- umfo.sub %*% lda.umfo$scaling

# create data frame with scores
umfo.CV <- data.frame(CVA.scores)
umfo.CV$morph <- umfo_csv$morph

umfo.cva.plot <- ggplot(umfo.CV, aes(x = LD1, y = LD2)) +
    geom_point(aes(color = morph, shape = morph), alpha = 1) +
    labs(x = "CV1", y = "CV2") + 
    coord_fixed(ratio = 1)

umfo.cva.plot

mu.u <- ddply(umfo.CV, "morph", summarise, umfo.mean=mean(LD1))
head(mu.u)

umfo.cva.density <- ggplot(umfo.CV, aes(x = LD1)) +
    geom_density(color = "black", aes(fill = morph), alpha = 0.3) +
    theme(panel.grid.major.x = element_blank()) +
    geom_vline(data = mu.u, aes(xintercept = umfo.mean, color = morph), 
               linetype = "dashed", size = 1) +
    labs(x = "CV1")
umfo.cva.density

#Chi
chi2 = qchisq(0.05, 2, lower.tail=FALSE)
chi2

CIregions.mean.and.pop <- umfo.CV %>% group_by(morph) %>%
    dplyr::summarize(CV1.mean = mean(LD1),
                     CV2.mean = mean(LD2),
                     mean.radii = sqrt(chi2/n()),
                     popn.radii = sqrt(chi2))
CIregions.mean.and.pop

#gg force
umfo.cva.plot2 <- umfo.cva.plot +
    geom_circle(data = CIregions.mean.and.pop, 
                mapping = aes(x0 = CV1.mean, y0 = CV2.mean, r = mean.radii,
                              color = morph, fill = morph), alpha =0.3,
                inherit.aes = FALSE) +
    geom_circle(data = CIregions.mean.and.pop,
                mapping = aes(x0 = CV1.mean, y0 = CV2.mean, r = popn.radii,
                              color = morph, fill = morph), alpha =0.1,
                linetype = "dashed",
                inherit.aes = FALSE)
umfo.cva.plot2

umfo.cva.plot2 <- umfo.cva.plot2 + geom_rug(aes(color = morph), sides = "b")

umfo.cva.plot2

arrange_U_CVA <- ggarrange(umfo.cva.plot2, umfo.cva.density, labels = "AUTO", 
                           font.label = list(face = "italic"), nrow = 2,
                           ncol = 1,
                           widths = c(5,1), heights = c(4,1.5))
arrange_U_CVA

#### PCA Umfolozia
res.pca <- prcomp(umfo_csv[,c(2:9)], center = TRUE, scale. = TRUE)
summary(res.pca)
str(res.pca)

ggbiplot(res.pca, ellipse = TRUE, ellipse.prob = 0.75, groups = umfo_csv$morph)

eig.val <- get_eigenvalue(res.pca)
eig.val

## Results for PCA Variables
res.var <- get_pca_var(res.pca)
res.var$coord          # Coordinates
res.var$contrib        # Contributions to the PCs
res.var$cos2           # Quality of representation 

#### Umfolozia boxplots ####
Umfolozia_TWEW <- 
    ggplot(data = data_umfo, aes(x = morph, y = TWEW, fill = morph)) + 
    geom_boxplot() +
    xlab("") + 
    ylab("Trackway external width (mm)") +
    theme(panel.grid.major.x = element_blank()) +
    theme(axis.title = element_text(size = 8)) +
    stat_summary(fun = mean, geom = "point", shape = 23, size = 4, 
                 colour = "black") +
    theme(legend.position = "none")

Umfolozia_ITEW <- 
    ggplot(data = data_umfo, aes(x = morph, y = ITEW, fill = morph)) + 
    geom_boxplot() +
    scale_x_discrete(limits = c("Ut")) +
    theme(panel.grid.major.x = element_blank()) +
    theme(axis.title = element_text(size = 8)) +
    xlab("") + 
    ylab("Internal track external width (mm)") +
    stat_summary(fun = mean, geom = "point", shape = 23, size = 4, 
                 colour = "black") +
    theme(legend.position = "none")

Umfolozia_TRW <- 
    ggplot(data = data_umfo, aes(x = morph, y = TRW, fill = morph)) + 
    geom_boxplot() +
    theme(panel.grid.major.x = element_blank()) +
    theme(axis.title = element_text(size = 8)) +
    xlab("") + 
    ylab("Track row width (mm)") +
    stat_summary(fun = mean, geom = "point", shape = 23, size = 4, 
                 colour = "black") +
    theme(legend.position = "none")

Umfolozia_SL <- 
    ggplot(data = data_umfo, aes(x = morph, y = SL, fill = morph)) + 
    geom_boxplot() +
    theme(panel.grid.major.x = element_blank()) +
    theme(axis.title = element_text(size = 8)) +
    scale_x_discrete(limits = c("Ut")) +
    xlab("") + 
    ylab("Internal track series lenght (mm)") +
    stat_summary(fun = mean, geom ="point", shape = 23, size = 4, 
                 colour = "black") +
    theme(legend.position="none")

Umfolozia_TAM <- 
    ggplot(data = data_umfo, aes(x = morph, y = TAM, fill = morph)) + 
    geom_boxplot() +
    theme(panel.grid.major.x = element_blank()) +
    theme(axis.title = element_text(size = 8)) +
    xlab("") + 
    ylab("Track row impression angle to midline (°)") +
    stat_summary(fun = mean, geom = "point", shape = 23, size = 4, 
                 colour = "black") +
    theme(legend.position="none")

Umfolozia_IW <- 
    ggplot(data = data_umfo, aes(x = morph, y = IW, fill = morph)) + 
    geom_boxplot() +
    theme(panel.grid.major.x = element_blank()) +
    theme(axis.title = element_text(size = 8)) +
    xlab("") + 
    ylab("Internal width (mm)") +
    stat_summary(fun = mean, geom = "point", shape = 23, size = 4, 
                 colour = "black") +
    theme(legend.position = "none")

Umfolozia_ITW <- 
    ggplot(data = data_umfo, aes(x = morph, y = ITW, fill = morph)) + 
    geom_boxplot() +
    theme(panel.grid.major.x = element_blank()) +
    theme(axis.title = element_text(size = 8)) +
    scale_x_discrete(limits = c("Ut")) +
    xlab("") + 
    ylab("Internal track width (mm)") +
    stat_summary(fun = mean, geom = "point", shape = 23, size = 4, 
                 colour = "black") +
    theme(legend.position = "none")

Umfolozia_ESL <- 
    ggplot(data = data_umfo, aes(x = morph, y = ESL, fill = morph)) + 
    geom_boxplot() +
    theme(panel.grid.major.x = element_blank()) +
    theme(axis.title = element_text(size = 8)) +
    xlab("") + 
    ylab("Series length (mm)") +
    stat_summary(fun = mean, geom = "point", shape = 23, size = 4, 
                 colour = "black") +
    theme(legend.position = "none")

Umfolozia_ETRW <- 
    ggplot(data = data_umfo, aes(x = morph, y = ETRW, fill = morph)) + 
    geom_boxplot() +
    theme(panel.grid.major.x = element_blank()) +
    theme(axis.title = element_text(size = 8)) +
    xlab("") + 
    ylab("Track row impression width (mm)") +
    stat_summary(fun = mean, geom = "point", shape = 23, size = 4, 
                 colour = "black") +
    theme(legend.position = "none")

Umfolozia_REIW <- 
    ggplot(data = data_umfo, aes(x = morph, y = REIW, fill = morph)) + 
    geom_boxplot() +
    theme(panel.grid.major.x = element_blank()) +
    theme(axis.title = element_text(size = 8)) +
    xlab("") + 
    ylab("Ratio of External / Internal width") +
    stat_summary(fun = mean, geom = "point", shape = 23, size = 4, 
                 colour = "black") +
    theme(legend.position = "none")

Umfolozia_RESL <- 
    ggplot(data = data_umfo, aes(x = morph, y = RESL, fill = morph)) + 
    geom_boxplot() +
    theme(panel.grid.major.x = element_blank()) +
    theme(axis.title = element_text(size = 8)) +
    scale_x_discrete(limits = c("Ut")) +
    xlab("") + 
    ylab("Ratio of External width / Internal series length") +
    stat_summary(fun = mean, geom = "point", shape = 23, size = 4, 
                 colour = "black") +
    theme(legend.position = "none")

Umfolozia_REELS <- 
    ggplot(data = data_umfo, aes(x = morph, y = REELS, fill = morph)) + 
    geom_boxplot() +
    theme(panel.grid.major.x = element_blank()) +
    theme(axis.title = element_text(size = 8)) +
    xlab("") + 
    ylab("Ratio of External width / Series length") +
    stat_summary(fun = mean, geom = "point", shape = 23, size = 4, 
                 colour = "black") +
    theme(legend.position = "none")

Umfolozia_RTWIT <- 
    ggplot(data = data_umfo, aes(x = morph, y = RTWIT, fill = morph)) + 
    geom_boxplot() +
    theme(panel.grid.major.x = element_blank()) +
    theme(axis.title = element_text(size = 8)) +
    scale_x_discrete(limits = c("Ut")) +
    xlab("") + 
    ylab("Ratio of External width / Internal track width") +
    stat_summary(fun = mean, geom = "point", shape = 23, size = 4, 
                 colour = "black") +
    theme(legend.position = "none")


arrange_umfolozia <-
    ggarrange(Umfolozia_TWEW, Umfolozia_TRW, Umfolozia_TAM, Umfolozia_IW, 
              Umfolozia_ESL, Umfolozia_ETRW, Umfolozia_REIW, Umfolozia_REELS,
              common.legend = TRUE, labels = "AUTO", 
              font.label = list(face = "italic"), nrow = 3,
              ncol = 3, widths = c(1,1,1), heights = c(1,1,1))
arrange_umfolozia

arrange_UT <- 
    ggarrange(Umfolozia_SL, Umfolozia_ITEW, Umfolozia_ITW, Umfolozia_RESL, 
              Umfolozia_RTWIT, labels = "AUTO", 
              font.label = list(face = "italic"), nrow = 1,
              ncol = 5, widths = c(1,1,1), heights = c(1,1))
arrange_UT


#### Warvichnium CVA and PCA ####
#### CVA Warvichnium
warv_pallete <- 
    define_palette(
        swatch = c("black", "#E74C3C", "#9B59B6", "#A46843"),
        gradient = c(lower = 'red', upper = 'green'),
        background = "#ECF0F1",
        text = c("#444444", "#444444"),
        line = c("#6e6e6e", "#6e6e6e"),
        gridline = "#BDC3C7")

ggthemr(warv_pallete)

lda.warv <- lda(morph ~ TWEW + TRW + SL + TAM + IW + REIW + RESL,
                data = warv_csv)

lda.warv$scaling

warv.sub <- warv_csv %>% dplyr::select(-morph) %>% # drop Species column 
   as.matrix # cast to matrix for calculations

# calculate CV scores
CVA.scores <- warv.sub %*% lda.warv$scaling

# create data frame with scores
warv.CV <- data.frame(CVA.scores)
warv.CV$morph <- warv_csv$morph

warv.cva.plot <- 
    ggplot(warv.CV, aes(x = LD1, y = LD2)) +
    geom_point(aes(color = morph, shape = morph), alpha = 1) +
    labs(x = "CV1", y = "CV2") +
    coord_fixed(ratio = 1)

mu.w <- ddply(warv.CV, "morph", summarise, warv.mean = mean(LD1))
head(mu.w)

warv.cva.density <- 
    ggplot(warv.CV, aes(x = LD1)) +
    geom_density(color = "black", aes(fill = morph), alpha = 0.3) +
    theme(panel.grid.major.x = element_blank()) +
    geom_vline(data = mu.w, 
               aes(xintercept = warv.mean, color = morph), 
               linetype = "dashed", size = 1) +
    labs(x = "CV1")
warv.cva.density

#Chi
chi2 = qchisq(0.05, 2, lower.tail = FALSE)
chi2

CIregions.mean.and.pop <- warv.CV %>% group_by(morph) %>%
    dplyr::summarize(CV1.mean = mean(LD1),
             CV2.mean = mean(LD2),
             mean.radii = sqrt(chi2/n()),
             popn.radii = sqrt(chi2))
CIregions.mean.and.pop

#gg force
warv.cva.plot2 <- warv.cva.plot +
    geom_circle(data = CIregions.mean.and.pop,
                mapping = aes(x0 = CV1.mean, y0 = CV2.mean, r = mean.radii,
                              color = morph, fill = morph), alpha =0.3,
                inherit.aes = FALSE) +
    geom_circle(data = CIregions.mean.and.pop,
                mapping = aes(x0 = CV1.mean, y0 = CV2.mean, r = popn.radii,
                              color = morph, fill = morph), alpha =0.1,
                linetype = "dashed",
                inherit.aes = FALSE)
warv.cva.plot2

warv.cva.plot2 <- warv.cva.plot2 + geom_rug(aes(color = morph), sides = "b")
warv.cva.plot2

arrange_W_CVA <-
    ggarrange(warv.cva.plot2, warv.cva.density,
              labels = "AUTO", font.label = list(face = "italic"), 
              nrow = 2, ncol = 1,
              widths = c(5,1), heights= c(4,1.5))
arrange_W_CVA

#### PCA Warvichnium
res.pca <- prcomp(warv_csv[,c(2:8)], center = TRUE, scale. = TRUE)
summary(res.pca)
str(res.pca)

ggbiplot(res.pca, ellipse=TRUE, ellipse.prob = 0.75, groups = warv_csv$morph)

eig.val <- get_eigenvalue(res.pca)
eig.val

## Results for PCA Variables
res.var <- get_pca_var(res.pca)
res.var$coord          # Coordinates
res.var$contrib        # Contributions to the PCs
res.var$cos2           # Quality of representation 

#### Warvichnium boxplots ####
Warvichnium_TWEW <- 
    ggplot(data = data_warv, aes(x = morph,y = TWEW, fill = morph)) +
    geom_boxplot() +
    scale_x_discrete(limits = c("Is", "Ip", "Wu")) +
    theme(panel.grid.major.x = element_blank()) +
    theme(axis.title = element_text(size = 8)) +
    xlab("") + 
    ylab("Trackway external width (mm)") +
    stat_summary(fun = mean, geom = "point", shape = 23, size = 4, 
                 colour = "black") +
    theme(legend.position = "none")

Warvichnium_MTEW <- 
    ggplot(data = data_warv, aes(x = morph, y = MTEW, fill = morph)) +
    geom_boxplot() +
    theme(panel.grid.major.x = element_blank()) +
    theme(axis.title = element_text(size = 8)) +
    scale_x_discrete(limits = c("Wu")) +
    xlab("") + 
    ylab("Medial track external width (mm)") +
    stat_summary(fun = mean, geom = "point", shape = 23, size = 4, 
                 colour = "black") +
    theme(legend.position = "none")

Warvichnium_ITEW <- 
    ggplot(data = data_warv, aes(x = morph, y = ITEW, fill = morph)) +
    geom_boxplot() +
    theme(panel.grid.major.x = element_blank()) +
    theme(axis.title = element_text(size = 8)) +
    scale_x_discrete(limits = c("Is")) +
    xlab("") + 
    ylab("Internal track external width (mm)") +
    stat_summary(fun = mean, geom = "point", shape = 23, size = 4, 
                 colour = "black") +
    theme(legend.position = "none")

Warvichnium_IIL <- 
    ggplot(data = data_warv, aes(x = morph, y = IIL, fill = morph)) +
    geom_boxplot() +
    theme(panel.grid.major.x = element_blank()) +
    theme(axis.title = element_text(size = 8)) +
    scale_x_discrete(limits = c("Wu")) +
    xlab("") + 
    ylab("Internal imprint lenght (mm)") +
    stat_summary(fun = mean, geom = "point", shape = 23, size = 4, 
                 colour = "black") +
    theme(legend.position = "none")

Warvichnium_MIAM <- 
    ggplot(data = data_warv, aes(x = morph, y = MIAM, fill = morph)) +
    geom_boxplot() +
    theme(panel.grid.major.x = element_blank()) +
    theme(axis.title = element_text(size = 8)) +
    scale_x_discrete(limits = c("Wu")) +
    xlab("") + 
    ylab("Middle imprint angle in relation to midline (°)") +
    stat_summary(fun = mean, geom = "point", shape = 23, size = 4, 
                 colour = "black") +
    theme(legend.position = "none")

Warvichnium_ERW <- 
    ggplot(data = data_warv, aes(x = morph, y = ERW, fill = morph)) +
    geom_boxplot() +
    theme(panel.grid.major.x = element_blank()) +
    theme(axis.title = element_text(size = 8)) +
    scale_x_discrete(limits = c("Wu")) +
    xlab("") +
    ylab("External track row width (mm)") +
    stat_summary(fun = mean, geom = "point", shape = 23, size = 4, 
                 colour = "black") +
    theme(legend.position = "none")

Warvichnium_MISL <- 
    ggplot(data = data_warv, aes(x = morph, y = MISL, fill = morph)) +
    geom_boxplot() +
    theme(panel.grid.major.x = element_blank()) +
    theme(axis.title = element_text(size = 8)) +
    scale_x_discrete(limits = c("Wu")) +
    xlab("") + 
    ylab("Medial track row series lenght (mm)") +
    stat_summary(fun = mean, geom = "point", shape = 23, size = 4, 
                 colour = "black") +
    theme(legend.position = "none")

Warvichnium_TRW <- 
    ggplot(data = data_warv, aes(x = morph, y = TRW, fill = morph)) +
    geom_boxplot() +
    theme(panel.grid.major.x = element_blank()) +
    theme(axis.title = element_text(size = 8)) +
    xlab("") + 
    ylab("Track row width (mm)") +
    stat_summary(fun = mean, geom = "point", shape = 23, size = 4, 
                 colour = "black") +
    theme(legend.position = "none")

Warvichnium_SL <- 
    ggplot(data = data_warv, aes(x = morph, y = SL, fill = morph)) +
    geom_boxplot() +
    theme(panel.grid.major.x = element_blank()) +
    theme(axis.title = element_text(size = 8)) +
    xlab("") + 
    ylab("Series length (mm)") +
    stat_summary(fun = mean, geom = "point", shape = 23, size = 4, 
                 colour = "black") +
    theme(legend.position = "none")

Warvichnium_TAM <- 
    ggplot(data = data_warv, aes(x = morph, y = TAM, fill = morph)) +
    geom_boxplot() +
    theme(panel.grid.major.x = element_blank()) +
    theme(axis.title = element_text(size = 8)) +
    xlab("") + 
    ylab("Track impression angle to midline (°)") +
    stat_summary(fun = mean, geom = "point", shape = 23, size = 4, 
                 colour = "black") +
    theme(legend.position = "none")

Warvichnium_IW <- 
    ggplot(data = data_warv, aes(x = morph, y = IW, fill = morph)) +
    geom_boxplot() +
    theme(panel.grid.major.x = element_blank()) +
    theme(axis.title = element_text(size = 8)) +
    xlab("") + 
    ylab("Internal width (mm)") +
    stat_summary(fun = mean, geom = "point", shape = 23, size = 4, 
                 colour = "black") +
    theme(legend.position = "none")

Warvichnium_ITW <- 
    ggplot(data = data_warv, aes(x = morph, y = ITW, fill = morph)) +
    geom_boxplot() +
    theme(panel.grid.major.x = element_blank()) +
    theme(axis.title = element_text(size = 8)) +
    scale_x_discrete(limits = c("Is")) +
    xlab("") + 
    ylab("Internal track width (mm)") +
    stat_summary(fun = mean, geom = "point", shape = 23, size = 4, 
                 colour = "black") +
    theme(legend.position = "none")

Warvichnium_ETRW <- 
    ggplot(data = data_warv, aes(x = morph, y = ETRW, fill = morph)) +
    geom_boxplot() +
    theme(panel.grid.major.x = element_blank()) +
    theme(axis.title = element_text(size = 8)) +
    scale_x_discrete(limits = c("Is")) +
    xlab("") + 
    ylab("External track row imprint width (mm)") +
    stat_summary(fun = mean, geom = "point", shape = 23, size = 4, 
                 colour = "black") +
    theme(legend.position = "none")

Warvichnium_REIW <- 
    ggplot(data = data_warv, aes(x = morph, y = REIW, fill = morph)) +
    geom_boxplot() +
    theme(panel.grid.major.x = element_blank()) +
    theme(axis.title = element_text(size = 8)) +
    xlab("") + 
    ylab("Ratio of External / Internal width") +
    stat_summary(fun = mean, geom = "point", shape = 23, size = 4, 
                 colour = "black") +
    theme(legend.position = "none")

Warvichnium_RESL <- 
    ggplot(data = data_warv, aes(x = morph, y = RESL, fill = morph)) +
    geom_boxplot() +
    theme(panel.grid.major.x = element_blank()) +
    theme(axis.title = element_text(size = 8)) +
    xlab("") + 
    ylab("Ratio of External width / Internal series length") +
    stat_summary(fun = mean, geom = "point", shape = 23, size = 4, 
                 colour = "black") +
    theme(legend.position = "none")

Warvichnium_RTWIT <- 
    ggplot(data = data_warv, aes(x = morph, y = RTWIT, fill = morph)) +
    geom_boxplot() +
    theme(panel.grid.major.x = element_blank()) +
    theme(axis.title = element_text(size = 8)) +
    scale_x_discrete(limits = c("Is")) +
    xlab("") + 
    ylab("Ratio of External width / Internal track width") +
    stat_summary(fun = mean, geom = "point", shape = 23, size = 4, 
                 colour = "black") +
    theme(legend.position = "none")

arrange_warvichnium <-
    ggarrange(Warvichnium_TWEW, Warvichnium_TRW, Warvichnium_SL, 
              Warvichnium_TAM, Warvichnium_IW, Warvichnium_REIW, 
              Warvichnium_RESL, 
              common.legend = TRUE, labels = "AUTO", 
              font.label = list(face = "italic"))
arrange_warvichnium

arrange_SA <-
    ggarrange(Warvichnium_MTEW, Warvichnium_IIL, Warvichnium_MIAM, 
              Warvichnium_ERW, Warvichnium_MISL,Warvichnium_ITEW, 
              Warvichnium_ITW, Warvichnium_ETRW, Warvichnium_RTWIT,
              labels = "AUTO", font.label = list(face = "italic"), 
              nrow = 2, ncol = 5, widths = c(1,1,1))
arrange_SA


#### Statistical tests ####
#### Wilcoxon-Mann-Whitney ####
kruskal.test(TWEW + TRW + TAM + IW + ESL + ETRW + REIW + REELS ~ morph, data = umfo_csv)

pairwise.wilcox.test(umfo_csv$TWEW, umfo_csv$morph, p.adjust.method = "BH")
pairwise.wilcox.test(umfo_csv$TRW, umfo_csv$morph, p.adjust.method = "BH")
pairwise.wilcox.test(umfo_csv$TAM, umfo_csv$morph, p.adjust.method = "BH")
pairwise.wilcox.test(umfo_csv$IW, umfo_csv$morph, p.adjust.method = "BH")
pairwise.wilcox.test(umfo_csv$ESL, umfo_csv$morph, p.adjust.method = "BH")
pairwise.wilcox.test(umfo_csv$ETRW, umfo_csv$morph, p.adjust.method = "BH")
pairwise.wilcox.test(umfo_csv$REIW, umfo_csv$morph, p.adjust.method = "BH")
pairwise.wilcox.test(umfo_csv$REELS, umfo_csv$morph, p.adjust.method = "BH")


kruskal.test(TWEW + TRW + SL + TAM + IW + REIW + RESL ~ morph, data = warv_csv)

pairwise.wilcox.test(warv_csv$TWEW, warv_csv$morph, p.adjust.method = "BH")
pairwise.wilcox.test(warv_csv$TRW, warv_csv$morph, p.adjust.method = "BH")
pairwise.wilcox.test(warv_csv$SL, warv_csv$morph, p.adjust.method = "BH")
pairwise.wilcox.test(warv_csv$TAM, warv_csv$morph, p.adjust.method = "BH")
pairwise.wilcox.test(warv_csv$IW, warv_csv$morph, p.adjust.method = "BH")
pairwise.wilcox.test(warv_csv$REIW, warv_csv$morph, p.adjust.method = "BH")
pairwise.wilcox.test(warv_csv$RESL, warv_csv$morph, p.adjust.method = "BH")

#### Paper Information ####

## A review of the glacial environment arthropod trace fossils Umfolozia 
##  and Warvichnium with description of new ichnotaxa

## Gabriel E. B. de Barros; Bernardo de C. P. e M. Peixoto; João H. D. Lima ; 
##  Nicholas J. Minter; Daniel Sedorko

## gbareabarros@gmail.com (G.E.B. de Barros, Corresponding author)