install.packages("demoniche", repos="http://R-Forge.R-project.org")

library(demoniche)

?demoniche_model
View(demoniche_model)


data(Hmontana)


noCC_nodispersal <- demoniche_model(modelname = "Hmontana", Niche = FALSE,
                                    + Dispersal = FALSE, repetitions = 2,
                                    + foldername = "noCC_nodispersal")

demoniche_setup

## create population grid
Populations_mine <- data.frame(patchID = seq(1:900),
                               X = rep(1:30, 30),
                               Y = rep(1:30, each = 30),
                               area = rep(1, 1))

library(popbio)
data(hudvrs)
data(hudsonia)
matrices_mine <- cbind(meanmatrix = as.vector(hudmxdef(hudvrs$mean)),
                      sapply(hudsonia, unlist))
colnames(matrices_mine) <- c("Reference_matrix", "Matrix_1", "Matrix_2", "Matrix_3", "Matrix_4")

matrices_var_mine <- matrix(0.01, ncol = 1, nrow = nrow(matrices_mine), dimnames = list(NULL, "sd"))

stages_mine <- colnames(hudsonia$A85)

density_individuals_mine <- 20000

proportion_initial_mine <- c(0.9818098089, 0.0006907668, 0.0069076675,
                             0.0036840893, 0.0057563896, 0.0011512779)

no_yrs_mine <- 4

Nichemap_mine <- data.frame(gridID = seq(1:900),
                            X = rep(1:30, 30),
                            Y = rep(1:30, each = 30),
                            period_2000 = rep(0.8, 900),
                            period_2010 = rep(0.7, 900),
                            period_2020 = rep(0.6, 900),
                            period_2030 = rep(0.5, 900))

niche_formulas <- as.formula(paste(paste(colnames(Nichemap_mine)[-c(1:3)],
                                         collapse="+"),"X+Y",sep="~"))
print(levelplot(niche_formulas, Nichemap_mine, col.regions=rev(heat.colors(100)),
                main = "Niche Values"))

## run minimal setup
demoniche_setup(modelname = "Hmontana_minimal",
                Populations = Populations_mine, 
                matrices_var = matrices_var_mine,
                matrices = matrices_mine,
                stages = stages_mine,
                proportion_initial = proportion_initial_mine,
                density_individuals = density_individuals_mine,
                no_yrs = no_yrs_mine,
                Nichemap = Nichemap_mine)

example_minimal <- demoniche_model(modelname = "Hmontana_minimal", Niche = TRUE,
                                   Dispersal = FALSE, repetitions = 1,
                                   foldername = "Hmontana_minimal")


View(Hmontana)
demoniche_model(modelname = "Hmontana",
                Niche = FALSE,
                Dispersal = FALSE, repetitions = 1,
                foldername = "Hmontana_minimal")



