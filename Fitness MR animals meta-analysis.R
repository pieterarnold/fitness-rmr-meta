#### Fitness and metabolic rate in animals meta-analysis script ####
#### Authors: Pieter A. Arnold, Steven Delean, Phill Cassey, Craig R. White ####


#### Load packages ####
library(nlme)
library(car)
library(lme4)
library(metafor)
library(ggplot2)
library(dplyr)


#### Load in colour palettes ####
# Palette URL: http://paletton.com/#uid=72Y0u0khy-e7yT7cSJTlxtFqRqX
# Primary color:
pal.green  <- c("#B6EDCC", "#82D9A5", "#57C182", "#38AA66", "#199A4D")
# Secondary color (1):
pal.blue   <- c("#B4D8E8", "#7EB5CD", "#5392AD", "#357895", "#1B6787")
# Secondary color (2):
pal.orange <- c("#FFE3C3", "#FFCF98", "#FFBD73", "#ECA14D", "#D78223")
# Complement color:
pal.red <- c("#FFCEC3", "#FFAC98", "#FF8E73", "#EC6C4D", "#D74523")


#### Load in and subset data ####

## Import data and subsetting
rawdata <- read.csv("MR_Fitness_Data_revised.csv")
str(rawdata)

rawdata <- within(rawdata, {
    Phylum <- as.factor(Phylum)
    Class <- as.factor(Class)
    Order <- as.factor(Order)
    Family <- as.factor(Family)
    Species <- as.factor(Species)
    FitnessComponent <- as.factor(FitnessComponent)
    FitnessClassification <- as.factor(FitnessClassification)
    MRmeasurement <- as.factor(MRmeasurement)
    LabField <- as.factor(LabField)
    MR.Temp <- as.numeric(MR.Temp)
    Mean.MR <- as.numeric(Mean.MR)
    Mean.Body.mass <- as.numeric(Mean.Body.mass)
    n.rep <- as.numeric(n.rep)
    p <- as.numeric(p)
    Reference <- as.factor(Reference)
    Journal <- as.factor(Journal)
})
str(rawdata)
summary(rawdata)

## Add columns and Subset data
rawdata$logMR <- log10(rawdata$Mean.MR)
rawdata$logM <- log10(rawdata$Mean.Body.mass)
rawdata$log.n.rep <- log10(rawdata$n.rep)
rawdata <- droplevels(subset(rawdata, !is.na(r)))
rawdata <- droplevels(subset(rawdata, !is.na(n.rep)))


## Create thermoregulation variable
unique(rawdata$Class)
rawdata$thermo <- as.factor(ifelse(as.character(rawdata$Class) %in% c("Aves", "Mammalia"), "Endotherm", "Ectotherm"))

#### Calculate transformed effect sizes and precision ####
meta.dat <- escalc(measure = "ZCOR", ri = r, ni = n.rep, data = rawdata, replace = FALSE)
summary(as.data.frame(meta.dat))


#### Subset data to remove data without body mass values ####
meta.dat.subMass <- droplevels(subset(meta.dat, !is.na(meta.dat$Mean.Body.mass)))
meta.dat.subMass$FitnessClassification <- as.factor(meta.dat.subMass$FitnessClassification)

## Create interaction thermoregulation by fitness
with(meta.dat.subMass, table(thermo, FitnessClassification))
meta.dat.subMass$fitness.endo <- with(meta.dat.subMass, interaction(FitnessClassification, thermo, drop = TRUE))

## Observation-level effect
meta.dat.subMass$obs <- seq(nrow(meta.dat.subMass))

# Remove the Dassis sea lion data row due to its extreme influence on the model selection
# due to extreme values of logM, logMR and Zr
meta.dat.subMass <- meta.dat.subMass[-112,]

## Some summaries
with(meta.dat.subMass, table(thermo, LabField))
# Labfield is poorly represented across most categories so we will not include it in the meta-analysis
with(meta.dat.subMass, table(FitnessClassification, LabField))
with(meta.dat.subMass, table(FitnessClassification, thermo))
lapply(split(meta.dat.subMass$logM, meta.dat.subMass$thermo), summary)


# Table S2: Summary
meta.dat.subMass %>%
    group_by(Class) %>%
    summarise(Count = n_distinct(Species), Minimum = min(n.rep, na.rm = TRUE), Maximum = max(n.rep, na.rm = TRUE), .groups = "drop")

(lightest.ecto <- min(meta.dat.subMass$Mean.Body.mass[meta.dat.subMass$thermo=="Ectotherm"]))
(heaviest.ecto <- max(meta.dat.subMass$Mean.Body.mass[meta.dat.subMass$thermo=="Ectotherm"]))

(lightest.endo <- min(meta.dat.subMass$Mean.Body.mass[meta.dat.subMass$thermo=="Endotherm"]))
(heaviest.endo <- max(meta.dat.subMass$Mean.Body.mass[meta.dat.subMass$thermo=="Endotherm"]))

# Table S3: summary
with(meta.dat.subMass, table(Class, FitnessClassification))


#### All and null meta-regression models ####

# with a priori determined combinations of moderators ####

# All moderator variables
meta.remv.allVars <- rma.mv(yi = yi, V = vi, mods = ~ FitnessClassification*thermo + logM*thermo, data = meta.dat.subMass, random = list(~1|Ref, ~1|obs))
summary(meta.remv.allVars)

# Null model
meta.remv.null <- rma.mv(yi = yi, V = vi, mods = ~ 1, data = meta.dat.subMass, random = list(~1|Ref, ~1|obs))
summary(meta.remv.null)


#### Main effects separately - meta-regression models ####

# Fitness class
meta.remv.fit <- rma.mv(yi = yi, V = vi, mods = ~ FitnessClassification, data = meta.dat.subMass, random = list(~1|Ref, ~1|obs))
summary(meta.remv.fit)

# Thermoregulatory strategy
meta.remv.thermo <- rma.mv(yi = yi, V = vi, mods = ~ thermo, data = meta.dat.subMass, random = list(~1|Ref, ~1|obs))
summary(meta.remv.thermo)

# Body mass
meta.remv.mass <- rma.mv(yi = yi, V = vi, mods = ~ logM, data = meta.dat.subMass, random = list(~1|Ref, ~1|obs))
summary(meta.remv.mass)


#### Interactions meta-regression models ####

# Fitness x thermo
meta.remv.fitTh <- rma.mv(yi = yi, V = vi, mods = ~ FitnessClassification*thermo, data = meta.dat.subMass, random = list(~1|Ref, ~1|obs))
summary(meta.remv.fitTh)

# Mass x thermo
meta.remv.massTh <- rma.mv(yi = yi, V = vi, mods = ~ logM*thermo, data = meta.dat.subMass, random = list(~1|Ref, ~1|obs))
summary(meta.remv.massTh)

## Fitness x mass
meta.remv.massFit <- rma.mv(yi = yi, V = vi, mods = ~ logM*FitnessClassification, data = meta.dat.subMass, random = list(~1|Ref, ~1|obs))
summary(meta.remv.massFit)

meta.remv.fitTh.pM <- rma.mv(yi = yi, V = vi, mods = ~ FitnessClassification*thermo + logM, data = meta.dat.subMass, random = list(~1|Ref, ~1|obs))
summary(meta.remv.fitTh.pM)

## Multiple interactions (same as all interactions model)
meta.remv.fitTh.massTh <- rma.mv(yi = yi, V = vi, mods = ~ FitnessClassification*thermo + logM*thermo, data = meta.dat.subMass, random = list(~1|Ref, ~1|obs))
summary(meta.remv.fitTh.massTh)


#### Extract AICc values from each model ####

## Extract AICc values and create summary table to rank models
aiccVals <- unlist(lapply(list(meta.remv.null, meta.remv.fit, meta.remv.thermo, meta.remv.mass, meta.remv.fitTh,
                               meta.remv.massTh, meta.remv.massFit, meta.remv.fitTh.massTh, meta.remv.fitTh.pM),
                          function(x) x$fit.stats[5, 1]))

aiccTab <- as.data.frame(cbind(aiccVals, aiccVals - min(aiccVals)))
names(aiccTab) <- c("AICc", "deltaAICc")
rownames(aiccTab) <- c("null", "fit", "thermo", "mass", "fitTh",
                       "massTh", "massFit", "fitTh.massTh", "fitTh.pM")

ze <- with(aiccTab, exp(-(aiccTab$deltaAICc/2)))
aiccTab$AICcWeight <- ze/sum(ze)
aiccTab$AICcWeightr <- round(aiccTab$AICcWeight, 3)
aiccTab <- aiccTab[order(aiccTab$deltaAICc), ]

#### Table S4 ####
aiccTab



#### Final model for interpretation and Table S5 ####

# Fitness x thermo
meta.remv.fitTh <- rma.mv(yi = yi, V = vi, mods = ~ FitnessClassification*thermo, data = meta.dat.subMass, random = list(~1|Ref, ~1|obs))
summary(meta.remv.fitTh)


#### Diagnostics on final model ####
# Plot standardised residuals
par(mfrow = c(2, 2), mar = c(4, 4, 2, 1) + 0.1)
plot(rstandard(meta.remv.fitTh)$z ~ fitted(meta.remv.fitTh))
qqnorm(rstandard(meta.remv.fitTh)$z)
qqline(rstandard(meta.remv.fitTh)$z)
par(mfrow = c(1, 1), mar = c(5, 4, 4, 2) + 0.1)

# Leverage
mod.hat <- hatvalues(meta.remv.fitTh)

# Cook's Distance
mod.cd <- cooks.distance(meta.remv.fitTh) # takes a long time to run

# Plot leverage vs. std.residual
plot(mod.hat, rstandard(meta.remv.fitTh)$z)
plot(mod.cd, rstandard(meta.remv.fitTh)$z)


#### Heterogeneity statistics for final model ####

# I^2 statistics for highest ranked moderator model
# Measurement error
Weight <- 1/meta.dat.subMass$vi
MV <- sum((Weight*(length(Weight)-1))/(sum(Weight)^2-sum(Weight^2)))
MV

# I2 for Study
(i2s <- 100*(meta.remv.fitTh$sigma2[1])/(meta.remv.fitTh$sigma2[1] + meta.remv.fitTh$sigma2[2] + MV))

# I2 for residuals (within Authors)
(i2ws <- 100*(meta.remv.fitTh$sigma2[2])/(meta.remv.fitTh$sigma2[1] + meta.remv.fitTh$sigma2[2] + MV))

# I2 for measurement error
(i2e <- 100*(MV)/(meta.remv.fitTh$sigma2[1] + meta.remv.fitTh$sigma2[2] + MV))


#### Heterogeneity statistics for null model ####

# I2 statistics for null moderator model
# Measurement error
Weight <- 1/meta.dat.subMass$vi
MV <- sum((Weight*(length(Weight)-1))/(sum(Weight)^2-sum(Weight^2)))
MV

# I2 for Study
(i2s_null <- 100*(meta.remv.null$sigma2[1])/(meta.remv.null$sigma2[1] + meta.remv.null$sigma2[2] + MV))

# I2 for residuals (within Authors)
(i2ws_null <- 100*(meta.remv.null$sigma2[2])/(meta.remv.null$sigma2[1] + meta.remv.null$sigma2[2] + MV))

# I2 for measurement error
(i2e_null <- 100*(MV)/(meta.remv.null$sigma2[1] + meta.remv.null$sigma2[2] + MV))

# (almost same as for above - so moderators account for variation at both between and within study level)


#### Summary statistics for overall effects (Zr) ####

# Summary statistics for overall effects (Zr): Across all data
mean.yi_all <- mean(meta.dat.subMass$yi)
SE.yi_all <- sqrt(var(meta.dat.subMass$yi))/sqrt(355)
(E.yi_all <- qnorm(0.975)*SE.yi_all)
mean.yi_all + c(-E.yi_all, 0, +E.yi_all)

# Summary statistics for overall effects (Zr): Only where Zr was reported as significant
mean.yi_sig <- mean(meta.dat.subMass$yi[meta.dat.subMass$sig=="1"])
SE.yi_sig <- sqrt(var(meta.dat.subMass$yi[meta.dat.subMass$sig=="1"]))/sqrt(length(meta.dat.subMass$yi[meta.dat.subMass$sig=="1"]))
(E.yi_sig <- qnorm(0.975)*SE.yi_sig)
mean.yi_sig + c(-E.yi_sig, 0, +E.yi_sig)

# Summary statistics for overall effects (Zr): Only where Zr was reported as non-significant
mean.yi_nonsig <- mean(meta.dat.subMass$yi[meta.dat.subMass$sig=="0"])
SE.yi_nonsig <- sqrt(var(meta.dat.subMass$yi[meta.dat.subMass$sig=="0"]))/sqrt(length(meta.dat.subMass$yi[meta.dat.subMass$sig=="0"]))
(E.yi_nonsig <- qnorm(0.975)*SE.yi_nonsig)
mean.yi_nonsig + c(-E.yi_nonsig, 0, +E.yi_nonsig)

# Summary statistics for overall effects (Zr): Excluding active MR
mean.yi_exAMR <- mean(meta.dat.subMass$yi[!meta.dat.subMass$FitnessClassification=="Active MR"])
SE.yi_exAMR <- sqrt(var(meta.dat.subMass$yi[!meta.dat.subMass$FitnessClassification=="Active MR"]))/
  sqrt(length(meta.dat.subMass$yi[!meta.dat.subMass$FitnessClassification=="Active MR"]))
(E.yi_exAMR <- qnorm(0.975)*SE.yi_exAMR)
mean.yi_exAMR + c(-E.yi_exAMR, 0, +E.yi_exAMR)

# Summary statistics for overall effects (Zr): Only of fitness-*related* traits (categories 1-6)
mean.yi_related <- mean(meta.dat.subMass$yi[meta.dat.subMass$FitnessClassification!="Survival" |
                                      meta.dat.subMass$FitnessClassification!="Reproduction"])
SE.yi_related <- sqrt(var(meta.dat.subMass$yi[meta.dat.subMass$FitnessClassification!="Survival" |
                                        meta.dat.subMass$FitnessClassification!="Reproduction"]))/
  sqrt(length(meta.dat.subMass$yi[meta.dat.subMass$FitnessClassification!="Survival" |
                                    meta.dat.subMass$FitnessClassification!="Reproduction"]))
(E.yi_related <- qnorm(0.975)*SE.yi_related)
mean.yi_related + c(-E.yi_related, 0, +E.yi_related)

# Summary statistics for overall effects (Zr): Only of 'true' fitness (categories 7-8)
mean.yi_fitness <- mean(meta.dat.subMass$yi[meta.dat.subMass$FitnessClassification=="Reproduction" |
                                      meta.dat.subMass$FitnessClassification=="Survival"])
SE.yi_fitness <- sqrt(var(meta.dat.subMass$yi[meta.dat.subMass$FitnessClassification=="Reproduction" |
                                        meta.dat.subMass$FitnessClassification=="Survival"]))/
  sqrt(length(meta.dat.subMass$yi[meta.dat.subMass$FitnessClassification=="Survival" |
                                    meta.dat.subMass$FitnessClassification=="Reproduction"]))
(E.yi_fitness <- qnorm(0.975)*SE.yi_fitness)
mean.yi_fitness + c(-E.yi_fitness, 0, +E.yi_fitness)

# Summary statistics for overall effects (Zr): Only of ectotherms only
mean.yi_ecto <- mean(meta.dat.subMass$yi[meta.dat.subMass$thermo=="Ectotherm"])
SE.yi_ecto <- sqrt(var(meta.dat.subMass$yi[meta.dat.subMass$thermo=="Ectotherm"]))/
  sqrt(length(meta.dat.subMass$yi[meta.dat.subMass$thermo=="Ectotherm"]))
E.yi_ecto <- qnorm(0.975)*SE.yi_ecto
mean.yi_ecto + c(-E.yi_ecto, 0, +E.yi_ecto)

# Summary statistics for overall effects (Zr): Only of endotherms only
mean.yi_endo <- mean(meta.dat.subMass$yi[meta.dat.subMass$thermo=="Endotherm"])
SE.yi_endo <- sqrt(var(meta.dat.subMass$yi[meta.dat.subMass$thermo=="Endotherm"]))/
  sqrt(length(meta.dat.subMass$yi[meta.dat.subMass$thermo=="Endotherm"]))
E.yi_endo <- qnorm(0.975)*SE.yi_endo
mean.yi_endo + c(-E.yi_endo, 0, +E.yi_endo)


#### Figure 3: Fitness x thermo interaction main plot ####

## Create new data for prediction of fitness*thermo interaction
newd <- with(meta.dat.subMass, expand.grid(levels(FitnessClassification), levels(thermo)))
names(newd) <- c("FitnessClassification", "thermo")

meta.remv.fitTh <- rma.mv(yi = yi, V = vi, mods = ~ FitnessClassification + thermo +
                            FitnessClassification:thermo,
                          data = meta.dat.subMass, random = list(~1|Ref, ~1|obs))

newd$logM <- median(meta.dat.subMass$logM[meta.dat.subMass$thermo == "Endotherm"], na.rm = TRUE)
newdV <- newd

## Coerce new data to the appropriate model matrix
newd <- as.data.frame(model.matrix(~ FitnessClassification * thermo, data = newd)[, -1])
newd

## Source adapted predict function
source("predictRMANew.R")

## Predict marginal means from model for new data values
preds.fitTh <- cbind(newdV, predict.rmaNew(meta.remv.fitTh, as.matrix(newd)))
preds.fitTh

levels(preds.fitTh$FitnessClassification)[length(levels(preds.fitTh$FitnessClassification)) + 1] <- ""
preds.fitTh1 <- preds.fitTh[1:6,]
blankrow1 <- data.frame(preds.fitTh[1,])
blankrow1[1,1] <- NA
blankrow1[1,2] <- "Ectotherm"
blankrow1[1,3:9] <- NA
preds.fitTh2 <- preds.fitTh[7:14,]
blankrow2 <- data.frame(preds.fitTh[1,])
blankrow2[1,1] <- NA
blankrow2[1,2] <- "Endotherm"
blankrow2[1,3:9] <- NA
preds.fitTh3 <- preds.fitTh[15:16,]
preds.fitTh_final <- rbind(preds.fitTh1, blankrow1, preds.fitTh2, blankrow2, preds.fitTh3)

fills <- c("steelblue2", "orangered")
pd <- position_dodge(0.6)
p <- ggplot(data = preds.fitTh_final, aes(x = FitnessClassification, y = pred, colour = thermo, group = thermo)) +
  geom_abline(intercept = 0, slope = 0, lty = 2) +
  geom_point(size = 3, position = pd) + expand_limits(y = c(-0.5, 1)) +
  geom_errorbar(aes(ymin = ci.lb, ymax = ci.ub), width = 0.25, position = pd) +
  xlab("                       Fitness-related traits                                               Fitness") +
  ylab(expression(paste(plain("Effect size ("), italic("Zr"), plain(")")))) +
  theme_classic(base_size = 12) + theme(legend.position = c(0.85,0.9), legend.title = element_blank()) +
  theme(axis.text.x = element_text(angle = 45, size = 11, vjust = 0.6, colour = "black"),
        axis.ticks.x = element_blank()) +
  scale_colour_manual(values = fills, drop = FALSE) +
  scale_x_discrete(limits = c("Active MR", "Aerobic capacity", "Boldness", "Dominance",
                             "Growth", "Movement/activity", "", "Reproduction", "Survival"))
p
pdf("Fig3_Fitnessthermo_R1.pdf", width = 7, height = 6)
p
dev.off()


#### Figure 4: Effect size vs sample size plot ####

pdf("Fig4_Samplesize_effectsize_R1.pdf", width = 7, height = 6)
par(mfrow = c(1,1), oma = c(1,1,1,1), mar = c(4,4,2,2), mgp = c(3,0.5,0))
plot(meta.dat.subMass$yi[meta.dat.subMass$sig=="0"] ~ meta.dat.subMass$n.rep[meta.dat.subMass$sig=="0"],
     data = meta.dat.subMass, ylim = c(-1.5, 3),
     pch = 16, col = rgb(col2rgb(pal.blue[5])[1,], col2rgb(pal.blue[5])[2,],
                         col2rgb(pal.blue[5])[3,], 0.5*255, maxColorValue = 255, names = NULL),
     bty = "l", log = "x", xlim = c(5, 1600), xlab = "", ylab = "")
points(meta.dat.subMass$yi[meta.dat.subMass$sig=="1"] ~ meta.dat.subMass$n.rep[meta.dat.subMass$sig=="1"],
       data = meta.dat.subMass, pch = 16,
       col = rgb(col2rgb(pal.orange[5])[1,], col2rgb(pal.orange[5])[2,],
                 col2rgb(pal.orange[5])[3,],0.5*255,maxColorValue = 255, names = NULL))

mtext(expression(paste(plain("Effect size ("), italic("Zr"), plain(")"))), line = 1.75, side = 2)
mtext(expression(paste(plain("Sample size ("), italic("n"), plain(")"))), line = 1.75, side = 1)

abline(0,0, lty = 2)

## Null for significant
meta.remv.null.sig <- rma.mv(yi = yi, V = vi, mods = ~ 1,
                             data = subset(meta.dat.subMass, meta.dat.subMass$sig=="1"),
                             random = list(~1|Ref, ~1|obs))
## Null for non-significant
meta.remv.null.nonsig <- rma.mv(yi = yi, V = vi, mods = ~ 1,
                                data = subset(meta.dat.subMass, meta.dat.subMass$sig=="0"),
                                random = list(~1|Ref, ~1|obs))

abline(h = meta.remv.null.nonsig$b, col = rgb(col2rgb(pal.blue[5])[1,], col2rgb(pal.blue[5])[2,],
                                              col2rgb(pal.blue[5])[3,], 0.75*255,maxColorValue = 255,
                                              names = NULL), lty = 1)
abline(h = meta.remv.null.nonsig$ci.lb, col = rgb(col2rgb(pal.blue[5])[1,], col2rgb(pal.blue[5])[2,],
                                                  col2rgb(pal.blue[5])[3,], 0.75*255,maxColorValue = 255,
                                                  names = NULL), lty = 2)
abline(h = meta.remv.null.nonsig$ci.ub, col = rgb(col2rgb(pal.blue[5])[1,], col2rgb(pal.blue[5])[2,],
                                                  col2rgb(pal.blue[5])[3,], 0.75*255,maxColorValue = 255,
                                                  names = NULL), lty = 2)

abline(h = meta.remv.null.sig$b, col = rgb(col2rgb(pal.orange[5])[1,], col2rgb(pal.orange[5])[2,],
                                           col2rgb(pal.orange[5])[3,], 0.75*255,maxColorValue = 255,
                                           names = NULL), lty = 1)
abline(h = meta.remv.null.sig$ci.lb, col = rgb(col2rgb(pal.orange[5])[1,], col2rgb(pal.orange[5])[2,],
                                               col2rgb(pal.orange[5])[3,], 0.75*255,maxColorValue = 255,
                                               names = NULL), lty = 2)
abline(h = meta.remv.null.sig$ci.ub, col = rgb(col2rgb(pal.orange[5])[1,], col2rgb(pal.orange[5])[2,],
                                               col2rgb(pal.orange[5])[3,], 0.75*255,maxColorValue = 255,
                                               names = NULL), lty = 2)

legend(400, 2, legend = c(expression(paste(italic("p"), plain(" < 0.05"))),
                          expression(paste(italic("p"), plain(" > 0.05")))), pch = 21,
       bty = "n", col = c(rgb(col2rgb(pal.orange[5])[1,], col2rgb(pal.orange[5])[2,],
                              col2rgb(pal.orange[5])[3,], 0.5*255,maxColorValue = 255, names = NULL),
                          rgb(col2rgb(pal.blue[5])[1,], col2rgb(pal.blue[5])[2,],
                              col2rgb(pal.blue[5])[3,], 0.5*255,maxColorValue = 255, names = NULL)),
       pt.bg = c(rgb(col2rgb(pal.orange[5])[1,], col2rgb(pal.orange[5])[2,],
                     col2rgb(pal.orange[5])[3,], 0.5*255,maxColorValue = 255, names = NULL),
                 rgb(col2rgb(pal.blue[5])[1,], col2rgb(pal.blue[5])[2,],
                     col2rgb(pal.blue[5])[3,], 0.5*255,maxColorValue = 255, names = NULL)),
       title = "Reported as", pt.cex = c(1,1))
dev.off()


#### Forest plot set up ####

## Forest plot from metafor
summary(meta.remv.null)

## Create data frame of best linear unbiased predictors (BLUP) of the random effects
studyRE.null.Ref <- ranef(meta.remv.null)[[1]]

## Add variables for Study, number of observations per study, thermo classification
studyRE.null.Ref$study <- as.character(meta.dat.subMass$Reference[match(row.names(studyRE.null.Ref),
                                                                        as.character(meta.dat.subMass$Ref))])

zt <- table(meta.dat.subMass$Ref)
studyRE.null.Ref$nInStudy <- as.numeric(zt[match(row.names(studyRE.null.Ref), names(zt))])
studyRE.null.Ref

studyRE.null.Ref$thermo <- meta.dat.subMass$thermo[match(row.names(studyRE.null.Ref),
                                                         as.character(meta.dat.subMass$Ref))]
summary(studyRE.null.Ref)

studyRE.null.Ref$Fit <- meta.dat.subMass$FitnessClassification[match(row.names(studyRE.null.Ref),
                                                                     as.character(meta.dat.subMass$Ref))]

## Intercept to add to BLUPs
int <- as.numeric(meta.remv.null$b)

## Source adapted forest plot function
source("forestNew.R")


#### Figure 1: Forest plot Ectotherms ####

pdf("Figure Forest Study Ectotherms_R1.pdf", width = 7, height = 9.5)
forestNew(studyRE.null.Ref$intrcpt + int, ci.lb = studyRE.null.Ref$pi.lb + int,
         ci.ub = studyRE.null.Ref$pi.ub + int, cex = 1.2, cex.lab = 1.2, cex.axis = 1, xlim = c(-3.5, 1.5),
         slab = studyRE.null.Ref$study, subset = studyRE.null.Ref$thermo == "Ectotherm",
         efac = 0, refline = as.numeric(meta.remv.null$b),
         refband = c(meta.remv.null$ci.lb, meta.remv.null$ci.ub), ilab = studyRE.null.Ref$nInStudy,
         ilab.xpos = -1.3, annotate = FALSE, pch = 21, bg = "steelblue2", main = "Ectotherms",
         xlab = expression(paste(plain("Effect size ("), italic("Zr"), plain(")"))),
         psize = log10(studyRE.null.Ref$nInStudy) + 0.55)
text(-3.25, 56, "Study", cex = 1.2)
text(-1.3, 56, expression(italic(n)), cex = 1.2)
dev.off()


#### Figure 2: Forest plot Endotherms ####

pdf("Figure Forest Study Endotherms_R1.pdf", width = 7, height = 9.5)
forestNew(studyRE.null.Ref$intrcpt + int, ci.lb = studyRE.null.Ref$pi.lb + int,
         ci.ub = studyRE.null.Ref$pi.ub + int, cex = 1.2, cex.lab = 1.2, cex.axis = 1, xlim = c(-3.5, 1.5),
         slab = studyRE.null.Ref$study, subset = studyRE.null.Ref$thermo == "Endotherm",
         efac = 0, refline = as.numeric(meta.remv.null$b),
         refband = c(meta.remv.null$ci.lb, meta.remv.null$ci.ub), ilab = studyRE.null.Ref$nInStudy,
         ilab.xpos = -1.3, annotate = FALSE, pch = 21, bg = "orangered", main = "Endotherms",
         xlab = expression(paste(plain("Effect size ("), italic("Zr"), plain(")"))),
         psize = log10(studyRE.null.Ref$nInStudy) + 0.55)
text(-3.25, 61, "Study", cex = 1.2)
text(-1.3, 61, expression(italic(n)), cex = 1.2)
dev.off()


#### Figure S2: Funnel plot to check asymmetry in the sample size distribution ####

meta.remv.fitTh.ecto <- rma.mv(yi = yi, V = vi, mods = ~ FitnessClassification,
                               data = subset(meta.dat.subMass, meta.dat.subMass$thermo=="Ectotherm"),
                               random = list(~1|Ref, ~1|obs))

meta.remv.fitTh.endo <- rma.mv(yi = yi, V = vi, mods = ~ FitnessClassification,
                               data = subset(meta.dat.subMass, meta.dat.subMass$thermo=="Endotherm"),
                               random = list(~1|Ref, ~1|obs))

# Moderator variable of ref and obs because rma.mv does not work for this
res <- meta.remv.fitTh
ranktest(res)

meta2 <- rma(yi = yi, vi = vi, mods = ~ Ref + obs, data = meta.dat.subMass)
regtest(meta2)

par(mfrow = c(2,2))
funnel(res, main = "Standard Error")
funnel(res, yaxis = "vi", main = "Sampling Variance")
funnel(res, yaxis = "seinv", main = "Inverse Standard Error")
funnel(res, yaxis = "vinv", main = "Inverse Sampling Variance")

par(mfrow = c(1,2))
funnel(res, level = c(90, 95, 99), shade = c("white", "gray", "darkgray"), refline = 0)
funnel(meta.remv.fitTh, level = c(90, 95, 99), shade = c("white", "gray", "darkgray"), refline = 0,
       xlab = expression(paste(plain("Effect size ("), italic("Zr"), plain(")"))))


#### Trim and fill ####
# trimfill only works on rma.uni models, so used only null univariate model
# these results should not be considered strong evidence; a visualisation of potential publication bias only
meta.dat.subMass.trim <- subset(meta.dat.subMass, !abs(meta.dat.subMass$yi) > 1.5)
nrow(meta.dat.subMass); nrow(meta.dat.subMass.trim)
meta.dat.subMass.trim.ecto <- subset(meta.dat.subMass.trim, meta.dat.subMass.trim$thermo=="Ectotherm")
meta.dat.subMass.trim.endo <- subset(meta.dat.subMass.trim, meta.dat.subMass.trim$thermo=="Endotherm")
nrow(meta.dat.subMass.trim.ecto); nrow(meta.dat.subMass.trim.endo)

# Overall
meta.trim <- rma.uni(yi = yi, vi = vi, data = meta.dat.subMass.trim)
summary(meta.trim)
taf <- trimfill(meta.trim)
taf
funnel(taf, legend = TRUE)

# Ectotherms
meta.trim.ecto <- rma.uni(yi = yi, vi = vi, data = meta.dat.subMass.trim.ecto)
summary(meta.trim.ecto)
taf.ecto <- trimfill(meta.trim.ecto)
taf.ecto
funnel(taf.ecto, legend = TRUE)

# Endotherms
meta.trim.endo <- rma.uni(yi = yi, vi = vi, data = meta.dat.subMass.trim.endo)
summary(meta.trim.endo)
taf.endo <- trimfill(meta.trim.endo)
taf.endo
funnel(taf.endo, legend = TRUE)

# Endotherm and ectotherms separately
pdf("FigS2_Funnel_diagnosticplots_R1.pdf", width = 8, height = 6)
par(mfrow = c(1,2))
res <- meta.remv.fitTh.endo
funnel(res, level = c(90, 95, 99), shade = c("white", "gray", "darkgray"), refline = 0,
       xlab = expression(paste(plain("Effect size ("), italic("Zr"), plain(")"))),
       xlim = c(-2.6,2.6), ylim = c(0, 0.8), pch = 21, col = "black", bg = "orangered")
mtext("a", font = 2, side = 3, padj = -1, adj = 0)
res <- meta.remv.fitTh.ecto
funnel(res, level = c(90, 95, 99), shade = c("white", "gray", "darkgray"), refline = 0,
       xlab = expression(paste(plain("Effect size ("), italic("Zr"), plain(")"))),
       xlim = c(-2.6,2.6), ylim = c(0, 0.8), pch = 21, col = "black", bg = "steelblue2")
mtext("b", font = 2, side = 3, padj = -1, adj = 0)
dev.off()

# Trim and fill plots
funnel(taf.endo, level = c(90, 95, 99), shade = c("white", "gray", "darkgray"), refline = 0,
       xlab = expression(paste(plain("Effect size ("), italic("Zr"), plain(")"))),
       xlim = c(-2.6,2.6), ylim = c(0, 0.8), pch = 21, col = "black", bg = "orangered")
funnel(taf.ecto, level = c(90, 95, 99), shade = c("white", "gray", "darkgray"), refline = 0,
       xlab = expression(paste(plain("Effect size ("), italic("Zr"), plain(")"))),
       xlim = c(-2.6,2.6), ylim = c(0, 0.8), pch = 21, col = "black", bg = "steelblue2")


#### Recalculating metrics for the rest of the supplementary figures ####

head(meta.dat.subMass)

# Frequency Calculations
frequency_counts <- meta.dat.subMass %>% group_by(Journal, JournalIF, sig) %>%
  summarize(abs_zr = mean(abs(yi)), n.in.journal = n()) %>% arrange(sig, JournalIF)
frequency_counts <- as.data.frame(frequency_counts)
sifdata <- frequency_counts
head(sifdata)

# Correlation
cor.test(sifdata$abs_zr, sifdata$JournalIF)


#### Figure S3: Zr ~ IF split by significance ####

pdf("FigS3_impactfactor_R1.pdf", width = 8, height = 7)
par(mfrow = c(1,1), oma = c(1,1,1,3), mar = c(4,4,2,6), mgp = c(3,0.5,0))
symbols(sifdata$JournalIF[sifdata$sig=="1"], sifdata$abs_zr[sifdata$sig=="1"],
        circles = sqrt(sifdata$n.in.journal[sifdata$sig=="1"]/pi), xlim = c(0.9,10), ylim = c(0,2.2),
        xlab = "", ylab = "", inches = 0.25, bty = "l",
        bg = rgb(col2rgb(pal.orange[5])[1,], col2rgb(pal.orange[5])[2,],
                 col2rgb(pal.orange[5])[3,], 0.4*255,maxColorValue = 255, names = NULL),
        fg = rgb(col2rgb(pal.orange[5])[1,], col2rgb(pal.orange[5])[2,],
                 col2rgb(pal.orange[5])[3,], 0.4*255,maxColorValue = 255, names = NULL))

symbols(sifdata$JournalIF[sifdata$sig=="0"], sifdata$abs_zr[sifdata$sig=="0"],
        circles = sqrt(sifdata$n.in.journal[sifdata$sig=="0"]/pi), inches = 0.4,
        bg = rgb(col2rgb(pal.blue[5])[1,], col2rgb(pal.blue[5])[2,],
                 col2rgb(pal.blue[5])[3,], 0.4*255, maxColorValue = 255, names = NULL),
        fg = rgb(col2rgb(pal.blue[5])[1,], col2rgb(pal.blue[5])[2,],
                 col2rgb(pal.blue[5])[3,], 0.4*255, maxColorValue = 255, names = NULL), add = TRUE)

mtext("Journal impact factor", line = 2, side = 1)
mtext(expression(paste(plain("Effect size (| "), italic("Zr"), plain(" |)"))), line = 2, side = 2)

csize <- sort(2.5*(sqrt(sifdata$n.in.journal/pi)), decreasing = FALSE)

legend(5, 2, legend = c("1", "15", "50"), title = expression(paste(italic("n"), plain(" published relationships"))),
       pch = 21, bty = "n", col = "black", pt.cex = c(csize[1], csize[61], csize[65]), xpd = TRUE,
       y.intersp = c(1,1.5,2), x.intersp = 3)

legend(6, 1, legend = c(expression(paste(italic("p"), plain(" < 0.05"))),
                          expression(paste(italic("p"), plain(" > 0.05")))), pch = 21,
       bty = "n", col = c(rgb(col2rgb(pal.orange[5])[1,], col2rgb(pal.orange[5])[2,],
                              col2rgb(pal.orange[5])[3,], 0.4*255,maxColorValue = 255, names = NULL),
                          rgb(col2rgb(pal.blue[5])[1,], col2rgb(pal.blue[5])[2,],
                              col2rgb(pal.blue[5])[3,], 0.4*255,maxColorValue = 255, names = NULL)),
       pt.bg = c(rgb(col2rgb(pal.orange[5])[1,], col2rgb(pal.orange[5])[2,],
                     col2rgb(pal.orange[5])[3,], 0.4*255,maxColorValue = 255, names = NULL),
                 rgb(col2rgb(pal.blue[5])[1,], col2rgb(pal.blue[5])[2,],
                     col2rgb(pal.blue[5])[3,], 0.4*255,maxColorValue = 255, names = NULL)),
       title = expression(paste(plain("Reported as"))), pt.cex = c(1,1), xpd = TRUE, y.intersp = 0.75)
dev.off()


#### Figure S4: Zr ~ Year split by significance ####

# Frequency Calculations
yearzrdata <- cbind(meta.dat.subMass$sig, abs(meta.dat.subMass$yi), meta.dat.subMass$Year)
yearzrdata <- as.data.frame(yearzrdata)
names(yearzrdata) <- c("sig", "abs_zr", "Year")
head(yearzrdata)

# Correlation
cor.test(yearzrdata$abs_zr, yearzrdata$Year)

pdf("FigS4_year_R1.pdf", width = 8, height = 7)
par(mfrow = c(1,1), oma = c(1,1,1,3), mar = c(4,4,2,1), mgp = c(3,0.5,0))
plot(abs_zr[sig=="1"] ~ Year[sig=="1"], data = yearzrdata, xlim = c(1980, 2018), ylim = c(0, 3),
     xlab = "", ylab = "", bty = "l", pch = 21,
     bg = rgb(col2rgb(pal.orange[5])[1,], col2rgb(pal.orange[5])[2,],
              col2rgb(pal.orange[5])[3,], 0.4*255,maxColorValue = 255, names = NULL),
     col = rgb(col2rgb(pal.orange[5])[1,], col2rgb(pal.orange[5])[2,],
               col2rgb(pal.orange[5])[3,], 0.4*255,maxColorValue = 255, names = NULL))
points(abs_zr[sig=="0"] ~ Year[sig=="0"], data = yearzrdata, xlim = c(1980, 2018), ylim = c(0, 3),
       xlab = "", ylab = "", bty = "l", pch = 21,
       bg = rgb(col2rgb(pal.blue[5])[1,], col2rgb(pal.blue[5])[2,],
                col2rgb(pal.blue[5])[3,], 0.4*255, maxColorValue = 255, names = NULL),
       col = rgb(col2rgb(pal.blue[5])[1,], col2rgb(pal.blue[5])[2,],
                 col2rgb(pal.blue[5])[3,], 0.4*255, maxColorValue = 255, names = NULL))
abline(lm(yearzrdata$abs_zr[yearzrdata$sig=="1"] ~ yearzrdata$Year[yearzrdata$sig=="1"]),
       col = rgb(col2rgb(pal.orange[5])[1,], col2rgb(pal.orange[5])[2,],
                 col2rgb(pal.orange[5])[3,], 0.5*255,maxColorValue = 255, names = NULL))
abline(lm(yearzrdata$abs_zr[yearzrdata$sig=="0"] ~ yearzrdata$Year[yearzrdata$sig=="0"]),
       col = rgb(col2rgb(pal.blue[5])[1,], col2rgb(pal.blue[5])[2,],
                 col2rgb(pal.blue[5])[3,], 0.5*255,maxColorValue = 255, names = NULL))
legend(2010, 3, legend = c(expression(paste(italic("p"), plain(" < 0.05"))),
                           expression(paste(italic("p"), plain(" > 0.05")))), pch = 21,
       bty = "n", col = c(rgb(col2rgb(pal.orange[5])[1,], col2rgb(pal.orange[5])[2,],
                              col2rgb(pal.orange[5])[3,], 0.4*255,maxColorValue = 255, names = NULL),
                          rgb(col2rgb(pal.blue[5])[1,], col2rgb(pal.blue[5])[2,],
                              col2rgb(pal.blue[5])[3,], 0.4*255,maxColorValue = 255, names = NULL)),
       pt.bg = c(rgb(col2rgb(pal.orange[5])[1,], col2rgb(pal.orange[5])[2,],
                     col2rgb(pal.orange[5])[3,], 0.4*255,maxColorValue = 255, names = NULL),
                 rgb(col2rgb(pal.blue[5])[1,], col2rgb(pal.blue[5])[2,],
                     col2rgb(pal.blue[5])[3,], 0.4*255,maxColorValue = 255, names = NULL)),
       title = expression(paste(plain("Reported as"))), pt.cex = c(1,1), xpd = TRUE, y.intersp = 1)

mtext("Publication year", line = 2, side = 1)
mtext(expression(paste(plain("Effect size (| "), italic("Zr"), plain(" |)"))), line = 2, side = 2)
dev.off()

####
