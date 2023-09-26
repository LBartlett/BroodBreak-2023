
#### BOA Study Analysis ####

# Analysis by Lewis Bartlett (lewis.bartlett@uga.edu)
# Study lead: Jennifer Berry (JBee@uga.edu)

# Read in libraries #

library(afex)
library(emmeans)

# Read in data

FullData <- read.csv('FullData.csv',
                     header = TRUE)

head(FullData)

FullData$PCM1 <- FullData$AW1Mites/FullData$AW1Bees
FullData$PCM2 <- FullData$AW2Mites/FullData$AW2Bees
FullData$PCM3 <- FullData$AW3Mites/FullData$AW3Bees
FullData$PCM4 <- FullData$AW4Mites/FullData$AW4Bees

FullData$DMT1 <- FullData$PCM2 - FullData$PCM1
FullData$DMT2 <- FullData$PCM3 - FullData$PCM2
FullData$DMT3 <- FullData$PCM4 - FullData$PCM3
FullData$DMFull <- FullData$PCM4 - FullData$PCM1

FullDeltaMod1.1 <- mixed(DMFull ~ BroodBreak + OADose + BroodBreak:OADose + (1|Apiary),
                         data = FullData)

nice(FullDeltaMod1.1)

emmeans(FullDeltaMod1.1, specs = 'BroodBreak', by = 'OADose')

emmeans(FullDeltaMod1.1, specs = 'OADose', by = 'BroodBreak')

FullData$OADF <- as.factor(FullData$OADose)

FullDeltaMod1.2 <- mixed(DMFull ~ BroodBreak + OADF + BroodBreak:OADF + (1|Apiary),
                         data = FullData)

nice(FullDeltaMod1.2)

emmeans(FullDeltaMod1.2, specs = 'BroodBreak', by = 'OADF')

emmeans(FullDeltaMod1.2, specs = 'OADF', by = 'BroodBreak')


FullDeltaMod2.1 <- mixed(DMFull ~ BroodBreak + OADose + (1|Apiary),
                         data = FullData)

nice(FullDeltaMod2.1)

summary(FullDeltaMod2.1)

FullDeltaMod2.2 <- mixed(DMFull ~ BroodBreak + OADF + (1|Apiary),
                         data = FullData)

nice(FullDeltaMod2.2)

summary(FullDeltaMod2.2)

emmeans(FullDeltaMod2.2, specs = 'BroodBreak', by = 'OADF')

emmeans(FullDeltaMod2.2, specs = 'OADF', by = 'BroodBreak')

FullData$OATreat <- FullData$OADose != 0

FullDeltaMod3.1 <- mixed(DMFull ~ BroodBreak + OATreat + BroodBreak:OATreat + (1|Apiary),
                         data = FullData)

nice(FullDeltaMod3.1)

emmeans(FullDeltaMod3.1, specs = 'BroodBreak', by = 'OATreat')

# Timeseries analysis

AWDataTime <- data.frame(Apiary = rep(FullData$Apiary, times = 4),
                         Colony = rep(FullData$Colony, times = 4),
                         BroodBreak = rep(FullData$BroodBreak, times = 4),
                         OADose = rep(FullData$OADose, times = 4),
                         OATreat = rep(FullData$OATreat, times = 4),
                         TimePoint = sort(rep(c(0,1,2,3), times = NROW(FullData))),
                         PCM.R = c(FullData$PCM1, FullData$PCM2, FullData$PCM3, FullData$PCM4)
)

AWDataTime$PCM <- ceiling(100*AWDataTime$PCM.R)

FullTempMod1.1 <- mixed(PCM ~ TimePoint + TimePoint:BroodBreak + TimePoint:OADose + TimePoint:OADose:BroodBreak  + (1|Apiary) + (1|Apiary:Colony),
                        family = 'poisson', method = 'LRT',
                        data = AWDataTime)
nice(FullTempMod1.1)


FullTempMod2.1 <- mixed(PCM ~ TimePoint + TimePoint:BroodBreak + TimePoint:as.factor(OADose) + TimePoint:as.factor(OADose):BroodBreak  + (1|Apiary) + (1|Apiary:Colony),
                        family = 'poisson', method = 'LRT',
                        data = AWDataTime)
nice(FullTempMod2.1)

emmeans(FullTempMod2.1, specs = c('BroodBreak','OADose'), by = 'TimePoint')


FullTempMod3.1 <- mixed(PCM ~ TimePoint + TimePoint:BroodBreak + TimePoint:OATreat + TimePoint:OATreat:BroodBreak  + (1|Apiary) + (1|Apiary:Colony),
                        family = 'poisson', method = 'LRT',
                        data = na.exclude(AWDataTime))
nice(FullTempMod3.1)

emmeans(FullTempMod3.1, specs = c('BroodBreak','OATreat'), by = 'TimePoint')

emtrends(FullTempMod3.1, specs = c('BroodBreak','OATreat'), var = 'TimePoint')

##############
# Of the above: end-point delta PCM models useful for 'real value' emmeans quoting
# temporal models for actual quoting significance
# both for saying 2 and 3 undifferentiable
# arguably both for showing 'brood break only or oa only very similar'

# dont at this point have anything about where we get to with however many treatments

##############

# Sticky screen data

FullData$DSS1 <- (FullData$SS2 - FullData$SS1)/FullData$SS1

FullDSSMod1.1 <- mixed(DSS1 ~ BroodBreak + OADose + BroodBreak:OADose + (1|Apiary),
                         data = FullData)

nice(FullDSSMod1.1)

emmeans(FullDSSMod1.1, specs = 'BroodBreak', by = 'OADose')

emmeans(FullDSSMod1.1, specs = 'OADose', by = 'BroodBreak')

FullDSSMod1.2 <- mixed(DSS1 ~ BroodBreak + OADF + BroodBreak:OADF + (1|Apiary),
                         data = FullData)

nice(FullDSSMod1.2)

emmeans(FullDSSMod1.2, specs = 'BroodBreak', by = 'OADF')

emmeans(FullDSSMod1.2, specs = 'OADF', by = 'BroodBreak')




FullDSSMod3.1 <- mixed(DSS1 ~ BroodBreak + OATreat + BroodBreak:OATreat + (1|Apiary),
                         data = FullData)

nice(FullDSSMod3.1)

emmeans(FullDSSMod3.1, specs = 'BroodBreak', by = 'OATreat')

emmeans(FullDSSMod3.1, specs = c('BroodBreak','OATreat'))

# now THERE'S an interaction effect if ever I saw one
# again, see earlier models for little difference in 2 vs 3
# final one for the 'big result'

################

# Think about what plots you now want
# time series [obviously] of AW
# start vs end PCM
# sticky screen result

FullData$TCode <- ((FullData$BroodBreak)*2)+FullData$OATreat
FullData$TCode[which(FullData$TCode == 0)] <- 'B0V0'
FullData$TCode[which(FullData$TCode == 1)] <- 'B0V1'
FullData$TCode[which(FullData$TCode == 2)] <- 'B1V0'
FullData$TCode[which(FullData$TCode == 3)] <- 'B1V1'

FullData$TCode <- factor(FullData$TCode, levels = c('B0V0','B1V0','B0V1','B1V1'))


# Quick transparency function for plotting
Transpa <- function(color, percent) {
  
  rgb.val <- col2rgb(color)
  
  t.col <- rgb(rgb.val[1], rgb.val[2], rgb.val[3],
               max = 255,
               alpha = (100-percent)*255/100)
  
  return(t.col)
  
}

# Set plotting colours
ColRef <- data.frame(Treatment = c('B0V0','B1V0','B0V1','B1V1'), Col =  c('pink3','blue2','red2','purple4'))

# Makes the plots we would like
par(mar=c(5,8,2,2))

# Change in mite drop (death metric)
boxplot(FullData$DSS1 ~ FullData$TCode, 
        main = NA, ylab = 'Change in Mite Drop', xlab = 'Treatment', 
        border = 'transparent', 
        cex.axis = 1.3, cex.lab = 1.5, outline = TRUE, lwd = 1.2,
        boxlty = 1, whisklty = 0, staplelty = 1, boxwex = 0.00, 
        col = sapply(X = as.character(ColRef$Col), FUN = Transpa, percent = 100))

abline(h = 0 , lty = 3)

stripchart(FullData$DSS1 ~ FullData$TCode,
           col = sapply(X = as.character(ColRef$Col[match(ColRef$Treatment, levels(FullData$TCode))]), FUN = Transpa, percent = 40),
           vertical = T, add = T, pch = 4, cex = 0.95, 
           method = 'jitter', lwd = 2)

# (this is not real analysis it is for plotting purposes only - STRIP THE WILD CONTROL OUTLIER OUT for this)
PlotMod <- glm(FullData$DSS1[-71] ~ FullData$TCode[-71],
               family = 'gaussian')

PlotCIs <- as.data.frame(emmeans(PlotMod, specs =c('TCode')))

for(L in 1:NROW(PlotCIs)){
  
  segments(x0 = L, x1 = L, y0 = PlotCIs$asymp.LCL[L], y1 = PlotCIs$asymp.UCL[L],
           col = as.character(ColRef$Col[which(as.character(ColRef$Treatment) == as.character(PlotCIs$TCode[L]))]),
           lwd = 3)
  
}

# Change in mite wash (start to end)

boxplot(FullData$DMFull ~ FullData$TCode, 
        main = NA, ylab = 'Change in PMI', xlab = 'Treatment', 
        border = 'transparent', 
        cex.axis = 1.3, cex.lab = 1.5, outline = TRUE, lwd = 1.2,
        boxlty = 1, whisklty = 0, staplelty = 1, boxwex = 0.00, 
        col = sapply(X = as.character(ColRef$Col), FUN = Transpa, percent = 100))

abline(h = 0 , lty = 3)

stripchart(FullData$DMFull ~ FullData$TCode,
           col = sapply(X = as.character(ColRef$Col[match(ColRef$Treatment, levels(FullData$TCode))]), FUN = Transpa, percent = 40),
           vertical = T, add = T, pch = 4, cex = 0.95, 
           method = 'jitter', lwd = 2)

PlotMod <- glm(FullData$DMFull ~ FullData$TCode,
               family = 'gaussian')

PlotCIs <- as.data.frame(emmeans(PlotMod, specs =c('TCode')))

for(L in 1:NROW(PlotCIs)){
  
  segments(x0 = L, x1 = L, y0 = PlotCIs$asymp.LCL[L], y1 = PlotCIs$asymp.UCL[L],
           col = as.character(ColRef$Col[which(as.character(ColRef$Treatment) == as.character(PlotCIs$TCode[L]))]),
           lwd = 3)
  
}

### Time series of alcohol washes

# Two parts - messy 'join the dots' vs smoothed model effects?

AWDataTime$UCol <- paste0(AWDataTime$Apiary, AWDataTime$Colony)

AWDataTime$TCode <- ((AWDataTime$BroodBreak)*2)+AWDataTime$OATreat
AWDataTime$TCode[which(AWDataTime$TCode == 0)] <- 'B0V0'
AWDataTime$TCode[which(AWDataTime$TCode == 1)] <- 'B0V1'
AWDataTime$TCode[which(AWDataTime$TCode == 2)] <- 'B1V0'
AWDataTime$TCode[which(AWDataTime$TCode == 3)] <- 'B1V1'

AWDataTime$TCode <- factor(AWDataTime$TCode, levels = c('B0V0','B1V0','B0V1','B1V1'))

par(mar=c(5,8,2,2))

# Make blank plot
plot(AWDataTime$PCM.R ~ AWDataTime$TimePoint, 
     type = 'n', 
     xaxp = c(0,3,3), 
     ylab = 'Percent Mite Intensity', 
     xlab ='Timepoint',
     cex.lab = 1.3)

# Add points showing each colony's mite intensity across time, colour coded by treatment

for(C in unique(AWDataTime$UCol)){
  
  CCol <- as.character(ColRef$Col[which(ColRef$Treatment == unique(AWDataTime$TCode[which(AWDataTime$UCol == C)]))])
  
  points(x = jitter(AWDataTime$TimePoint[which(AWDataTime$UCol == C)], factor = 0.25),
         y = AWDataTime$PCM.R[which(AWDataTime$UCol == C)],
         col = Transpa(CCol, 40), type = 'b', lwd = 0.5, pch = 4, lty = 3, cex = 2)
  
}


# Add big fucking averages

for(M in unique(AWDataTime$TCode)){
  
  CCol <- as.character(ColRef$Col[which(ColRef$Treatment == M)])
  
  XP <- c(0,1,2,3)
  YP <- c(mean(AWDataTime$PCM.R[which(AWDataTime$TimePoint == 0 & AWDataTime$TCode == M)], na.rm = T),
          mean(AWDataTime$PCM.R[which(AWDataTime$TimePoint == 1 & AWDataTime$TCode == M)], na.rm = T),
          mean(AWDataTime$PCM.R[which(AWDataTime$TimePoint == 2 & AWDataTime$TCode == M)], na.rm = T),
          mean(AWDataTime$PCM.R[which(AWDataTime$TimePoint == 3 & AWDataTime$TCode == M)], na.rm = T)          )
  
  points(x = XP,
         y = YP,
         col = Transpa(CCol, 15), type = 'b', lwd = 4, pch = 45, lty = 1, cex = 8)
  
}





