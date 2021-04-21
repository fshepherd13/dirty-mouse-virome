# load CRAN packages
library(tidyverse)
library(seqinr)
library(ggrepel)
library(readxl)

# load github package for codon pair bias analysis
install.packages("remotes")
remotes::install_github("alex-sbu/CPBias")
library(CPBias)

# load Biocmanager packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("Biostrings")
library(Biostrings)


# load fasta files for narnavirus genomes
reepicheep <- read.fasta("Reepicheep.fasta")
matryoshka_1 <- read.fasta("Matryoshka1.fasta")


# codon usage
  
## seqinr::uco -- caluclating Codon Usage Bias metrics for narnaviruses
Reepicheep_uco <- uco(reepicheep$Reepicheep, frame = 0, as.data.frame = TRUE) %>%
  arrange(codon)
Matryoshka_1_uco <- uco(matryoshka_1$MN698829.1, as.data.frame = TRUE) %>% 
  arrange(codon)

## read in host codon usage tables from HIVE Biochemistry
### mouse
Mouse_Codons <- read.csv("mouseCodons.csv") %>% arrange(Codon)
Mouse_Codons$Codon <- tolower(Mouse_Codons$Codon)
Mouse_Codons$AA <- Reepicheep_uco$AA
Mouse_Codons$AA <- as.factor(Mouse_Codons$AA)

### cryptosporidium
Crypto_Codons <- read.csv("cryptoCodons.csv") %>% arrange(Codon)
Crypto_Codons$Codon <- tolower(Crypto_Codons$Codon)
Crypto_Codons$AA <- Reepicheep_uco$AA
Crypto_Codons$AA <- as.factor(Crypto_Codons$AA)

### plasmodium
Plasmo_Codons <- read.csv("plasmoCodons.csv") %>% arrange(Codon)
Plasmo_Codons$Codon <- tolower(Plasmo_Codons$Codon)
Plasmo_Codons$AA <- Reepicheep_uco$AA
Plasmo_Codons$AA <- as.factor(Plasmo_Codons$AA)

### human
Human_Codons <- read.csv("humanCodons.csv") %>%
  arrange(Codon)
Human_Codons$Codon <- tolower(Human_Codons$Codon)
Human_Codons$AA <- Reepicheep_uco$AA
Human_Codons$AA <- as.factor(Human_Codons$AA)

## calculate the rscu for the hosts

### mouse
## create a table with the # of codons that code for an AA and the occurrences of all codons that code for the same AA
Synonymous_Codons_Mouse <- Mouse_Codons %>% 
  group_by(AA) %>% 
  summarize(AA_count = n(), Codon_Count = sum(Count))

## rscu calculation
Mouse_Codons$RSCU <- ((Mouse_Codons$Count * Synonymous_Codons_Mouse$AA_count[match(Mouse_Codons$AA, Synonymous_Codons_Mouse$AA)])/Synonymous_Codons_Mouse$Codon_Count[match(Mouse_Codons$AA, Synonymous_Codons_Mouse$AA)])


### Cryptosporidium
Synonymous_Codons_Crypto <- Crypto_Codons %>% 
  group_by(AA) %>% 
  summarize(AA_count = n(), Codon_Count = sum(Count))

Crypto_Codons$RSCU <- ((Crypto_Codons$Count * Synonymous_Codons_Crypto$AA_count[match(Crypto_Codons$AA, Synonymous_Codons_Crypto$AA)])/Synonymous_Codons_Crypto$Codon_Count[match(Crypto_Codons$AA, Synonymous_Codons_Crypto$AA)])

### Plasmodium
Synonymous_Codons_Plasmo <- Plasmo_Codons %>% 
  group_by(AA) %>% 
  summarize(AA_count = n(), Codon_Count = sum(Count))

Plasmo_Codons$RSCU <- ((Plasmo_Codons$Count * Synonymous_Codons_Plasmo$AA_count[match(Plasmo_Codons$AA, Synonymous_Codons_Plasmo$AA)])/Synonymous_Codons_Plasmo$Codon_Count[match(Plasmo_Codons$AA, Synonymous_Codons_Plasmo$AA)])

### human
Synonymous_Codons_Human <- Human_Codons %>% 
  group_by(AA) %>% 
  summarize(AA_count = n(), Codon_Count = sum(Count))

Human_Codons$RSCU <- ((Human_Codons$Count * Synonymous_Codons_Human$AA_count[match(Human_Codons$AA, Synonymous_Codons_Human$AA)])/Synonymous_Codons_Human$Codon_Count[match(Human_Codons$AA, Synonymous_Codons_Human$AA)])

## create a dataframe that contains the rscu values for the narnaviruses and hosts
MDS_Data_Frame <- data.frame(Reepicheep_uco$codon, Reepicheep_uco$RSCU, Mouse_Codons$RSCU, Crypto_Codons$RSCU, Matryoshka_1_uco$RSCU, Plasmo_Codons$RSCU, Human_Codons$RSCU) %>% t()

colnames(MDS_Data_Frame) <- MDS_Data_Frame[1,]

MDS_Data_Frame <- MDS_Data_Frame[-1,]

## MDS
d <- dist(MDS_Data_Frame) 
fit <- cmdscale(d,eig=TRUE, k=2) 
fit 

Fit <- data.frame(fit$points)

## names for the MDS plot
row.names(Fit) <- c("Reepicheep", "Mouse", "Cryptosporidium", "Matryoshka 1", "Plasmodium", "Human")

## plot
ggplot(Fit, aes(x = X1, y = X2)) +
  geom_point() +
  ggtitle("Relative Synonymous Codon Usage Metric MDS") +
  xlab("Coordinate 1") +
  ylab("Coordinate 2") +
  geom_label_repel(aes(label = row.names(Fit)),
                   box.padding   = 0.35, 
                   point.padding = 0.5,
                   segment.color = 'grey50') +
  theme_bw() + 
  theme(panel.background = element_blank(), 
        panel.grid.major = element_blank(),  #remove major-grid labels
        panel.grid.minor = element_blank(),  #remove minor-grid labels
        plot.background = element_blank(), 
        plot.title = element_text(hjust = 0.5)) #center the title


# dinucleotide bias
  
## seqinr::rho to calculate narnavirus dinucleotide bias
  
### reepicheep rdrp
Dinuc_Reepicheep <- rho(reepicheep[[1]])

#matryoshka1
Dinuc_Matryoshka1 <- rho(matryoshka_1[[1]])

## hosts

## read in host dinucleotide tables from HIVE Biochemistry
Dinuc_Mouse_Freq <- read.csv("mouseDinuc.csv")
Dinuc_Human_Freq <- read.csv("humanDinuc.csv")
Dinuc_Crypto_Freq <- read.csv("cryptoDinuc.csv")
Dinuc_Plasmo_Freq <- read.csv("plasmoDinuc.csv")

## rho calculation

### mouse
## calculating the freqency of each nucleotide
Mouse_FreqA <- (Dinuc_Mouse_Freq$AA *2 + Dinuc_Mouse_Freq$AC + Dinuc_Mouse_Freq$AG + Dinuc_Mouse_Freq$AT + Dinuc_Mouse_Freq$CA + Dinuc_Mouse_Freq$GA + Dinuc_Mouse_Freq$TA)/(Dinuc_Mouse_Freq$X.DINUCLEOTIDES * 2)

Mouse_FreqC <- (Dinuc_Mouse_Freq$CC *2 + Dinuc_Mouse_Freq$CA + Dinuc_Mouse_Freq$CG + Dinuc_Mouse_Freq$CT + Dinuc_Mouse_Freq$AC + Dinuc_Mouse_Freq$GC + Dinuc_Mouse_Freq$TC)/(Dinuc_Mouse_Freq$X.DINUCLEOTIDES * 2)

Mouse_FreqG <- (Dinuc_Mouse_Freq$GG *2 + Dinuc_Mouse_Freq$GA + Dinuc_Mouse_Freq$GC + Dinuc_Mouse_Freq$GT + Dinuc_Mouse_Freq$CG + Dinuc_Mouse_Freq$AG + Dinuc_Mouse_Freq$TG)/(Dinuc_Mouse_Freq$X.DINUCLEOTIDES * 2)

Mouse_FreqT <- (Dinuc_Mouse_Freq$TT *2 + Dinuc_Mouse_Freq$TC + Dinuc_Mouse_Freq$TG + Dinuc_Mouse_Freq$TA + Dinuc_Mouse_Freq$CT + Dinuc_Mouse_Freq$GT + Dinuc_Mouse_Freq$AT)/(Dinuc_Mouse_Freq$X.DINUCLEOTIDES * 2)

## calculating the rho scores
Dinuc_Mouse <- c()

Dinuc_Mouse[1] <- (Dinuc_Mouse_Freq$AA/Dinuc_Mouse_Freq$X.DINUCLEOTIDES)/(Mouse_FreqA*Mouse_FreqA)

Dinuc_Mouse[2] <- (Dinuc_Mouse_Freq$AC/Dinuc_Mouse_Freq$X.DINUCLEOTIDES)/(Mouse_FreqA*Mouse_FreqC)

Dinuc_Mouse[3] <- (Dinuc_Mouse_Freq$AG/Dinuc_Mouse_Freq$X.DINUCLEOTIDES)/(Mouse_FreqA*Mouse_FreqG)

Dinuc_Mouse[4] <- (Dinuc_Mouse_Freq$AT/Dinuc_Mouse_Freq$X.DINUCLEOTIDES)/(Mouse_FreqA*Mouse_FreqT)

Dinuc_Mouse[6] <- (Dinuc_Mouse_Freq$CC/Dinuc_Mouse_Freq$X.DINUCLEOTIDES)/(Mouse_FreqC*Mouse_FreqC)

Dinuc_Mouse[5] <- (Dinuc_Mouse_Freq$CA/Dinuc_Mouse_Freq$X.DINUCLEOTIDES)/(Mouse_FreqC*Mouse_FreqA)

Dinuc_Mouse[7] <- (Dinuc_Mouse_Freq$CG/Dinuc_Mouse_Freq$X.DINUCLEOTIDES)/(Mouse_FreqC*Mouse_FreqG)

Dinuc_Mouse[8] <- (Dinuc_Mouse_Freq$CT/Dinuc_Mouse_Freq$X.DINUCLEOTIDES)/(Mouse_FreqC*Mouse_FreqT)

Dinuc_Mouse[11] <- (Dinuc_Mouse_Freq$GG/Dinuc_Mouse_Freq$X.DINUCLEOTIDES)/(Mouse_FreqG*Mouse_FreqG)

Dinuc_Mouse[9] <- (Dinuc_Mouse_Freq$GA/Dinuc_Mouse_Freq$X.DINUCLEOTIDES)/(Mouse_FreqG*Mouse_FreqA)

Dinuc_Mouse[10] <- (Dinuc_Mouse_Freq$GC/Dinuc_Mouse_Freq$X.DINUCLEOTIDES)/(Mouse_FreqG*Mouse_FreqC)

Dinuc_Mouse[12] <- (Dinuc_Mouse_Freq$GT/Dinuc_Mouse_Freq$X.DINUCLEOTIDES)/(Mouse_FreqG*Mouse_FreqT)

Dinuc_Mouse[16] <- (Dinuc_Mouse_Freq$TT/Dinuc_Mouse_Freq$X.DINUCLEOTIDES)/(Mouse_FreqT*Mouse_FreqT)

Dinuc_Mouse[13] <- (Dinuc_Mouse_Freq$TA/Dinuc_Mouse_Freq$X.DINUCLEOTIDES)/(Mouse_FreqT*Mouse_FreqA)

Dinuc_Mouse[14] <- (Dinuc_Mouse_Freq$TC/Dinuc_Mouse_Freq$X.DINUCLEOTIDES)/(Mouse_FreqT*Mouse_FreqC)

Dinuc_Mouse[15] <- (Dinuc_Mouse_Freq$TG/Dinuc_Mouse_Freq$X.DINUCLEOTIDES)/(Mouse_FreqT*Mouse_FreqG)

### human
## nucleotide frequencies
Human_FreqA <- (Dinuc_Human_Freq$AA *2 + Dinuc_Human_Freq$AC + Dinuc_Human_Freq$AG + Dinuc_Human_Freq$AT + Dinuc_Human_Freq$CA + Dinuc_Human_Freq$GA + Dinuc_Human_Freq$TA)/(Dinuc_Human_Freq$X.DINUCLEOTIDES * 2)

Human_FreqC <- (Dinuc_Human_Freq$CC *2 + Dinuc_Human_Freq$CA + Dinuc_Human_Freq$CG + Dinuc_Human_Freq$CT + Dinuc_Human_Freq$AC + Dinuc_Human_Freq$GC + Dinuc_Human_Freq$TC)/(Dinuc_Human_Freq$X.DINUCLEOTIDES * 2)

Human_FreqG <- (Dinuc_Human_Freq$GG *2 + Dinuc_Human_Freq$GA + Dinuc_Human_Freq$GC + Dinuc_Human_Freq$GT + Dinuc_Human_Freq$CG + Dinuc_Human_Freq$AG + Dinuc_Human_Freq$TG)/(Dinuc_Human_Freq$X.DINUCLEOTIDES * 2)

Human_FreqT <- (Dinuc_Human_Freq$TT *2 + Dinuc_Human_Freq$TC + Dinuc_Human_Freq$TG + Dinuc_Human_Freq$TA + Dinuc_Human_Freq$CT + Dinuc_Human_Freq$GT + Dinuc_Human_Freq$AT)/(Dinuc_Human_Freq$X.DINUCLEOTIDES * 2)

## rho scores
Dinuc_Human <- c()

Dinuc_Human[1] <- (Dinuc_Human_Freq$AA/Dinuc_Human_Freq$X.DINUCLEOTIDES)/(Human_FreqA*Human_FreqA)

Dinuc_Human[2] <- (Dinuc_Human_Freq$AC/Dinuc_Human_Freq$X.DINUCLEOTIDES)/(Human_FreqA*Human_FreqC)

Dinuc_Human[3] <- (Dinuc_Human_Freq$AG/Dinuc_Human_Freq$X.DINUCLEOTIDES)/(Human_FreqA*Human_FreqG)

Dinuc_Human[4] <- (Dinuc_Human_Freq$AT/Dinuc_Human_Freq$X.DINUCLEOTIDES)/(Human_FreqA*Human_FreqT)

Dinuc_Human[6] <- (Dinuc_Human_Freq$CC/Dinuc_Human_Freq$X.DINUCLEOTIDES)/(Human_FreqC*Human_FreqC)

Dinuc_Human[5] <- (Dinuc_Human_Freq$CA/Dinuc_Human_Freq$X.DINUCLEOTIDES)/(Human_FreqC*Human_FreqA)

Dinuc_Human[7] <- (Dinuc_Human_Freq$CG/Dinuc_Human_Freq$X.DINUCLEOTIDES)/(Human_FreqC*Human_FreqG)

Dinuc_Human[8] <- (Dinuc_Human_Freq$CT/Dinuc_Human_Freq$X.DINUCLEOTIDES)/(Human_FreqC*Human_FreqT)

Dinuc_Human[11] <- (Dinuc_Human_Freq$GG/Dinuc_Human_Freq$X.DINUCLEOTIDES)/(Human_FreqG*Human_FreqG)

Dinuc_Human[9] <- (Dinuc_Human_Freq$GA/Dinuc_Human_Freq$X.DINUCLEOTIDES)/(Human_FreqG*Human_FreqA)

Dinuc_Human[10] <- (Dinuc_Human_Freq$GC/Dinuc_Human_Freq$X.DINUCLEOTIDES)/(Human_FreqG*Human_FreqC)

Dinuc_Human[12] <- (Dinuc_Human_Freq$GT/Dinuc_Human_Freq$X.DINUCLEOTIDES)/(Human_FreqG*Human_FreqT)

Dinuc_Human[16] <- (Dinuc_Human_Freq$TT/Dinuc_Human_Freq$X.DINUCLEOTIDES)/(Human_FreqT*Human_FreqT)

Dinuc_Human[13] <- (Dinuc_Human_Freq$TA/Dinuc_Human_Freq$X.DINUCLEOTIDES)/(Human_FreqT*Human_FreqA)

Dinuc_Human[14] <- (Dinuc_Human_Freq$TC/Dinuc_Human_Freq$X.DINUCLEOTIDES)/(Human_FreqT*Human_FreqC)

Dinuc_Human[15] <- (Dinuc_Human_Freq$TG/Dinuc_Human_Freq$X.DINUCLEOTIDES)/(Human_FreqT*Human_FreqG)

### crypto
## nucleotide frequencies
Crypto_FreqA <- (Dinuc_Crypto_Freq$AA *2 + Dinuc_Crypto_Freq$AC + Dinuc_Crypto_Freq$AG + Dinuc_Crypto_Freq$AT + Dinuc_Crypto_Freq$CA + Dinuc_Crypto_Freq$GA + Dinuc_Crypto_Freq$TA)/(Dinuc_Crypto_Freq$X.DINUCLEOTIDES * 2)

Crypto_FreqC <- (Dinuc_Crypto_Freq$CC *2 + Dinuc_Crypto_Freq$CA + Dinuc_Crypto_Freq$CG + Dinuc_Crypto_Freq$CT + Dinuc_Crypto_Freq$AC + Dinuc_Crypto_Freq$GC + Dinuc_Crypto_Freq$TC)/(Dinuc_Crypto_Freq$X.DINUCLEOTIDES * 2)

Crypto_FreqG <- (Dinuc_Crypto_Freq$GG *2 + Dinuc_Crypto_Freq$GA + Dinuc_Crypto_Freq$GC + Dinuc_Crypto_Freq$GT + Dinuc_Crypto_Freq$CG + Dinuc_Crypto_Freq$AG + Dinuc_Crypto_Freq$TG)/(Dinuc_Crypto_Freq$X.DINUCLEOTIDES * 2)

Crypto_FreqT <- (Dinuc_Crypto_Freq$TT *2 + Dinuc_Crypto_Freq$TC + Dinuc_Crypto_Freq$TG + Dinuc_Crypto_Freq$TA + Dinuc_Crypto_Freq$CT + Dinuc_Crypto_Freq$GT + Dinuc_Crypto_Freq$AT)/(Dinuc_Crypto_Freq$X.DINUCLEOTIDES * 2)

## rho scores
Dinuc_Crypto <- c()

Dinuc_Crypto[1] <- (Dinuc_Crypto_Freq$AA/Dinuc_Crypto_Freq$X.DINUCLEOTIDES)/(Crypto_FreqA*Crypto_FreqA)

Dinuc_Crypto[2] <- (Dinuc_Crypto_Freq$AC/Dinuc_Crypto_Freq$X.DINUCLEOTIDES)/(Crypto_FreqA*Crypto_FreqC)

Dinuc_Crypto[3] <- (Dinuc_Crypto_Freq$AG/Dinuc_Crypto_Freq$X.DINUCLEOTIDES)/(Crypto_FreqA*Crypto_FreqG)

Dinuc_Crypto[4] <- (Dinuc_Crypto_Freq$AT/Dinuc_Crypto_Freq$X.DINUCLEOTIDES)/(Crypto_FreqA*Crypto_FreqT)

Dinuc_Crypto[6] <- (Dinuc_Crypto_Freq$CC/Dinuc_Crypto_Freq$X.DINUCLEOTIDES)/(Crypto_FreqC*Crypto_FreqC)

Dinuc_Crypto[5] <- (Dinuc_Crypto_Freq$CA/Dinuc_Crypto_Freq$X.DINUCLEOTIDES)/(Crypto_FreqC*Crypto_FreqA)

Dinuc_Crypto[7] <- (Dinuc_Crypto_Freq$CG/Dinuc_Crypto_Freq$X.DINUCLEOTIDES)/(Crypto_FreqC*Crypto_FreqG)

Dinuc_Crypto[8] <- (Dinuc_Crypto_Freq$CT/Dinuc_Crypto_Freq$X.DINUCLEOTIDES)/(Crypto_FreqC*Crypto_FreqT)

Dinuc_Crypto[11] <- (Dinuc_Crypto_Freq$GG/Dinuc_Crypto_Freq$X.DINUCLEOTIDES)/(Crypto_FreqG*Crypto_FreqG)

Dinuc_Crypto[9] <- (Dinuc_Crypto_Freq$GA/Dinuc_Crypto_Freq$X.DINUCLEOTIDES)/(Crypto_FreqG*Crypto_FreqA)

Dinuc_Crypto[10] <- (Dinuc_Crypto_Freq$GC/Dinuc_Crypto_Freq$X.DINUCLEOTIDES)/(Crypto_FreqG*Crypto_FreqC)

Dinuc_Crypto[12] <- (Dinuc_Crypto_Freq$GT/Dinuc_Crypto_Freq$X.DINUCLEOTIDES)/(Crypto_FreqG*Crypto_FreqT)

Dinuc_Crypto[16] <- (Dinuc_Crypto_Freq$TT/Dinuc_Crypto_Freq$X.DINUCLEOTIDES)/(Crypto_FreqT*Crypto_FreqT)

Dinuc_Crypto[13] <- (Dinuc_Crypto_Freq$TA/Dinuc_Crypto_Freq$X.DINUCLEOTIDES)/(Crypto_FreqT*Crypto_FreqA)

Dinuc_Crypto[14] <- (Dinuc_Crypto_Freq$TC/Dinuc_Crypto_Freq$X.DINUCLEOTIDES)/(Crypto_FreqT*Crypto_FreqC)

Dinuc_Crypto[15] <- (Dinuc_Crypto_Freq$TG/Dinuc_Crypto_Freq$X.DINUCLEOTIDES)/(Crypto_FreqT*Crypto_FreqG)

### plasmo
## nucleotide frequencies
Plasmo_FreqA <- (Dinuc_Plasmo_Freq$AA *2 + Dinuc_Plasmo_Freq$AC + Dinuc_Plasmo_Freq$AG + Dinuc_Plasmo_Freq$AT + Dinuc_Plasmo_Freq$CA + Dinuc_Plasmo_Freq$GA + Dinuc_Plasmo_Freq$TA)/(Dinuc_Plasmo_Freq$X.DINUCLEOTIDES * 2)

Plasmo_FreqC <- (Dinuc_Plasmo_Freq$CC *2 + Dinuc_Plasmo_Freq$CA + Dinuc_Plasmo_Freq$CG + Dinuc_Plasmo_Freq$CT + Dinuc_Plasmo_Freq$AC + Dinuc_Plasmo_Freq$GC + Dinuc_Plasmo_Freq$TC)/(Dinuc_Plasmo_Freq$X.DINUCLEOTIDES * 2)

Plasmo_FreqG <- (Dinuc_Plasmo_Freq$GG *2 + Dinuc_Plasmo_Freq$GA + Dinuc_Plasmo_Freq$GC + Dinuc_Plasmo_Freq$GT + Dinuc_Plasmo_Freq$CG + Dinuc_Plasmo_Freq$AG + Dinuc_Plasmo_Freq$TG)/(Dinuc_Plasmo_Freq$X.DINUCLEOTIDES * 2)

Plasmo_FreqT <- (Dinuc_Plasmo_Freq$TT *2 + Dinuc_Plasmo_Freq$TC + Dinuc_Plasmo_Freq$TG + Dinuc_Plasmo_Freq$TA + Dinuc_Plasmo_Freq$CT + Dinuc_Plasmo_Freq$GT + Dinuc_Plasmo_Freq$AT)/(Dinuc_Plasmo_Freq$X.DINUCLEOTIDES * 2)

## rho scores
Dinuc_Plasmo <- c()

Dinuc_Plasmo[1] <- (Dinuc_Plasmo_Freq$AA/Dinuc_Plasmo_Freq$X.DINUCLEOTIDES)/(Plasmo_FreqA*Plasmo_FreqA)

Dinuc_Plasmo[2] <- (Dinuc_Plasmo_Freq$AC/Dinuc_Plasmo_Freq$X.DINUCLEOTIDES)/(Plasmo_FreqA*Plasmo_FreqC)

Dinuc_Plasmo[3] <- (Dinuc_Plasmo_Freq$AG/Dinuc_Plasmo_Freq$X.DINUCLEOTIDES)/(Plasmo_FreqA*Plasmo_FreqG)

Dinuc_Plasmo[4] <- (Dinuc_Plasmo_Freq$AT/Dinuc_Plasmo_Freq$X.DINUCLEOTIDES)/(Plasmo_FreqA*Plasmo_FreqT)

Dinuc_Plasmo[6] <- (Dinuc_Plasmo_Freq$CC/Dinuc_Plasmo_Freq$X.DINUCLEOTIDES)/(Plasmo_FreqC*Plasmo_FreqC)

Dinuc_Plasmo[5] <- (Dinuc_Plasmo_Freq$CA/Dinuc_Plasmo_Freq$X.DINUCLEOTIDES)/(Plasmo_FreqC*Plasmo_FreqA)

Dinuc_Plasmo[7] <- (Dinuc_Plasmo_Freq$CG/Dinuc_Plasmo_Freq$X.DINUCLEOTIDES)/(Plasmo_FreqC*Plasmo_FreqG)

Dinuc_Plasmo[8] <- (Dinuc_Plasmo_Freq$CT/Dinuc_Plasmo_Freq$X.DINUCLEOTIDES)/(Plasmo_FreqC*Plasmo_FreqT)

Dinuc_Plasmo[11] <- (Dinuc_Plasmo_Freq$GG/Dinuc_Plasmo_Freq$X.DINUCLEOTIDES)/(Plasmo_FreqG*Plasmo_FreqG)

Dinuc_Plasmo[9] <- (Dinuc_Plasmo_Freq$GA/Dinuc_Plasmo_Freq$X.DINUCLEOTIDES)/(Plasmo_FreqG*Plasmo_FreqA)

Dinuc_Plasmo[10] <- (Dinuc_Plasmo_Freq$GC/Dinuc_Plasmo_Freq$X.DINUCLEOTIDES)/(Plasmo_FreqG*Plasmo_FreqC)

Dinuc_Plasmo[12] <- (Dinuc_Plasmo_Freq$GT/Dinuc_Plasmo_Freq$X.DINUCLEOTIDES)/(Plasmo_FreqG*Plasmo_FreqT)

Dinuc_Plasmo[16] <- (Dinuc_Plasmo_Freq$TT/Dinuc_Plasmo_Freq$X.DINUCLEOTIDES)/(Plasmo_FreqT*Plasmo_FreqT)

Dinuc_Plasmo[13] <- (Dinuc_Plasmo_Freq$TA/Dinuc_Plasmo_Freq$X.DINUCLEOTIDES)/(Plasmo_FreqT*Plasmo_FreqA)

Dinuc_Plasmo[14] <- (Dinuc_Plasmo_Freq$TC/Dinuc_Plasmo_Freq$X.DINUCLEOTIDES)/(Plasmo_FreqT*Plasmo_FreqC)

Dinuc_Plasmo[15] <- (Dinuc_Plasmo_Freq$TG/Dinuc_Plasmo_Freq$X.DINUCLEOTIDES)/(Plasmo_FreqT*Plasmo_FreqG)

## create a dataframe with virus/host as rows, dinuc as columns, and rho scores as values

Dinuc_MDS_Dataframe <- Dinuc_Reepicheep %>%
  cbind(Dinuc_Human) %>%
  cbind(Dinuc_Mouse) %>%
  cbind(Dinuc_Crypto) %>%
  cbind(Dinuc_Plasmo) %>%
  cbind(Dinuc_Matryoshka1)

colnames(Dinuc_MDS_Dataframe)[1] <- "Dinuc_Reepicheep"

Dinuc_MDS_Dataframe <- t(Dinuc_MDS_Dataframe)

## mds
d2 <- dist(Dinuc_MDS_Dataframe)
fit2 <- cmdscale(d2 ,eig=TRUE, k=2)
fit2

Fit2 <- data.frame(fit2$points)

## set names for the MDS plot
row.names(Fit2) <- c("Reepicheep", "Human", "Mouse", "Cryptosporidium", "Plasmodium", "Matryoshka 1")

## plot
ggplot(Fit2 , aes(x = X1, y = X2)) +
  geom_point() +
  ggtitle("Dinucleotide Rho Scores Metric MDS") +
  xlab("Coordinate 1") +
  xlim(c(-.6, .6)) +
  ylim(c(-.5, .65)) +
  ylab("Coordinate 2") +
  geom_label_repel(aes(label = row.names(Fit2)),
                   box.padding   = 0.35, 
                   point.padding = 0.5,
                   segment.color = 'grey50') +
  theme_bw() + 
  theme(panel.background = element_blank(), 
        panel.grid.major = element_blank(),  #remove major-grid labels
        panel.grid.minor = element_blank(),  #remove minor-grid labels
        plot.background = element_blank(), 
        plot.title = element_text(hjust = 0.5)) #center the title


# codon pair bias 
  
## calculate CPS for the narnaviruses using CPBias::CPBtable
Reepicheep_CPB <- CPBtable("reepicheep.fasta")
Reepicheep_CPB <- Reepicheep_CPB$complete.CPBtable %>%
  arrange(codonpair)

Matryoshka_1_CPB <- CPBtable("Matryoshka1.fasta") 
Matryoshka_1_CPB <- Matryoshka_1_CPB$complete.CPBtable %>%
  arrange(codonpair)

## load in file of values for mouse and humans which have more complex genomes
Host_CPS <- read_excel("hostCodonPairScores.xlsx")

## arrange codon pairs so they're in the same order as in the other CPB tables
Host_CPS$`Codon pair` <- str_remove(Host_CPS$`Codon pair`, "-") 
Host_CPS$`Codon pair` <- tolower(Host_CPS$`Codon pair`)

### Human
Human_CPB <- Host_CPS[-1, c(1,2,3)] %>%
  arrange(`Codon pair`)
Human_CPB$CPS <- as.numeric(Human_CPB$CPS)

### Mouse
Mouse_CPB <- Host_CPS[-1, c(1,2,5)] %>%
  arrange(`Codon pair`)
Mouse_CPB$...5 <- as.numeric(Mouse_CPB$...5)

## for plasmo and crypto i have the codon pair COUNTS from HIVE Biochemistry

### Plasmo

## load in the file and get it so its just codon pairs
Plasmo_Codon_Pair <- read.csv("plasmoCodonPair.csv")
Plasmo_Codon_Pair <- t(Plasmo_Codon_Pair)
Plasmo_Codon_Pair <- Plasmo_Codon_Pair[-c(1:10),]
Plasmo_Codon_Pair <- data.frame(Plasmo_Codon_Pair)
Plasmo_Codon_Pair[,1] <- as.numeric(Plasmo_Codon_Pair[,1])

## split up the codon pairs
codon_pairs <- rownames(Plasmo_Codon_Pair)
split_codon_pairs <- strsplit(codon_pairs, "(?<=.{3})", perl = TRUE)

## put them in their own columns
for(i in 1:length(split_codon_pairs)){
  Plasmo_Codon_Pair[i,2] <- split_codon_pairs[[i]][1]
}

for(i in 1:length(split_codon_pairs)){
  Plasmo_Codon_Pair[i,3] <- split_codon_pairs[[i]][2]
}

## translate the codon pairs into their corresponding AAs
Plasmo_Codon_Pair[,4] <- GENETIC_CODE[Plasmo_Codon_Pair[,2]] 
AA_list_factors <- as.factor(Plasmo_Codon_Pair[,4])
Plasmo_Codon_Pair[,5] <- GENETIC_CODE[Plasmo_Codon_Pair[,3]]

colnames(Plasmo_Codon_Pair) <- c("frequency", "codon1", "codon2", "AA1", "AA2")

## setting this up so I can make a table of codons later
Plasmo_Codon_Pair[,2] <- as.factor(Plasmo_Codon_Pair[,2]) 

## make columns that combine the codon/AA pairs for later
Plasmo_Codon_Pair <- Plasmo_Codon_Pair %>% 
  unite("Codon_pair", codon1:codon2, sep = "", remove = FALSE) %>%
  unite("AA_pair", AA1:AA2, sep = "", remove = FALSE)

## make a table of all AA combinations
AA_factor_grid <- expand.grid(levels(AA_list_factors), levels(AA_list_factors))
AA_factor_grid <- AA_factor_grid %>% 
  unite("AA_pair", Var1:Var2, sep = "", remove = FALSE)
AA_factor_grid[,4] <- 0

## count the different AA combinations
for(i in 1:441){
  for(j in 1:4096){
    if(AA_factor_grid[i,2] == Plasmo_Codon_Pair[j,6] & AA_factor_grid[i,3] == Plasmo_Codon_Pair[j,7]){
      AA_factor_grid[i,4] <- AA_factor_grid[i,4] + Plasmo_Codon_Pair[j,1]
    }
  }
}

## make a table of the different AAs
AA_factor <- as.data.frame(levels(AA_list_factors))
AA_factor[,2] <-0

## count the occurrence of each AA
for(i in 1:21){
  for(j in 1:4096){
    if(AA_factor[i,1] == Plasmo_Codon_Pair[j,6]){
      AA_factor[i,2] <- AA_factor[i,2] + Plasmo_Codon_Pair[j,1]
    }
    if(AA_factor[i,1] == Plasmo_Codon_Pair[j,7]){
      AA_factor[i,2] <- AA_factor[i,2] + Plasmo_Codon_Pair[j,1]
    }
  }
}

## make a table of the different codons
Codon_factor <- as.data.frame(levels(Plasmo_Codon_Pair[,3]))
Codon_factor[,2] <- 0

## count the occurence of each codon
for(i in 1:64){
  for(j in 1:4096){
    if(Codon_factor[i,1] == Plasmo_Codon_Pair[j,3]){
      Codon_factor[i,2] <- Codon_factor[i,2] + Plasmo_Codon_Pair[j,1]
    }
    if(Codon_factor[i,1] == Plasmo_Codon_Pair[j,4]){
      Codon_factor[i,2] <- Codon_factor[i,2] + Plasmo_Codon_Pair[j,1]
    }
  }
}

## calculate the CPS 
Plasmo_Codon_Pair[,8] <- 0

for(i in 1:4096){
  CPS <- (Plasmo_Codon_Pair[i,1] * AA_factor[match(Plasmo_Codon_Pair[i,6], AA_factor[,1]),2] * AA_factor[match(Plasmo_Codon_Pair[i,7], AA_factor[,1]),2])/ (AA_factor_grid[match(Plasmo_Codon_Pair[i,5], AA_factor_grid[,1]), 4] * Codon_factor[match(Plasmo_Codon_Pair[i,3], Codon_factor[,1]),2] * Codon_factor[match(Plasmo_Codon_Pair[i,4], Codon_factor[,1]),2] )
  Plasmo_Codon_Pair[i,8] <- CPS
}

Plasmo_Codon_Pair[,9] <- log(Plasmo_Codon_Pair[,8])


### Crypto

## load in the file and get it so its just codon pairs
Crypto_Codon_Pair <- read.csv("cryptoCodonPair.csv")
Crypto_Codon_Pair <- t(Crypto_Codon_Pair)
Crypto_Codon_Pair <- Crypto_Codon_Pair[-c(1:10),]
Crypto_Codon_Pair <- data.frame(Crypto_Codon_Pair)
Crypto_Codon_Pair[,1] <- as.numeric(Crypto_Codon_Pair[,1])


## split up the codon pairs
codon_pairs <- rownames(Crypto_Codon_Pair)
split_codon_pairs <- strsplit(codon_pairs, "(?<=.{3})", perl = TRUE)

## put them in their own columns
for(i in 1:length(split_codon_pairs)){
  Crypto_Codon_Pair[i,2] <- split_codon_pairs[[i]][1]
}

for(i in 1:length(split_codon_pairs)){
  Crypto_Codon_Pair[i,3] <- split_codon_pairs[[i]][2]
}

## translate the codon pairs into their corresponding AAs
Crypto_Codon_Pair[,4] <- GENETIC_CODE[Crypto_Codon_Pair[,2]] 
AA_list_factors_2 <- as.factor(Crypto_Codon_Pair[,4])
Crypto_Codon_Pair[,5] <- GENETIC_CODE[Crypto_Codon_Pair[,3]]

colnames(Crypto_Codon_Pair) <- c("frequency", "codon1", "codon2", "AA1", "AA2")

## setting this up so I can make a table of codons later
Crypto_Codon_Pair[,2] <- as.factor(Crypto_Codon_Pair[,2]) 

## make columns that combine the codon/AA pairs for later
Crypto_Codon_Pair <- Crypto_Codon_Pair %>% 
  unite("Codon_pair", codon1:codon2, sep = "", remove = FALSE) %>%
  unite("AA_pair", AA1:AA2, sep = "", remove = FALSE)

## make a table of all AA combinations
AA_factor_2_grid <- expand.grid(levels(AA_list_factors_2), levels(AA_list_factors_2))
AA_factor_2_grid <- AA_factor_2_grid %>% 
  unite("AA_pair", Var1:Var2, sep = "", remove = FALSE)
AA_factor_2_grid[,4] <- 0

## count the different AA combinations
for(i in 1:441){
  for(j in 1:4096){
    if(AA_factor_2_grid[i,2] == Crypto_Codon_Pair[j,6] & AA_factor_2_grid[i,3] == Crypto_Codon_Pair[j,7]){
      AA_factor_2_grid[i,4] <- AA_factor_2_grid[i,4] + Crypto_Codon_Pair[j,1]
    }
  }
}

## make a table of the different AAs
AA_factor_2 <- as.data.frame(levels(AA_list_factors_2))
AA_factor_2[,2] <-0

## count the occurrence of each AA
for(i in 1:21){
  for(j in 1:4096){
    if(AA_factor_2[i,1] == Crypto_Codon_Pair[j,6]){
      AA_factor_2[i,2] <- AA_factor_2[i,2] + Crypto_Codon_Pair[j,1]
    }
    if(AA_factor_2[i,1] == Crypto_Codon_Pair[j,7]){
      AA_factor_2[i,2] <- AA_factor_2[i,2] + Crypto_Codon_Pair[j,1]
    }
  }
}

## make a table of the different codons
Codon_factor_2 <- as.data.frame(levels(Crypto_Codon_Pair[,3]))
Codon_factor_2[,2] <- 0

## count the occurence of each codon
for(i in 1:64){
  for(j in 1:4096){
    if(Codon_factor_2[i,1] == Crypto_Codon_Pair[j,3]){
      Codon_factor_2[i,2] <- Codon_factor_2[i,2] + Crypto_Codon_Pair[j,1]
    }
    if(Codon_factor_2[i,1] == Crypto_Codon_Pair[j,4]){
      Codon_factor_2[i,2] <- Codon_factor_2[i,2] + Crypto_Codon_Pair[j,1]
    }
  }
}

## calculate the CPS 
Crypto_Codon_Pair[,8] <- 0

for(i in 1:4096){
  CPS <- (Crypto_Codon_Pair[i,1] * AA_factor_2[match(Crypto_Codon_Pair[i,6], AA_factor_2[,1]),2] * AA_factor_2[match(Crypto_Codon_Pair[i,7], AA_factor_2[,1]),2])/ (AA_factor_2_grid[match(Crypto_Codon_Pair[i,5], AA_factor_2_grid[,1]), 4] * Codon_factor_2[match(Crypto_Codon_Pair[i,3], Codon_factor_2[,1]),2] * Codon_factor_2[match(Crypto_Codon_Pair[i,4], Codon_factor_2[,1]),2] )
  Crypto_Codon_Pair[i,8] <- CPS
}

Crypto_Codon_Pair[,9] <- log(Crypto_Codon_Pair[,8])

## get all of the data frames in the same order
## removing start/stop codons from crypto and plasmo
Crypto_Codon_Pair <- arrange(Crypto_Codon_Pair, Codon_pair)
Crypto_Codon_Pair <- Crypto_Codon_Pair[!grepl('\\*', Crypto_Codon_Pair$AA_pair), ]
Plasmo_Codon_Pair <- arrange(Plasmo_Codon_Pair, Codon_pair)
Plasmo_Codon_Pair <- Plasmo_Codon_Pair[!grepl('\\*', Plasmo_Codon_Pair$AA_pair), ]

## create a data frame with the CPS scores for the narnaviruses and hosts
MDS_CPS_Data_Frame <- data.frame(Reepicheep_CPB[,15], Mouse_CPB[,3], Plasmo_Codon_Pair[,9], Matryoshka_1_CPB[,15], Crypto_Codon_Pair[,9], Human_CPB[,3]) %>%
  t()

colnames(MDS_CPS_Data_Frame) <- Reepicheep_CPB$codonpair

row.names(MDS_CPS_Data_Frame) <- c("Reepicheep", "Mouse", "Plasmodium", "Matryoshka 1", "Cryptosporidium", "Human")


## make a copy of the data frame for other analysis
MDS_CPS_Data_Frame2 <- MDS_CPS_Data_Frame

## un-natural log transforming the CPS scores (many of the codon pairs were not present in narnavirus genomes)
MDS_CPS_Data_Frame[is.na(MDS_CPS_Data_Frame)] <- 0
MDS_CPS_Data_Frame[MDS_CPS_Data_Frame == -Inf] <- 0

MDS_CPS_Data_Frame <- apply(MDS_CPS_Data_Frame, MARGIN = 2, FUN = exp)

## transform all the 1s back to 0s 
MDS_CPS_Data_Frame[MDS_CPS_Data_Frame == 1] <- 0

## MDS
d6 <- dist(MDS_CPS_Data_Frame) 
fit6 <- cmdscale(d6,eig=TRUE, k=2) 
fit6 

Fit6 <- data.frame(fit6$points)

#names for the MDS plot
row.names(Fit6) <- c("Reepicheep", "Mouse", "Plasmodium", "Matryoshka 1", "Cryptosporidium", "Human")

# plot solution
ggplot(Fit6, aes(x = X1, y = X2)) +
  geom_point() +
  ggtitle("Untransformed Codon Pair Score Metric MDS") +
  xlab("Coordinate 1") +
  ylab("Coordinate 2") +
  xlim(c(-300, 200)) +
  ylim(c(-150, 150)) +
  geom_label_repel(aes(label = row.names(Fit6)),
                   box.padding   = 0.35, 
                   point.padding = 0.5,
                   segment.color = 'grey50') +
  theme_bw() + 
  theme(panel.background = element_blank(), 
        panel.grid.major = element_blank(),  #remove major-grid labels
        panel.grid.minor = element_blank(),  #remove minor-grid labels
        plot.background = element_blank(), 
        plot.title = element_text(hjust = 0.5)) #center the title


## remove all columns with NA from the copy dataframe
MDS_CPS_Data_Frame2[MDS_CPS_Data_Frame2 == -Inf] <- NA
MDS_CPS_Data_Frame2[is.nan(MDS_CPS_Data_Frame2)] <- NA

na_cols <- apply(MDS_CPS_Data_Frame2,2,function(x)!any(is.na(x)))
MDS_CPS_Data_Frame2 <- MDS_CPS_Data_Frame2[ ,na_cols]

## keep the log transformation

## MDS
d7 <- dist(MDS_CPS_Data_Frame2) 
fit7 <- cmdscale(d7,eig=TRUE, k=2) 
fit7

Fit7 <- data.frame(fit7$points)

#names for the MDS plot
row.names(Fit7) <- c("Reepicheep", "Mouse", "Plasmodium", "Matryoshka 1", "Cryptosporidium", "Human")

# plot solution
ggplot(Fit7, aes(x = X1, y = X2)) +
  geom_point() +
  ggtitle("Codon Pair Score (complete cases) Metric MDS") +
  xlab("Coordinate 1") +
  ylab("Coordinate 2") +
  geom_label_repel(aes(label = row.names(Fit7)),
                   box.padding   = 0.35, 
                   point.padding = 0.5,
                   segment.color = 'grey50') +
  theme_bw() + 
  theme(panel.background = element_blank(), 
        panel.grid.major = element_blank(),  #remove major-grid labels
        panel.grid.minor = element_blank(),  #remove minor-grid labels
        plot.background = element_blank(), 
        plot.title = element_text(hjust = 0.5)) #center the title
