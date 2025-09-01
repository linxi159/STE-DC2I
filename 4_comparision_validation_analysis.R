
#######comparison validate analysis #######

#####[1]Loading validated and computationally predicted cancer driver genes (CDGs)
cdg_CGCv100_canonical_CRC <- read.csv("CGC_v100_GRCh38_CRC_DriverGege_ordered.csv", header = TRUE, sep = ",")
cdg_NCG7.1_canonical <- read.csv("NCG7.1_canonical_drivers.csv", header = TRUE, sep = ",")
cdg_NCG7.1_candidate_CRC <- read.csv("NCG7.1_candidate_drivers_CRC.csv", header = TRUE, sep = ",")


#####[2]The number of cancer driver genes (CDGs) validated by other methods and computationally predicted#######
#(1)IntOGen: 7 methods
#dndscv  cbase  mutpanning
#oncodriveclustl  hotmaps  smregions
#oncodrivefml
IntOGen_computed_drivers_CRC <- read.csv("./IntOGen_comparasion/IntOGen_computed_drivers_CRC.csv", header = TRUE, sep = ",")
IntOGen_dndscv <- read.csv("./IntOGen_comparasion/IntOGen_computed_drivers_CRC_1_dndscv.csv", header = TRUE, sep = ",")
IntOGen_cbase <- read.csv("./IntOGen_comparasion/IntOGen_computed_drivers_CRC_2_cbase.csv", header = TRUE, sep = ",")
IntOGen_mutpanning <- read.csv("./IntOGen_comparasion/IntOGen_computed_drivers_CRC_3_mutpanning.csv", header = TRUE, sep = ",")
IntOGen_oncodriveclustl <- read.csv("./IntOGen_comparasion/IntOGen_computed_drivers_CRC_4_oncodriveclustl.csv", header = TRUE, sep = ",")
IntOGen_hotmaps <- read.csv("./IntOGen_comparasion/IntOGen_computed_drivers_CRC_5_hotmaps.csv", header = TRUE, sep = ",")
IntOGen_smregions <- read.csv("./IntOGen_comparasion/IntOGen_computed_drivers_CRC_6_smregions.csv", header = TRUE, sep = ",")
IntOGen_oncodrivefml <- read.csv("./IntOGen_comparasion/IntOGen_computed_drivers_CRC_7_oncodrivefml.csv", header = TRUE, sep = ",")

#(2)DriverDBv4: 8 methods
#comet,dawnrank,DIABLO,DriverSubNet
#MOFA, OPLSDA, SDGCCA, SGCCA
DriverDBv4_comet <- read.csv("./DriverDBv4_comparasion/DriverDBv4_drivers_CRC_1_comet.csv", header = TRUE, sep = ",")
DriverDBv4_dawnrank <- read.csv("./DriverDBv4_comparasion/DriverDBv4_drivers_CRC_2_dawnrank.csv", header = TRUE, sep = ",")
DriverDBv4_DIABLO <- read.csv("./DriverDBv4_comparasion/DriverDBv4_drivers_CRC_3_DIABLO.csv", header = TRUE, sep = ",")
DriverDBv4_DriverSubNet <- read.csv("./DriverDBv4_comparasion/DriverDBv4_drivers_CRC_4_DriverSubNet.csv", header = TRUE, sep = ",")
DriverDBv4_MOFA <- read.csv("./DriverDBv4_comparasion/DriverDBv4_drivers_CRC_5_MOFA.csv", header = TRUE, sep = ",")
DriverDBv4_OPLSDA <- read.csv("./DriverDBv4_comparasion/DriverDBv4_drivers_CRC_6_OPLSDA.csv", header = TRUE, sep = ",")
DriverDBv4_SDGCCA <- read.csv("./DriverDBv4_comparasion/DriverDBv4_drivers_CRC_7_SDGCCA.csv", header = TRUE, sep = ",")
DriverDBv4_SGCCA <- read.csv("./DriverDBv4_comparasion/DriverDBv4_drivers_CRC_8_SGCCA.csv", header = TRUE, sep = ",")










