# Data processing
library("dplyr")
library("tibble")

# Load fermentation data
ferm.p = read.delim( "raw-data/AvgKO_Ferm_P.txt" )

# Molecule names should be char
ferm.p$Molecule.Name = as.character( ferm.p$Molecule.Name )

# Distinguish the two glyceric acid rows
ferm.p = mutate( ferm.p, Molecule.Name = replace( Molecule.Name, Molecule.Name == "Glyceric acid", paste0( "Glyceric acid ", 1:2 ) ) )

# Save
save( ferm.p, file = "clean-data/fermentation.RData" )

# Split out proteins
ferm.p.proteins = ferm.p %>% filter( Molecule.Type == "Protein" ) %>% select( -(Molecule.Type) )

# Drop type column in ferm.p
ferm.p = ferm.p %>% select( -Molecule.Type )

# Cols to names
ferm.p = column_to_rownames( ferm.p, "Molecule.Name" )
ferm.p.proteins = column_to_rownames( ferm.p.proteins, "Molecule.Name" )

# Save to clean data
save( ferm.p, ferm.p.proteins, file = "clean-data/ferm-for-nem.RData" )

# Load respiratory data
#resp.p = read.delim( "raw-data/AvgKO_Resp_P.txt" )

# Save to clean data

# Load resp-RDR data

# Save to clean data