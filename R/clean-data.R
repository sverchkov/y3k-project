# Data processing
library("dplyr")
library("tibble")

# Load fermentation data
ferm.p = read.delim( "raw-data/AvgKO_Ferm_P.txt" )

# Molecule names should be char
ferm.p$Molecule.Name = as.character( ferm.p$Molecule.Name )

# Cleanup: Drop empty column and distinguish the two glyceric acid rows
ferm.p = ferm.p %>%
  select( -X ) %>%
  mutate( Molecule.Name = replace( Molecule.Name, Molecule.Name == "Glyceric acid", paste0( "Glyceric acid ", 1:2 ) ) )

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

print("Fermentation data prepped.")

####
####

# Load respiratory data
resp.p = read.delim( "raw-data/AvgKO_Resp_P.txt" )

# Molecule names should be char
resp.p$Molecule.Name = as.character( resp.p$Molecule.Name )

# Cleanup: drop empty column
resp.p = resp.p %>% select( -X )

# Save to clean data
save( resp.p, file = "clean-data/respiration.RData" )

# Split out proteins
resp.p.proteins = resp.p %>% filter( Molecule.Type == "Protein" ) %>% select( -(Molecule.Type) )

# Drop type column in resp.p
resp.p = resp.p %>% select( -Molecule.Type )

# Cols to names
resp.p = column_to_rownames( resp.p, "Molecule.Name" )
resp.p.proteins = column_to_rownames( resp.p.proteins, "Molecule.Name" )

# Save to clean data
save( resp.p, resp.p.proteins, file = "clean-data/resp-for-nem.RData" )

print("Respiration data prepped.")

####
####

# Load resp-RDR data
resp.rdr.p = read.delim( "raw-data/AvgKO_Resp-RDR_P.txt" )

# Molecule names should be char
resp.rdr.p$Molecule.Name = as.character( resp.rdr.p$Molecule.Name )

# Cleanup: drop empty column
resp.rdr.p = resp.rdr.p %>% select( -X )

# Save to clean data
save( resp.rdr.p, file = "clean-data/respiration-rdr.RData" )

# Split out proteins
resp.rdr.p.proteins = resp.rdr.p %>% filter( Molecule.Type == "Protein" ) %>% select( -(Molecule.Type) )

# Drop type column in resp.rdr.p
resp.rdr.p = resp.rdr.p %>% select( -Molecule.Type )

# Cols to names
resp.rdr.p = column_to_rownames( resp.rdr.p, "Molecule.Name" )
resp.rdr.p.proteins = column_to_rownames( resp.rdr.p.proteins, "Molecule.Name" )

# Save to clean data
save( resp.rdr.p, resp.rdr.p.proteins, file = "clean-data/resp-rdr-for-nem.RData" )

print("Respiration-RDR data prepped.")

####
####
