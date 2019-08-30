library(devtools)
install_github("gonzaparra/frustratometeR")

library(frustratometeR)
library(bio3d)

# Note: 1mj0 no anda

PdbFile="/home/gonzalo/frustratometeR/1n0r.pdb"

ResultsDir="/home/gonzalo/Desktop/"
Modes="configurational"

# Calculate frustration for a givenfile
Pdb=calculate_frustration(PdbFile=PdbFile, Modes = "configurational", ResultsDir = ResultsDir)

# Calculate frustration for a structure to be downloaded from the PDB
Pdb=calculate_frustration(PdbID="1mj0", Modes = "configurational", ResultsDir = ResultsDir)

# Calculate frustration for a particular chain in a structure to be downloaded from the PDB
Pdb=calculate_frustration(PdbID="1ikn", Chain="D", Modes = "configurational", ResultsDir = ResultsDir)


plot_5Andens(Pdb, chain=NULL)
plot_5Adens_proportions(Pdb, chain=NULL)
plot_contact_map(Pdb, chain="A")

view_frustration_pymol(Pdb)
