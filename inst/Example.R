library(devtools)
install_github("gonzaparra/frustratometeR")

library(frustratometeR)
library(bio3d)

PdbFile="/home/gonzalo/frustratometeR/1a0i_A.pdb"

ResultsDir="/home/gonzalo/Desktop/"
Modes="configurational"

Pdb=calculate_frustration(PdbFile=PdbFile, Modes = "configurational", ResultsDir = ResultsDir)

plot_5Andens(Pdb, chain=NULL)
plot_5Adens_proportions(Pdb, chain=NULL)
plot_contact_map(Pdb, chain="A")

view_frustration_pymol(Pdb)
