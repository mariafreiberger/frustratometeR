library(devtools)
install_github("gonzaparra/frustratometeR")

library(frustratometeR)
library(bio3d)

PdbFile="/home/gonzalo/Desktop/frustratometer2-master/1n0r.pdb"

ResultsDir="/home/gonzalo/Desktop/"
Modes="configurational"
seqdist=12

Electrostatics_K=3.1


Pdb=calculate_frustration(PdbFile=PdbFile, Modes = "configurational", Electrostatics_K = NULL, ResultsDir = ResultsDir, seqdist=12)

plot_5Andens(Pdb, chain=NULL)
plot_5Adens_proportions(Pdb, chain=NULL)
plot_contact_map(Pdb, chain="A")

view_frustration_pymol(Pdb)
