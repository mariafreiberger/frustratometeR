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


create_object <-function(PdbFile=NULL, PdbID=NULL, Chain=NULL)
{
  if(!is.null(PdbFile))
  {
    Pdb <- read.pdb(PdbFile)
    PdbBase <- basename.pdb(PdbFile)
    Pdb[["PdbBase"]] <- PdbBase

  }else if(!is.null(PdbID))
  {
    if(is.null(Chain))
    {
      boolsplit=F
    }else{
      boolsplit=T
    }
    tempfolder <-tempfile()
    PdbAux<-get.pdb(PdbID, split = T, path = tempfolder)
    if(is.null(Chain))
    {
      Pdb <- read.pdb(paste(tempfolder,"/",PdbID,".pdb", sep=""))
    }else{
      Pdb <- read.pdb(paste(tempfolder,"/split_chain/",PdbID, "_", Chain, ".pdb", sep=""))
    }
    Pdb[["PdbBase"]] <- PdbID
  }
  return(Pdb)
}


check=create_object(PdbID = "1n0r")

