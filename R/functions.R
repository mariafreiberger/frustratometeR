#' Pdb Equivalences
#'
#' Internal function that produces auxiliar files to run the frustratometer pipeline.
#'
#' @param Pdb Frustration object
#' @return Pdb Frustration object with backbone completed
#'
pdb_equivalences <- function(Pdb)
{
  existing_res <- unique(cbind(Pdb$atom$chain[which(Pdb$atom$type=="ATOM")], Pdb$atom$resno[which(Pdb$atom$type=="ATOM")], Pdb$atom$resid[which(Pdb$atom$type=="ATOM")]))
  equivalences <- cbind(existing_res[,1], seq(1:length(existing_res[,1])), existing_res[,2], existing_res[,3])
  write.table(equivalences, file = paste(Pdb$JobDir, "/", Pdb$PdbBase, ".pdb_equivalences.txt", sep=""), quote = F, col.names = F, row.names = F, sep="\t")

  Pdb[["equivalences"]] <- equivalences
  return(Pdb)
}

#' Check Backbone Complete
#'
#' Checks the backbone of a given protein structure to be processed by the frustratometer pipeline
#'
#' @param Pdb Frustration object
#' @return Pdb Frustration object with backbone completed

check_backbone_complete <- function(Pdb)
{
  # Select backbone atoms
  backpdb <- atom.select(Pdb, elety=c("N", "CA", "C", "O", "CB"), value=TRUE)
  # Only for those atoms with ATOM coords
  atom_res <- backpdb$atom[backpdb$atom$type=="ATOM",]
  # Check there are no missing atoms in the backbone
  Complete <- length(which(atom_res$elety=="CA"))==dim(Pdb$equivalences)[1] & length(which(atom_res$elety=="O"))==dim(Pdb$equivalences)[1] & length(which(atom_res$elety=="CB")) + length(which(Pdb$equivalences[,4]=="GLY"))==dim(Pdb$equivalences)[1]

  # print(paste("CA", length(which(atom_res$elety=="CA"))))
  # print(paste("O", length(which(atom_res$elety=="O"))))
  # print(paste("CB", length(which(atom_res$elety=="CB"))))

  # print(paste("C", length(which(atom_res$elety=="C"))))
  # print(paste("N", length(which(atom_res$elety=="N"))))

  return(Complete)
}

#' Complete Backbone
#'
#' Completes the backbone of a given protein structure to be processed by the frustratometer pipeline
#'
#' @param Pdb Frustration object
#' @return Pdb Frustration object with backbone completed
complete_backbone <- function(Pdb)
{
  if(!check_backbone_complete(Pdb))
  {
    system(paste("python ", Pdb$scriptsDir, "MissingAtoms.py ",  Pdb$JobDir, Pdb$PdbBase, ".pdb", sep=""))
    system(paste("mv ", Pdb$JobDir, Pdb$PdbBase, ".pdb_completed", " ", Pdb$JobDir, Pdb$PdbBase, ".pdb", sep=""))
  }
}

#' Calculate Frustration
#'
#' Calculates local energetic frustration for a protein structure
#'
#' @param PdbFile File containing the protein structure. The full path to the file is needed
#' @param PdbID File containing the protein structure. The full path to the file is needed
#' @param Chain File containing the protein structure. The full path to the file is needed
#' @param Modes Local frustration index to be calculated (configurational, mutational, singleresidue). Default=configurational
#' @param Electrostatics_K K constant to use in the electrostatics mode. Default: NULL (no electrostatics is considered).
#' @param seqdist Sequence at which contacts are considered to interact.
#' @param ResultsDir Path to the folder where results will be stored.
#' @return Pdb Frustration object
#'
#' @export

calculate_frustration <- function(PdbFile=NULL, PdbID=NULL, Chain=NULL, Electrostatics_K=NULL, seqdist=12, Modes="configurational", ResultsDir)
{

  if(is.null(PdbFile))
  {
    tempfolder <-tempfile()
    if(is.null(Chain))
    {
      boolsplit=F
    }else{
      boolsplit=T
    }
    PdbAux<-get.pdb(PdbID, split = boolsplit, path = tempfolder)
    if(is.null(Chain))
    {
      PdbFile=paste(tempfolder,"/",PdbID,".pdb", sep="")
    }else{
      PdbFile=paste(tempfolder,"/split_chain/",PdbID, "_", Chain, ".pdb", sep="")
    }
  }

  Pdb <- read.pdb(PdbFile, ATOM.only=T)

  PdbBase <- basename.pdb(PdbFile)

  JobDir=paste(ResultsDir, PdbBase, ".done/", sep="")

  Pdb[["PdbBase"]] <- PdbBase
  Pdb[["JobDir"]] <- JobDir
  Pdb[["scriptsDir"]] <- paste(find.package("frustratometeR"), "/Scripts/", sep="")
  Pdb[["mode"]] <- Modes

  #Creates JobDir
  system(paste("mkdir ", Pdb$JobDir, sep=""))
  system(paste("cp ", PdbFile, " ", Pdb$JobDir, sep=""))

  # Save equivalences
  Pdb <- pdb_equivalences(Pdb)
  complete_backbone(Pdb)

  # PdbBase <- basename.pdb(PdbFile)

  print("Preparing files..")
  #Prepare the PDB file to get awsem input files, create the workdir and move neccessary files to it.
  system(paste("cd ", Pdb$JobDir, "; pwd; ", Pdb$scriptsDir, "AWSEMFiles/AWSEMTools/PdbCoords2Lammps.sh ", Pdb$PdbBase, " ", Pdb$PdbBase, " ", Pdb$scriptsDir, sep=""))
  system(paste("cp ", Pdb$scriptsDir, "AWSEMFiles/*.dat* ", Pdb$JobDir, sep=""))

  print("Setting options...")
  #Modify the .in file to run a single step - Modify the fix_backbone file to change the mode and set options

  system(paste("cd ", JobDir, "; sed -i 's/run		10000/run		0/g' ", Pdb$PdbBase, ".in; sed -i 's/mutational/", Pdb$mode, "/g' fix_backbone_coeff.data",  sep=""))

  if(!is.null(Electrostatics_K))
  {
    print("Setting electrostatics...")
    system(paste("cd ", Pdb$JobDir, "; sed -i 's/\\[DebyeHuckel\\]-/\\[DebyeHuckel\\]/g' fix_backbone_coeff.data; sed -i 's/4.15 4.15 4.15/", Electrostatics_K, " ", Electrostatics_K, " ", Electrostatics_K, "/g' fix_backbone_coeff.data;", sep=""))
    print("Setting electrostatics...")
    system(paste("cd ", Pdb$JobDir, "; python ", Pdb$scriptsDir, "Pdb2Gro.py ", Pdb$PdbBase, ".pdb ", Pdb$PdbBase, ".pdb.gro; perl ", Pdb$scriptsDir, "GenerateChargeFile.pl ", Pdb$PdbBase, ".pdb.gro > ", JobDir, "charge_on_residues.dat", sep=""))
  }

  print("calculating...")
  system(paste("cp ", Pdb$scriptsDir, "lmp_serial_", seqdist, " ", Pdb$JobDir, "; cd ", Pdb$JobDir, "; ./lmp_serial_", seqdist, " < ", Pdb$PdbBase, ".in", sep=""))

  if(Pdb$mode == "configurational" | Pdb$mode == "mutational")
  {
    system(paste("perl ", Pdb$scriptsDir, "5Adens.pl ", Pdb$PdbBase, ".pdb ", gsub(".$", "", Pdb$JobDir), " ", Pdb$mode, sep=""))
  }
  system(paste("perl ", Pdb$scriptsDir, "RenumFiles.pl ", Pdb$PdbBase, " ", Pdb$JobDir, " ", Pdb$mode, sep="" ))

  system(paste("perl ", Pdb$scriptsDir, "GenerateVisualizations.pl ", Pdb$PdbBase, "_", Pdb$mode, ".pdb_auxiliar ", Pdb$PdbBase, " ", gsub(".$", "", Pdb$JobDir), " ", Pdb$mode, sep=""))
  system(paste("cp ", Pdb$scriptsDir, "draw_links.py ", Pdb$JobDir, sep=""))

  return(Pdb)
}
