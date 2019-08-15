# dependencies:
# numpy
# bio3d

library(bio3d)

pdb_equivalences <- function(Pdb, JobDir, PdbBase)
{
  existing_res <- unique(cbind(Pdb$atom$chain[which(Pdb$atom$type=="ATOM")], Pdb$atom$resno[which(Pdb$atom$type=="ATOM")], Pdb$atom$resid[which(Pdb$atom$type=="ATOM")]))
  equivalences <- cbind(existing_res[,1], seq(1:length(existing_res[,1])), existing_res[,2], existing_res[,3])
  write.table(equivalences, file = paste(JobDir, "/", PdbBase, ".pdb_equivalences.txt", sep=""), quote = F, col.names = F, row.names = F, sep="\t")
  
  Pdb[["equivalences"]] <- equivalences
  return(Pdb)
}

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

complete_backbone <- function(Pdb, scriptsDir, JobDir, PdbBase)
{
  if(!check_backbone_complete(Pdb))
  {
    system(paste("python ", scriptsDir, "MissingAtoms.py ",  JobDir, PdbBase, ".pdb", sep=""))
    system(paste("mv ", JobDir, PdbBase, ".pdb_completed", " ", JobDir, PdbBase, ".pdb", sep=""))
  }
}

calculate_frustration <- function(PdbFile=PdbFile, Electrostatics_K=3.1, seqdist=12, Modes="configurational", scriptsDir=scriptsDir, ResultsDir="/home/gonzalo/Desktop/")
{

  Pdb <- read.pdb(PdbFile)
  
  PdbBase <- basename.pdb(PdbFile)
  
  JobDir=paste(ResultsDir, PdbBase, ".done/", sep="")
  
  Pdb[["PdbBase"]] <- PdbBase
  Pdb[["JobDir"]] <- JobDir
  Pdb[["scriptsDir"]] <- scriptsDir
  Pdb[["mode"]] <- Modes
  
  #Creates JobDir
  system(paste("mkdir ", Pdb$JobDir, sep=""))
  system(paste("cp ", PdbFile, " ", Pdb$JobDir, sep=""))
  
  # Save equivalences
  Pdb <- pdb_equivalences(Pdb, Pdb$JobDir, Pdb$PdbBase)
  complete_backbone(Pdb, scriptsDir, Pdb$JobDir, Pdb$PdbBase)
  
  # PdbBase <- basename.pdb(PdbFile)
  
  print("Preparing files..")
  #Prepare the PDB file to get awsem input files, create the workdir and move neccessary files to it.
  system(paste("cd ", Pdb$JobDir, "; pwd; ", scriptsDir, "AWSEMFiles/AWSEMTools/PdbCoords2Lammps.sh ", Pdb$PdbBase, " ", Pdb$PdbBase, " ", scriptsDir, sep=""))
  system(paste("cp ", scriptsDir, "AWSEMFiles/*.dat* ", Pdb$JobDir, sep=""))
  
  print("Setting options...")
  #Modify the .in file to run a single step - Modify the fix_backbone file to change the mode and set options
  
  system(paste("cd ", JobDir, "; sed -i 's/run		10000/run		0/g' ", Pdb$PdbBase, ".in; sed -i 's/mutational/", Pdb$mode, "/g' fix_backbone_coeff.data",  sep=""))
  
  if(!is.null(Electrostatics_K))
  {
    print("Setting electrostatics...")
    system(paste("cd ", Pdb$JobDir, "; sed -i 's/\\[DebyeHuckel\\]-/\\[DebyeHuckel\\]/g' fix_backbone_coeff.data; sed -i 's/4.15 4.15 4.15/", Electrostatics_K, " ", Electrostatics_K, " ", Electrostatics_K, "/g' fix_backbone_coeff.data;", sep=""))
    print("Setting electrostatics...")
    system(paste("cd ", Pdb$JobDir, "; python ", scriptsDir, "Pdb2Gro.py ", Pdb$PdbBase, ".pdb ", Pdb$PdbBase, ".pdb.gro; perl ", scriptsDir, "GenerateChargeFile.pl ", Pdb$PdbBase, ".pdb.gro > ", JobDir, "charge_on_residues.dat", sep=""))
  }
  
  print("calculating...")
  system(paste("cp ", scriptsDir, "lmp_serial_", seqdist, " ", Pdb$JobDir, "; cd ", Pdb$JobDir, "; ./lmp_serial_", seqdist, " < ", Pdb$PdbBase, ".in", sep=""))

  if(Pdb$mode == "configurational" | Pdb$mode == "mutational")
  {
    system(paste("perl ", scriptsDir, "5Adens.pl ", Pdb$PdbBase, ".pdb ", gsub(".$", "", Pdb$JobDir), " ", Pdb$mode, sep=""))
  }
  system(paste("perl ", scriptsDir, "RenumFiles.pl ", Pdb$PdbBase, " ", Pdb$JobDir, " ", Pdb$mode, sep="" ))
  
  system(paste("perl ", scriptsDir, "GenerateVisualizations.pl ", Pdb$PdbBase, "_", Pdb$mode, ".pdb_auxiliar ", Pdb$PdbBase, " ", gsub(".$", "", Pdb$JobDir), " ", Pdb$mode, sep=""))
  system(paste("cp ", scriptsDir, "draw_links.py ", Pdb$JobDir, sep=""))
  
  return(Pdb)
}

plot_5Andens <- function(Pdb, chain=NULL)
{
  JobID=Pdb$PdbBase;
  Dir=Pdb$JobDir;
  mode=Pdb$mode;
  
  AdensTable=read.table(file=paste(Dir,"/", JobID,".pdb_", Pdb$mode, "_5adens", sep=""))
  
  
  Positions=as.numeric(AdensTable[,1])
  Chains=AdensTable[,2]
  MaximallyFrst= as.numeric(AdensTable[,4])
  NeutrallyFrst=as.numeric(AdensTable[,5])
  MinimallyFrst=as.numeric(AdensTable[,6])
  Total=as.numeric(AdensTable[,3])
  
  PositionsTotal=seq(from=1, to=length(Positions), by=1)
  
  if(is.null(chain))
  {
    #Primer Grafico
    # png(filename = paste(Dir, "/", JobID, "_", mode, ".png_5Adens", ".png", sep=""), width = 540, height = 420)
    Maximum=max(c(max(MaximallyFrst),max(MinimallyFrst),max(NeutrallyFrst), max(Total) ))
    par(mar=c(5,5,3,2), xpd=TRUE)
    plot(PositionsTotal, Total, type ="l", col="black", lwd=1, ylab="Local frustration density (5A sphere)", xlab="Position", ylim=(c(0, Maximum)), cex.axis=1.5, cex.lab=1.5)
    lines(PositionsTotal, NeutrallyFrst, type ="l", col="gray", lwd=2, ylab="Local frustration density (5A sphere)", xlab="Position")
    lines(PositionsTotal, MaximallyFrst, type ="l", col="red", lwd=2, ylab="Local frustration density (5A sphere)", xlab="Position")
    lines(PositionsTotal, MinimallyFrst, type ="l", col="green", lwd=2, ylab="Local frustration density (5A sphere)", xlab="Position")
    legend(x="top",  inset=c(-0.5,-0.12), legend=c("minimally frustrated", "neutral", "highly frustrated", "total"), pch=c(15, 15, 15, 15), col=c("green", "gray", "red", "black"), horiz = T, x.intersp=0.5, xjust=0,yjust=0, bty = "n")
    box(lwd=2)
  }else{
    # png(filename = paste(Dir, "/", JobID, "_", mode, "_5Adens_chain", i, ".png", sep=""), width = 540, height = 420)
    Maximum=max(c(max(MaximallyFrst[which(Chains==chain)]),max(MinimallyFrst[which(Chains==chain)]),max(NeutrallyFrst[which(Chains==chain)]), max(Total[which(Chains==chain)]) ))
    par(mar=c(5,5,3,2), xpd=TRUE)
    plot(Positions[which(Chains==chain)], Total[which(Chains==chain)], type ="l", col="black", lwd=1, ylab="Local frustration density (5A sphere)", xlab="Position", ylim=(c(0, Maximum)), cex.axis=1.5, cex.lab=1.5)
    lines(Positions[which(Chains==chain)], NeutrallyFrst[which(Chains==chain)], type ="l", col="gray", lwd=2, ylab="Local frustration density (5A sphere)", xlab="Position")
    lines(Positions[which(Chains==chain)], MaximallyFrst[which(Chains==chain)], type ="l", col="red", lwd=2, ylab="Local frustration density (5A sphere)", xlab="Position")
    lines(Positions[which(Chains==chain)], MinimallyFrst[which(Chains==chain)], type ="l", col="green", lwd=2, ylab="Local frustration density (5A sphere)", xlab="Position")
    legend(x="top",  inset=c(-0.5,-0.12), legend=c("minimally frustrated", "neutral", "highly frustrated", "total"), pch=c(15, 15, 15, 15), col=c("green", "gray", "red", "black"), horiz = T, x.intersp=0.5, xjust=0,yjust=0, bty = "n")
    box(lwd=2)
    # dev.off() 
  }
}

plot_5Adens_proportions <- function(Pdb, chain=NULL)
{
  JobID=Pdb$PdbBase;
  Dir=Pdb$JobDir;
  mode=Pdb$mode;
  
  AdensTable=read.table(file=paste(Dir,"/", JobID,".pdb_", Pdb$mode, "_5adens", sep=""))
  
  Positions=as.numeric(AdensTable[,1])
  Chains=AdensTable[,2]
  MaximallyFrst= as.numeric(AdensTable[,4])
  NeutrallyFrst=as.numeric(AdensTable[,5])
  MinimallyFrst=as.numeric(AdensTable[,6])
  Total=as.numeric(AdensTable[,3])
  
  PositionsTotal=seq(from=1, to=length(Positions), by=1)
  
  MinimallyFrst=as.numeric(AdensTable[,9])
  NeutrallyFrst=as.numeric(AdensTable[,8])
  MaximallyFrst=as.numeric(AdensTable[,7])
  
  if(is.null(chain))
  {
    #Segundo Grafico total
    # png(filename = paste(Dir, "/", JobID, "_", mode, "_5Adens_around.png", sep=""), width = 540, height = 420)
    par(mar=c(5,5,3,2), xpd=TRUE)
    par(mgp=c(2,0.5,0))
    barplot(rbind(MinimallyFrst, NeutrallyFrst, MaximallyFrst), col=c("green", "gray", "red"), axis.lty=1, xlab="Position", ylab="Density arround 5A sphere (%)", cex.lab=1.5, border = NA, space = c(0), xaxt = "n")
    AuxAxis=seq(from=0, to=length(PositionsTotal))
    axis(side = 1, at = AuxAxis , labels = (AuxAxis+min(PositionsTotal)),  tick = FALSE)
    box(lwd=2)
    legend(x="top",inset=c(-0.5,-0.12), legend=c("minimally frustrated", "neutral", "highly frustrated"), pch=c(15, 15, 15), col=c("green", "gray", "red"), horiz = T, bty = "n")
    # dev.off()
  }
  else{
    #Segundo grafico por chain
    # png(filename = paste(Dir, "/", JobID, "_", mode, "_5Adens_around_chain", chain, ".png", sep=""), width = 540, height = 420)
    par(mar=c(5,5,3,2), xpd=TRUE)
    par(mgp=c(2,0.5,0))
    barplot(rbind(MinimallyFrst[which(Chains==chain)], NeutrallyFrst[which(Chains==chain)], MaximallyFrst[which(Chains==chain)]), col=c("green", "gray", "red"), axis.lty=1, xlab="Position", ylab="Density arround 5A sphere (%)", cex.lab=1.5, border = NA, space = c(0), xaxt = "n")
    AuxAxis=seq(from=0, to=length(Positions[which(Chains==chain)]))
    axis(side = 1, at = AuxAxis , labels = (AuxAxis+min(Positions[which(Chains==chain)])),  tick = FALSE)
    box(lwd=2)
    legend(x="top",inset=c(-0.5,-0.12), legend=c("minimally frustrated", "neutral", "highly frustrated"), pch=c(15, 15, 15), col=c("green", "gray", "red"), horiz = T, bty = "n")
    # dev.off()
  }
}

plot_contact_map <-function(Pdb, chain=NULL)
{
  JobID=Pdb$PdbBase;
  Dir=Pdb$JobDir;
  mode=Pdb$mode;
  
  AdensTable=read.table(file=paste(Dir,"/", JobID,".pdb_", Pdb$mode, "_5adens", sep=""))
  
  Positions=as.numeric(AdensTable[,1])
  Chains=AdensTable[,2]
  MaximallyFrst= as.numeric(AdensTable[,4])
  NeutrallyFrst=as.numeric(AdensTable[,5])
  MinimallyFrst=as.numeric(AdensTable[,6])
  Total=as.numeric(AdensTable[,3])
  
  PositionsTotal=seq(from=1, to=length(Positions), by=1)
  
  datos<-read.table(file=paste(Dir,"/", JobID, ".pdb_", Pdb$mode ,sep=""),stringsAsFactors = F)
  
  chains<-sort(unique(c(datos$V3,datos$V4)))
  positions<-matrix(ncol=3,nrow=length(chains))
  auxPosVec<-c()
  for(i in seq_along(chains)){
    positions[i,1:2]<-range(c(datos$V1[which(datos$V3==chains[i])],datos$V2[which(datos$V4==chains[i])]))
    positions[i,3]<-positions[i,2]-positions[i,1]+1
    auxPosVec<-c(auxPosVec,positions[i,1]:positions[i,2])
  }
  
  #Renumero posiciones
  datos$pos1 <- NA
  datos$pos2 <- NA
  for(i in seq_along(chains)){
    if(i==1){bias <- 0}else{bias <- sum(positions[1:(i-1),3])}
    idx <- which(datos$V3==chains[i] )
    datos$pos1[idx] <- datos$V1[idx]-positions[i,1]+bias+1
    idx <- which(datos$V4==chains[i])
    datos$pos2[idx] <- datos$V2[idx]-positions[i,1]+bias+1
  }
  
  posNEW<-matrix(ncol=3,nrow=length(chains))
  for(i in seq_along(chains)){
    posNEW[i,1:2]<-range(c(datos$pos1[which(datos$V3==chains[i])],datos$pos2[which(datos$V4==chains[i])]))
    posNEW[i,3]<-posNEW[i,2]-posNEW[i,1]+1
  }
  
  total.positions<-sum(apply(positions,1,function(x){x[2]-x[1]+1}))
  matrz <- matrix(NA,ncol=total.positions,nrow=total.positions)
  for(i in 1:nrow(datos)){
    matrz[datos$pos1[i],datos$pos2[i]]<-datos$V12[i]
    matrz[datos$pos2[i],datos$pos1[i]]<-datos$V12[i]
    
  }
  # png(filename = paste(Dir, "/", JobID, "_", mode, "_map.png", sep=""),width = 600,height = 540)
  cotaM <- 4
  n.breaks=50
  vecvalores <- seq(-cotaM,cotaM,length.out = 2*n.breaks)
  # vecvalores <- c(vecvalores[1:n.breaks],0,vecvalores[(n.breaks+1):(2*n.breaks)])
  colores.paleta<-colorRampPalette(c("red","grey","green"))
  colores<-colores.paleta(length(vecvalores)-1)
  par(mar=c(5,5,5,5),cex.lab=1.5,cex.axis=1.5)
  layout(mat = matrix(1:2,ncol=2),widths = c(5,1))
  image(x=1:total.positions,y=1:total.positions,matrz,axes=F,xlab="Residue i",ylab="Residue j",col=colores,asp=1,breaks = vecvalores)
  
  if(length(chains)>1){
    abline(v=cumsum(posNEW[,3])+.5,lty=3,col="darkgrey",lwd=1.5)
    abline(h=cumsum(posNEW[,3])+.5,lty=3,col="darkgrey",lwd=1.5)
    axis(side=3,at=apply(posNEW[,1:2],1,mean),tck=0,las=1,lwd=2,labels = chains,tick = 0)
    axis(side=4,at=apply(posNEW[,1:2],1,mean),tck=0,las=1,lwd=2,labels = chains,tick = 0)
    mtext(text = "Chain",line=3,side=3,cex=1.5)
    mtext(text = "Chain",line=3,side=4,cex=1.5)
  }
  box(lwd=2)
  
  ticksPos <- axTicks(1)[which(axTicks(1)!=0)]
  axis(side=1,at=ticksPos,labels = auxPosVec[ticksPos],tcl=.5,lwd=2)
  axis(side=2,at=ticksPos,labels = auxPosVec[ticksPos],tcl=.5,lwd=2)
  axis(side=3,at=ticksPos,labels = NA,tcl=.5,lwd=2)
  axis(side=4,at=ticksPos,labels = NA,tcl=.5,lwd=2)
  
  par(mar=c(5,1,5,5))
  image(y=vecvalores,z=matrix(vecvalores,nrow=1),col=colores,axes=F,ylab="")
  axis(side=4,at=seq(-3,3,by=1),tcl=0.3,las=1,lwd=2)
  axis(side=4,at=c(-4,4),tcl=0,las=1,lwd=2)
  mtext(text = "Local Configurational Frustration Index",line=3,side=4,cex=1.5)
  axis(side=2,at=seq(-3,3,by=1),tcl=0.3,labels = NA,lwd=2)
  
  box(lwd=2)
  
  # dev.off()
}

view_frustration_pymol <- function(pdb)
{
  system(paste("cd ", Pdb$JobDir,  "; pymol ", Pdb$JobDir, Pdb$PdbBase, ".pdb_", Pdb$mode, ".pml", sep=""))
}

############# Input ####################################################
################################################

PdbFile="/home/gonzalo/Desktop/frustratometer2-master/1n0r.pdb"

ResultsDir="/home/gonzalo/Desktop/"
Modes="configurational"
seqdist=12
scriptsDir="/home/gonzalo/Desktop/frustratometer2-master/Scripts/"
Electrostatics_K=3.1

########################################################################
# Calculate Frustration
Pdb=calculate_frustration(PdbFile=PdbFile, Modes = "configurational", Electrostatics_K = NULL, scriptsDir = scriptsDir, ResultsDir = ResultsDir, seqdist=12)

plot_5Andens(Pdb, chain=NULL)
plot_5Adens_proportions(Pdb, chain=NULL)
plot_contact_map(Pdb, chain="A")

view_frustration_pymol(Pdb)