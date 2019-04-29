read_ligand_pdbqt <- function(route){
  tmp <- readLines(route)
  library(stringr)
  retain <- str_detect(tmp, "ATOM") #judge whether there is "ATOM" in string
  tmp <- tmp[retain]
  tmp <- data.frame(tmp, stringsAsFactors = F )
  extract_info_line <- function(x){
    t <- unlist(str_split(x, " "))
    t <- t[-which(t=="")]
    return(t)
  }
  tmp_df <- t(apply(tmp, 1, extract_info_line))
  tmp_df <- data.frame(tmp_df, stringsAsFactors = F )
  tmp_df[, c(2, 6:11)] <- sapply(tmp_df[, c(2, 6:11)], as.numeric)
  colnames(tmp_df) <- c("atom", "atom_id", "atom_name", "","", "x", "y", "z", "", "", "partial_charge", "atom_type")
  return(tmp_df)
}

ligand <- read_ligand_pdbqt("file:///D:/jupyter notebook/code/ligand_ref.pdbqt")
#save(ligand, file = "ligand.Rdata")

read_receptor_pdbqt <- function(route){
  tmp <- readLines(route)
  library(stringr)
  retain <- str_detect(tmp, "ATOM") #judge whether there is "ATOM" in string
  tmp <- tmp[retain]
  tmp <- data.frame(tmp, stringsAsFactors = F )
  extract_info_line <- function(x){
    t <- unlist(str_split(x, " "))
    t <- t[-which(t=="")]
    return(t)
  }
  tmp_df <- t(apply(tmp, 1, extract_info_line))
  tmp_df <- data.frame(tmp_df, stringsAsFactors = F )
  tmp_df[, c(2, 6:11)] <- sapply(tmp_df[, c(2, 6:11)], as.numeric)
  colnames(tmp_df) <- c("atom", "atom_id", "atom_name", "res_name","", "x", "y", "z", "", "", "partial_charge", "atom_type")
  return(tmp_df)
}

receptor <- read_receptor_pdbqt("file:///D:/jupyter notebook/code/receptor.pdbqt")
#save(receptor, file = "receptor.Rdata")
cal_distance <- function(atom1, atom2){
  d <- sqrt( (atom1$x - atom2$x)^2 + 
             (atom1$y - atom2$y)^2 + 
             (atom1$z - atom2$z)^2 )
  return(d)
}

#Define box parameters
CUBE_LENGTH = 24
GRID_SIZE = 0.5

##Define information channels:
RECEPTOR_ATOM_TYPE <- c('C', # Non H-bonding Aliphatic Carbon
                     'A', # Non H-bonding Aromatic Carbon
                     'N', # Non H-bonding Nitrogen
                     'NA',# Acceptor 1 H-bond Nitrogen
                     'NS',# Acceptor S Spherical Nitrogen
                     'OA',# Acceptor 2 H-bonds Oxygen
                     'OS',# Acceptor S Spherical Oxygen
                     'CA',# Non H-bonding Calcium
                     'FE',# Non H-bonding Iron
                     'MG',# Non H-bonding Magnesium
                     'ZN',# Non H-bonding Zinc
                     'P' ,# Non H-bonding Phosphorus
                     'S' ,# Non H-bonding Sulphur
                     'SA',# Acceptor 2 H-bonds Sulphur
                     'HD'
                     )
LIGAND_ATOM_TYPE <- c('C', # Non H-bonding Aliphatic Carbon
                      'A', # Non H-bonding Aromatic Carbon
                      'N', # Non H-bonding Nitrogen
                      'NA',# Acceptor 1 H-bond Nitrogen
                      'NS',# Acceptor S Spherical Nitrogen
                      'OA',# Acceptor 2 H-bonds Oxygen
                      'OS',# Acceptor S Spherical Oxygen
                      'F', # Non H-bonding Fluorine
                      'Cl',# Non H-bonding Chlorine
                      'Br',# Non H-bonding Bromine
                      'I', # Non H-bonding Iodine
                      'P' ,# Non H-bonding Phosphorus
                      'S' ,# Non H-bonding Sulphur
                      'SA',# Acceptor 2 H-bonds Sulphur
                      'HD'
                    )

grid_cube_center <- function(ref_ligand_dict){
  n <- length(ref_ligand_dict[,1])
  x <- y <- z <- NULL
  ##Center of the grid cube is geometry center
  Xc <- sum(ref_ligand_dict$x)/n
  Yc <- sum(ref_ligand_dict$y)/n
  Zc <- sum(ref_ligand_dict$z)/n
  cube_center = c(Xc,Yc,Zc)
  return(cube_center)
}

in_cube <- function(atom_center, cube_center, cube_length){
  #Judge if x is in the range of cube
  if(atom_center[1] < cube_center[1] - 0.5*cube_length) return(FALSE)
  else if(atom_center[1] > cube_center[1] + 0.5*cube_length) return(FALSE)
  #Judge if y is in the range of cube
  else if(atom_center[2] < cube_center[2] - 0.5*cube_length) return(FALSE)
  else if(atom_center[2] > cube_center[2] + 0.5*cube_length) return(FALSE)
  #Judge if z is in the range of cube
  else if(atom_center[3] < cube_center[3] - 0.5*cube_length) return(FALSE)
  else if(atom_center[3] > cube_center[3] + 0.5*cube_length) return(FALSE)
  else return(TRUE)
}

feature_extract <- function(receptor_dict,ref_ligand_dict,ligand_dict){
  #Step 1. Calculating cube center
  cube_center <- grid_cube_center(ref_ligand_dict)
  #Step 2. Calculating grid dimension, m
  m <- floor(CUBE_LENGTH/GRID_SIZE)
  #Step 3. Calculating channel size, n
  nr <- length(RECEPTOR_ATOM_TYPE)
  nl <- length(LIGAND_ATOM_TYPE)
  n <- nr+nl
  #Step 4. Define the feature matrix shape=m*m*m*n
  f_m <- array(0,c(m,m,m,n))
  
  ligand_n <- 0
  #ligand_dict <- ligand[1,]
  for(i in 1:length(ligand_dict[,1])){
    tmp_line <- ligand_dict[i, ]
    atom_center <- c(tmp_line$x, tmp_line$y, tmp_line$z)
    atom_type <- tmp_line$atom_type
    if(in_cube(atom_center, cube_center, CUBE_LENGTH)){#only consider atom is inside cube
      #judge which grid the atom is
      i <- floor(floor(atom_center[1]-cube_center[1]+0.5*CUBE_LENGTH)/GRID_SIZE)
      j <- floor(floor(atom_center[2]-cube_center[2]+0.5*CUBE_LENGTH)/GRID_SIZE)
      k <- floor(floor(atom_center[3]-cube_center[3]+0.5*CUBE_LENGTH)/GRID_SIZE)
      #Judge which atom type
      if(!(atom_type %in% LIGAND_ATOM_TYPE)){
        cat(atom_type, "is not in the defined ligand atom type list\n")
        next
      }
      l <- nr + which(atom_type==LIGAND_ATOM_TYPE)
      #Assign the value in feature matrix
      f_m[i,j,k,l] <- 1
      ligand_n <- ligand_n + 1
    }
  }
  cat("Found", ligand_n, "ligand atoms\n")

  #Step 6. Annotate the receptor atom channel
  receptor_n <- 0
  for(i in 1:length(receptor_dict[,1])){
    tmp_line <- receptor_dict[i, ]
    atom_center <- c(tmp_line$x, tmp_line$y, tmp_line$z)
    atom_type <- tmp_line$atom_type
    if(in_cube(atom_center, cube_center, CUBE_LENGTH)){#only consider atom is inside cube
      #judge which grid the atom is
      i <- floor(floor(atom_center[1]-cube_center[1]+0.5*CUBE_LENGTH)/GRID_SIZE)
      j <- floor(floor(atom_center[2]-cube_center[2]+0.5*CUBE_LENGTH)/GRID_SIZE)
      k <- floor(floor(atom_center[3]-cube_center[3]+0.5*CUBE_LENGTH)/GRID_SIZE)
      #Judge which atom type
      if(!(atom_type %in% LIGAND_ATOM_TYPE)){
        cat(atom_type, "is not in the defined ligand atom type list\n")
        next
      }
      l <- which(atom_type==RECEPTOR_ATOM_TYPE)
      #Assign the value in feature matrix
      f_m[i,j,k,l] <- 1
      receptor_n <- receptor_n + 1
    }
  }
  cat("Found", receptor_n, "receptor atoms")
  return(f_m)
}

feature_extract_grid <- function(receptor_file, ref_ligand_file, ligand_file){
  receptor_dict <- read_receptor_pdbqt(receptor_file)
  ref_ligand_dict <- read_ligand_pdbqt(ref_ligand_file)
  ligand_dict <- read_ligand_pdbqt(ligand_file)
  feature_dict <- feature_extract(receptor_dict, ref_ligand_dict, ligand_dict)
  return(feature_dict)
}

receptor_file <- "file:///D:/jupyter notebook/code/receptor.pdbqt"
ref_ligand_file <- "file:///D:/jupyter notebook/code/ligand_ref.pdbqt"
ligand_file <- "file:///D:/jupyter notebook/code/ligand_ref.pdbqt"
res <- feature_extract_grid(receptor_file, ref_ligand_file, ligand_file)

