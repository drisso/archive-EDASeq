############################################################
# 
# author: Ludwig Geistlinger
# date: 2015-03-10 13:32:37
# 
# descr: get.gene.length.and.gc.content
# 
############################################################

#
# @input:
#   - id: one or more gene IDs (ensembl or entrez)
#   - org: organism three letter code, e.g. 'hsa' for 'Homo sapiens'
#   - mode: 1. biomart (supports all ensembl organisms, but might be a time-consuming)
#           2. org.db (based on BioC annotation, which is much faster but only 
#                    for organisms with a respective TxDb, BSgenome, and OrgDb package)
getGeneLengthAndGCContent <- function(id, org, mode=c("biomart", "org.db"))
{
    id.type <- .autoDetectGeneIdType(id[1])
    if(is.na(id.type))
        stop("Only ENTREZ or ENSEMBL gene IDs are supported")
    
    mode <- match.arg(mode)

    # (a) based on BioC annotation utilities:
    #       (0) OrgDb: map between identifiers (if necessary)
    #       (1) TxDB: get genomic coordinates of genes
    #       (2) BSgenome: get sequences of genomic coordinates
    #
    if(mode=="org.db")
    {
        txdb.pkg <- .org2pkg(org, type="TxDb")
        bsgen.pkg <- .org2pkg(org, type="BSgenome") 
        require(txdb.pkg, character.only=TRUE)
        require(bsgen.pkg, character.only=TRUE)

        txdb.spl <- unlist(strsplit(txdb.pkg, "\\."))
        txdb.id.type <- txdb.spl[length(txdb.spl)]
        if(txdb.id.type == "ensGene") txdb.id.type <- "ensembl"
        else if(txdb.id.type == "knownGene") txdb.id.type <- "entrez"
        else stop(paste("TxDb does not use ENSEMBL or ENTREZ gene IDs"))

        # (0) map ensembl <-> entrez, 
        # if given id.type is entrez, but Txdb uses ensembl (or vice versa)
        if(id.type != txdb.id.type)
        {
            orgdb.pkg <- .org2pkg(org)
            require(orgdb.pkg, character.only=TRUE)
            orgdb.pkg <- get(orgdb.pkg)
            if(id.type == "entrez") id.map <- select(orgdb.pkg, 
                keys=id, columns="ENSEMBL", keytype="ENTREZID") 
            else id.map <- select(orgdb.pkg, 
                keys=id, columns= "ENTREZID", keytype="ENSEMBL")
            id <- id.map[!duplicated(id.map[,1]), 2]
        }
        
        # (1) get genomic coordinates
        txdb.pkg <- get(txdb.pkg)
        coords <- genes(txdb.pkg, vals=list(gene_id=id))
        
        # (2) get sequences
        bsgen.pkg <- get(bsgen.pkg)
        seqs <- getSeq(bsgen.pkg, coords) 
        if(length(seqs) > 1) seqs <- seqs[match(id, names(seqs))]
        len <- width(seqs)
    }
    # (b) based on BioMart
    #
    #
    else
    {
        require(biomaRt)
        require(Biostrings)
        id.type <- paste0(id.type, ifelse(id.type=="entrez", "gene", "_gene_id"))
        
        # setting mart
        message("Connecting to BioMart ...")
        ensembl <- useMart("ensembl")
        ds <- listDatasets(ensembl)[,"dataset"]
        ds <- grep(paste0("^", org), ds, value=TRUE)
        if(length(ds) == 0) 
            stop(paste("Mart not found for:", org))
        else if(length(ds) > 1)
        {
            message("Found several marts")
            sapply(ds, function(d)
                message(paste(which(ds==d), d, sep=": ")))
            n <- readline(paste0("Choose mart (1-", length(ds),") : "))
            ds <- ds[as.integer(n)]
        }

        ensembl <- useDataset(ds, mart=ensembl)

        message( paste0( "Downloading sequence", 
            ifelse(length(id) > 1, "s", ""), " ..."))
        if(length(id) > 100) message("This may take a few minutes ...")

        seqs <- getSequence(id=id, 
            type=id.type, seqType="gene_exon_intron", mart=ensembl)
        seqs <- seqs[!duplicated(seqs[,2]),]
        seqs <- seqs[match(id, seqs[,2]),1]
        seqs <- sapply(seqs, DNAString, USE.NAMES=FALSE)
        len <- sapply(seqs, length)
    }

    gc.cont <- sapply(seqs, function(s) 
        sum(alphabetFrequency(s, as.prob=TRUE)[c("C","G")])) 

    res <- cbind(len, gc.cont)
    colnames(res) <- c("length", "gc")
    rownames(res) <- id
    return(res)
}

.getOrgIdType <- function(org)
{
    it <- "eg"
    if(org == "At") it <- "tair"
    else if(org == "Pf") it <- "plasmo"
    else if(org == "Sc") it <- "sgd"
    return(it)
}   

.supportedOrganisms <- function() sub(".db0$", "", available.db0pkgs())

.availableOrgPkgs <- function(type=c("OrgDb", "TxDb", "BSgenome"))
{
    type <- type[1]
    pkgs <- available.packages(BIOC.ANNO.URL)[, "Package"]
    
    org.string <- "^org.[A-z][a-z]+.[a-z]+.db$"
    if(type == "TxDb") 
        org.string <- "^TxDb.[A-Z][a-z]+.UCSC.[a-z]{2}[A-Za-z]*[0-9]{1,3}.[a-z]{3,5}Gene$"
    else if(type == "BSgenome") 
        org.string <- "^BSgenome.[A-Z][a-z]+.UCSC.[a-z]{2}[A-Za-z]*[0-9]{1,3}$"
    org.pkgs <- grep(org.string, pkgs, value=TRUE)
    names(org.pkgs) <- NULL 
    return(org.pkgs)
}

.org2pkg <- function(org, type=c("OrgDb", "TxDb", "BSgenome"))
{
   BIOC.ANNO.URL <- "http://bioconductor.org/packages/release/data/annotation/src/contrib"

    SPECIES <- rbind(
        c("anopheles", "Anopheles gambiae", "Ag", "aga", "anoGam", "7165"),
        c("arabidopsis", "Arabidopsis thaliana", "At", "ath", NA, "3702"),
        c("bovine", "Bos taurus", "Bt", "bta", "bosTau", "9913"),
        c("canine", "Canis familiaris", "Cf", "cfa", "canFam", "9615"),
        c("chicken", "Gallus gallus", "Gg", "gga", "galGal", "9031"), 
        c("chimp", "Pan troglodytes", "Pt", "ptr", "PanTro", "9598"),
        c("ecoliK12", "Escherichia coli K12", "EcK12", "eco", NA, "562,83333,511145"), 
        c("ecoliSakai", "Escherichia coli Sakai", "EcSakai", "ecs", NA, "83334"),
        c("fly", "Drosophila melanogaster", "Dm", "dme", "dm", "7227"),
        c("human", "Homo sapiens", "Hs", "hsa", "hg", "9606"),
        c("malaria", "Plasmodium falciparum", "Pf", "pfa", NA, "5833"),
        c("mouse", "Mus musculus", "Mm", "mmu", "mm", "10090"),
        c("pig", "Sus scrofa", "Ss", "ssc", "susScr", "9823"),
        c("rat", "Rattus norvegicus", "Rn", "rno", "rn", "10116"), 
        c("rhesus", "Macaca mulatta", "Mmu", "mcc", "rheMac", "9544"),  
        c("worm", "Caenorhabditis elegans", "Ce", "cel", "ce", "6239"),
        c("xenopus", "Xenopus laevis", "Xl", "xla", "NA", "8355"),
        c("yeast", "Saccharomyces cerevisiae", "Sc", "sce", "sacCer", "4932,559292"),
        c("zebrafish", "Danio rerio", "Dr", "dre", "danRer", "7955")
    )
    colnames(SPECIES) <- c("common", "tax", "bioc", "kegg", "ucsc", "ncbi")

    type <- match.arg(type)
    ind <- apply(SPECIES, 1, function(r) org %in% r)
    if(any(ind)) i <- which(ind)[1]
    else stop(paste0("unrecognized organism ID \'", org, "\'"))
    bioc.id <- SPECIES[i, "bioc"]
    ucsc.id <- SPECIES[i, "ucsc"]

    if(type %in% c("TxDb", "BSgenome"))
    {
        pkg.string <- paste0("^", type, ".", bioc.id, "[a-z]+.UCSC.", ucsc.id)
        pkg <- grep(pkg.string, .availableOrgPkgs(type), value=TRUE)
        if(length(pkg) == 0)
            stop(paste("No corresponding", type, "package for", org))
        else if(length(pkg) > 1)
        {
            message("Found several genome assemblies")
            sapply(pkg, function(p) 
                message(paste(which(pkg==p), p, sep=": ")))
            n <- readline(paste0("Choose assembly (1-", length(pkg),") : "))
            pkg <- pkg[as.integer(n)]

            #message("Found several genome assemblies")
            #message(paste("Using latest:", pkg))
            #ver <- sapply(pkg, 
            #    function(p)
            #    {
            #        spl <- unlist(strsplit(p, "\\."))
            #        ind <- length(spl)
            #        if(type == "TxDb") ind <- ind - 1
            #        ass <- spl[ind]
            #        ver <- sub("^[a-zA-Z]+", "", ass)
            #        return(as.integer(ver))
            #    })
            #pkg <- pkg[which.max(ver)]
        }
    }
    else
    {
        id.type <- .getOrgIdType(bioc.id)
        pkg <- paste("org", bioc.id, id.type, "db", sep=".")
    }
    return(pkg)
}

.autoDetectGeneIdType <- function(id)
{
    type <- NA
    if(length(grep("^ENS[A-Z]{0,3}G[0-9]+", id))) type <- "ensembl"
    else if(length(grep("^[0-9]+$", id))) type <- "entrez"
    return(type)
}

