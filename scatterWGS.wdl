workflow ScatterWGS {

  Array[String] symbols
  Int radius

  scatter (symbol in symbols) {
    call doVariantWorkflow { input:symbol=symbol,radius=radius }
  }
}

task doVariantWorkflow {
  String symbol
  Int radius
  command {
    R -e "BiocManager::install('variants', version = '3.9', update=TRUE, ask=FALSE); \
		library('variants'); \
		file <- system.file('vcf', 'NA06985_17.vcf.gz', package = 'cgdv17'); \
        genesym <- '${symbol}'; \
		geneid <- select(org.Hs.eg.db, keys=genesym, keytype='SYMBOL', \
		         columns='ENTREZID'); \
		txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene; \
		seqlevelsStyle(txdb) = 'NCBI'; \
		txdb <- keepSeqlevels(txdb, '17'); \
		txbygene = transcriptsBy(txdb, 'gene'); \
		gnrng <- unlist(range(txbygene[geneid[['ENTREZID']]]), use.names=FALSE); \
		names(gnrng) <- geneid[['SYMBOL']]; \
		param <- ScanVcfParam(which = gnrng+${radius}, info = 'DP', geno = c('GT', 'cPd')); \
		vcf <- readVcf(file, 'hg19', param); \
		ans = locateVariants(vcf, txdb, AllVariants()); \
		table(mcols(ans)[['LOCATION']]); \
        names(ans) = make.names(names(ans),unique=TRUE); \
                write.csv(as.data.frame(ans), '${symbol}.csv');"               
  }
  output {
  	File out= "${symbol}.csv"
  }
  	
  runtime {
    docker: "waldronlab/bioconductor_devel"
    bootDiskSizeGb:50
    }
}
