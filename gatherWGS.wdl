workflow ScatterWGS {

  Array[String] symbols
  Int radius

  scatter (symbol in symbols) {
    call doVariantWorkflow { input:symbol=symbol,radius=radius }
  }
  call gatherCSVs{input:files=doVariantWorkflow.out}
}

task doVariantWorkflow {
  String symbol
  Int radius
  command {
    R -e "library('variants'); \
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
        ans = data.frame(ans); \
        ans = ans[,1:12]; \
                write.csv(ans, 'out.csv');"               
  }
  output {
  	File out= "out.csv"
  }
  	
  runtime {
    docker: "waldronlab/bioconductor_devel"
    bootDiskSizeGb:50
    }
}

task gatherCSVs {

    Array[File] files

    command <<<
        R -e "files=strsplit('${sep= ' ' files}',split=' ')[[1]]; \
              res=lapply(files, function(x){read.csv(x,header=TRUE,row.names=1)}); \
              save(res,file='res.rda');"             
    >>>
    output { File out = "res.rda" }
    
    runtime {
    docker: "reshg/devel_variants"
    bootDiskSizeGb:50
    }
}
