
apply_QTWAS <- function(tissue.name, snp.col, snp.col.name, special.end, phecode_long, add_savemodel, add_GWAS, add_savepval){
  
  add_org = getwd()
  
  interval.mat = matrix(c(0.05,0.35,0.25,0.55,0.45,0.75,0.65,0.95), nrow = 4, ncol=2, byrow = T)
  interval.name = apply(interval.mat, 1, function(x){paste0("(", x[1], ",", x[2], ")")})
  
  
  # read database file
  driver <- dbDriver('SQLite')
  conn <- dbConnect(drv = driver, paste0(add_savemodel,'QTWAS_gtex_v8_', tissue.name, '.db'))
  mytable.beta <- dbReadTable(conn,"beta")
  mytable.cov_mat <- dbReadTable(conn,"cov_mat")
  dbDisconnect(conn)
  
  gene.list = unique(mytable.beta$gene)
  ngene = length(gene.list)
  
  istart = 1
  all_pval_mat = matrix(NA, ncol = 11, nrow = ngene)
  rownames(all_pval_mat) = gene.list
  colnames(all_pval_mat) = c("qr", "qr(0.05,0.35)", "qr(0.25,0.55)", "qr(0.45,0.75)", "qr(0.65,0.95)",  "n_used", "n_model", "Z_Q1", "Z_Q2", "Z_Q3", "Z_Q4")
  
  ##### start the loop
  for(igene in istart:ngene){
    gene.id = gene.list[igene]
    # get rsid and extract snps from gwas data
    beta.set = which(mytable.beta$gene == gene.id)
    cov.set = which(mytable.cov_mat$gene == gene.id)
    
    # extract gwas snp
    setwd(add_GWAS)
    rsid.set = unique(mytable.beta$rsid[beta.set])
    bs = paste0(" '$", snp.col, "==", "\"" ,snp.col.name,  "\"" )
    
    for(i in 1:length(rsid.set)){
      bs  = paste0(bs, " || $", snp.col, "==",  "\"", rsid.set[i], "\"")
    }
    bs = paste0(bs, "' ", phecode_long, ".txt > gwas_temp_",tissue.name, "_",phecode,"_", gene.id, ".txt")
    
    system2('awk', args = bs)
    gwas.table = read.table(paste0("gwas_temp_",tissue.name, "_", phecode, "_",gene.id, ".txt"), header = TRUE)
    system(paste0("rm gwas_temp_",tissue.name, "_",phecode, "_", gene.id, ".txt"))
    
    if(nrow(gwas.table) == 0){
      all_pval_mat[igene, 6:7] = c(0, length(rsid.set))
      #system(paste0("rm gwas_temp_", gene.id, ".txt"))
      
      setwd(add_savepval)
      if(is.null(select.genes)){
        write.table(all_pval_mat,quote = FALSE, file = paste0(tissue.name, "_chr", gene.chr, "_",special.end, ".txt"))
      }
      cat("-------------------------------------------------------")
      print(paste0("Gene ",  igene,"/",ngene, " DONE!"))
      cat("-------------------------------------------------------")
      setwd(add_org)
      next
    }
    
    # check gwas flipped or not
    final.set = intersect(gwas.table$SNP, rsid.set)
    gwas.table.final = gwas.table[match(final.set, gwas.table$SNP),]
    indT1 = which(gwas.table.final$A1 == T) # check if any T has been read as the logic variable TRUE
    if(length(indT1) > 0){
      gwas.table.final$A1[indT1] = "T"
    }
    indT2 = which(gwas.table.final$A2 == T)
    if(length(indT2) > 0){
      gwas.table.final$A2[indT2] = "T"
    }
    
    beta.final = (mytable.beta[beta.set, ])[which(!is.na(match(mytable.beta$rsid[beta.set], final.set)) == TRUE), ]
    covX.final = (mytable.cov_mat[cov.set, ])[which(!is.na(match(mytable.cov_mat$rsid1[cov.set], final.set)) == TRUE & !is.na(match(mytable.cov_mat$rsid2[cov.set], final.set)) == TRUE), ]
    
    # reset the covX matrix
    preset.covX = matrix(0, length(final.set), length(final.set))
    uptri.loc = upper.tri(preset.covX, diag = T)
    preset.covX[uptri.loc] = covX.final$value
    preset.covX = preset.covX + t(preset.covX)
    diag(preset.covX) = diag(preset.covX)/2
    covX = preset.covX
    sdX =sqrt(diag(covX))
    
    for(ii in 1: length(final.set)){
      gwas.ref = gwas.table.final$A1[ii]
      gwas.alt = gwas.table.final$A2[ii]
      beta.ref = beta.final$ref[match(final.set[ii], beta.final$rsid)]
      beta.alt = beta.final$alt[match(final.set[ii], beta.final$rsid)]
      
      if(gwas.ref == beta.alt & gwas.alt == beta.ref){
        gwas.table.final$Z[ii] = -1*gwas.table.final$Z[ii]
      }
      
    }
    gwas.z = gwas.table.final$Z
    
    # reconstruct beta_A and covX
    interval.index = match(beta.final$interval, interval.name)
    pval_qr = runif(4)
    zscore_qr = rep(NA, 4)
    
    for(jj in unique(interval.index)){
      beta_qr = beta.final$weight[interval.index == jj]
      rsid_temp = beta.final$rsid[interval.index == jj]
      qr_model_id = match(rsid_temp, final.set)
      if(length(qr_model_id) > 1){
        gwas_id = match(rsid_temp, gwas.table.final$SNP)
        Z = beta_qr%*%diag(sdX[qr_model_id])%*%gwas.z[gwas_id]
        Sigma_z = t(beta_qr)%*%covX[qr_model_id, qr_model_id]%*%beta_qr
        pval_qr[jj] = 2*(pnorm( abs(as.numeric(Z/sqrt(Sigma_z))) , lower.tail=F))
        zscore_qr[jj] = as.numeric(Z/sqrt(Sigma_z))
        
      }
      if(length(qr_model_id) == 1){
        gwas_id = match(rsid_temp, gwas.table.final$SNP)
        Z = beta_qr*sdX[qr_model_id]*gwas.z[gwas_id]
        Sigma_z = t(beta_qr)*covX[qr_model_id, qr_model_id]*beta_qr
        pval_qr[jj] = 2*(pnorm( abs(as.numeric(Z/sqrt(Sigma_z))) , lower.tail=F))
        zscore_qr[jj] = as.numeric(Z/sqrt(Sigma_z))
        
      }
      
    }
    p_combo =Get.cauchy(pval_qr)
    all_pval_mat[igene, ] = c(p_combo, pval_qr, length(final.set), length(rsid.set), zscore_qr)
    print(c(p_combo, pval_qr))
    
    setwd(add_savepval)
    if(is.null(select.genes)){
      write.table(all_pval_mat,quote = FALSE, file = paste0(tissue.name, "_chr", gene.chr, "_",special.end, ".txt"))
    }
    cat("-------------------------------------------------------")
    print(paste0("Gene ",   igene,"/",ngene, " DONE!"))
    cat("-------------------------------------------------------")
    
    setwd(add_org)
  }
  print(paste0(tissue.name, " is all DONE!"))
  
  return(all_pval_mat)
}


Get.cauchy<-function(p){
  p[p>0.99]<-0.99
  is.small<-(p<1e-16) & !is.na(p)
  is.regular<-(p>=1e-16) & !is.na(p)
  temp<-rep(NA,length(p))
  temp[is.small]<-1/p[is.small]/pi
  temp[is.regular]<-as.numeric(tan((0.5-p[is.regular])*pi))
  
  cct.stat<-mean(temp,na.rm=T)
  if(is.na(cct.stat)){return(NA)}
  if(cct.stat>1e+15){return((1/cct.stat)/pi)}else{
    return(1-pcauchy(cct.stat))
  }
}


