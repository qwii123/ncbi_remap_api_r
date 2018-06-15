main = function(
  mode=c("batches","asm-asm","asm-rsg","rsg-asm"),
  verbose=c(0,1),
  data_dir="",
  data_name="",
  out_mode=c("download","pass"),
  out_dir="",
  out_name="",
  from_acc="",
  dest_acc="",
  min_cov=0.5,
  max_exp=2.0,
  allow_dupes=c("on","off"),
  merge=c("on","off"),
  in_format=c("guess"),
  out_format=c("guess")
  ){
  # Libraries ----
  require(httr)
  require(jsonlite)
  require(rvest)

  # Variables ----
  URL.BASE = "https://www.ncbi.nlm.nih.gov:443/projects/genome/tools/remap/"
  CGI.SERVICE = "remapservice.cgi"
  CGI.INFO = "remapinfo.cgi"
  URL.CHECK = "https://www.ncbi.nlm.nih.gov/genome/tools/remap"
  URL.INFO = paste0(URL.BASE,CGI.INFO)
  URL.SERVICE = paste0(URL.BASE,CGI.SERVICE)
  verbose = as.integer(verbose)

  # Acquired variables ----
  FILE.DATA = file.path(data_dir,data_name)
  FILE.OUT = file.path(out_dir,out_name)
  if(FILE.OUT=="/"){
    FILE.OUT = file.path(data_dir,paste0("remapped_",data_name))
  }
  ID.BATCH = GetAsmBatchId(URL.INFO,from_acc,dest_acc)

  # Print Batches ----
  if(mode == "batches"){
    return(PrintBatches(URL.INFO,verbose))
  }

  # Main Submission ----
  MainSubmission(mode,verbose,URL.SERVICE,URL.CHECK,FILE.DATA,out_mode,FILE.OUT,ID.BATCH,from_acc,dest_acc,min_cov,max_exp,allow_dupes,merge,in_format,out_format)
}

MainSubmission = function(mode,verbose,URL.SERVICE,URL.CHECK,FILE.DATA,out_mode,FILE.OUT,ID.BATCH,acc_one,acc_two,min_cov,max_exp,allow_dupes,merge,in_format,out_format){
  # Variables ----
  # Setting from_arg and dest_arg according to mode
  from_arg = "source-assembly"
  dest_arg = "target-assembly"

  if((mode == "rsg-asm") | (mode == "asm-rsg")){
    from_arg = "datamapfrom"
    dest_arg = "datamapto"
  }

  # Miscellaneous variables
  with_refseqgenes = "on"
  without_refseqgenes = "on"
  gbenchout = "off"

  # API Request ----
  request = httr::POST(URL.SERVICE,
    add_headers("Content-Type" = "multipart/form-data"),
    body=list(
      "mode" = mode,
      "batch" = ID.BATCH,
      from_arg = acc_one,
      dest_arg = acc_two,
      "allow_dupes" = allow_dupes,
      "merge" = merge,
      "min_cov" = min_cov,
      "max_exp" = max_exp,
      "in_format" = in_format,
      "out_format" = out_format,
      "with_refseqgenes" = with_refseqgenes,
      "without_refseqgenes" = without_refseqgenes,
      "gbenchout" = gbenchout,
      "api" = "true",
      annot = upload_file(FILE.DATA)))

  # Response ----
  JSID = strsplit(request[["url"]],split="/")[[1]][4] # API returns from incorrect URL.
  check.URL = paste0(URL.CHECK,JSID)
  check.request = httr::GET(check.URL)
  message(paste0("Checking results at ",check.URL))
  if(!(status_code(check.request)==200)){
    stop("Error: Check your parameters.")
  }

  # RecheckSubmission ----
  download.list = RecheckSubmission(URL.CHECK)

  # out_mode "download" outputs downloads
  if(out_mode=="download"){
    download.file(url=download.list[3],destfile=FILE.OUT)
    download.file(url=download.list[2],destfile=FILE.OUT)
  }
  # out_mode "pass" outputs data as object
  if(out_mode=="pass"){
    return(download.list)
  }
}

RecheckSubmission = function(URL.CHECK){
  repeat{
    xml.doc = read_html(URL.CHECK)
    xml.title = xml.doc %>% html_nodes("title") %>% html_text()
    if(xml.title=="NCBI Remap: Results"){
      break
    }else{
      message("Refreshing...")
      Sys.sleep(3)
    }
  }
  xml.download.list = xml.doc %>% html_nodes("a.downloadExcel") %>% rvest::html_attr(name="href")
  return(xml.download.list)
}

PrintBatches = function(URL.INFO,verbose){
  if(verbose == 1){
    print(paste0("URL: ",URL.INFO,"\n"))
  }

  request = httr::GET(url=URL.INFO)

  if(verbose == 1){
    message("Batches Verbose:")
    message(paste0("header: ",header(request)))
    message(paste0("content: ",content(request)))
  }
  if(!((status_code(request)>=200)&(status_code(request)<300))){
    stop("Error: Could not retrieve alignment batch list.")
  }

  RMP.Modes = httr::content(request)
  RMP.Taxlist = RMP.Modes[[1]][[1]]
  RMP.Assemblylist = RMP.Modes[[2]][[1]]
  RMP.Batchlist = RMP.Modes[[3]][[1]]
  RMP.Refseqgenelist = RMP.Modes[[4]][[1]]
  RMP.Altlocilist = RMP.Modes[[5]][[1]]

  tax.hr = integer()
  for(tax in RMP.Taxlist){
    tax.id = tax["taxid"]
    tax.hr = unlist(c(tax.hr,tax.id))
  }

  asm.hr = character()
  for(asm in RMP.Assemblylist){
    asm.id = asm["accession"]
    asm.hr = unlist(c(asm.hr,asm.id))
  }


  df.batch = data.frame(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA)
  colnames(df.batch) = c("#batch_id","query_species","query_name","query_ucsc","query_acc","target_species","target_name","target_ucsc","target_acc","alignment_date")
  for(batch in RMP.Batchlist){
    batch.id = batch["batch_id"][[1]]
    query.acc = batch["query_acc"][[1]]
    target.acc = batch["target_acc"][[1]]
    datetime = batch["create_date"][[1]]
    date = strsplit(datetime,split=" ")[[1]][1]
    dtime = strsplit(datetime,split=" ")[[1]][2]

    query.idx = which(asm.hr == query.acc)
    target.idx = which(asm.hr == target.acc)
    query = RMP.Assemblylist[query.idx][[1]]
    target = RMP.Assemblylist[target.idx][[1]]

    qtax.idx = which(tax.hr == query["taxid"])
    ttax.idx = which(tax.hr == target["taxid"])
    qtax = RMP.Taxlist[qtax.idx][[1]]
    ttax = RMP.Taxlist[ttax.idx][[1]]

    df.loop = data.frame(batch.id,qtax["species"],query["name"],query["ucsc"],query["accession"],ttax["species"],target["name"],target["ucsc"],target["accession"],date)
    colnames(df.loop) = c("#batch_id","query_species","query_name","query_ucsc","query_acc","target_species","target_name","target_ucsc","target_acc","alignment_date")
    df.batch = rbind(df.batch,df.loop)
  }
  df.batch = df.batch[-1,]
  rownames(df.batch) = 1:nrow(df.batch)

  df.refseqgene = data.frame(NA,NA,NA,NA,NA,NA,NA,NA,NA)
  colnames(df.refseqgene) = c("#batch_id","species","name","ucsc","accession","species","RefSeqGene","","RefSeqGene")
  for(assem in RMP.Refseqgenelist){
    tax = RMP.Taxlist[[which(tax.hr==assem["taxid"])]]
    df.loop = data.frame(0,tax["species"],assem["name"],assem["ucsc"],assem["accession"],tax["species"],"RefSeqGene","","RefSeqGene")
    colnames(df.loop) = c("#batch_id","species","name","ucsc","accession","species","RefSeqGene","","RefSeqGene")
    df.refseqgene=rbind(df.refseqgene,df.loop)
  }
  df.refseqgene = df.refseqgene[-1,]
  rownames(df.refseqgene) = 1:nrow(df.refseqgene)

  df.altloci = data.frame(NA,NA,NA,NA,NA,NA,NA,NA,NA)
  colnames(df.altloci) = c("#batch_id","species","name","ucsc","accession","species","AltLoci","","AltLoci")
  for(assem in RMP.Altlocilist){
    tax = RMP.Taxlist[[which(tax.hr==assem["taxid"])]]
    df.loop = data.frame(0,tax["species"],assem["name"],assem["ucsc"],assem["accession"],tax["species"],"AltLoci","","AltLoci")
    colnames(df.loop) = c("#batch_id","species","name","ucsc","accession","species","AltLoci","","AltLoci")
    df.altloci=rbind(df.altloci,df.loop)
  }
  df.altloci = df.altloci[-1,]
  rownames(df.altloci) = 1:nrow(df.altloci)


  return(df.batch)
}

GetAsmBatchId = function(URL.INFO,acc_one,acc_two){
  request = httr::GET(URL.INFO)

  RMP.Modes = httr::content(request)
  RMP.Taxlist = RMP.Modes[[1]][[1]]
  RMP.Assemblylist = RMP.Modes[[2]][[1]]
  RMP.Batchlist = RMP.Modes[[3]][[1]]
  RMP.Refseqgenelist = RMP.Modes[[4]][[1]]
  RMP.Altlocilist = RMP.Modes[[5]][[1]]

  for(batch in RMP.Batchlist){
    batch.id = batch["batch_id"]
    query.acc = batch["query_acc"]
    target.acc = batch["target_acc"]

    if(((query.acc == acc_one) & (target.acc == acc_two)) | ((query.acc == acc_two) & (target.acc == acc_one))){
      return(batch.id[[1]])
    }
  }
  return(0)
}
