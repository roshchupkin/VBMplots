

SNPs<<-c("APOE4","PICALM_rs10792832_G",  "CELF1_rs10838725_C",   "CD2AP_rs10948363_G" ,  "SORL1_rs11218343_T" ,  "EPHA1_rs11771145_G"  , "ZCWPW1_rs1476679_T"  ,
         "FERMT2_rs17125944_C" , "MEF2C_rs190982_A"  ,   "NME8_rs2718058_A" ,    "PTK2B_rs28834970_C" ,  "INPP5D_rs35349669_T" , "ABCA7_rs4147929_A" ,  
         "CR1_rs6656401_A" ,     "BIN1_rs6733839_T" ,    "CASS4_rs7274581_T",    "CLU_rs9331896_T" ,     "MS4A6A_rs983392_A" ,   "SLC24A4_rs10498633_G")


Hammer_color_74<<-rainbow(74)
Hammer_color_83<<-rainbow(83)

#################################################################################
HammerAtlas<-function(exclude=NULL) {
  
  
  if (is.null(exclude)) {  
    t<-read.table(file.path(getwd(),'/data/Hammers_mith_atlases_structure_names_r83.txt'), sep='\t')
    
    code<-factor(t[,1])
    
    name<-t[,2]
    
    result<-data.frame(name, code, Hammer_color_83)
    
    colnames(result)<-c('names','code', 'colors')
    
    return (result)
    
    
                        }
  
  if (exclude=='VBM_Hammer'){
    
    
    exclude_regions<-c(17,18,19,44,45,46,47,48,49)
    
    index<-vector()
    
    t<-read.table(file.path(getwd(),'/data/Hammers_mith_atlases_structure_names_r83.txt'), sep='\t')
    
    for (i in exclude_regions) {index<-c(index,which(t[,1]==i))}
    
    code<-factor(t[-index,1])
    
    
    name<-t[-index,2]
    
    result<-data.frame(name, code, Hammer_color_74)
    
    colnames(result)<-c('names','code', 'colors')
    
    return (result)
    
    
    
                      }
  
}
###############################################################################


Read_p_values_per_region<-function(region_code){
  
  
  
  Directory<-getwd()
  
  work_dir<-vector()
  
  table_p_values<-list()
  
  
  for (i in 1:length(SNPs)){
    
    setwd(file.path(Directory , 'data' ,SNPs[i]) )
    
    infile <- paste("data",paste(region_code,'_p_values.bin',sep=''),sep='')
    con <- file(infile, "rb")
    dim <- readBin(con, "integer", 2)
    Mat <- matrix( readBin(con, "numeric", prod(dim)), dim[1], dim[2])
    close(con)
    
    table_p_values[[ SNPs[i] ]]<-Mat[4,]   
    
  }
  
  
  setwd(Directory)
  
  brain_regions<-vector()
  
  p_values<-vector()
  
  snps_names<-vector()
  
  snps_codes<-vector()
  
  
  for (i in 1:length(SNPs)) {
    
    p_values<-c(p_values, table_p_values[[ SNPs[i] ]])
    
    snps_codes<-c(snps_codes, rep( i , length(table_p_values[[ SNPs[i] ]] ) )  )
    
    snps_names<-c(snps_names, rep( SNPs[i] , length(table_p_values[[ SNPs[i] ]] ) ) )
    
  }
  
  
  data_frame_p_values<-data.frame(p_values, snps_names )
  
  
  names(data_frame_p_values)<-c('p_values', 'SNPs')
  
  
  ##convert p values to -logp
  
  
  data_frame_p_values[['p_values']]=-log10(data_frame_p_values[['p_values']])
  
  
  ##generate a new column with a series of numbers
  data_frame_p_values[['order']]=1:nrow(data_frame_p_values)
  
  data_frame_p_values[['SNPs']]<-factor(data_frame_p_values[['SNPs']])
  
  data_frame_p_values[['codes']]<-snps_codes
  
  return (data_frame_p_values)

  
}

#################################################################################

ggvis_plot_region<-function(data_frame, threshold , delete, name) {
  
  
  if (threshold<2) {threshold=2}
  if (threshold>4) {threshold=4}
  
  df=data_frame[data_frame["p_values"]>threshold,]

  df_new=subset(df)
  
  if (length(delete)!=0){
    
    for (i in delete) {df_new=subset(df_new, codes!=i)}
    
    
                        }
  
  
  rn<-as.vector(df_new$SNPs)
  
  df_new$SNPs<-rn
  
  df_new %>%
    
    
    ggvis(x=~order, y=~p_values, fill=~SNPs) %>%
    add_legend( c("p_values", "fill"))%>%
    add_axis("x", orient = "top", ticks = 0, title = name,
             properties = axis_props(
               axis = list(stroke = "white"),
               labels = list(fontSize = 0)))%>%
    add_axis("y", title = "-log10(p)") %>%
    add_tooltip(function(x) {
      if(is.null(x)) return(NULL)
      row<-x[,-1]
      paste0(names(row), ": ", format(row), collapse = "<br />")
    }, "hover" )%>%
    layer_points()
  
}


#########################################################
Read_p_values<-function(SNP, atlas, study_name=NULL){
  
  atlas<-1
  
  Directory<-getwd()
  
  work_dir<-vector()
  codes<-vector()
  
  work_dir<-'/results/'   
  
  Hammer<-HammerAtlas(exclude='VBM_Hammer')
  
  codes<-as.numeric(as.vector(Hammer[,2]))
  regions<- as.vector(Hammer[,1]) 
  colors<-as.vector(Hammer[,3])
  
  if (study_name=='RS_VBM'){
    setwd(file.path(Directory , "data",SNPs[as.numeric(SNP)]) )
    
                            }
  
  table_p_values<-list()
  
  
  for (i in 1:length(codes) ) {      
    
    infile <- paste("data",paste(codes[i],'_p_values.bin',sep=''),sep='')     
    con <- file(infile, "rb")
    dim <- readBin(con, "integer", 2)
    Mat <- matrix( readBin(con, "numeric", prod(dim)), dim[1], dim[2])
    close(con)    
    table_p_values[[ as.character(codes[i]) ]]<-Mat[4,] 
    
                                }
  
  setwd(Directory)
  
  
  brain_regions<-vector()
  
  p_values<-vector()
  
  region_codes<-vector()
  
  region_colors<-vector()
  
  
  
  for (i in 1:length(codes)) {
    
    p_values<-c(p_values, table_p_values[[as.character(codes[i])]])
    brain_regions<-c(brain_regions, rep( regions[i] , length(table_p_values[[as.character(codes[i])]] ) )  )
    region_codes<-c(region_codes, rep( codes[i] , length(table_p_values[[as.character(codes[i])]] ) ) )
    region_colors<-c(region_colors, rep( colors[i] , length(table_p_values[[as.character(codes[i])]] ) ) )
                              }
  
  data_frame_p_values<-data.frame(-log10(p_values), factor(brain_regions),region_colors, region_codes )
  
  
  names(data_frame_p_values)<-c('p_values', 'Region_name', 'colors', 'codes')
  
  
  data_frame_p_values[['order']]<-1:nrow(data_frame_p_values)
  
  
  return (data_frame_p_values)
  
  
}


##################################################################################
ggvis_plot<-function(data_frame, min_c, max_c, threshold , delete, name, study_name) {
  
  Hammer_codes<-as.numeric(as.vector(HammerAtlas(exclude='VBM_Hammer')[,2]))
  
  
  if (threshold<2) {threshold=2}
  if (threshold>4) {threshold=4}
  
  
  df<-data_frame[data_frame["p_values"]>threshold,]

  df_new<<-subset(df, codes>=Hammer_codes[min_c] & codes<=Hammer_codes[max_c])
  
  
  if (length(delete)!=0){
    for (i in delete) {df_new=subset(df_new, codes!=i)}
                        }
  
  rn<-factor(df_new$Region_name, levels=unique(as.vector(df_new$Region_name)))
  
  df_new$Region_name<-rn

  df_new %>%
    ggvis(x=~order, y=~p_values, fill=~Region_name, key:=~order) %>%
    add_axis("x", title = "Voxels")%>%
    add_axis("x", orient = "top", ticks = 0, title = name,
             properties = axis_props(
               axis = list(stroke = "white"),
               labels = list(fontSize = 0)))%>%
    add_axis("y", title = "-log10(p)") %>%
    add_tooltip(function(x) {
      if(is.null(x)) return(NULL)
      row <- df_new[df_new$order == x$order,1:2]
      paste0(names(row), ": ", format(row), collapse = "<br />")
    }, "hover")%>%
    layer_points()
  
  
  
}






