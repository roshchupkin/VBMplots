
library(oro.nifti)
library(dplyr)
result_path<<-paste(getwd(),'/data/results/',sep='')
data_frame_path<<-paste(getwd(),'/data_frames/',sep='')
maps<<-dir(result_path)
atlas_path<<-paste(getwd(),'/data/atlas/Atlas.nii.gz',sep='')


#################################################################################
Atlas<-function() {
  
    t<-read.table(file.path(data_frame_path,'Atlas.txt'), sep='\t')
    
    code<-factor(t[,1])
    
    name<-t[,2]
    
    result<-data.frame(name, code, rainbow(length(t[,1])))
    
    colnames(result)<-c('Region_name','codes', 'color')
    
    result[['Region_name']]<-as.character(result[['Region_name']])
    result[['codes']]<-as.numeric(result[['codes']])
    
    return (result)

  
}
###############################################################################
Read_atlas_image<-function(atlas_name=NULL){
  
  if (is.null(atlas_name)) {
    Atlas_image<-readNIfTI(file.path(atlas_path))
  }
  
  return(Atlas_image)
  
}

###############################################################################
Get_data<-function(name=NULL,df_atlas=NULL, atlas_name=NULL){
  
  Atlas_image<-Read_atlas_image()
  
  data_frame<-list()
  
  file_name<-dir(paste(result_path,name,sep = ''))
  I<-readNIfTI(file.path(result_path,name,file_name))
  Atlas_image<-as.numeric(as.vector(Atlas_image))
  
  I<-I[Atlas_image!=0]
  Atlas_image<-Atlas_image[Atlas_image!=0]
  
  data_frame[['data_values']]<-I
  data_frame[['codes']]<-Atlas_image
  
  t1<-Sys.time()
  result<-merge(data_frame,df_atlas, by.y=c('codes'), by.x=c('codes') )
  print (Sys.time()-t1)

  return(result)
  
  
}
####################################################################################
Read_data<-function(map_index, atlas=NULL,data_type='p-value', atlas_name=NULL){
  
  atlas<-Atlas()
  
  
  t1<-Sys.time()
  data_frame<-Get_data(maps[as.numeric(map_index)], df_atlas = atlas )
  
  
  data_frame[['order']]=1:nrow(data_frame)
  data_frame[['Region_name']]<-factor(data_frame[['Region_name']])

  if (data_type=='p-value')  {
    index=which(data_frame[['data_values']]!=0)
    data_frame[['data_values']][index] = -log10(data_frame[['data_values']][index])
  }
  
  print(Sys.time() - t1)
  
  return (data_frame)
  
  
}
##############################################################################
Get_region_coordinates<-function(region_code=NULL, atlas_name=NULL){
  
  Atlas_image<-Read_atlas_image()
  
  index<-which(Atlas_image==region_code)
  
  return(index)  
  
}
##############################################################################
Get_region_data<-function(name=NULL,region_code=NULL, atlas_name=NULL, index=NULL){
  
  print (c(name,region_code))
  
  if (is.null(index)) {
    Atlas_image<-Read_atlas_image()
    index<-which(Atlas_image==region_code)
                      }

  file_name<-dir(paste(result_path,name,sep = ''))
  I=readNIfTI(file.path(result_path,name,file_name))
  t1<-Sys.time()
  d<-I[index]
  print (paste("Time to read one region data...",as.character(Sys.time()- t1),sep='' ))
  return(d)
  
}
###############################################################################
Read_region_data<-function(region_code, data_type='p-value', atlas_name=NULL){
  
  t1<-Sys.time()
  
  data_frame<-list()
  data_frame[['data_values']]<-vector()
  data_frame[['names']]<-vector()
  data_frame[['codes']]<-vector()
  
  index<-Get_region_coordinates(region_code = region_code)
  
  for (i in 1:length(maps)){
    
    d<-Get_region_data(name=maps[i],region_code=region_code, index = index ) 
    data_frame[['data_values']]<-c( data_frame[['data_values']], d )
    data_frame[['names']]<-c(data_frame[['names']], rep( maps[i], length(d) ) )
    data_frame[['codes']]<-c(data_frame[['codes']],rep( i, length(d)  ) )
    
  }
  data_frame<-data.frame(data_frame)
  
  data_frame[['order']]=1:nrow(data_frame)
  data_frame[['names']]<-factor(data_frame[['names']])
  
  data_frame<-data.frame(data_frame)
  
  if (data_type=='p-value')  {
    index=which(data_frame[['data_values']]!=0)
    data_frame[['data_values']][index] = -log10(data_frame[['data_values']][index])
    }
  
  print (paste("Time to read regions data...",as.character(Sys.time()- t1),sep='' ))
  return (data_frame)
  
}

#################################################################################

ggvis_plot_region<-function(data_frame, threshold , delete, name) {
  
  
  if (threshold<2) {threshold=2}
  if (threshold>4) {threshold=4}
  
  df=data_frame[data_frame["data_values"]>threshold,]

  df_new=subset(df)
  
  if (length(delete)!=0){
    for (i in delete) {df_new=subset(df_new, codes!=i)}
                        }
  
  
  rn<-as.vector(df_new$SNPs)
  df_new$SNPs<-rn

  df_new %>%
    ggvis(x=~order, y=~data_values, fill=~names) %>%
    add_legend( c("values", "fill"))%>%
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
    layer_points()%>%
    set_options(width = 1000, height = 600)
  
}
##################################################################################
ggvis_plot<-function(data_frame, min_c, max_c, threshold , delete, name) {
  
  atlas_codes<-as.numeric(as.vector(Atlas()[,2]))
  
  print(dim(data_frame))
  
  if (threshold<2) {threshold=2}
  if (threshold>4) {threshold=4}
  
  
  df<-data_frame[data_frame["data_values"]>threshold,]

  df_new<<-subset(df, codes>=atlas_codes[min_c] & codes<=atlas_codes[max_c])
  
  print(dim(df_new))
  
  if (length(delete)!=0){
    for (i in delete) {df_new=subset(df_new, codes!=i)}
                        }
  
  rn<-factor(df_new$Region_name, levels=unique(as.vector(df_new$Region_name)))
  
  df_new$Region_name<-rn

  df_new %>%
    ggvis(x=~order, y=~data_values, fill=~Region_name, key:=~order) %>%
    add_axis("x", title = "Voxels")%>%
    add_axis("x", orient = "top", ticks = 0, title = name,
             properties = axis_props(
               axis = list(stroke = "white"),
               labels = list(fontSize = 0)))%>%
    add_axis("y", title = "-log10(p)") %>%
    add_tooltip(function(x) {
      if(is.null(x)) return(NULL)
      row <- df_new[df_new$order == x$order,2:3]
      paste0(names(row), ": ", format(row), collapse = "<br />")
    }, "hover")%>%
    layer_points()%>%
    set_options(width = 1000, height = 600)
  
  
  
}






