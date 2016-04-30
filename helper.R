
library(oro.nifti)

result_path<<-paste(getwd(),'/data/results/',sep='')
data_frame_path<<-paste(getwd(),'/data_frames/',sep='')
maps<<-dir(result_path)
atlas_path<<-''


#################################################################################
Atlas<-function() {
  
    t<-read.table(file.path(data_frame_path,'Atlas.txt'), sep='\t')
    
    code<-factor(t[,1])
    
    name<-t[,2]
    
    result<-data.frame(name, code, rainbow(length(t[,1])))
    
    colnames(result)<-c('names','code', 'colors')
    
    return (result)

  
}
###############################################################################
Get_region_data<-function(name=NULL,region_code=NULL, atlas_name=NULL){
  
  if (is.null(atlas_name)) {
    Atlas_image<-readNIfTI(file.path(atlas_path))
                            }
  
  file_name<-dir(paste(result_path,name,sep = ''))
  I<-readNIfTI(file.path(result_path,name,file_name))
  d<-I[Atlas_image==region_code]
  return(d)
  
}
###############################################################################
Read_region_data<-function(region_code,data_type='p-value', atlas_name=NULL){
  

  data_frame<-list()
  data_frame[['data_values']]<-vector()
  data_frame[['names']]<-vector()
  data_frame[['codes']]<-vector()
  
  
  for (i in 1:length(maps)){
    
    d<-Get_region_data(name=maps[i],region_code=region_code ) 
    data_frame[['data_values']]<-c( data_frame[['data_values']], d )
    data_frame[['names']]<-rep( maps[i], length(d) )
    data_frame[['codes']]<-rep( i, length(d)  )
    
                            }
  data_frame[['order']]=1:nrow(data_frame)
  data_frame[['names']]<-factor(data_frame[['names']])
  
  data_frame<-data.frame(data_frame)
  
  if (data_type=='p-value')  {data_frame[['data_values']] = -log10(data_frame[['data_values']])}

  return (data_frame_data_values)
  
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
Read_data<-function(map, atlas=NULL,data_type='p-value', atlas_name=NULL){
  
  atlas<-Atlas()

  data_frame<-list()
  data_frame[['data_values']]<-vector()
  data_frame[['Region_names']]<-vector()
  data_frame[['codes']]<-vector()
  data_frame[['color']]<-vector()
  
  
  for (i in 1:length(atlas[,1]) ) {      
    
    d<-Get_region_data(name=map,region_code=atlas[i,2]) 
    data_frame[['data_values']]<-c( data_frame[['data_values']], d )
    data_frame[['Region_name']]<-rep(atlas[i,1], length(d) )
    data_frame[['codes']]<-rep( atlas[i,2], length(d)  )
    data_frame[['color']]<-rep( atlas[i,3], length(d)  )
    
                                }
  
  
  data_frame[['order']]=1:nrow(data_frame)
  data_frame[['Region_name']]<-factor(data_frame[['Region_name']])
  
  data_frame<-data.frame(data_frame)
  
  if (data_type=='p-value')  {data_frame[['data_values']] = -log10(data_frame[['data_values']])}
 
  return (data_frame)
  
  
}
##################################################################################
ggvis_plot<-function(data_frame, min_c, max_c, threshold , delete, name) {
  
  atlas_codes<-as.numeric(as.vector(Atlas()[,2]))
  
  
  if (threshold<2) {threshold=2}
  if (threshold>4) {threshold=4}
  
  
  df<-data_frame[data_frame["data_values"]>threshold,]

  df_new<<-subset(df, codes>=atlas_codes[min_c] & codes<=atlas_codes[max_c])
  
  
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






