# server.R

library(ggvis)
library(dplyr)

source("helper.R")

shinyServer(function(input, output) {
  
  
  data_frame1<-reactive( { Read_data(input$select) } )
  
  
  observe({

    roi<-input$range
    
    th<-input$num
    
    del<-input$checkGroup
    
    ggvis_plot(data_frame1(),roi[1], roi[2] , th, del,
               "Brain Manhattan plot. 
               Every dot is a single voxel from gray matter.Right corner: resize plot; mouse: tooltips" 
                )%>% bind_shiny("plot")
    
  })
  
  
  data_frame2<-reactive( { Read_region_data(input$select2 ) } )
  
  
  vis2<-reactive({
    
    th2<-input$num2
    
    del2<-input$checkGroup2
    
    ggvis_plot_region(data_frame2(),th2, del2, 
    "Brain Manhattan plot. Effect from all maps shown for the same voxels from selected brain regions." )
    
  })
  
  
  vis2 %>% bind_shiny("plot2")
  
  
  
  })