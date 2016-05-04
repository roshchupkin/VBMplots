# ui.R

source("helper.R")
library(ggvis)
library(dplyr)

check<-list()

names<-as.vector(Atlas()[,1])
codes<-as.numeric(as.vector(Atlas()[,2]))

for (i in 1:length(names)) {
  check[[paste(names[i], as.character(i),sep=" - ") ]]<-codes[i]
                            }



choose<-list()
for (i in 1:length(maps)) {choose[[maps[i]]]<-i}

shinyUI(navbarPage(
  
  title="VBM",
  
  tabPanel("Manhattan plot",
           
           sidebarLayout(
             
             sidebarPanel(
               
               numericInput("num", 
                            label = h3("Set value threshold"), 
                            value = 3),
               
               helpText("Choose Brain Regions to display on plot"),
               
               sliderInput("range", 
                           label = "Range of interest:",
                           min = 1, max = 72, value = c(1, 72)),               
               
               
               
               selectInput("select", label = h3("Select map"), 
                           choices = choose, selected = 1),
               
               
               checkboxGroupInput("checkGroup",
                                  
                                  label = h3("Delete Brain Regions"), 
                                  choices = check #, selected=5
                                  
                                  )
               
               
               
               
             ),
             
             mainPanel(
               
              ggvisOutput("plot")
             )
             
             
           )
           
  ),
  
  
  
  tabPanel("Brain Region Manhattan plot",
           
           sidebarLayout(
             
             sidebarPanel(
               
               numericInput("num2", 
                            label = h3("Set p value threshold"), 
                            value = 2),
               
               selectInput("select2", label = h3("Select Brain Region"), 
                           choices = check, selected = 1),
               
               checkboxGroupInput("checkGroup2",
                                  label = h3("Delete map"), 
                                  choices = choose #, selected=5
                                  )
               
                        ),
             
             mainPanel(
               ggvisOutput("plot2")    
                      )
             
                    )
           
  )
  
  
  ) )













