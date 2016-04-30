# ui.R

source("helper.R")
library(ggvis)
library(dplyr)

check<-list()

names<-as.vector(HammerAtlas(exclude='VBM_Hammer')[,1])
codes<-as.numeric(as.vector(HammerAtlas(exclude='VBM_Hammer')[,2]))

for (i in 1:length(names)) {
  check[[paste(names[i], as.character(i),sep=" - ") ]]<-codes[i]
                            }




snp<-c("APOE4","PICALM_rs10792832_G",  "CELF1_rs10838725_C", 
       "CD2AP_rs10948363_G" ,  "SORL1_rs11218343_T" ,  "EPHA1_rs11771145_G"  , "ZCWPW1_rs1476679_T"  ,
       "FERMT2_rs17125944_C" , "MEF2C_rs190982_A"  ,   "NME8_rs2718058_A" ,    "PTK2B_rs28834970_C" ,  "INPP5D_rs35349669_T" , "ABCA7_rs4147929_A" ,  
       "CR1_rs6656401_A" ,     "BIN1_rs6733839_T" ,    "CASS4_rs7274581_T",    "CLU_rs9331896_T" ,     "MS4A6A_rs983392_A" ,   "SLC24A4_rs10498633_G")





choose<-list()
for (i in 1:length(snp)) {choose[[snp[i]]]<-i}

shinyUI(navbarPage(
  
  title="VBM",
  
  tabPanel("Manhattan plot",
           
           sidebarLayout(
             
             sidebarPanel(
               
               numericInput("num", 
                            label = h3("Set p value threshold"), 
                            value = 3),
               
               helpText("Choose Brain Regions to display on plot"),
               
               sliderInput("range", 
                           label = "Range of interest:",
                           min = 1, max = 72, value = c(1, 72)),               
               
               
               
               selectInput("select", label = h3("Select SNPs"), 
                           choices = choose, selected = 3),
               
               
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
                                  label = h3("Delete SNPs"), 
                                  choices = choose #, selected=5
                                  )
               
                        ),
             
             mainPanel(
               ggvisOutput("plot2")    
                      )
             
                    )
           
  )
  
  
  ) )













