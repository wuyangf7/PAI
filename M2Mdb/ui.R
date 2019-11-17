library(shiny)
library(ggplot2)  # for the diamonds dataset

shinyUI(pageWithSidebar(
  headerPanel('Methylome to Methylome database (M2Mdb)'),
 
  sidebarPanel(
    #selectInput("Trait",
    #            label = "Choose a trait",
    #            choices = c("BMI","HDL", "UC","CAD", "Height","RA","T2D","WHRadjBMI","CD","LDL","SCZ",
    #                       "TG"),
    #           selected = "BMI"),
    helpText('This database is designed to query the promoter-anchored chromatin interactions identified in the Wu et al. (In preparation) study. '),
    selectInput("Chr",
                label = "Chromosome",
                choices = c("All","1", "2","3", "4","5","6","7","8","9","10",
                            "11", "12","13", "14","15","16","17","18","19","20",
                            "21","22"),
                selected = "All"),
    textInput("meprobeid", label = ("Exposure DNA methylation probe"), 
              value ="", placeholder="e.g. cg03123370"),
    textInput("eprobeid", label = ("Outcome DNA methylation probe"), 
              value ="", placeholder="e.g. cg00008971"),
    textInput("geneid", label = ("Gene"), 
              value ="", placeholder="e.g. CCDC163P"),
    textInput("psmr", label = ("SMR p-value threshold"), 
              value ="",placeholder="e.g. 1.57e-9"),
    textInput("phet", label = ("HEIDI p-value threshold"), 
              value ="",placeholder="e.g. 0.01"),
    submitButton("Submit"),
    tags$p(),
    tags$p(tags$strong("Citation:"),"Wu Y, Qi T, Wang H, Zhang F, Zheng Z, Deary IJ, McRae AF, Wray NR, Zeng J & Yang J (2018) Promoter-anchored chromatin interactions predicted from genetic analysis of epigenomic data. In 
preparation."),
    tags$p(tags$strong("Credits:"),tags$a(href = "mailto:y.wu2@uq.edu.au", "Yang Wu"), ", ",tags$a(href = "https://researchers.uq.edu.au/researcher/15871", "Ting Qi"),
           "and ", tags$a(href = "http://researchers.uq.edu.au/researcher/2713", "Jian Yang"), " at the Program in ",tags$a(href = "http://cnsgenomics.com/", "Complex Trait Genomics"),
           ", The University of Queensland."),
    tags$p(tags$strong("Troubleshooting:"), tags$a(href = "mailto:jian.yang@imb.uq.edu.au", "Jian Yang")),
    width = 3
  ),
  mainPanel(
    tags$head(includeScript("google-analytics.js")),
    tabsetPanel(
      tabPanel('Methylome to Methylome',
               DT::dataTableOutput("mymapping"))
      ,
      br(),
      tags$p(),
      tags$p(tags$strong(helpText('Note:')))
      ,
      helpText('Exposure DNAm: DNA methylation probe used as exposure')
      ,
      helpText('Exposure DNAm gene: nearest gene to the exposure DNAm')
      ,
      helpText('Outcome DNAm: DNA methylation probe used as outcome')
      ,
      helpText('Outcome DNAm gene: nearest gene to the outcome DNAm')
      ,
      helpText('pSMR: p-value from the SMR test')
      ,
      helpText('pHEIDI: p-value from the HEIDI test')
      
      #,
      #tabPanel('Plot',
      #         htmlOutput("text1"),
      #        plotOutput("myplot",width = "100%",height = "800px"))
      #,
      #tabPanel('About',
      #          htmlOutput("about"))
    )
  )
))
