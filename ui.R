#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(shinyjs)
library(shinythemes)

# Define UI for application
shinyUI(
  navbarPage("MG",
             #shinythemes::themeSelector(),
             tabPanel("Home",
                      fluidRow(
                        column(1),
                        column(10,
                               h3("Introduction"),
                               p("Wound healing assay is a simple and cost-effective in vitro assay for assessing therapeutic impacts on cell migration.  Its major limitation is the possible confoundment by other cellular phenotypes, causing misinterpretation of the experimental outcome.  Here, we developed a novel methodology that enables the quantification of therapeutic impacts on different cellular phenotypes. Our technique requires only the concurrent measurement of cell death together with the relative wound density changes and is easily adapted for high-throughput drug screening.  By using context-dependent wound closure time, our approach can minimize influences of extrinsic noises such as cell seeding density, ensuring that more robust inference of therapeutic effect on different phenotypic outcomes."),
                               img(src='MG_home.png',style="display: block; margin-left: auto; margin-right: auto; width: 90%;", align = "center"),
                               h3("Web development"),
                               p("Somchai, P., Phongkitkarun, K. and Sampattavanich, S.")
                        ),
                        column(1)
                      )
             ),
             
             tabPanel("About MG parameters",
                      fluidPage(
                        fluidRow(
                          column(1),
                          column(10,
                                 titlePanel(h3("General Information")),
                                 tabsetPanel(
                                   tabPanel("Introduction",
                                            p(),
                                            p("Wound healing assay has been used ubiquitously for the assessment of cell migration in biomedical research but was also known to be confounded by complex influences of multiple cellular phenotypes.  We have demonstrated in this study how to modify the standard wound healing assay for accurate scoring of different phenotypic contributions.  Since the conventional wound healing assay cannot be used for scoring the contributions of different phenotypes, past efforts in identifying potential inhibitors of cell migration often misinterpret phenotypic influences by different compounds.  By concurrent monitoring of wound closure change together with cell death quantification, we showed that therapeutic effects on cellular phenotypes can be scored based on the established relationship of wound closure rate and changes of cell death at different drug concentrations.  Fractional contribution of cell death and cell migration can then be approximated from the slope of the fitted relationship."),
                                            img(src='MG_introduction.png',style="display: block; margin-left: auto; margin-right: auto; width: 90%;", align = "center"),
                                            p("Fig. 1: The established relationship between cumulative changes of caspase signal and relative wound density (RWD) changes at different drug concentrations (left). Scoring of therapeutic impacts on different phenotypics at increasing drug concentrations. (right)."),
                                            h3("Funding Sources:"),
                                            p("This work was funded by NSTDA (P-15-50208) and partially supported by Faculty of Medicine Siriraj Hospital, Mahidol University")
                                   ),
                                   tabPanel("About MG parameter",
                                            withMathJax(),
                                            h4("Calculation of cell migration/death scores:"),
                                            p("To separate the contribution of cell migration and cell death from the observed delayed in wound closure, we generated the relationship of the cumulative changes of cell death signal and the relative wound density."),
                                            p("For calculation of “wound closure change” (\\({WC}_{xi}\\)), we determined the AUC of the RWD from each drug treatment condition (\\({AUC(RWD)}_{xi}\\)) and normalized this value by that from the DMSO control group (\\({AUC(RWD)}_{ctrl}\\)), all with MMC pre-treatment."),
                                            p("$${WC}_{xi}=\\frac{({AUC(RWD)}_{ctrl}-{AUC(RWD)}_{xi})}{{AUC(RWD)}_{ctrl}}$$"),
                                            p("It is important to note that some drug cannot completely inhibit wound closure while others can inhibit wound closure almost entirely. Thus, we define another parameter called drug-independent wound closure change (\\({WC}_{ind}\\)) and can summarize the relationship between the two types of wound closure changes from the following equation:"),
                                            p("$$1={WC}_{xi}+{WC}_{ind}$$ , where"),
                                            p("$${WC}_{ind}=\\frac{{AUC(RWD)}_{xi}}{{AUC(RWD)}_{ctrl}}$$"),
                                            p("For “Death change” (\\({DC}_{xi}\\)), we determined the area under curve from caspase 3/7 activity (\\({AUC(D)}_{xi}\\)) at each drug treatment condition and normalized this value by that of the DMSO control (\\({AUC(D)}_{ctrl}\\)), both with MMC pre-treatment."),
                                            p("$${DC}_{xi}=\\frac{({AUC(D)}_{xi}-{AUC(D)}_{ctrl})}{{AUC(D)}_{max}}$$"),
                                            p("The data obtained from (\\({WC}_{xi}\\)) and (\\({DC}_{xi}\\)) were used to determine the contribution of cell death and cell migration from the drug-induced wound closure change.  We attempted to capture this relationship by fitting with 1st-order rate equation.  If the fitting quality is inadequate (i.e. \\(r^2\\lt 0.6\\)), a vertical straight line is used for fitting instead.  If the fitting with straight line is favorable, we concluded that the observed wound change is dependent only on one cellular phenotype."),
                                            p("When the relationship is explainable by 1st order rate equation, we approximated the fractional contribution of cell migration and cell death by calculating the arctan of the slope (\\(\\theta_{xi}\\) in radian).  The summation of the migration score (\\({Migration}_{xi}\\)) and death score (\\({Death}_{xi}\\)) should be equal to the drug-induced wound closure change (\\({WC}_{xi}\\)).  Mathematically, such relationship can be defined as follows:"),
                                            p("$${WC}_{xi}={Migration}_{xi}+{Death}_{xi}$$ , where"),
                                            p("$${Death}_{xi}={WC}_{xi}\\times\\left(\\frac{\\theta_{xi}}{2}\\right)$$")      
                                   ),
                                   tabPanel("How to use our tool",
                                            h4("Formatting input files"),
                                            p("Input files may be either comma-separated files (.csv). Example dataset can be seen in the \"Online MG Calculator\" section."),
                                            h4("Instructions"),
                                            p("To calculate migration score and death score using wound healing assay, the calculator requires measurement of cell death (ex. From caspase 3/7 activity dye) together with Mitomycin-C pre-treatment, before adding the desired drug.  Users must provide a data file in which each row represents a separate treatment condition and the columns specify the keys (variables) that define the treatment condition (e.g. Cell_Name, Density, Replicate, Drug_Name, Dose, Mitomycin_C, Time, Relative_Wound_Density (or RWD%) and Relative_Caspase). Interactive analysis and visualization tools are provided. Detailed instructions can be found below."),
                                            h4("Step 1:"),
                                            h5("Load the data file containing RWD% and Caspase 3/7 signal for treated and control cells. MGcalculator accepts comma-separated (.csv) input files."),
                                            p("Click the Load Example button and then click ‘Load Example’ to view a sample data file in the Data Tables tab."),
                                            h4("Step 2:"),
                                            h5("Select the variables (ex. Cell Line, Cell Density and Drug)."),
                                            p("Select key variables for establishing the relationship between ‘Wound closure change’ or WC and ‘Death change’ or DC. By default, all variables are sorted by an alphabetical."),
                                            h4("Step 3:"),
                                            h5("Select the method for calculating the phenotypic scores by clicking ‘Define analytical endpoint’"),
                                            p("User may choose different endpoint analytical methods for calculating migration and death scores.  In our application, the options include 1) Defined Time (e.g. 72hrs) 2) Percent Slope (e.g. time when the slope reach 2%) and 3) Percent of RWD (e.g. time when RWD reach  50%). ‘Percent Slope’ is the default option.")
                                   ),
                                   tabPanel("Authors",
                                            p(),
                                            p("This website was designed and built by the following members of the Siriraj Initialtive in Systems Pharmacology (SISP), Faculty of Medicine Siriraj Hospital, Mahidol University and Department of Biomedical Engineering, Faculty of Engineering, Mahidol University"),
                                            tags$li("Somchai, P."),
                                            tags$li("Phongkitkarun, K."),
                                            tags$li("Kueanjinda, P."),
                                            tags$li("Jamnongsong, S."),
                                            tags$li("Vaeteewoottacharn, K."),
                                            tags$li("Luvira, V."),
                                            tags$li("Okada, S."),
                                            tags$li("Jirawatnotai, S."),
                                            tags$li("Sampattavanich, S.")
                                   ))),
                          column(1)
                        )
                      )
             ),
             
             tabPanel("Online MG Calculator",
                      fluidPage(
                        tags$head(tags$style(HTML(".shiny-output-error-validation {color: red;}"))),
                        sidebarPanel(
                          useShinyjs(), # Define shinyjs function for show/hidden ui
                          helpText("Create migration-death relationship from would healing assay."),
                          actionButton("loadexample","Load Example"),
                          fileInput("dataFile",h5("Import data file"),
                                    accept = c(
                                      "text/csv",
                                      "text/comma-separated-values,text/plain",
                                      ".csv")),
                          selectInput("selectCell",h5("Cell Line"), 
                                      choices = "", selected = ""),
                          selectInput("selectDensity",h5("Cell Density"),
                                      choices = "", selected = ""),
                          selectInput("selectDrug",h5("Drug"), 
                                      choices = "", selected = ""),
                          radioButtons("selectTime", h5("Define Half Time"),
                                       choices = list("Method 1 : Define Time" = 1,
                                                      "Method 2 : Percent Slope" = 2,
                                                      "Method 3 : Percent RWD" = 3),selected = 2),
                          uiOutput("manualTime"),
                          uiOutput("manualSlope"),
                          uiOutput("manualRWD")
                        ),
                        # Show a plot of the generated distribution
                        mainPanel(
                          tabsetPanel(
                            tabPanel("Data Tables",
                                     tableOutput("contents")),
                            tabPanel("RWD-Time Control",
                                     plotOutput("graphRWDT", height="auto",width = "100%"),
                                     textOutput("textRWD.halftime"),
                                     textOutput("textRWD.rsquare")),
                            tabPanel("WC-DC",
                                     plotOutput("graphWCDC", height="auto",width = "100%"),
                                     textOutput("textWCDC.rsquare")),
                            tabPanel("Affected Score",
                                     plotOutput("graphAS", height="auto",width = "100%")),
                            tabPanel("Summary",
                                     tableOutput("outputs"))
                          )
                        )
                      )        
             ),
             tabPanel("Support",
                      fluidPage(
                        fluidRow(
                          column(1),
                          column(10,
                                 titlePanel(h3("Technical Support")),
                                 p(),
                                 h5("For further information, please contact sisyspharm@gmail.com")
                                 ),
                          column(1)
                        )))
  )
)

