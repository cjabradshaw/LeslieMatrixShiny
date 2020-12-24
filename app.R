# Leslie matrix projection R Shiny app
# Corey Bradshaw
# Flinders University

## remove everything
rm(list = ls())

# load libraries
library(shiny)
library(shinybusy)
library(ggplot2)
library(tsoutliers)

## call functions
source(file.path("./functions/", "createLmatrix.R"), local=T)
source(file.path("./functions/", "maxLambda.R"), local=T)
source(file.path("./functions/", "maxR.R"), local=T)
source(file.path("./functions/", "Gval.R"), local=T)
source(file.path("./functions/", "Rval.R"), local=T)
source(file.path("./functions/", "StableStageDist.R"), local=T)
source(file.path("./functions/", "projLmatrix.R"), local=T)
source(file.path("./functions/", "projStochMatrix.R"), local=T)
source(file.path("./functions/", "projStochPulse.R"), local=T)
source(file.path("./functions/", "projStochPress.R"), local=T)
source(file.path("./functions/", "StochMVP.R"), local=T)
source(file.path("./functions/", "LogPowFunc.R"), local=T)
source(file.path("./functions/", "predictNLS.R"), local=T)
source(file.path("./functions/", "estBetaParams.R"), local=T)
source(file.path("./functions/", "setBackgroundColor.R"), local=T)

ui <- fluidPage(
  
  # title of app
  titlePanel("Leslie-matrix projection of a population"),
  
  tags$style("@import url(https://use.fontawesome.com/releases/v5.7.2/css/all.css);"),
  
  setBackgroundColor(
    color = "#7fda5d",
    gradient = c("linear"),
    direction = c("bottom")
  ),

  wellPanel(style = "background: #d7f9da",
    tags$img(height = 150, src = "LmatrixLifeCycleGraph.png", style="float:right"),
    tags$p(style="font-family:Avenir", tags$i(class="fab fa-r-project", title="R Project"),"Shiny App by", tags$a(href="https://globalecologyflinders.com/people/#CJAB", "Corey Bradshaw "),
    tags$a(href = "mailto:corey.bradshaw@flinders.edu.au","(",tags$i(class="far fa-envelope"),"e-mail"),";",
    tags$a(href = "https://github.com/cjabradshaw", tags$i(class="fab fa-github",title="Github"),"Github)")),
    tags$h4(style="font-family:Avenir", "Preamble"),
    tags$p(style="font-family:Avenir", "This app projects a user-defined Leslie (age-classified) matrix to examine population changes through time. This fully
           customisable app includes a density-feedback function on survival relative to desired initial population size and carrying capacity,
           stochastic projections with user-defined variance for survival, a generationally scaled catastrophe function, a single 'pulse' disturbance
           function, and a 'press' disturbance function with a user-defined time application."),
    tags$p(style="font-family:Avenir", "A detailed instructions",tags$i(class="fas fa-directions"), "tab (tab I) is included for guidance, but a brief sequence description is included below.
           User-defined settings in each tab are carried over to subsequent tabs."),
    tags$ol(type="A", tags$li(tags$p(style="font-family:Avenir",tags$strong("SET-UP"),tags$i(class="fas fa-pencil-ruler"),": set matrix dimensions (longevity), age",tags$em("x"),"-specific survival
                           (",tags$em("s", tags$sub("x")),") and fertility (",tags$em("f",tags$sub("x")), ") probabilities,
                           offspring sex ratio, % variance around survival/fertility probabilities, and whether lifespan is abrupt or diffuse.")),
            tags$li(tags$p(style="font-family:Avenir",tags$strong("MATRIX PROPERTIES"),tags$i(class="fas fa-table"),": shows Leslie matrix according to settings in tab A (or for a previously defined matrix that you can upload),
                           as well as the dominant eigen value",tags$em("λ"), "instantaneous rate of population change", tags$em("r"),"generation length",
                           tags$em("G"), ", and reproduction number R0 (number of ♀ offspring/adult ♀).")),
            tags$li(tags$p(style="font-family:Avenir",tags$strong("DENSITY FEEDBACK"),tags$i(class="fas fa-exchange-alt"),": set initial population size and carrying capacity", tags$em("K"),
                           ", as well as the three coefficients (", tags$em("a"),",",tags$em("b"),",",tags$em("c"),") from a logistic power function to define the relationship between a survival modifier",
                           tags$em("S"),tags$sub("mod"),"and population size.")),
            tags$img(height = 100, src = "Amatrix.png", style="float:right"),
            tags$li(tags$p(style="font-family:Avenir",tags$strong("PROJECT"),tags$i(class="fas fa-chart-line"),": deterministic projection of the population, setting the number of years
                           (or generations) to project the population, initial population size, and whether to invoke the density-feedback function
                           set in the previous tab.")),
            tags$li(tags$p(style="font-family:Avenir",tags$strong("STOCHASTIC"),tags$i(class="fas fa-bolt"),": stochastic projection of the population based on previous settings
                           (including the % variances set in the first tab); the user can set the number of iterations to repeat the stochastic 
                           resampling, the quasi-extinction threshold (population size below which it is considered functionally extinct), and whether to
                           invoke a generationally scaled catatastrophic mortality probability (the magnitude and variance of which can be set by the user).")),
            tags$li(tags$p(style="font-family:Avenir",tags$strong("SINGLE PULSE"),tags$i(class="fas fa-level-down-alt"),": a single 'pulse' disturbance, where the user can set the disturbance
                           to be either a proportion of the total population that is removed, or a fixed number of individuals removed, at the time (year)
                           the user wishes to invoke the pulse.")),
            tags$img(height = 100, src = "GEL Logo Kaurna transparent.png", style="float:right"),
            tags$li(tags$p(style="font-family:Avenir",tags$strong("PRESS"),tags$i(class="fas fa-compress-arrows-alt"),": a press disturbance, where the user can set the disturbance
                           to be either a proportion of the total population that is removed, or a fixed number of individuals removed, during the interval
                           over which the user wishes to invoke the press.")),
            tags$li(tags$p(style="font-family:Avenir",tags$strong("MVP"),tags$i(class="fas fa-search-minus"),": calculate the minimum viable population
                           size according to the parameters set in previous tabs."))
            
    ),
        tags$p(style="font-family:Avenir", "This", tags$i(class="fab fa-github"), "Github ",
           tags$a(href = "https://github.com/cjabradshaw/LeslieMatrixShiny", "repository"),
           "provides all the 'under-the-bonnet'",tags$i(class="fab fa-r-project"),"code for the app."),
  ),
  
  tabsetPanel(id="tabs",
              tabPanel(value="tab1", title=tags$strong("A. SET-UP", tags$i(class="fas fa-pencil-ruler")),
                wellPanel(style = "background: #d7f9da",                  
                    tags$h3(style="font-family:Avenir", tags$i(class="fas fa-pencil-ruler"), "set base parameters for matrix set-up"),
                    fluidRow(
                      column(2,
                      selectInput("agemax", label=tags$p("1. ",tags$i(class='fas fa-hourglass-end'), "max age (years,",tags$em("X"),")"),
                                  choices = seq(1,150,1), selected=(5))),
                      column(2,
                             numericInput(inputId = "survmax", label=tags$p("2. ", tags$i(class='fas fa-heart'), "adult", tags$em("s"),tags$sub("max")),
                                          value=0.9, min=0, max=1, step=0.01)),
                      column(2,
                             numericInput(inputId = "fertmax", label=tags$p("3. ",tags$i(class='fas fa-egg'), "max offspring/♀"), value=(2))),
                      column(2,
                             numericInput(inputId = "primiparity", label=tags$p("4. ",tags$i(class='fas fa-venus'), "age 1st breeding",tags$em("α")), value=(1))),
                      column(2,
                             sliderInput(inputId = "sexratio", label=tags$p("5. ",tags$i(class='fas fa-venus-mars'), "offspring sex ratio (% ♀)"),
                                         min=0, max=100, value=(50), round=F, ticks=F, step=1)),
                      column(2,
                             radioButtons("longevAbr", label=tags$p("6. ",tags$i(class='fas fa-skull-crossbones'), "abrupt end to lifespan?"), inline=T,
                                          choiceNames = list((icon("fas fa-thumbs-down")), (icon("fas fa-thumbs-up"))), choiceValues = list("no","abrupt")))
                    ) # end fluidRow
                ), # end wellPanel
                
                  mainPanel(
                    fluidRow(
                      column(3,
                             actionButton("SFfill", label=tags$p(style="font-family:Avenir", "set",tags$em("S"),"/",tags$em("F"),"vectors"),
                                          icon=shiny::icon("fas fa-list-ol"))),
                      tags$br()
                    ),
                    
                    fluidRow(
                      tags$br(),
                      column(3,
                             wellPanel(style = "background: #d7f9da",
                                       tags$p(style="font-family:Avenir", tags$strong("7.",tags$i(class='fas fa-heart'), "set", tags$em("s"),
                                       tags$sub(tags$em("x")))),uiOutput("Ssliders"))),
                      column(3,
                             wellPanel(style = "background: #d7f9da",
                                       tags$p(style="font-family:Avenir", tags$strong("8.",tags$i(class='fas fa-wave-square'), "set SD(",tags$em("s"),
                                       tags$sub(tags$em("x"))),") (%)"),uiOutput("SSDsliders"))),
                      column(3,
                             wellPanel(style = "background: #d7f9da",
                                       tags$p(style="font-family:Avenir", tags$strong("9.",tags$i(class='fas fa-egg'), "set", tags$em("f"),
                                                                                      tags$sub(tags$em("x")))),uiOutput("Fsliders"))),
                      column(3,
                             wellPanel(style = "background: #d7f9da",
                                       tags$p(style="font-family:Avenir", tags$strong("10.",tags$i(class='fas fa-wave-square'), "set SD(",tags$em("f"),
                                                                                      tags$sub(tags$em("x"))),") (%)"),uiOutput("FSDsliders"))))
                  ) # close mainPanel
              ), # end tab 1
              
              tabPanel(value="tab2", title=tags$strong("B. MATRIX PROPERTIES",tags$i(class="fas fa-table")),
                       sidebarLayout(
                         sidebarPanel(style = "background: #d7f9da",
                                      
                           tags$h3(tags$p(style="font-family:Avenir",tags$i(class="fas fa-table"), "matrix properties")),
                                      textOutput("maxlambda"),
                                      textOutput("Rmax"),
                                      textOutput("RV"),
                                      textOutput("gen"),
                                      tags$br(),
                           tags$h3(tags$p(style="font-family:Avenir", "stable age distribution")),
                                      plotOutput("SSDplot")
                         ), # end sidebarPanel

                       mainPanel(
                                wellPanel(style = "background: #d7f9da",
                                         fluidRow(
                                           column(3,
                                              actionButton("makeMatrix", label="generate/update matrix",icon=shiny::icon("fas fa-calculator"))),
                                           column(1,
                                                  tags$p(style="font-family:Avenir","or")),
                                           column(3,
                                                  fileInput("uploadMatrix", label=tags$p(tags$i(class='fas fa-upload'),"upload existing matrix"),
                                                            multiple=F, buttonLabel="choose file", placeholder="no file selected")),
                                           column(1,
                                                  tags$p(style="font-family:Avenir","and")),
                                           column(3,
                                                  fileInput("uploadSDs", label=tags$p(tags$i(class='fas fa-upload'),"upload S & F SDs"),
                                                            multiple=F, buttonLabel="choose file", placeholder="no file selected")),
                                           ),
                                         
                                tags$h3(tags$p(style="font-family:Avenir", "Leslie matrix:")),
                                 add_busy_spinner(spin="fading-circle", color="#17ca3a", timeout=500, position="bottom-right", height = 250, width = 250),
                                 tableOutput("matrix"),
                                 tags$h3(tags$p(style="font-family:Avenir", "survival & fertility standard deviations (%):")),
                                 tableOutput("SFSDs"),

                                fluidRow(
                                  column(3,
                                         downloadButton('downloadMat', 'download matrix',icon = shiny::icon("download"))),
                                  column(3,
                                         downloadButton('downloadSDs', 'download S & F SDs',icon = shiny::icon("download")))),
                                
                                ), # end wellPanel
                       ) # close mainPanel
                      ) # end sidebar Layout
                       
              ), # end tab2
              
              tabPanel(value="tab3", title=tags$strong("C. DENSITY FEEDBACK",tags$i(class="fas fa-exchange-alt")),
                       sidebarLayout(
                         sidebarPanel(style = "background: #d7f9da",
                           tags$h3(tags$p(style="font-family:Avenir",tags$i(class="fas fa-exchange-alt"), "density feedback")),
                           numericInput(inputId = "N1", label=tags$p("1.", tags$i(class='fas fa-play'), "initial", tags$em("N"), "(♀+♂)"), value=1000),
                           numericInput(inputId = "carCap", label=tags$p("2.",tags$i(class='fas fa-mountain'), "carrying capacity (♀+♂)", tags$em("K")), value=2000),

                             sliderInput(inputId = "DFa", label=tags$p("3.",tags$i(class='fab fa-modx'), tags$em("a")),
                                         min=0.0, max=2, step=0.01, value=1),
                           
                             sliderInput(inputId = "DFb", label=tags$p("4.",tags$i(class='fab fa-modx'), tags$em("b")),
                                         min=0, max=30000, value=1000),
                           
                             sliderInput(inputId = "DFc", label=tags$p("5.",tags$i(class='fab fa-modx'), tags$em("c")),
                                         min=0, max=20, step=0.1, value=8)
                           
                         ), # end sidebarPanel

                         mainPanel(
                           wellPanel(style = "background: #d7f9da",
                                     actionButton("DFcalc", label="calculate density-feedback function",icon=shiny::icon("fas fa-exchange-alt")),
                                     tags$hr(),
                                     add_busy_spinner(spin="fading-circle", color="#17ca3a", timeout=500, position="bottom-right", height = 250, width = 250),
                                     plotOutput("DDrelPlot"),
                                     tags$h3(tags$p(style="font-family:Avenir", "updated",tags$em("r"),"at",tags$em("K"))),
                                     textOutput("RmaxK"),
                                     
                           ), # end wellPanel
                         ) # close mainPanel
                       ) # end sidebar Layout
                       
              ), # end tab3
              
              
              tabPanel(value="tab4", title=tags$strong("D. PROJECT",tags$i(class="fas fa-chart-line")),
                       sidebarLayout(
                         sidebarPanel(style = "background: #d7f9da",
                           tags$h3(tags$p(style="font-family:Avenir", tags$i(class="fas fa-chart-line"),"projection parameters")),
                           radioButtons("yrsOrGen", label=tags$p("1.",tags$i(class='fas fa-user-clock'), "project years or generations?"), inline=T,
                                        choiceNames = list(icon("fas fa-calendar-alt"), icon("fas fa-seedling")),
                                        choiceValues = list("years","gens")),
                           tags$hr(),
                           conditionalPanel(
                             condition = "input.yrsOrGen == 'years'",
                             numericInput(inputId = "yrsFutureProj", label=tags$p("2.",tags$i(class='fas fa-calendar-alt'), "years to project into the future"), value=(10)),
                           ),
                           conditionalPanel(
                             condition = "input.yrsOrGen == 'gens'",
                             numericInput(inputId = "gensFutureProj", label=tags$p("2.",tags$i(class='fas fa-seedling'), "generations to project into the future"), value=(40)),
                           ),
                           numericInput(inputId = "Nstart", label=tags$p("3.",tags$i(class='fas fa-braille'), "initial", tags$em("N"), "(♀+♂)"), value=(1000)),
                           tags$hr(),
                           radioButtons("DFinvoke", label=tags$p("4.",tags$i(class='fas fa-exchange-alt'), "invoke density feedback?"), inline=T,
                                        choiceNames = list((icon("fas fa-thumbs-down")), (icon("fas fa-thumbs-up"))), choiceValues = list("no","yes")),
                           tags$br(),
                           actionButton("projectMatrix", label="project deterministic population",icon=shiny::icon("fas fa-chart-line"))
                         ), # end sidebarPanel
                         
                         mainPanel(
                           wellPanel(style = "background: #d7f9da",
                                     add_busy_spinner(spin="fading-circle", color="#17ca3a", timeout=500, position="bottom-right", height = 250, width = 250),
                                     plotOutput("detProjPlot"),
                                     tags$br(),
                                     downloadButton('downloadDetN', 'download projection',icon = shiny::icon("download"))
                           ), # end wellPanel
                         ) # close mainPanel
                       ) # end sidebar Layout
                       
              ), # end tab4

              tabPanel(value="tab5", title=tags$strong("E. STOCHASTIC",tags$i(class="fas fa-bolt")),
                       sidebarLayout(
                         sidebarPanel(style = "background: #d7f9da",
                           tags$h3(tags$p(style="font-family:Avenir",tags$i(class="fas fa-bolt"), "stochastic projection")),
                           radioButtons("iter", label=tags$p("1.",tags$i(class='fas fa-dice'), "number of iterations"), inline=T,
                                        choiceNames = list("100","1000","10000","100000"), choiceValues = list(100,1000,10000,10000)),
                           numericInput(inputId = "Qthresh", label=tags$p("2.",tags$i(class='fas fa-skull'), "quasi-extinction threshold (min ♀+♂)"), value=(100)),
                           radioButtons("CatInvoke", label=tags$p("3.",tags$i(class='fas fa-poo-storm'), "invoke catastrophes?"), inline=T,
                                        choiceNames = list((icon("fas fa-thumbs-down")), (icon("fas fa-thumbs-up"))), choiceValues = list("no","yes")),
                           conditionalPanel(
                             condition = "input.CatInvoke == 'yes'",
                                sliderInput(inputId = "catMag", label=tags$p("4.",tags$i(class='fas fa-arrows-alt-v'), "catastrophe magnitude"),
                                       min=0, max=100, value=(50), round=F, ticks=F, step=1),
                                sliderInput(inputId = "catMagSD", label=tags$p("5.",tags$i(class='fas fa-wave-square'), "catastrophe magnitude",
                                         tags$em("σ"),"%"), min=0, max=50, value=(10), round=F, ticks=F, step=1)
                             ),
                           tags$hr(),
                           tags$br(),
                           actionButton("projectStoch", label="project/update",icon=shiny::icon("fas fa-bolt"))
                           
                         ), # end sidebarPanel
                         
                         mainPanel(
                           wellPanel(style = "background: #d7f9da",
                                     tags$h3(tags$p(style="font-family:Avenir", "a. population size")),
                                     #add_busy_gif(src = "expand.gif", timeout=500, position='full-page', height = 150, width = 150), # add busy gif
                                     add_busy_spinner(spin="fading-circle", color="#17ca3a", timeout=500, position="bottom-right", height = 250, width = 250),
                                     plotOutput("projectStochPlot"),
                                     tags$br(),
                                     downloadButton('downloadStochN', 'download projection',icon = shiny::icon("download")),
                                     tags$hr(),
                                     textOutput("PrExt"),
                                     textOutput("minPop"),
                                     textOutput("rMn"),
                           ), # end wellPanel
                           
                           wellPanel(style = "background: #d7f9da",
                                     tags$h3(tags$p(style="font-family:Avenir", "b. population rate of change (r)")),
                                     plotOutput("plotRplot")
                           ) # end wellPanel
                           
                         ) # close mainPanel
                       ) # end sidebar Layout
                       
              ), # end tab5
              
              tabPanel(value="tab6", title=tags$strong("F. SINGLE PULSE",tags$i(class="fas fa-level-down-alt")),
                       sidebarLayout(
                         sidebarPanel(style = "background: #d7f9da",
                           tags$h3(tags$p(style="font-family:Avenir", tags$i(class="fas fa-level-down-alt"),"pulse disturbance")),
                           radioButtons("percOrFix", label=tags$p("1.",tags$i(class='fas fa-user-cog'), "percentage or fixed?"), inline=T,
                                        choiceNames = list(icon("fas fa-percentage"), icon("fas fa-screwdriver")),
                                        choiceValues = list("perc","fix")),
                           tags$hr(),
                           conditionalPanel(
                             condition = "input.percOrFix == 'perc'",
                             sliderInput(inputId = "percPulse", label=tags$p("2. select offtake", tags$i(class='fas fa-percentage')),
                                         min=1, max=100, value=50, round=F, ticks=F, step=1)
                           ),
                           conditionalPanel(
                             condition = "input.percOrFix == 'fix'",
                             numericInput(inputId = "fixPulse", label=tags$p("2.",tags$i(class='fas fa-screwdriver'),"select # individuals for offtake"), value=500)
                           ),
                           radioButtons("setOrRand", label=tags$p("3.",tags$i(class='fas fa-user-clock'), "random or set timing of single pulse disturbance?"), inline=T,
                                        choiceNames = list(tags$p(style="font-family:Avenir","random",icon("fas fa-random")),
                                                           tags$p(style="font-family:Avenir","set", icon("fas fa-wrench"))),
                                        choiceValues = list("random","set")),
                           conditionalPanel(
                             condition = "input.setOrRand == 'set'",
                             numericInput(inputId = "setOnePulse", label=tags$p("4.",tags$i(class='fas fa-wrench'),"select year of single pulse"), value=5)
                           ),
                           tags$hr(),
                           tags$br(),
                           actionButton("projectPulse", label="project/update",icon=shiny::icon("fas fa-level-down-alt"))
                           
                         ), # end sidebarPanel
                         
                         mainPanel(
                           wellPanel(style = "background: #d7f9da",
                                     tags$h3(tags$p(style="font-family:Avenir", "population size")),
                                     add_busy_spinner(spin="fading-circle", color="#17ca3a", timeout=500, position="bottom-right", height = 250, width = 250),
                                     plotOutput("projectPulsePlot"),
                                     tags$br(),
                                     downloadButton('downloadPulseN', 'download projection',icon = shiny::icon("download")),
                                     tags$hr(),
                                     textOutput("PrExtPulse"),
                                     textOutput("minPopPulse"),
                                     textOutput("rMnPulse")
                           ), # end wellPanel
                           
                         ) # close mainPanel
                       ) # end sidebar Layout
                       
              ), # end tab6
              
              tabPanel(value="tab7", title=tags$strong("G. PRESS",tags$i(class="fas fa-compress-arrows-alt")),
                       sidebarLayout(
                         sidebarPanel(style = "background: #d7f9da",
                                      tags$h3(tags$p(style="font-family:Avenir",tags$i(class="fas fa-compress-arrows-alt"), "press disturbance")),
                                      radioButtons("percOrFixPress", label=tags$p("1.",tags$i(class='fas fa-user-cog'), "percentage or fixed?"), inline=T,
                                                   choiceNames = list(icon("fas fa-percentage"), icon("fas fa-screwdriver")),
                                                   choiceValues = list("percPress","fixPress")),
                                      tags$hr(),
                                      conditionalPanel(
                                        condition = "input.percOrFixPress == 'percPress'",
                                        sliderInput(inputId = "percentPress", label=tags$p("2. select offtake", tags$i(class='fas fa-percentage'),"per year"),
                                                    min=1, max=100, value=50, round=F, ticks=F, step=1)
                                      ),
                                      conditionalPanel(
                                        condition = "input.percOrFixPress == 'fixPress'",
                                        numericInput(inputId = "fixedPress", label=tags$p("2.",tags$i(class='fas fa-screwdriver'),"select # individuals (♀+♂) for offtake per year"), value=500)
                                      ),
                                      radioButtons("fullOrSubsetPress", label=tags$p("3.",tags$i(class='fas fa-chart-pie'), "full or subset timing of press disturbance?"), inline=T,
                                                   choiceNames = list(tags$p(style="font-family:Avenir","full",icon("far fa-circle")),
                                                                      tags$p(style="font-family:Avenir","subset", icon("fas fa-pizza-slice"))),
                                                   choiceValues = list("fullInt","subsetInt")),
                                      conditionalPanel(
                                        condition = "input.fullOrSubsetPress == 'subsetInt'",
                                        sliderInput(inputId = "setSubsetPress", label=tags$p("4.",tags$i(class='fas fa-pizza-slice'),"select % period of
                                                                                              projection interval to which press is applied"), value=c(20,90),
                                                    min=1, max=100, step=1, ticks=F, round=F)
                                      ),
                                      tags$hr(),
                                      tags$br(),
                                      actionButton("projectPress", label="project/update",icon=shiny::icon("fas fa-compress-arrows-alt"))
                                      
                         ), # end sidebarPanel
                         
                         mainPanel(
                           wellPanel(style = "background: #d7f9da",
                                     tags$h3(tags$p(style="font-family:Avenir", "population size")),
                                     add_busy_spinner(spin="fading-circle", color="#17ca3a", timeout=500, position="bottom-right", height = 250, width = 250),
                                     plotOutput("projectPressPlot"),
                                     tags$br(),
                                     downloadButton('downloadPressN', 'download projection',icon = shiny::icon("download")),
                                     tags$hr(),
                                     textOutput("PrExtPress"),
                                     textOutput("minPopPress"),
                                     textOutput("rMnPress")
                           ), # end wellPanel
                           
                         ) # close mainPanel
                       ) # end sidebar Layout
                       
              ), # end tab7
              
              tabPanel(value="tab8", title=tags$strong("H. MVP",tags$i(class="fas fa-search-minus")),
                       sidebarLayout(
                         sidebarPanel(style = "background: #d7f9da",
                                      tags$h3(tags$p(style="font-family:Avenir",tags$i(class="fas fa-search-minus"), "minimum viable population size",
                                                     tags$em("N"),tags$sub("MVP"))),
                                      radioButtons("genOrYrs", label=tags$p("1.",tags$i(class='fas fa-user-clock'), "project generations or years?"), inline=T,
                                                   choiceNames = list(icon("fas fa-seedling"), icon("fas fa-calendar-alt")),
                                                   choiceValues = list("gen","yr")),
                                      tags$hr(),
                                      conditionalPanel(
                                        condition = "input.genOrYrs == 'gen'",
                                        numericInput(inputId = "genProj", label=tags$p("2.",tags$i(class='fas fa-seedling'), "generations to project into the future"), value=(40)),
                                      ),
                                      conditionalPanel(
                                        condition = "input.genOrYrs == 'yr'",
                                        numericInput(inputId = "yrProj", label=tags$p("2.",tags$i(class='fas fa-calendar-alt'), "years to project into the future"), value=(100)),
                                      ),
                                      sliderInput(inputId = "persistPr", label=tags$p("3.",tags$i(class='fas fa-dice-six'),"select persistence probability"),
                                                  value=0.99, min=0.800, max=0.999, step=0.001, ticks=F, round=F),
                                      numericInput(inputId = "iterMVP", label=tags$p("4.",tags$i(class='fas fa-dice'), "MVP iterations"), value=(100)),
                                      numericInput(inputId = "Nhigh", label=tags$p("5.",tags$i(class='fas fa-angle-double-up'), "upper initial", tags$em("N"), "(♀+♂)"), value=(1000)),
                                      numericInput(inputId = "Nlo", label=tags$p("6.",tags$i(class='fas fa-angle-double-down'), "lower initial", tags$em("N"), "(♀+♂)"), value=(10)),
                                      numericInput(inputId = "Nstep", label=tags$p("7.",tags$i(class='fas fa-walking'), "step size", tags$em("N")), value=(10)),

                         ), # end sidebarPanel
                         
                         mainPanel(
                           wellPanel(style = "background: #d7f9da",
                                     actionButton("calcMVP", label=tags$p("calculate",tags$em("N"),tags$sub("MVP")),icon=shiny::icon("fas fa-search-minus")),
                                     tags$h3(tags$p(style="font-family:Avenir", "change in persistence probability")),
                                     add_busy_spinner(spin="fading-circle", color="#17ca3a", timeout=500, position="top-right", height = 250, width = 250),
                                     plotOutput("MVPPlot"),
                                     tags$hr(),
                                     textOutput("MVPest"),
                                     textOutput("MVPstepChange"),
                                     textOutput("MVPasymptote")
                           ), # end wellPanel
                           
                         ) # close mainPanel
                       ) # end sidebar Layout
                       
              ), # end tab8
              
              
              tabPanel(value="tab9", title=tags$strong("I. INSTRUCTIONS",tags$i(class="fas fa-directions")), style = "background: #d7f9da",
                       wellPanel(style = "background: #d7f9da",
                                 tags$h3(style="font-family:Avenir",tags$i(class="fas fa-directions"),"Detailed instructions and notes"),
                                 tags$a(href="https://flinders.edu.au/", tags$img(height = 100, src = "F_V_CMYK.png", style="float:right",title="Flinders University")),
                                 tags$h4(style="font-family:Avenir",tags$em("Preamble")),
                                 tags$p(style="font-family:Avenir", "I designed this app for ecologists who might be daunted by the prospect of coding
                                        their own Leslie matrix models, not to mention the mathematics underlying them. There is no need to become a 
                                        modeller to use this app, but a reasonable grounding in the basic components is a good idea all the same."),
                                 tags$p(style="font-family:Avenir","My secondary intention is that once you are able to customise the app to your purposes, you might feel more
                                        capable of designing your own code from scratch (relying heavily on the functions I have provided on",
                                        tags$i(class="fab fa-github"), tags$a(href="https://github.com/cjabradshaw/LeslieMatrixShiny","Github"),"). For the more adventurous, I can recommend the following sources to find
                                        out more:"),
                                 tags$ul(tags$li(tags$p(style="font-family:Avenir", tags$i(class="fas fa-book"),"Hal Caswell's 2006 (second edition) book", tags$a(href="https://global.oup.com/academic/product/matrix-population-models-9780878931217?cc=au&lang=en&#", tags$em("Matrix Population Models:
                                                 Construction, Analysis, and Interpretation")), "(view", tags$a(href="https://www.whoi.edu/cms/files/mpm2e_tableofcontents_116984.pdf",
                                                                                                                "table of contents)"))),
                                         tags$li(tags$p(style="font-family:Avenir",tags$i(class="fab fa-github"), "Kevin Shoemaker's", tags$a(href="https://kevintshoemaker.github.io/NRES-470/LECTURE7.html",
                                                                                                                tags$em("R Matrix Population Models")), "course on Github")),
                                         tags$li(tags$p(style="font-family:Avenir",tags$i(class="fas fa-database"), "The", tags$a(href="https://compadredb.wordpress.com/2015/10/05/introducing-the-comadre-animal-matrix-database/",
                                                                                                  tags$em("COMADRE")), "(animal) and",tags$a(href="https://compadredb.wordpress.com/", tags$em("COMPADRE")), "(plant) matrix databases")),
                                 ), # end ul
                                 tags$br(),
                                 tags$h4(style="font-family:Avenir",tags$em("What is a Leslie matrix?")),
                                 tags$a(href="https://epicaustralia.org.au/", tags$img(height = 150, src = "CABAHlogo.png",
                                                                                       style="float:right", title="ARC Centre of Excellence for Australian Biodiversity and Heritage")),
                                 tags$p(style="font-family:Avenir","A 'Leslie' matrix (named after Patrick H. Leslie of Oxford in 1945) is used to 'project'
                                        an age-classified population. All matrices used in population biology are what we call 'transition' matrices,
                                        because you set the probability of transitioning from one cell to another. In the special case of Leslie matrices,
                                        the transitions are pegged to movement between discrete age classes, so there can only be one direction (i.e., an
                                        individual transitions from age", tags$em("x"), "to", tags$em("x"), "+ 1, but can never return to age", tags$em("x"),
                                        "). You therefore need to set the transition probabilities related to the aging process (age-specific survival),
                                        and to the production of new individuals (age-specific fertility). Defining these two vectors from birth to
                                        maximum age (i.e., the species' longevity) are the essential components of the Leslie matrix. Note that this app
                                        cannot be used for species with more complex life-history transitions where regression is possible, as is the
                                        case in many stage-classified life cycles. As for the timing of the discrete age classes, I generally assume this will
                                        be equivalent to 'years' (i.e., one breeding event to another). While this assumption can be relaxed in that you can
                                        set the intervals to a shorter/longer duration, all plots are expressed in terms of 'years' (number of projection units)."), 
                                 tags$br(),
                                 tags$h4(style="font-family:Avenir",tags$em("Instructions for using this app")),
                                 tags$p(style="font-family:Avenir","I have set up this app so that the settings in each tab are carried over to subsequent
                                        tabs, meaning that changes in a specific tab's settings will be reactively inherited by others farther down the
                                        chain. I recommend that you follow the order prescribed below before getting fancy and jumping around between tabs.
                                        Another important consideration is that while for many parameters the values are set according to both sexes, the 
                                        projections themselves only consider the female element of the population. All relevant parameters are halved automatically
                                        to estimate the total number of females per category assuming a 50:50 adult sex ratio. The graphs produced always
                                        show only the female component of the population. Finally, you have the capacity to save and subsequently upload
                                        all parameters set within the matrix and the survival/fertility standard deviations (see tab B)."),
                                 
                                 tags$hr(),
                                 
                                 tags$ol(type = "A", tags$li(tags$p(style="font-family:Avenir", tags$strong("SET-UP"), tags$i(class="fas fa-pencil-ruler")),
                                         tags$ol(tags$li(tags$p(style="font-family:Avenir","You first need to define how long your species will live — this is the
                                         maximum age. The matrix will have the same number of dimensions vertically and horizontally defined as maximum
                                         age + 1, because we start the process off at birth (age 0).")),
                                         tags$li(tags$p(style="font-family:Avenir","Next you can set the highest survival probability among all age classes.
                                         While this is not absolutely necessary, it assists with constructing the change in survival probabilities with
                                         age (see item A7 below).")),
                                         tags$li(tags$p(style="font-family:Avenir","Set the maximum number of total offspring (i.e., both daughters and
                                         sons) in the most fertile age class. Again, this is not absolutely necessary at this stage, but it helps defining
                                         the age-specific fertility vector later (see item A9 below).")),
                                         tags$li(tags$p(style="font-family:Avenir","Set the age of first breeding (primiparity, or α). While you can make
                                         adjustments to the fertility vector in item A9, setting this now makes things easier later.")),
                                         tags$li(tags$p(style="font-family:Avenir","Sometimes the sex ratio at birth is skewed to one gender or another,
                                         although these are usually special cases (e.g., temperature-dependent sex determination in some reptiles, or 
                                         biased maternal investment in one sex following variation in resource availability). For most species, there is usually an
                                         average of a 50:50 sex ratio at birth.")),
                                         tags$li(tags$p(style="font-family:Avenir","This item can be a bit tricky to understand (or measure in real
                                         populations). An 'abrupt' end to an organism's lifespan can occur, such as the catastrophic mortality experienced
                                         by dasyurid marsupials (including devils, quolls, and antechinus) at the end of life, or in truly semelparous
                                         species. In many other species, some rare individuals can exceed life expectancy by chance alone. The general
                                         rule of thumb is that if too many individuals are showing up in the 'final' age class, your longevity estimate
                                         might be too low, or your final-age survival probability might be too high (you can view the stable age structure in tab B).")),
                                         tags$li(tags$p(style="font-family:Avenir","Once you have set the six previous parameters, click the 'set",tags$em("S"),"/",
                                         tags$em("F"), "vectors' button to set up the survival and fertility vectors, and their respective standard deviations (%SD; 
                                         more on these below in items A8 and A10). The number of slider bars that appear in items A7-10 is set as maximum age (item A1)
                                         + 1. You can then adjust each age-specific survival probability up to maximum longevity.")),
                                         tags$li(tags$p(style="font-family:Avenir","Here you can set the standard deviation of survival (set as a
                                         percentage of that age's survival probability) for all ages. These parameters are only invoked in the stochastic projections
                                         (tabs E—G), but if you intend to make stochastic projections, it is a good idea to set these up here.")),
                                         tags$li(tags$p(style="font-family:Avenir","Just like you set the age-specific survival probabilities in item A7, here you
                                         can set the age-specific fertilities (total number of daughters and sons produced per mother).")),
                                         tags$li(tags$p(style="font-family:Avenir","As for item A8, here you can set the %SD for age-specific fertilities
                                         for use in the stochastic projections in tabs E—G.")),
                                         ) # end numbered ol
                                 ), # end A lis
                                 
                                 tags$hr(),
                                 
                                 tags$li(tags$p(style="font-family:Avenir", tags$strong("MATRIX PROPERTIES"), tags$i(class="fas fa-table")),
                                         tags$p(style="font-family:Avenir", "Once you are happy with the parameters set in tab A, click on 
                                         tab B, and then click the 'generate/update matrix' button in the main panel (alternatively, you can upload a previously
                                         saved matrix and the associated standard deviations for survival/fertility). This action will calculate the matrix
                                         and display it in the main panel. If you instead choose to upload a previous matrix file (and survival/fertility standard
                                         deviations file) — these can be saved using the 'download' buttons at the very bottom of the main panel — you will 
                                         be able to view the uploaded matrix and standard deviations in the main panel. Either option results in the
                                         calculation of some basic properties of the matrix displayed in the side panel:"),
                                         tags$ol(type="i",tags$li(tags$p(style="font-family:Avenir","This is the dominant eigen value (",tags$em("λ"),") of
                                         the matrix, which essentially means it is the matrix's capacity to produce either an increasing (",tags$em("λ"),
                                         "> 1), stable (",tags$em("λ"), "= 1), or declining (",tags$em("λ"), "< 1) population. For example, if",tags$em("λ"),
                                         "= 1.1, then a population can be expected to grow from one time interval to the next by 1.1×.")),
                                         tags$li(tags$p(style="font-family:Avenir","Parameter",tags$em("r"), "is simply the natural logarithm of",tags$em("λ"),
                                         ", which equates to the instantaneous rate of population increase. When",tags$em("r"),"> 0, the population grows,
                                         when",tags$em("r"),"= 0, the population is stable, and when",tags$em("r"),"< 0, the population declines.")),
                                         tags$li(tags$p(style="font-family:Avenir","Net reproductive rate",tags$em("R"),tags$sub("0"), "is the mean number of
                                         offspring by which a newborn individual will be replaced by the end of its life (i.e., rate of population increase
                                         from one generation to the next).")),
                                         tags$li(tags$p(style="font-family:Avenir","Generation time", tags$em("G"), "is the mean time required for the population
                                         to increase by a factor of",tags$em("R"),tags$sub("0"),".")),
                                         ) # end numbered ol
                                    ), # end B li
                                 
                                 tags$hr(),
                                 
                                 tags$li(tags$p(style="font-family:Avenir", tags$strong("DENSITY FEEDBACK"), tags$i(class="fas fa-exchange-alt")),
                                         tags$p(style="font-family:Avenir", "Without some sort of restriction on growth, any population with", 
                                         tags$em("r"), "> 0 will increase exponentially forever. Clearly this is impossible, so some upper limit must be
                                         imposed. There are many ways to do this, with the most basic being a simple cap on the total allowable population
                                         size. Imposing such a cap is mathematically simple, but biologically unrealistic (like all the individuals in a box
                                         suddenly running out of room, with every additional individual dying instantaneously). In reality, populations are more
                                         complex, so some sort of feedback mechanism typically slows at least one of the demographic rates. While
                                         biologically more realistic, measuring where density feedback actually operates can be challenging. For this reason,
                                         I have chosen to implement a middle-of-the-road solution by calculating a compensatory feedback modifier (",
                                         tags$em("S"),tags$sub("mod"), ") that reduces the survival vector as the population approaches carrying
                                         capacity", tags$em("K"),". The function in question is a basic, 3-parameter logistic-power relationship of the
                                         form:", tags$em("S"),tags$sub("mod"),"=",tags$em("a"),"/(1+(",tags$em("N"),"/",tags$em("b"),")",tags$sup("c"),")",
                                         ", where",tags$em("N"),"= total population size at any given time. The five parameters to set below should be adjusted
                                         such that the function produces",tags$em("S"),tags$sub("mod"),"declining from 1 to some value < 1, such that the
                                         value of", tags$em("r"),"at", tags$em("K"),"≈ 1."),
                                         tags$ol(tags$li(tags$p(style="font-family:Avenir","Set the initial population size. Note that this parameter is
                                         specific to the density-feedback tab; you will have the opportunity to provide the 'real' start population for 
                                         the projections in subsequent tabs.")),
                                         tags$li(tags$p(style="font-family:Avenir","Set the long-term average environmental carrying capacity", tags$em("K"),
                                         ". As the population approaches this value, survival will decline according to",tags$em("S"),tags$sub("mod"),", meaning
                                         that the population rate of change will also slow.")),
                                         tags$li(tags$p(style="font-family:Avenir","The", tags$em("a"), "parameter of the logistic-power function adjusts
                                         the initial value of",tags$em("S"),tags$sub("mod"), ". Generally this should remain at or near 1.")),
                                         tags$li(tags$p(style="font-family:Avenir","The", tags$em("b"), "parameter of the logistic-power function changes
                                         how quickly",tags$em("S"),tags$sub("mod"), "approaches",tags$em("K"),".")),
                                         tags$li(tags$p(style="font-family:Avenir","The", tags$em("c"), "parameter of the logistic-power function the shape
                                         of the relationship between",tags$em("S"),tags$sub("mod"), "and",tags$em("K"),".")),
                                          
                                         ) # end numbered ol
                                 ), # end C li
                                 
                                 tags$hr(),
                                 
                                 tags$li(tags$p(style="font-family:Avenir", tags$strong("PROJECT"), tags$i(class="fas fa-chart-line")),
                                         tags$p(style="font-family:Avenir", "This is finally where you get to view some of the fruits of your parameter-setting
                                         labour. This process takes the deterministic matrix from tabs A—B, and projects an initial population forward either taking
                                         density feedback into account (relationship from tab C) or not."),
                                         tags$ol(tags$li(tags$p(style="font-family:Avenir","Your first choice is to set how far into the future your
                                         projection will run according to either years or generations. If you choose 'generations', then the projection window inherits the 
                                         deterministic value of", tags$em("G"), "from tab B")),
                                         tags$li(tags$p(style="font-family:Avenir","Choose the value in years or generations accordingly.")),
                                         tags$li(tags$p(style="font-family:Avenir","Next, choose your initial population size (both females and males). As mentioned
                                         above, this can be different from the initial value set in tab C, but it should be comparable so that the density-feedback
                                         relationship remains meaningful to the population range projected here.")),
                                         tags$li(tags$p(style="font-family:Avenir","Decide whether to invoke the density-feedback relationship established in tab C.")),
                                        
                                         ), # end numbered ol
                                         tags$p(style="font-family:Avenir","Next, click the 'project deterministic population' button to view the projection graph
                                         in the main panel. You can then adjust the parameters to move", tags$em("r"), "at",tags$em("K"), "to ~ 0 (note that you need to
                                         reclick the 'project deterministic population' button to update the", tags$em("r"), "at",tags$em("K"), "value, even if the graph
                                         updates automatically).")
                                 ), # end D li
                                 
                                 tags$hr(),
                                 
                                 tags$li(tags$p(style="font-family:Avenir", tags$strong("STOCHASTIC"), tags$i(class="fas fa-bolt")),
                                         tags$p(style="font-family:Avenir", "Up to here the neat, deterministic projection of your population belies the uncertainty
                                         in such mathematical wizardry. A stochastic projection incorporates reasonable uncertainty in both the parameters used to
                                         construct the matrix, as well as random, unforeseen 'catastrophic' die-offs caused by a variety of natural disturbances. The
                                         app allows you to incorporate this uncertainty by resampling your base survival and fertility vectors, using a",tags$em("β"),"distribution
                                         for the former and a Gaussian distribution for the latter. This is where the % standard deviations for each age-specific survival
                                         probability and fertility set in tab A come into play. There are three (or five, if you choose to invoke catastrophes)
                                         parameters to set in this tab:"),
                                         tags$ol(tags$li(tags$p(style="font-family:Avenir","The number of iterations to resample the demographic values. The more iterations
                                         you chooose, the better your confidence intervals will be estimated at the expense of longer computation time. I advise starting
                                         with the lowest number of iterations first to check out how the model behaves, and then increasing this later once you are happy with
                                         set-up.")),
                                         tags$li(tags$p(style="font-family:Avenir","Set what is known as the 'quasi'-extinction threshold. Because small populations tend
                                         to be susceptible to extinction for reasons that are typically different from those that caused them to become small in the 
                                         first place, we generally invoke a threshold >> 1 below which we deem the population to be 'functionally' extinct (i.e., high
                                         chance of going extinct anyway).", tags$a(href="http://doi.org/10.1016/j.biocon.2013.12.036", "Theory and data"), "suggest that this can be as high has 100 individuals, but you can adjust
                                         the threshold according to your own views. Note that all subsequent estimates of extinction probability Pr(ext) are base on this
                                         threshold.")),
                                         tags$li(tags$p(style="font-family:Avenir","Among vertebrates, there is", tags$a(href="http://doi.org/10.1017/S1367943003003147",
                                         "good evidence"), "to suggest that 'catastrophic' die-offs tend to happen with a predictable probability per generation. In this case,
                                         I've used the mean value of 0.14 per generation. Choose whether you want to invoke this type of event occurring in the subsequent
                                         stochastic projections.")),
                                         tags$li(tags$p(style="font-family:Avenir","If you choose to invoke the catastrophe function, you will have the option of defining
                                         what a 'catastrophic' die-off means in terms of average magnitude of increased mortality. This defaults to 50% mortality, which is the
                                         definition used in the",tags$a(href="http://doi.org/10.1017/S1367943003003147", "aforementioned paper"),
                                         " (note that this is",tags$em("β"),"-resampled per iteration according to the standard deviation % set below).")),
                                         tags$li(tags$p(style="font-family:Avenir","Set the standard deviation % for resampling the catastrophic die-off magnitude set in
                                         the previous slider.")),
                                                 
                                         ), # end numbered ol
                                         tags$p(style="font-family:Avenir","Next, click the 'project/update' button to view the graphs projecting population size",
                                         tags$em("N"), "(main panel graph a) and the instantaneous rate of population increase", tags$em("r"), " (main panel graph b).
                                         Both these graphs show the 95% confidence intervals of the temporal trend in grey. There are also three emergent values
                                         provided by the simulation shown below panel figure a (note that you need to reclick the 'project/update' button to
                                         update the following values even if the graph updates automatically):"),
                                         tags$ol(type="i",tags$li(tags$p(style="font-family:Avenir","Pr(Ext) is the probability of extinction (i.e., number of iterations where
                                         the total population size fell below the quasi-extinction threshold).")),
                                         tags$li(tags$p(style="font-family:Avenir","min N is minimum population size achieved during the projection window averaged across all
                                         iterations. The 95% confidence interval (range) of this value is also provied.")),
                                         tags$li(tags$p(style="font-family:Avenir","mean r is the mean population rate of increase from one time step to the next throughout the
                                         projection interval, averaged (+ 95% confidence interval range).")),

                                         ), # end numbered ol
                                         
                                 ), # end E li
                                 
                                 tags$hr(),
                                 
                                 tags$li(tags$p(style="font-family:Avenir", tags$strong("SINGLE PULSE"), tags$i(class="fas fa-level-down-alt")),
                                         tags$p(style="font-family:Avenir", "A 'pulse' disturbance is an acute perturbation that happens abruptly. I implemented this
                                         function to test the effect of one-off disturbance events like a harvest or known catastrophe (in addition or independent of the 
                                         'normal' probability of a generic catastrophes implemented in the previous tab). You have the choice to set a specific timing
                                         for the pulse disturbance, or let it happen at any time randomly during the projection window (and differently in every
                                         stochastic iteration). You can also choose to set the perturbation as a percentage of the current population size, or as a
                                         fixed number of individuals (spread across the entire age range according to the stable age distribution from panel B):"),
                                         tags$ol(tags$li(tags$p(style="font-family:Avenir","First set the disturbance pulse as either a percentage of the current population,
                                         or as a fixed number of indiviuals.")),
                                         tags$li(tags$p(style="font-family:Avenir","Set either the percentage or fixed value of the disturbance pulse.")),
                                         tags$li(tags$p(style="font-family:Avenir","Allow the disturbance pulse to happen randomly, or at a set time in during the projection
                                         window.")),
                                         tags$li(tags$p(style="font-family:Avenir","If you set the timing, when should this disturbance pulse occur?")),
                                                 
                                         ), # end numbered ol
                                         tags$p(style="font-family:Avenir","Next, click the 'project/update' button to view the graph projecting population size",
                                                tags$em("N"), "(main panel graph). The graph includes the 95% confidence intervals of the temporal trend in grey. There are
                                                also three emergent values provided by the simulation shown below the main figure panel (note that you need to reclick the 'project/update' button to
                                         update the following values even if the graph updates automatically):"),
                                         tags$ol(type="i",tags$li(tags$p(style="font-family:Avenir","Pr(Ext) is the probability of extinction (i.e., number of iterations where
                                         the total population size fell below the quasi-extinction threshold).")),
                                                 tags$li(tags$p(style="font-family:Avenir","min N is minimum population size achieved during the projection window averaged across all
                                         iterations. The 95% confidence interval (range) of this value is also provied.")),
                                                 tags$li(tags$p(style="font-family:Avenir","mean r is the mean population rate of increase from one time step to the next throughout the
                                         projection interval, averaged (+ 95% confidence interval range).")),
                                                 
                                         ), # end numbered ol
                                 ), # end F li
                                 
                                 tags$hr(),
                                 
                                 tags$li(tags$p(style="font-family:Avenir", tags$strong("PRESS"), tags$i(class="fas fa-compress-arrows-alt")),
                                         tags$p(style="font-family:Avenir", "A 'press' disturbance is a perturbation that happens over a longer time frame
                                         than a pulse disturbance. I implemented this function to test the effect of a sustained disturbance event like an annual
                                         harvest (in addition or independent of the 'normal' probability of a generic catastrophes implemented in tab E). You have
                                         the choice to set a specific interval for the press disturbance, or let it happen throughout the entire projection
                                         window. You can also choose to set the perturbation as a percentage of the current population size, or as a
                                         fixed number of individuals (spread across the entire age range according to the stable age distribution from panel B):"),
                                         tags$ol(tags$li(tags$p(style="font-family:Avenir","First set the disturbance press as either a percentage of the current population,
                                         or as a fixed number of indiviuals.")),
                                         tags$li(tags$p(style="font-family:Avenir","Set either the percentage or fixed value of the disturbance press.")),
                                         tags$li(tags$p(style="font-family:Avenir","Allow the disturbance press to occur over the entire projection window,
                                         or only during a set interval (percentage) of the projection window.")),
                                         tags$li(tags$p(style="font-family:Avenir","If you set the interval, between what percentages of the window should the press
                                         disturbance occur?")),
                                         
                                         ), # end numbered ol
                                         tags$p(style="font-family:Avenir","Next, click the 'project/update' button to view the graph projecting population size",
                                                tags$em("N"), "(main panel graph). The graph includes the 95% confidence intervals of the temporal trend in grey. There are
                                                also three emergent values provided by the simulation shown below the main figure panelå (note that you need to reclick the 'project/update' button to
                                         update the following values even if the graph updates automatically):"),
                                         tags$ol(type="i",tags$li(tags$p(style="font-family:Avenir","Pr(Ext) is the probability of extinction (i.e., number of iterations where
                                         the total population size fell below the quasi-extinction threshold).")),
                                                 tags$li(tags$p(style="font-family:Avenir","min N is minimum population size achieved during the projection window averaged across all
                                         iterations. The 95% confidence interval (range) of this value is also provied.")),
                                                 tags$li(tags$p(style="font-family:Avenir","mean r is the mean population rate of increase from one time step to the next throughout the
                                         projection interval, averaged (+ 95% confidence interval range).")),
                                                 
                                         ), # end numbered ol
                                 ), # end G li
                                 
                                 tags$hr(),
                                 
                                 tags$li(tags$p(style="font-family:Avenir", tags$strong("MVP"), tags$i(class="fas fa-search-minus")),
                                         tags$p(style="font-family:Avenir", "A minimum viable population (MVP) is a population of size", tags$em("N"),
                                         "at time 0 that has a specified probabiliy of persisting time",tags$em("t"),"into the future. Clearly, the value of",
                                         tags$em("N"),tags$sub("MVP"),"therefore depends not only on the demographic parameters of the population in question,
                                         it also depends on the choice of",tags$em("persistence probability"),"as well as the time to project the population
                                         into the future. These latter parameters can be somewhat arbitrary choices, but convention generally defaults to a 99%
                                         probability of persisting",tags$a(href="https://onlinelibrary.wiley.com/doi/full/10.1111/j.1461-0248.2006.00883.x", "40 generations"),
                                         "40 generations (preferred) or 100 years",tags$a(href="https://portals.iucn.org/library/node/10315","(IUCN Red List definition)"),
                                         ". Any compensatory feedback imposed on the projections also complicates the issue insofar as the user assumes a constant,
                                         long-term carrying capacity to which the population generally returns after a perturbation. Nonetheless, it is a useful concept
                                         fully stochastic demographic projections of populations demonstrate generally that estimates of", tags$em("N"),tags$sub("MVP"),"tend
                                         to align with the",tags$a(href="http://dx.doi.org/10.1016/j.biocon.2013.12.036", "'100/1000 rule'"),"(formerly known as the '50/500' rule) derived from genetic and evolutionary 
                                         considerations. For the calculation of MVP, the process invokes the characteristics defined in earlier tabs; however, the user can set
                                         the following parameters in addition:"),
                                         tags$ol(tags$li(tags$p(style="font-family:Avenir","Set whether to project for number of generations or years.")),
                                                 tags$li(tags$p(style="font-family:Avenir","Depending on the previous choice, set number of generations (default = 40) or
                                                                number of years (default = 100).")),
                                                 tags$li(tags$p(style="font-family:Avenir","Select the persistence probability (default = 0.99)")),
                                                 tags$li(tags$p(style="font-family:Avenir","Set the number of iterations for each initial population-size
                                                                run. Note that large values will slow down the calculation considerably, so I recommend
                                                                starting with low values first before settling on the desired parameters.")),
                                                 tags$li(tags$p(style="font-family:Avenir","Set the largest initial population size from which progressively
                                                                smaller values will be assessed.")),
                                                 tags$li(tags$p(style="font-family:Avenir","Set the minimum initial population size to test.")),
                                                 tags$li(tags$p(style="font-family:Avenir","Set the population step size that will be used to define the intervals
                                                                tested between the maximum and minimum inital population sizes set above. Note that a higher number
                                                                of intervals will slow down the calculation.")),
                                                 
                                         ), # end numbered ol
                                         tags$p(style="font-family:Avenir","Next, click the ' calculate", tags$em("N"),tags$sub("MVP"),"' button to begin the
                                                calculation. Once complete, a graph will appear showing the relationship between initial population size
                                                and persistence probability over the projection window set above. Below the graph, three outputs will appear:"),
                                         tags$ol(type="i",tags$li(tags$p(style="font-family:Avenir","Provided the simulations had at least one initial population
                                         size that resulted in the persistence probability being achieved or exceeded, the minimum viable population size (number
                                         of female individuals) will be displayed here.")),
                                                 tags$li(tags$p(style="font-family:Avenir","A 'breakpoint' is displayed showing the first precipitous reduction in
                                                                persistence probability as initial population size declines.")),
                                         tags$a(href="https://github.com/cjabradshaw/LeslieMatrixShiny/blob/main/LICENSE", tags$img(height = 50, src = "GNU GPL3.png", style="float:right", title="GNU General Public Licence v3.0")),
                                                 
                                                 tags$li(tags$p(style="font-family:Avenir","I have also included a logistic power function of the form: 
                                                                Pr(Persistence) =", tags$em("a"), "/ (1 + (", tags$em("N"), "/", tags$em("b"),")",
                                                                tags$sup(tags$em("c")), "), where parameters", tags$em("a, b, c"), "are constants calculated
                                                                using a non-linear optimisation function, and", tags$em("N"), "is the number of initial adult
                                                                females in the population.", "This logistic function (given as one of the outputs) provides
                                                                an 'asymptotic' estimate of minimum viable population size, with a corresponding confidence
                                                                interval derived from a multinomial resampling procedure (light red-shaded region in the
                                                                corresponding graph). Note that if the parameters are not chosen carefully, this function
                                                                cannot be fit and the asymptotic MVP size will not be displayed."))
                                                 
                                                 
                                         ), # end numbered ol
                                 ), # end H li
                                 
                                 ), # end upper-case letters ol
                                 
                       ) # end wellPanel
              ) # end tab9
              
  ) # end tabset
  
) # close fluidPage

server <- function(input, output, session) {

  observeEvent(input$SFfill, {

    if(input$tabs == "tab1"){
      
    output$Ssliders <- renderUI({
      ageclasses <- as.integer(input$agemax) + 1
      lapply(1:ageclasses, function(i) {
        sliderInput(inputId = paste0("S", (i-1)), label = paste("age", (i-1)),
                    min = 0, max = 1, value=ifelse((i-1) == 0, 0.5*input$survmax, input$survmax), round=F, ticks=F, step = 0.01)
      })
    })

    output$SSDsliders <- renderUI({
      ageclasses <- as.integer(input$agemax) + 1
      lapply(1:ageclasses, function(i) {
        sliderInput(inputId = paste0("SSD", (i-1)), label = paste("age", (i-1)),
                    min = 0, max = 100, value=10, round=F, ticks=F, step = 1)
      })
    })

    output$Fsliders <- renderUI({
      ageclasses <- as.integer(input$agemax) + 1
      lapply(1:ageclasses, function(i) {
        sliderInput(inputId = paste0("F", (i-1)), label = paste("age", (i-1)),
                    min = 0, max = input$fertmax, value= ifelse((i-1) < input$primiparity, 0, ifelse((i-1) == input$primiparity, 0.5*input$fertmax, input$fertmax)), round=F, ticks=F, step = 0.1)
      })
    })
    
    output$FSDsliders <- renderUI({
      ageclasses <- as.integer(input$agemax) + 1
      lapply(1:ageclasses, function(i) {
        sliderInput(inputId = paste0("FSD", (i-1)), label = paste("age", (i-1)),
                    min = 0, max = 100, value=10, round=F, ticks=F, step = 1)
      })
    })
    } # end tab1 if
    
  })
  
    
      observeEvent(input$makeMatrix, {
        
        if(input$tabs == "tab2"){
          
      inputsRctv <- reactiveValues()
      observe({
        inputsRctv$am <- as.numeric(input$agemax)
        inputsRctv$sr <- as.numeric(input$sexratio)
      })
      
      Srctv <- reactive({
        ageclasses <- as.integer(input$agemax) + 1
        Svec <<- as.numeric(unlist(lapply(1:ageclasses, function(i) {
          input[[paste0("S", (i-1))]]
        })))
      }) # end Srctv

      Frctv <- reactive({
        ageclasses <- as.integer(input$agemax) + 1
        Fvec <<- as.numeric(unlist(lapply(1:ageclasses, function(i) {
            input[[paste0("F", (i-1))]]
        })))
      }) # end Frctv
      
      DemRdat <<- isolate({
        dat <- data.frame(Srctv(), Frctv())
      })
      
      observe({
        output$matrix <- renderTable(bordered=T,colnames=F,striped=T, {
          Dmat <<- createLmatrix(age.max=inputsRctv$am, Svec=DemRdat[,1], Fvec=(inputsRctv$sr/100) * DemRdat[,2],
                                 finalStage=as.character(input$longevAbr)) # deterministic matrix
        })
      }) # end observe
      
      output$downloadMat <- downloadHandler(
        filename = function() {
          paste("matrixOut", "csv", sep = ".")
        },
        
        content = function(file) {
          sep <- ","
          
          write.table(Dmat, file, sep=sep, row.names = F, col.names = F)
        }
      )
      
      S_SDrctv <- reactive({
        ageclasses <- as.integer(input$agemax) + 1
        S_SDvec <<- as.numeric(unlist(lapply(1:ageclasses, function(i) {
          input[[paste0("SSD", (i-1))]]
        })))
      }) # end Srctv
      
      F_SDrctv <- reactive({
        ageclasses <- as.integer(input$agemax) + 1
        F_SDvec <<- as.numeric(unlist(lapply(1:ageclasses, function(i) {
          input[[paste0("FSD", (i-1))]]
        })))
      }) # end Frctv
      
      DemRSDdat <<- isolate({
        SDdat <- data.frame(S_SDrctv(), F_SDrctv())
      })
      
      observe({
        output$SFSDs <- renderTable(bordered=T,rownames=F,colnames=F,striped=T, {
          SDdat
        })
      })
      
      output$downloadSDs <- downloadHandler(
        filename = function() {
          paste("SFsd", "csv", sep = ".")
        },
        
        content = function(file2) {
          sep <- ","
          
          write.table(DemRSDdat, file2, sep=sep, row.names = F, col.names = F)
        }
      )
      
      observe({
        output$maxlambda <- renderText( {
          maxl <<- maxLambda(Dmat)
          paste("i. max λ = ", round(maxl, 4), "(λ < 1 → N↓; λ > 1 → N↑)")
        })
      }) # end observe
      observe({
        output$Rmax <- renderText( {
          maxr <<- maxR(Dmat)
          paste("ii. max r = ", round(maxr, 4), "(r < 0 → N↓; r > 0 → N↑)")
        })
      }) # end observe
      observe({
        output$gen <- renderText( {
          G <<- Gval(Dmat, inputsRctv$am+1)
          paste("iii. G = ", round(G, 4), "years")
        })
      }) # end observe
      observe({
        output$RV <- renderText( {
          R <<- Rval(Dmat, inputsRctv$am+1)
          paste("iv. R0 = ", round(R, 4), "(lifetime ♀ offspring/♀)")
        })
      }) # end observe
      
      observe({
        output$SSDplot <- renderPlot( {
          ssdDat <<- data.frame(0:(dim(Dmat)[1]-1), round(StableStageDist(Dmat), 3))
          colnames(ssdDat) <- c("age","relProp")
          
          Ctheme = theme(
            axis.title.x = element_text(size = 16),
            axis.text.x = element_text(size = 14),
            axis.title.y = element_text(size = 16),
            axis.text.y = element_text(size = 14))
          
          ggplot(ssdDat, aes(x=age, y=relProp)) +
            geom_point() +
            geom_path() +
            labs(x="age (years)", y="relative proportion in population") +
            Ctheme
            
        })
      }) # end observe
      
      } # end tab 2 if
      
    }) # end Events
  
      
      observeEvent(input$uploadMatrix, {
        
        if(input$tabs == "tab2"){
          
        observe({
          output$matrix <- renderTable(bordered=T,rownames=F,colnames=F,striped=T, {
            file_to_read = input$uploadMatrix
            if(is.null(file_to_read)){
              return()
            }
            read.table(file_to_read$datapath, sep=",", header=F)
          }) # end output table1
      
          Dmat2 <<- eventReactive(input$uploadMatrix, {
              as.matrix(read.table(input$uploadMatrix$datapath, sep=",", header = F))
            })
        })
        
          observe({
            output$SFSDs <- renderTable(bordered=T,rownames=F,colnames=F,striped=T, {
              file_to_read = input$uploadSDs
              if(is.null(file_to_read)){
                return()
              }
              read.table(file_to_read$datapath, sep=",", header=F)
            }) # end output table2
            
            SFSD <<- eventReactive(input$uploadSDs, {
              (read.table(input$uploadSDs$datapath, sep=",", header = F))
            })
          })
          
        observe({
          output$maxlambda <- renderText( {
            maxl <<- maxLambda(Dmat2())
            paste("i. max λ = ", round(maxl, 4), "(λ < 1 → N↓; λ > 1 → N↑)")
          })
        }) # end observe
        observe({
          output$Rmax <- renderText( {
            maxr <<- maxR(Dmat2())
            paste("ii. max r = ", round(maxr, 4), "(r < 0 → N↓; r > 0 → N↑)")
          })
        }) # end observe
        observe({
          output$gen <- renderText( {
            G <<- Gval(Dmat2(), dim(Dmat2())[1])
            paste("iii. G = ", round(G, 4), "years")
          })
        }) # end observe
        observe({
          output$RV <- renderText( {
            R <<- Rval(Dmat2(), dim(Dmat2())[1])
            paste("iv. R0 = ", round(R, 4), "(lifetime ♀ offspring/♀)")
          })
        }) # end observe
        
        observe({
          output$SSDplot <- renderPlot( {
            ssdDat <<- data.frame(0:(dim(Dmat2())[1]-1), round(StableStageDist(Dmat2()), 3))
            colnames(ssdDat) <- c("age","relProp")
            
            Ctheme = theme(
              axis.title.x = element_text(size = 16),
              axis.text.x = element_text(size = 14),
              axis.title.y = element_text(size = 16),
              axis.text.y = element_text(size = 14))
            
            ggplot(ssdDat, aes(x=age, y=relProp)) +
              geom_point() +
              geom_path() +
              labs(x="age (years)", y="relative proportion in population") +
              Ctheme
            
          })
        }) # end observe

      } # end tab 2 if
        
      }) # end Event

        
  observeEvent(input$DFcalc, {
    
    if(input$tabs == "tab3"){
      
      observe({
        
        output$DDrelPlot <- renderPlot({
          
          NKvec <- seq((input$N1/2), (input$carCap/2), 1)
          CompOut <<- data.frame(NKvec, as.numeric(input$DFa)/(1+(as.numeric(NKvec)/as.numeric(input$DFb))^as.numeric(input$DFc)))
          colnames(CompOut) <- c("Kcont", "predSmod")
          
          Ctheme = theme(
            axis.title.x = element_text(size = 16),
            axis.text.x = element_text(size = 14),
            axis.title.y = element_text(size = 16),
            axis.text.y = element_text(size = 14))
          
          ggplot(CompOut, aes(x=Kcont, y=predSmod)) +
            geom_path() +
            labs(x="N ♀", y=expression("S"[mod])) +
            Ctheme
          
        })
      }) # end observe

      DM2 <- reactive({
        if (exists("Dmat2")) {
          DM <<- Dmat2()
        }
      })

      Dmat <<- isolate({
        if(exists("Dmat2")) {
          datmat <- as.matrix(DM2())
        } else {
          datmat <- Dmat
        }
      })
  
      inputsRctv <- reactiveValues()
      observe({
        inputsRctv$am <- ifelse(exists("Dmat2")==T, dim(Dmat)[1] - 1, as.numeric(input$agemax))
        inputsRctv$sr <- as.numeric(input$sexratio)
      })
      
      Srctv <- reactive({
        ageclasses <- ifelse(exists("Dmat2")==T, dim(Dmat)[1], as.integer(input$agemax) + 1)
        if (exists("Dmat2")==T) {
          Svec <<- c(as.numeric(diag(Dmat[2:dim(Dmat)[1], ])), as.numeric(Dmat[dim(Dmat)[1],dim(Dmat)[1]]))
        } else {
        Svec <<- as.numeric(unlist(lapply(1:ageclasses, function(i) {
          input[[paste0("S", (i-1))]]
        })))}
      }) # end Srctv
      
      Frctv <- reactive({
        ageclasses <- ifelse(exists("Dmat2")==T, dim(Dmat)[1], as.integer(input$agemax) + 1)
        if (exists("Dmat2")==T) {
          Fvec <<- as.numeric(Dmat[1,])
        } else {
          Fvec <<- as.numeric(unlist(lapply(1:ageclasses, function(i) {
          input[[paste0("F", (i-1))]]
        })))}
      }) # end Frctv
      
      DemRdat <<- isolate({
        dat <- data.frame(Srctv(), Frctv())
      })
      
      observe({
        SmodK <<- as.numeric(input$DFa)/(1+(as.numeric(input$carCap/2)/as.numeric(input$DFb))^as.numeric(input$DFc))
        SvecUpdate <<- SmodK * DemRdat[,1]
        
        if (exists("Dmat2")==T) {
          SR <<- 1
          
          if (as.numeric(Dmat[dim(Dmat)[1],dim(Dmat)[1]]) == 0) {
            finalStageChar <<- "abrupt"
          }
          if (as.numeric(Dmat[dim(Dmat)[1],dim(Dmat)[1]]) > 0) {
            finalStageChar <<- "no"  
          }
        } else {
          finalStageChar <<- as.character(input$longevAbr)
          SR <<- inputsRctv$sr/100
        }
          
        DmatUpdate <<- createLmatrix(age.max=(dim(Dmat)[1] - 1), Svec=SvecUpdate, Fvec=(SR * DemRdat[,2]),
                               finalStage=finalStageChar)
      })
      
      reactiveVal({
        output$RmaxK <- renderText({
          Kmaxr <<- maxR(DmatUpdate)
          paste(" = ", round(Kmaxr, 4), "(target: 0)")
        })
      })  
      
    } # end tab3 if
    
  })
  
  observeEvent(input$projectMatrix, {
    
    if(input$tabs == "tab4"){
      
      observe({
        if (input$DFinvoke == "no" & input$yrsOrGen == "years") {
          
          output$detProjPlot <- renderPlot({
            NprojOut <<- projLmatrix(x=Dmat, initN=(as.numeric(input$Nstart)/2), projYrs=input$yrsFutureProj)
            
            Ctheme = theme(
              axis.title.x = element_text(size = 16),
              axis.text.x = element_text(size = 14),
              axis.title.y = element_text(size = 16),
              axis.text.y = element_text(size = 14))
              
             ggplot(NprojOut, aes(x=yrs, y=npred)) +
               geom_path(linetype=2) +
               labs(x="years into future", y="N (♀ only)") +
               Ctheme
          })
        } # end if
        
        if (input$DFinvoke == "yes" & input$yrsOrGen == "years") {
            
          output$detProjPlot <- renderPlot({
            NprojOut <<- projLmatrix(x=Dmat, initN=input$Nstart/2, projYrs=input$yrsFutureProj, DFparams = c(input$DFa, input$DFb, input$DFc))
            
            Ctheme = theme(
              axis.title.x = element_text(size = 16),
              axis.text.x = element_text(size = 14),
              axis.title.y = element_text(size = 16),
              axis.text.y = element_text(size = 14))
            
            ggplot(NprojOut, aes(x=yrs, y=npred)) +
              geom_path(linetype=2) +
              labs(x="years into future", y="N (♀ only)") +
              Ctheme
          })
          
        } # end if

        if (input$DFinvoke == "no" & input$yrsOrGen == "gens") {
          
          output$detProjPlot <- renderPlot({
            NprojOut <<- projLmatrix(x=Dmat, initN=(as.numeric(input$Nstart)/2), projGens=input$gensFutureProj)
            
            Ctheme = theme(
              axis.title.x = element_text(size = 16),
              axis.text.x = element_text(size = 14),
              axis.title.y = element_text(size = 16),
              axis.text.y = element_text(size = 14))
            
            ggplot(NprojOut, aes(x=yrs, y=npred)) +
              geom_path(linetype=2) +
              labs(x="years into future", y="N (♀ only)") +
              Ctheme
          })
        } # end if
        
        if (input$DFinvoke == "yes" & input$yrsOrGen == "gens") {
          
          output$detProjPlot <- renderPlot({
            NprojOut <<- projLmatrix(x=Dmat, initN=input$Nstart/2, projGens=input$gensFutureProj, DFparams = c(input$DFa, input$DFb, input$DFc))
            
            Ctheme = theme(
              axis.title.x = element_text(size = 16),
              axis.text.x = element_text(size = 14),
              axis.title.y = element_text(size = 16),
              axis.text.y = element_text(size = 14))
            
            ggplot(NprojOut, aes(x=yrs, y=npred)) +
              geom_path(linetype=2) +
              labs(x="years into future", y="N (♀ only)") +
              Ctheme
          })
          
        } # end if
        
        output$downloadDetN <- downloadHandler(
          filename = function() {
            paste("deterministicNproj", "csv", sep = ".")
          },
          
          content = function(file) {
            sep <- ","
            NprojDL <<- NprojOut
            NprojDL[,2] <- round(NprojOut[,2], 0)
            write.table(NprojDL, file, sep=sep, row.names = F, col.names = T)
          }
        )
        
        
      }) # end observe

    } # end tab 4 if
    
  }) # end tab Events
  
  observeEvent(input$projectStoch, {
    
    if(input$tabs == "tab5"){
      
      if (exists("SFSD")==T) {
        
        SD2 <- reactive({
            SD <<- SFSD()
        })
        
        DemRSDdat <<- isolate({
            SDdat <- SD2()
        })
        
      } else {
      
        S_SDrctv <- reactive({
          ageclasses <- as.integer(input$agemax) + 1
          S_SDvec <<- as.numeric(unlist(lapply(1:ageclasses, function(i) {
            input[[paste0("SSD", (i-1))]]
          })))
        }) # end Srctv
        
        F_SDrctv <- reactive({
          ageclasses <- as.integer(input$agemax) + 1
          F_SDvec <<- as.numeric(unlist(lapply(1:ageclasses, function(i) {
            input[[paste0("FSD", (i-1))]]
          })))
        }) # end Frctv
        
        DemRSDdat <<- isolate({
            SDdat <- data.frame(S_SDrctv(), F_SDrctv())
        })
      } # end else
      
      observe({
        ###################
        ## project by years
        ###################
        if (input$CatInvoke == "no" & input$DFinvoke == "no" & input$yrsOrGen == "years") {
          
          
          output$projectStochPlot <- renderPlot({
            NprojStochNoDF <<- projStochMatrix(x=Dmat, initN=input$Nstart/2, projYrs=input$yrsFutureProj, iter=as.numeric(input$iter),
                                               S_SD=DemRSDdat[,1], F_SD=DemRSDdat[,2], Qthresh=input$Qthresh/2)
            
            Ctheme = theme(
              axis.title.x = element_text(size = 16),
              axis.text.x = element_text(size = 14),
              axis.title.y = element_text(size = 16),
              axis.text.y = element_text(size = 14))
            
            ggplot(NprojStochNoDF$Nrange, aes(x=yrvec, y=Nmd)) +
              geom_ribbon(aes(ymin = Nlo, ymax = Nup), fill = "grey70", alpha=0.4) +
              geom_path(linetype=2) +
              labs(x="years into future", y="N (♀ only)") +
              Ctheme
          })
          
          output$downloadStochN <- downloadHandler(
            filename = function() {
              paste("stochasticNproj", "csv", sep = ".")
            },
            
            content = function(file) {
              sep <- ","
              stochNprojDL <<- NprojStochNoDF$Nrange
              stochNprojDL[,2] <- round(NprojStochNoDF$Nrange[,2], 0)
              stochNprojDL[,3] <- round(NprojStochNoDF$Nrange[,3], 0)
              stochNprojDL[,4] <- round(NprojStochNoDF$Nrange[,4], 0)
              write.table(stochNprojDL, file, sep=sep, row.names = F, col.names = T)
            }
          )
          
          reactiveVal({
            output$PrExt <- renderText({
              paste("i. Pr(Ext) = ", round(NprojStochNoDF$PrQext, 4))
            })
          }) # end observe
          
          reactiveVal({
            output$minPop <- renderText({
              paste("ii. min N = ", round(NprojStochNoDF$minMnN, 0), "(range: ", round(NprojStochNoDF$minLoN, 0), " to ", round(NprojStochNoDF$minUpN, 0),")")
            })
          }) # end observe
          
          reactiveVal({
            output$rMn <- renderText({
              paste("iii. mean r = ", round(NprojStochNoDF$mnRmn, 2), "(range: ", round(NprojStochNoDF$mnRlo, 2), " to ", round(NprojStochNoDF$mnRup, 2),")")
            })
          }) # end observe
          
          output$plotRplot <- renderPlot({
            Ctheme = theme(
              axis.title.x = element_text(size = 16),
              axis.text.x = element_text(size = 14),
              axis.title.y = element_text(size = 16),
              axis.text.y = element_text(size = 14))
            
            ggplot(NprojStochNoDF$rRange, aes(x=yrvec, y=rMn)) +
              geom_ribbon(aes(ymin = rLo, ymax = rUp), fill = "grey70", alpha=0.4) +
              geom_path(linetype=2) +
              labs(x="years into future", y="r") +
              Ctheme
          })
          
        } # end if
        
        if (input$CatInvoke == "yes" & input$DFinvoke == "no" & input$yrsOrGen == "years") {
          
         
          output$projectStochPlot <- renderPlot({
            
            NprojStochCatNoDF <<- projStochMatrix(x=Dmat, initN=input$Nstart/2, projYrs=input$yrsFutureProj, iter=as.numeric(input$iter),
                                                  S_SD=DemRSDdat[,1], F_SD=DemRSDdat[,2], Qthresh=input$Qthresh/2,
                                                  CatParams = c(as.numeric(input$catMag), as.numeric(input$catMagSD)))
            Ctheme = theme(
              axis.title.x = element_text(size = 16),
              axis.text.x = element_text(size = 14),
              axis.title.y = element_text(size = 16),
              axis.text.y = element_text(size = 14))
            
            ggplot(NprojStochCatNoDF$Nrange, aes(x=yrvec, y=Nmd)) +
              geom_ribbon(aes(ymin = Nlo, ymax = Nup), fill = "grey70", alpha=0.4) +
              geom_path(linetype=2) +
              labs(x="years into future", y="N (♀ only)") +
              Ctheme
          })
          
          output$downloadStochN <- downloadHandler(
            filename = function() {
              paste("stochasticNproj", "csv", sep = ".")
            },
            
            content = function(file) {
              sep <- ","
              stochNprojDL <<- NprojStochCatNoDF$Nrange
              stochNprojDL[,2] <- round(NprojStochCatNoDF$Nrange[,2], 0)
              stochNprojDL[,3] <- round(NprojStochCatNoDF$Nrange[,3], 0)
              stochNprojDL[,4] <- round(NprojStochCatNoDF$Nrange[,4], 0)
              write.table(stochNprojDL, file, sep=sep, row.names = F, col.names = T)
            }
          )
          
          reactiveVal({
            output$PrExt <- renderText({
              paste("i. Pr(Ext) = ", round(NprojStochCatNoDF$PrQext, 4))
            })
          }) # end observe
          
          reactiveVal({
            output$minPop <- renderText({
              paste("ii. min N = ", round(NprojStochCatNoDF$minMnN, 0), "(range: ", round(NprojStochCatNoDF$minLoN, 0), " to ", round(NprojStochCatNoDF$minUpN, 0),")")
            })
          }) # end observe
          
          reactiveVal({
            output$rMn <- renderText({
              paste("iii. mean r = ", round(NprojStochCatNoDF$mnRmn, 2), "(range: ", round(NprojStochCatNoDF$mnRlo, 2), " to ", round(NprojStochCatNoDF$mnRup, 2),")")
            })
          }) # end observe
          
          output$plotRplot <- renderPlot({
            Ctheme = theme(
              axis.title.x = element_text(size = 16),
              axis.text.x = element_text(size = 14),
              axis.title.y = element_text(size = 16),
              axis.text.y = element_text(size = 14))
            
            ggplot(NprojStochCatNoDF$rRange, aes(x=yrvec, y=rMn)) +
              geom_ribbon(aes(ymin = rLo, ymax = rUp), fill = "grey70", alpha=0.4) +
              geom_path(linetype=2) +
              labs(x="years into future", y="r") +
              Ctheme
          })
          
        } # end if
        
        if (input$CatInvoke == "no" & input$DFinvoke == "yes" & input$yrsOrGen == "years") {
          

          output$projectStochPlot <- renderPlot({
            NprojStoch <<- projStochMatrix(x=Dmat, initN=input$Nstart/2, projYrs=input$yrsFutureProj, iter=as.numeric(input$iter), 
                                           S_SD=DemRSDdat[,1], F_SD=DemRSDdat[,2], Qthresh=input$Qthresh/2,
                                           DFparams = c(input$DFa, input$DFb, input$DFc))
            
            Ctheme = theme(
              axis.title.x = element_text(size = 16),
              axis.text.x = element_text(size = 14),
              axis.title.y = element_text(size = 16),
              axis.text.y = element_text(size = 14))
            
            ggplot(NprojStoch$Nrange, aes(x=yrvec, y=Nmd)) +
              geom_ribbon(aes(ymin = Nlo, ymax = Nup), fill = "grey70", alpha=0.4) +
              geom_path(linetype=2) +
              labs(x="years into future", y="N (♀ only)") +
              Ctheme
          })
          
          output$downloadStochN <- downloadHandler(
            filename = function() {
              paste("stochasticNproj", "csv", sep = ".")
            },
            
            content = function(file) {
              sep <- ","
              stochNprojDL <<- NprojStoch$Nrange
              stochNprojDL[,2] <- round(NprojStoch$Nrange[,2], 0)
              stochNprojDL[,3] <- round(NprojStoch$Nrange[,3], 0)
              stochNprojDL[,4] <- round(NprojStoch$Nrange[,4], 0)
              write.table(stochNprojDL, file, sep=sep, row.names = F, col.names = T)
            }
          )
          
          reactiveVal({
            output$PrExt <- renderText({
              paste("i. Pr(Ext) = ", round(NprojStoch$PrQext, 4))
            })
          }) # end observe
          
          reactiveVal({
            output$minPop <- renderText({
              paste("ii. min N = ", round(NprojStoch$minMnN, 0), "(range: ", round(NprojStoch$minLoN, 0), " to ", round(NprojStoch$minUpN, 0),")")
            })
          }) # end observe
          
          reactiveVal({
            output$rMn <- renderText({
              paste("iii. mean r = ", round(NprojStoch$mnRmn, 2), "(range: ", round(NprojStoch$mnRlo, 2), " to ", round(NprojStoch$mnRup, 2),")")
            })
          }) # end observe
          
          output$plotRplot <- renderPlot({
            Ctheme = theme(
              axis.title.x = element_text(size = 16),
              axis.text.x = element_text(size = 14),
              axis.title.y = element_text(size = 16),
              axis.text.y = element_text(size = 14))
            
            ggplot(NprojStoch$rRange, aes(x=yrvec, y=rMn)) +
              geom_ribbon(aes(ymin = rLo, ymax = rUp), fill = "grey70", alpha=0.4) +
              geom_path(linetype=2) +
              labs(x="years into future", y="r") +
              Ctheme
          })
          
        } # end if
        
        if (input$CatInvoke == "yes" & input$DFinvoke == "yes" & input$yrsOrGen == "years") {
          

          output$projectStochPlot <- renderPlot({
            
            NprojStochCat <<- projStochMatrix(x=Dmat, initN=input$Nstart/2, projYrs=input$yrsFutureProj, iter=as.numeric(input$iter),
                                              S_SD=DemRSDdat[,1], F_SD=DemRSDdat[,2], Qthresh=input$Qthresh/2,
                                              DFparams = c(input$DFa, input$DFb, input$DFc), 
                                              CatParams = c(as.numeric(input$catMag), as.numeric(input$catMagSD)))
            
            Ctheme = theme(
              axis.title.x = element_text(size = 16),
              axis.text.x = element_text(size = 14),
              axis.title.y = element_text(size = 16),
              axis.text.y = element_text(size = 14))

            ggplot(NprojStochCat$Nrange, aes(x=yrvec, y=Nmd)) +
              geom_ribbon(aes(ymin = Nlo, ymax = Nup), fill = "grey70", alpha=0.4) +
              geom_path(linetype=2) +
              labs(x="years into future", y="N (♀ only)") +
              Ctheme
          })
          
          output$downloadStochN <- downloadHandler(
            filename = function() {
              paste("stochasticNproj", "csv", sep = ".")
            },
            
            content = function(file) {
              sep <- ","
              stochNprojDL <<- NprojStochCat$Nrange
              stochNprojDL[,2] <- round(NprojStochCat$Nrange[,2], 0)
              stochNprojDL[,3] <- round(NprojStochCat$Nrange[,3], 0)
              stochNprojDL[,4] <- round(NprojStochCat$Nrange[,4], 0)
              write.table(stochNprojDL, file, sep=sep, row.names = F, col.names = T)
            }
          )
          
          reactiveVal({
            output$PrExt <- renderText({
              paste("i. Pr(Ext) = ", round(NprojStochCat$PrQext, 4))
            })
          }) # end observe

          reactiveVal({
            output$minPop <- renderText({
              paste("ii. min N = ", round(NprojStochCat$minMnN, 0), "(range: ", round(NprojStochCat$minLoN, 0), " to ", round(NprojStochCat$minUpN, 0),")")
            })
          }) # end observe

          reactiveVal({
            output$rMn <- renderText({
              paste("iii. mean r = ", round(NprojStochCat$mnRmn, 2), "(range: ", round(NprojStochCat$mnRlo, 2), " to ", round(NprojStochCat$mnRup, 2),")")
            })
          }) # end observe
          
          output$plotRplot <- renderPlot({
            Ctheme = theme(
              axis.title.x = element_text(size = 16),
              axis.text.x = element_text(size = 14),
              axis.title.y = element_text(size = 16),
              axis.text.y = element_text(size = 14))
            
            ggplot(NprojStochCat$rRange, aes(x=yrvec, y=rMn)) +
              geom_ribbon(aes(ymin = rLo, ymax = rUp), fill = "grey70", alpha=0.4) +
              geom_path(linetype=2) +
              labs(x="years into future", y="r") +
              Ctheme
          })
          
          
      } # end if

      #########################  
      ## project by generations
      #########################  
        if (input$CatInvoke == "no" & input$DFinvoke == "no" & input$yrsOrGen == "gens") {
          
         
          output$projectStochPlot <- renderPlot({
            NprojStochNoDF <<- projStochMatrix(x=Dmat, initN=input$Nstart/2, projGens=input$gensFutureProj, iter=as.numeric(input$iter),
                                               S_SD=DemRSDdat[,1], F_SD=DemRSDdat[,2], Qthresh=input$Qthresh/2)
            
            Ctheme = theme(
              axis.title.x = element_text(size = 16),
              axis.text.x = element_text(size = 14),
              axis.title.y = element_text(size = 16),
              axis.text.y = element_text(size = 14))
            
            ggplot(NprojStochNoDF$Nrange, aes(x=yrvec, y=Nmd)) +
              geom_ribbon(aes(ymin = Nlo, ymax = Nup), fill = "grey70", alpha=0.4) +
              geom_path(linetype=2) +
              labs(x="years into future", y="N (♀ only)") +
              Ctheme
          })
          
          output$downloadStochN <- downloadHandler(
            filename = function() {
              paste("stochasticNproj", "csv", sep = ".")
            },
            
            content = function(file) {
              sep <- ","
              stochNprojDL <<- NprojStochNoDF$Nrange
              stochNprojDL[,2] <- round(NprojStochNoDF$Nrange[,2], 0)
              stochNprojDL[,3] <- round(NprojStochNoDF$Nrange[,3], 0)
              stochNprojDL[,4] <- round(NprojStochNoDF$Nrange[,4], 0)
              write.table(stochNprojDL, file, sep=sep, row.names = F, col.names = T)
            }
          )
          
          reactiveVal({
            output$PrExt <- renderText({
              paste("i. Pr(Ext) = ", round(NprojStochNoDF$PrQext, 4))
            })
          }) # end observe
          
          reactiveVal({
            output$minPop <- renderText({
              paste("ii. min N = ", round(NprojStochNoDF$minMnN, 0), "(range: ", round(NprojStochNoDF$minLoN, 0), " to ", round(NprojStochNoDF$minUpN, 0),")")
            })
          }) # end observe
          
          reactiveVal({
            output$rMn <- renderText({
              paste("iii. mean r = ", round(NprojStochNoDF$mnRmn, 2), "(range: ", round(NprojStochNoDF$mnRlo, 2), " to ", round(NprojStochNoDF$mnRup, 2),")")
            })
          }) # end observe
          
          output$plotRplot <- renderPlot({
            Ctheme = theme(
              axis.title.x = element_text(size = 16),
              axis.text.x = element_text(size = 14),
              axis.title.y = element_text(size = 16),
              axis.text.y = element_text(size = 14))
            
            ggplot(NprojStochNoDF$rRange, aes(x=yrvec, y=rMn)) +
              geom_ribbon(aes(ymin = rLo, ymax = rUp), fill = "grey70", alpha=0.4) +
              geom_path(linetype=2) +
              labs(x="years into future", y="r") +
              Ctheme
          })
          
        } # end if
        
        if (input$CatInvoke == "yes" & input$DFinvoke == "no" & input$yrsOrGen == "gens") {
          

          output$projectStochPlot <- renderPlot({
            
            NprojStochCatNoDF <<- projStochMatrix(x=Dmat, initN=input$Nstart/2, projGens=input$gensFutureProj, iter=as.numeric(input$iter),
                                                  S_SD=DemRSDdat[,1], F_SD=DemRSDdat[,2], Qthresh=input$Qthresh/2,
                                                  CatParams = c(as.numeric(input$catMag), as.numeric(input$catMagSD)))
            Ctheme = theme(
              axis.title.x = element_text(size = 16),
              axis.text.x = element_text(size = 14),
              axis.title.y = element_text(size = 16),
              axis.text.y = element_text(size = 14))
            
            ggplot(NprojStochCatNoDF$Nrange, aes(x=yrvec, y=Nmd)) +
              geom_ribbon(aes(ymin = Nlo, ymax = Nup), fill = "grey70", alpha=0.4) +
              geom_path(linetype=2) +
              labs(x="years into future", y="N (♀ only)") +
              Ctheme
          })
          
          output$downloadStochN <- downloadHandler(
            filename = function() {
              paste("stochasticNproj", "csv", sep = ".")
            },
            
            content = function(file) {
              sep <- ","
              stochNprojDL <<- NprojStochCatNoDF$Nrange
              stochNprojDL[,2] <- round(NprojStochCatNoDF$Nrange[,2], 0)
              stochNprojDL[,3] <- round(NprojStochCatNoDF$Nrange[,3], 0)
              stochNprojDL[,4] <- round(NprojStochCatNoDF$Nrange[,4], 0)
              write.table(stochNprojDL, file, sep=sep, row.names = F, col.names = T)
            }
          )
          
          reactiveVal({
            output$PrExt <- renderText({
              paste("i. Pr(Ext) = ", round(NprojStochCatNoDF$PrQext, 4))
            })
          }) # end observe
          
          reactiveVal({
            output$minPop <- renderText({
              paste("ii. min N = ", round(NprojStochCatNoDF$minMnN, 0), "(range: ", round(NprojStochCatNoDF$minLoN, 0), " to ", round(NprojStochCatNoDF$minUpN, 0),")")
            })
          }) # end observe
          
          reactiveVal({
            output$rMn <- renderText({
              paste("iii. mean r = ", round(NprojStochCatNoDF$mnRmn, 2), "(range: ", round(NprojStochCatNoDF$mnRlo, 2), " to ", round(NprojStochCatNoDF$mnRup, 2),")")
            })
          }) # end observe
          
          output$plotRplot <- renderPlot({
            Ctheme = theme(
              axis.title.x = element_text(size = 16),
              axis.text.x = element_text(size = 14),
              axis.title.y = element_text(size = 16),
              axis.text.y = element_text(size = 14))
            
            ggplot(NprojStochCatNoDF$rRange, aes(x=yrvec, y=rMn)) +
              geom_ribbon(aes(ymin = rLo, ymax = rUp), fill = "grey70", alpha=0.4) +
              geom_path(linetype=2) +
              labs(x="years into future", y="r") +
              Ctheme
          })
          
        } # end if
        
        if (input$CatInvoke == "no" & input$DFinvoke == "yes" & input$yrsOrGen == "gens") {
          

          output$projectStochPlot <- renderPlot({
            NprojStoch <<- projStochMatrix(x=Dmat, initN=input$Nstart/2, projGens=input$gensFutureProj, iter=as.numeric(input$iter), 
                                           S_SD=DemRSDdat[,1], F_SD=DemRSDdat[,2], Qthresh=input$Qthresh/2,
                                           DFparams = c(input$DFa, input$DFb, input$DFc))
            
            Ctheme = theme(
              axis.title.x = element_text(size = 16),
              axis.text.x = element_text(size = 14),
              axis.title.y = element_text(size = 16),
              axis.text.y = element_text(size = 14))
            
            ggplot(NprojStoch$Nrange, aes(x=yrvec, y=Nmd)) +
              geom_ribbon(aes(ymin = Nlo, ymax = Nup), fill = "grey70", alpha=0.4) +
              geom_path(linetype=2) +
              labs(x="years into future", y="N (♀ only)") +
              Ctheme
          })
          
          output$downloadStochN <- downloadHandler(
            filename = function() {
              paste("stochasticNproj", "csv", sep = ".")
            },
            
            content = function(file) {
              sep <- ","
              stochNprojDL <<- NprojStoch$Nrange
              stochNprojDL[,2] <- round(NprojStoch$Nrange[,2], 0)
              stochNprojDL[,3] <- round(NprojStoch$Nrange[,3], 0)
              stochNprojDL[,4] <- round(NprojStoch$Nrange[,4], 0)
              write.table(stochNprojDL, file, sep=sep, row.names = F, col.names = T)
            }
          )
          
          reactiveVal({
            output$PrExt <- renderText({
              paste("i. Pr(Ext) = ", round(NprojStoch$PrQext, 4))
            })
          }) # end observe
          
          reactiveVal({
            output$minPop <- renderText({
              paste("ii. min N = ", round(NprojStoch$minMnN, 0), "(range: ", round(NprojStoch$minLoN, 0), " to ", round(NprojStoch$minUpN, 0),")")
            })
          }) # end observe
          
          reactiveVal({
            output$rMn <- renderText({
              paste("iii. mean r = ", round(NprojStoch$mnRmn, 2), "(range: ", round(NprojStoch$mnRlo, 2), " to ", round(NprojStoch$mnRup, 2),")")
            })
          }) # end observe
          
          output$plotRplot <- renderPlot({
            Ctheme = theme(
              axis.title.x = element_text(size = 16),
              axis.text.x = element_text(size = 14),
              axis.title.y = element_text(size = 16),
              axis.text.y = element_text(size = 14))
            
            ggplot(NprojStoch$rRange, aes(x=yrvec, y=rMn)) +
              geom_ribbon(aes(ymin = rLo, ymax = rUp), fill = "grey70", alpha=0.4) +
              geom_path(linetype=2) +
              labs(x="years into future", y="r") +
              Ctheme
          })
          
        } # end if
        
        if (input$CatInvoke == "yes" & input$DFinvoke == "yes" & input$yrsOrGen == "gens") {
          

          output$projectStochPlot <- renderPlot({
            
            NprojStochCat <<- projStochMatrix(x=Dmat, initN=input$Nstart/2, projGens=input$gensFutureProj, iter=as.numeric(input$iter),
                                              S_SD=DemRSDdat[,1], F_SD=DemRSDdat[,2], Qthresh=input$Qthresh/2,
                                              DFparams = c(input$DFa, input$DFb, input$DFc), 
                                              CatParams = c(as.numeric(input$catMag), as.numeric(input$catMagSD)))
            
            Ctheme = theme(
              axis.title.x = element_text(size = 16),
              axis.text.x = element_text(size = 14),
              axis.title.y = element_text(size = 16),
              axis.text.y = element_text(size = 14))
            
            ggplot(NprojStochCat$Nrange, aes(x=yrvec, y=Nmd)) +
              geom_ribbon(aes(ymin = Nlo, ymax = Nup), fill = "grey70", alpha=0.4) +
              geom_path(linetype=2) +
              labs(x="years into future", y="N (♀ only)") +
              Ctheme
          })
          
          output$downloadStochN <- downloadHandler(
            filename = function() {
              paste("stochasticNproj", "csv", sep = ".")
            },
            
            content = function(file) {
              sep <- ","
              stochNprojDL <<- NprojStochCat$Nrange
              stochNprojDL[,2] <- round(NprojStochCat$Nrange[,2], 0)
              stochNprojDL[,3] <- round(NprojStochCat$Nrange[,3], 0)
              stochNprojDL[,4] <- round(NprojStochCat$Nrange[,4], 0)
              write.table(stochNprojDL, file, sep=sep, row.names = F, col.names = T)
            }
          )
          
          reactiveVal({
            output$PrExt <- renderText({
              paste("i. Pr(Ext) = ", round(NprojStochCat$PrQext, 4))
            })
          }) # end observe
          
          reactiveVal({
            output$minPop <- renderText({
              paste("ii. min N = ", round(NprojStochCat$minMnN, 0), "(range: ", round(NprojStochCat$minLoN, 0), " to ", round(NprojStochCat$minUpN, 0),")")
            })
          }) # end observe
          
          reactiveVal({
            output$rMn <- renderText({
              paste("iii. mean r = ", round(NprojStochCat$mnRmn, 2), "(range: ", round(NprojStochCat$mnRlo, 2), " to ", round(NprojStochCat$mnRup, 2),")")
            })
          }) # end observe
          
          output$plotRplot <- renderPlot({
            Ctheme = theme(
              axis.title.x = element_text(size = 16),
              axis.text.x = element_text(size = 14),
              axis.title.y = element_text(size = 16),
              axis.text.y = element_text(size = 14))
            
            ggplot(NprojStochCat$rRange, aes(x=yrvec, y=rMn)) +
              geom_ribbon(aes(ymin = rLo, ymax = rUp), fill = "grey70", alpha=0.4) +
              geom_path(linetype=2) +
              labs(x="years into future", y="r") +
              Ctheme
          })
          
        } # end if
        
      }) # end observe
      
    } # end tab 5 if
    
  }) # end tab Events
  
  
  observeEvent(input$projectPulse, {
    
    if(input$tabs == "tab6"){
      
      observe({
        ###################
        ## project by years
        ###################
        if (input$percOrFix == "perc" & input$setOrRand == "random" & input$yrsOrGen == "years") {
          
          output$projectPulsePlot <- renderPlot({
            
            NprojStochPulsePrecRand <<- projStochPulse(x=Dmat, initN=input$Nstart/2, projYrs=input$yrsFutureProj, iter=as.numeric(input$iter),
                                              S_SD=DemRSDdat[,1], F_SD=DemRSDdat[,2], Qthresh=input$Qthresh/2,
                                              DFparams = c(input$DFa, input$DFb, input$DFc), 
                                              CatParams = c(as.numeric(input$catMag), as.numeric(input$catMagSD)),
                                              PercOff=as.numeric(input$percPulse))
            
            Ctheme = theme(
              axis.title.x = element_text(size = 16),
              axis.text.x = element_text(size = 14),
              axis.title.y = element_text(size = 16),
              axis.text.y = element_text(size = 14))
            
            ggplot(NprojStochPulsePrecRand$Nrange, aes(x=yrvec, y=Nmd)) +
              geom_ribbon(aes(ymin = Nlo, ymax = Nup), fill = "grey70", alpha=0.4) +
              geom_path(linetype=2) +
              labs(x="years into future", y="N (♀ only)") +
              Ctheme
          })
          
          output$downloadPulseN <- downloadHandler(
            filename = function() {
              paste("pulseNproj", "csv", sep = ".")
            },
            
            content = function(file) {
              sep <- ","
              pulseNprojDL <<- NprojStochPulsePrecRand$Nrange
              pulseNprojDL[,2] <- round(NprojStochPulsePrecRand$Nrange[,2], 0)
              pulseNprojDL[,3] <- round(NprojStochPulsePrecRand$Nrange[,3], 0)
              pulseNprojDL[,4] <- round(NprojStochPulsePrecRand$Nrange[,4], 0)
              write.table(pulseNprojDL, file, sep=sep, row.names = F, col.names = T)
            }
          )
          
          reactiveVal({
            output$PrExtPulse <- renderText({
              paste("i. Pr(Ext) = ", round(NprojStochPulsePrecRand$PrQext, 4))
            })
          }) # end observe
          
          reactiveVal({
            output$minPopPulse <- renderText({
              paste("ii. min N = ", round(NprojStochPulsePrecRand$minMnN, 0), "(range: ", round(NprojStochPulsePrecRand$minLoN, 0), " to ", round(NprojStochPulsePrecRand$minUpN, 0),")")
            })
          }) # end observe
          
          reactiveVal({
            output$rMnPulse <- renderText({
              paste("iii. mean r = ", round(NprojStochPulsePrecRand$mnRmn, 2), "(range: ", round(NprojStochPulsePrecRand$mnRlo, 2), " to ", round(NprojStochPulsePrecRand$mnRup, 2),")")
            })
          }) # end observe
          
        } # end if
        
        if (input$percOrFix == "perc" & input$setOrRand == "set" & input$yrsOrGen == "years") {
          
          output$projectPulsePlot <- renderPlot({
            
            NprojStochPulsePrecSet <<- projStochPulse(x=Dmat, initN=input$Nstart/2, projYrs=input$yrsFutureProj, iter=as.numeric(input$iter),
                                                       S_SD=DemRSDdat[,1], F_SD=DemRSDdat[,2], Qthresh=input$Qthresh/2,
                                                       DFparams = c(input$DFa, input$DFb, input$DFc), 
                                                       CatParams = c(as.numeric(input$catMag), as.numeric(input$catMagSD)),
                                                       PercOff=as.numeric(input$percPulse), TimeOff=as.numeric(input$setOnePulse))
            
            Ctheme = theme(
              axis.title.x = element_text(size = 16),
              axis.text.x = element_text(size = 14),
              axis.title.y = element_text(size = 16),
              axis.text.y = element_text(size = 14))
            
            ggplot(NprojStochPulsePrecSet$Nrange, aes(x=yrvec, y=Nmd)) +
              geom_ribbon(aes(ymin = Nlo, ymax = Nup), fill = "grey70", alpha=0.4) +
              geom_path(linetype=2) +
              geom_vline(xintercept=as.numeric(input$setOnePulse), linetype=3, color="black", size=0.5) +
              labs(x="years into future", y="N (♀ only)") +
              Ctheme
          })
          
          output$downloadPulseN <- downloadHandler(
            filename = function() {
              paste("pulseNproj", "csv", sep = ".")
            },
            
            content = function(file) {
              sep <- ","
              pulseNprojDL <<- NprojStochPulsePrecSet$Nrange
              pulseNprojDL[,2] <- round(NprojStochPulsePrecSet$Nrange[,2], 0)
              pulseNprojDL[,3] <- round(NprojStochPulsePrecSet$Nrange[,3], 0)
              pulseNprojDL[,4] <- round(NprojStochPulsePrecSet$Nrange[,4], 0)
              write.table(pulseNprojDL, file, sep=sep, row.names = F, col.names = T)
            }
          )
          
          reactiveVal({
            output$PrExtPulse <- renderText({
              paste("i. Pr(Ext) = ", round(NprojStochPulsePrecSet$PrQext, 4))
            })
          }) # end observe
          
          reactiveVal({
            output$minPopPulse <- renderText({
              paste("ii. min N = ", round(NprojStochPulsePrecSet$minMnN, 0), "(range: ", round(NprojStochPulsePrecSet$minLoN, 0), " to ", round(NprojStochPulsePrecSet$minUpN, 0),")")
            })
          }) # end observe
          
          reactiveVal({
            output$rMnPulse <- renderText({
              paste("iii. mean r = ", round(NprojStochPulsePrecSet$mnRmn, 2), "(range: ", round(NprojStochPulsePrecSet$mnRlo, 2), " to ", round(NprojStochPulsePrecSet$mnRup, 2),")")
            })
          }) # end observe
          
        } # end if
        
        
        if (input$percOrFix == "fix" & input$setOrRand == "random" & input$yrsOrGen == "years") {
          
          output$projectPulsePlot <- renderPlot({
            
            NprojStochPulseFixRand <<- projStochPulse(x=Dmat, initN=input$Nstart/2, projYrs=input$yrsFutureProj, iter=as.numeric(input$iter),
                                                       S_SD=DemRSDdat[,1], F_SD=DemRSDdat[,2], Qthresh=input$Qthresh/2,
                                                       DFparams = c(input$DFa, input$DFb, input$DFc), 
                                                       CatParams = c(as.numeric(input$catMag), as.numeric(input$catMagSD)),
                                                       FixOff=as.numeric(input$fixPulse))
            
            Ctheme = theme(
              axis.title.x = element_text(size = 16),
              axis.text.x = element_text(size = 14),
              axis.title.y = element_text(size = 16),
              axis.text.y = element_text(size = 14))
            
            ggplot(NprojStochPulseFixRand$Nrange, aes(x=yrvec, y=Nmd)) +
              geom_ribbon(aes(ymin = Nlo, ymax = Nup), fill = "grey70", alpha=0.4) +
              geom_path(linetype=2) +
              labs(x="years into future", y="N (♀ only)") +
              Ctheme
          })
          
          output$downloadPulseN <- downloadHandler(
            filename = function() {
              paste("pulseNproj", "csv", sep = ".")
            },
            
            content = function(file) {
              sep <- ","
              pulseNprojDL <<- NprojStochPulseFixRand$Nrange
              pulseNprojDL[,2] <- round(NprojStochPulseFixRand$Nrange[,2], 0)
              pulseNprojDL[,3] <- round(NprojStochPulseFixRand$Nrange[,3], 0)
              pulseNprojDL[,4] <- round(NprojStochPulseFixRand$Nrange[,4], 0)
              write.table(pulseNprojDL, file, sep=sep, row.names = F, col.names = T)
            }
          )
          
          reactiveVal({
            output$PrExtPulse <- renderText({
              paste("i. Pr(Ext) = ", round(NprojStochPulseFixRand$PrQext, 4))
            })
          }) # end observe
          
          reactiveVal({
            output$minPopPulse <- renderText({
              paste("ii. min N = ", round(NprojStochPulseFixRand$minMnN, 0), "(range: ", round(NprojStochPulseFixRand$minLoN, 0), " to ", round(NprojStochPulseFixRand$minUpN, 0),")")
            })
          }) # end observe
          
          reactiveVal({
            output$rMnPulse <- renderText({
              paste("iii. mean r = ", round(NprojStochPulseFixRand$mnRmn, 2), "(range: ", round(NprojStochPulseFixRand$mnRlo, 2), " to ", round(NprojStochPulseFixRand$mnRup, 2),")")
            })
          }) # end observe
          
        } # end if
        
        if (input$percOrFix == "fix" & input$setOrRand == "set" & input$yrsOrGen == "years") {
          
          output$projectPulsePlot <- renderPlot({
            
            NprojStochPulseFixSet <<- projStochPulse(x=Dmat, initN=input$Nstart/2, projYrs=input$yrsFutureProj, iter=as.numeric(input$iter),
                                                      S_SD=DemRSDdat[,1], F_SD=DemRSDdat[,2], Qthresh=input$Qthresh/2,
                                                      DFparams = c(input$DFa, input$DFb, input$DFc), 
                                                      CatParams = c(as.numeric(input$catMag), as.numeric(input$catMagSD)),
                                                      FixOff=as.numeric(input$fixPulse), TimeOff=as.numeric(input$setOnePulse))
            
            Ctheme = theme(
              axis.title.x = element_text(size = 16),
              axis.text.x = element_text(size = 14),
              axis.title.y = element_text(size = 16),
              axis.text.y = element_text(size = 14))
            
            ggplot(NprojStochPulseFixSet$Nrange, aes(x=yrvec, y=Nmd)) +
              geom_ribbon(aes(ymin = Nlo, ymax = Nup), fill = "grey70", alpha=0.4) +
              geom_path(linetype=2) +
              geom_vline(xintercept=as.numeric(input$setOnePulse), linetype=3, color="black", size=0.5) +
              labs(x="years into future", y="N (♀ only)") +
              Ctheme
          })
          
          output$downloadPulseN <- downloadHandler(
            filename = function() {
              paste("pulseNproj", "csv", sep = ".")
            },
            
            content = function(file) {
              sep <- ","
              pulseNprojDL <<- NprojStochPulseFixSet$Nrange
              pulseNprojDL[,2] <- round(NprojStochPulseFixSet$Nrange[,2], 0)
              pulseNprojDL[,3] <- round(NprojStochPulseFixSet$Nrange[,3], 0)
              pulseNprojDL[,4] <- round(NprojStochPulseFixSet$Nrange[,4], 0)
              write.table(pulseNprojDL, file, sep=sep, row.names = F, col.names = T)
            }
          )
          
          reactiveVal({
            output$PrExtPulse <- renderText({
              paste("i. Pr(Ext) = ", round(NprojStochPulseFixSet$PrQext, 4))
            })
          }) # end observe
          
          reactiveVal({
            output$minPopPulse <- renderText({
              paste("ii. min N = ", round(NprojStochPulseFixSet$minMnN, 0), "(range: ", round(NprojStochPulseFixSet$minLoN, 0), " to ", round(NprojStochPulseFixSet$minUpN, 0),")")
            })
          }) # end observe
          
          reactiveVal({
            output$rMnPulse <- renderText({
              paste("iii. mean r = ", round(NprojStochPulseFixSet$mnRmn, 2), "(range: ", round(NprojStochPulseFixSet$mnRlo, 2), " to ", round(NprojStochPulseFixSet$mnRup, 2),")")
            })
          }) # end observe
          
        } # end if
        
        #########################
        ## project by generations
        #########################
        if (input$percOrFix == "perc" & input$setOrRand == "random" & input$yrsOrGen == "gens") {
          
          output$projectPulsePlot <- renderPlot({
            
            NprojStochPulsePrecRand <<- projStochPulse(x=Dmat, initN=input$Nstart/2, projGens=input$gensFutureProj, iter=as.numeric(input$iter),
                                                       S_SD=DemRSDdat[,1], F_SD=DemRSDdat[,2], Qthresh=input$Qthresh/2,
                                                       DFparams = c(input$DFa, input$DFb, input$DFc), 
                                                       CatParams = c(as.numeric(input$catMag), as.numeric(input$catMagSD)),
                                                       PercOff=as.numeric(input$percPulse))
            
            Ctheme = theme(
              axis.title.x = element_text(size = 16),
              axis.text.x = element_text(size = 14),
              axis.title.y = element_text(size = 16),
              axis.text.y = element_text(size = 14))
            
            ggplot(NprojStochPulsePrecRand$Nrange, aes(x=yrvec, y=Nmd)) +
              geom_ribbon(aes(ymin = Nlo, ymax = Nup), fill = "grey70", alpha=0.4) +
              geom_path(linetype=2) +
              labs(x="years into future", y="N (♀ only)") +
              Ctheme
          })
          
          output$downloadPulseN <- downloadHandler(
            filename = function() {
              paste("pulseNproj", "csv", sep = ".")
            },
            
            content = function(file) {
              sep <- ","
              pulseNprojDL <<- NprojStochPulsePrecRand$Nrange
              pulseNprojDL[,2] <- round(NprojStochPulsePrecRand$Nrange[,2], 0)
              pulseNprojDL[,3] <- round(NprojStochPulsePrecRand$Nrange[,3], 0)
              pulseNprojDL[,4] <- round(NprojStochPulsePrecRand$Nrange[,4], 0)
              write.table(pulseNprojDL, file, sep=sep, row.names = F, col.names = T)
            }
          )
          
          reactiveVal({
            output$PrExtPulse <- renderText({
              paste("i. Pr(Ext) = ", round(NprojStochPulsePrecRand$PrQext, 4))
            })
          }) # end observe
          
          reactiveVal({
            output$minPopPulse <- renderText({
              paste("ii. min N = ", round(NprojStochPulsePrecRand$minMnN, 0), "(range: ", round(NprojStochPulsePrecRand$minLoN, 0), " to ", round(NprojStochPulsePrecRand$minUpN, 0),")")
            })
          }) # end observe
          
          reactiveVal({
            output$rMnPulse <- renderText({
              paste("iii. mean r = ", round(NprojStochPulsePrecRand$mnRmn, 2), "(range: ", round(NprojStochPulsePrecRand$mnRlo, 2), " to ", round(NprojStochPulsePrecRand$mnRup, 2),")")
            })
          }) # end observe
          
        } # end if
        
        if (input$percOrFix == "perc" & input$setOrRand == "set" & input$yrsOrGen == "gens") {
          
          output$projectPulsePlot <- renderPlot({
            
            NprojStochPulsePrecSet <<- projStochPulse(x=Dmat, initN=input$Nstart/2, projGens=input$gensFutureProj, iter=as.numeric(input$iter),
                                                      S_SD=DemRSDdat[,1], F_SD=DemRSDdat[,2], Qthresh=input$Qthresh/2,
                                                      DFparams = c(input$DFa, input$DFb, input$DFc), 
                                                      CatParams = c(as.numeric(input$catMag), as.numeric(input$catMagSD)),
                                                      PercOff=as.numeric(input$percPulse), TimeOff=as.numeric(input$setOnePulse))
            
            Ctheme = theme(
              axis.title.x = element_text(size = 16),
              axis.text.x = element_text(size = 14),
              axis.title.y = element_text(size = 16),
              axis.text.y = element_text(size = 14))
            
            ggplot(NprojStochPulsePrecSet$Nrange, aes(x=yrvec, y=Nmd)) +
              geom_ribbon(aes(ymin = Nlo, ymax = Nup), fill = "grey70", alpha=0.4) +
              geom_path(linetype=2) +
              geom_vline(xintercept=as.numeric(input$setOnePulse), linetype=3, color="black", size=0.5) +
              labs(x="years into future", y="N (♀ only)") +
              Ctheme
          })
          
          output$downloadPulseN <- downloadHandler(
            filename = function() {
              paste("pulseNproj", "csv", sep = ".")
            },
            
            content = function(file) {
              sep <- ","
              pulseNprojDL <<- NprojStochPulsePrecSet$Nrange
              pulseNprojDL[,2] <- round(NprojStochPulsePrecSet$Nrange[,2], 0)
              pulseNprojDL[,3] <- round(NprojStochPulsePrecSet$Nrange[,3], 0)
              pulseNprojDL[,4] <- round(NprojStochPulsePrecSet$Nrange[,4], 0)
              write.table(pulseNprojDL, file, sep=sep, row.names = F, col.names = T)
            }
          )
          
          reactiveVal({
            output$PrExtPulse <- renderText({
              paste("i. Pr(Ext) = ", round(NprojStochPulsePrecSet$PrQext, 4))
            })
          }) # end observe
          
          reactiveVal({
            output$minPopPulse <- renderText({
              paste("ii. min N = ", round(NprojStochPulsePrecSet$minMnN, 0), "(range: ", round(NprojStochPulsePrecSet$minLoN, 0), " to ", round(NprojStochPulsePrecSet$minUpN, 0),")")
            })
          }) # end observe
          
          reactiveVal({
            output$rMnPulse <- renderText({
              paste("iii. mean r = ", round(NprojStochPulsePrecSet$mnRmn, 2), "(range: ", round(NprojStochPulsePrecSet$mnRlo, 2), " to ", round(NprojStochPulsePrecSet$mnRup, 2),")")
            })
          }) # end observe
          
        } # end if
        
        
        if (input$percOrFix == "fix" & input$setOrRand == "random" & input$yrsOrGen == "gens") {
          
          output$projectPulsePlot <- renderPlot({
            
            NprojStochPulseFixRand <<- projStochPulse(x=Dmat, initN=input$Nstart/2, projGens=input$gensFutureProj, iter=as.numeric(input$iter),
                                                      S_SD=DemRSDdat[,1], F_SD=DemRSDdat[,2], Qthresh=input$Qthresh/2,
                                                      DFparams = c(input$DFa, input$DFb, input$DFc), 
                                                      CatParams = c(as.numeric(input$catMag), as.numeric(input$catMagSD)),
                                                      FixOff=as.numeric(input$fixPulse))
            
            Ctheme = theme(
              axis.title.x = element_text(size = 16),
              axis.text.x = element_text(size = 14),
              axis.title.y = element_text(size = 16),
              axis.text.y = element_text(size = 14))
            
            ggplot(NprojStochPulseFixRand$Nrange, aes(x=yrvec, y=Nmd)) +
              geom_ribbon(aes(ymin = Nlo, ymax = Nup), fill = "grey70", alpha=0.4) +
              geom_path(linetype=2) +
              labs(x="years into future", y="N (♀ only)") +
              Ctheme
          })
          
          output$downloadPulseN <- downloadHandler(
            filename = function() {
              paste("pulseNproj", "csv", sep = ".")
            },
            
            content = function(file) {
              sep <- ","
              pulseNprojDL <<- NprojStochPulseFixRand$Nrange
              pulseNprojDL[,2] <- round(NprojStochPulseFixRand$Nrange[,2], 0)
              pulseNprojDL[,3] <- round(NprojStochPulseFixRand$Nrange[,3], 0)
              pulseNprojDL[,4] <- round(NprojStochPulseFixRand$Nrange[,4], 0)
              write.table(pulseNprojDL, file, sep=sep, row.names = F, col.names = T)
            }
          )
          
          reactiveVal({
            output$PrExtPulse <- renderText({
              paste("i. Pr(Ext) = ", round(NprojStochPulseFixRand$PrQext, 4))
            })
          }) # end observe
          
          reactiveVal({
            output$minPopPulse <- renderText({
              paste("ii. min N = ", round(NprojStochPulseFixRand$minMnN, 0), "(range: ", round(NprojStochPulseFixRand$minLoN, 0), " to ", round(NprojStochPulseFixRand$minUpN, 0),")")
            })
          }) # end observe
          
          reactiveVal({
            output$rMnPulse <- renderText({
              paste("iii. mean r = ", round(NprojStochPulseFixRand$mnRmn, 2), "(range: ", round(NprojStochPulseFixRand$mnRlo, 2), " to ", round(NprojStochPulseFixRand$mnRup, 2),")")
            })
          }) # end observe
          
        } # end if
        
        if (input$percOrFix == "fix" & input$setOrRand == "set" & input$yrsOrGen == "gens") {
          
          output$projectPulsePlot <- renderPlot({
            
            NprojStochPulseFixSet <<- projStochPulse(x=Dmat, initN=input$Nstart/2, projGens=input$gensFutureProj, iter=as.numeric(input$iter),
                                                     S_SD=DemRSDdat[,1], F_SD=DemRSDdat[,2], Qthresh=input$Qthresh/2,
                                                     DFparams = c(input$DFa, input$DFb, input$DFc), 
                                                     CatParams = c(as.numeric(input$catMag), as.numeric(input$catMagSD)),
                                                     FixOff=as.numeric(input$fixPulse), TimeOff=as.numeric(input$setOnePulse))
            
            Ctheme = theme(
              axis.title.x = element_text(size = 16),
              axis.text.x = element_text(size = 14),
              axis.title.y = element_text(size = 16),
              axis.text.y = element_text(size = 14))
            
            ggplot(NprojStochPulseFixSet$Nrange, aes(x=yrvec, y=Nmd)) +
              geom_ribbon(aes(ymin = Nlo, ymax = Nup), fill = "grey70", alpha=0.4) +
              geom_path(linetype=2) +
              geom_vline(xintercept=as.numeric(input$setOnePulse), linetype=3, color="black", size=0.5) +
              labs(x="years into future", y="N (♀ only)") +
              Ctheme
          })
          
          output$downloadPulseN <- downloadHandler(
            filename = function() {
              paste("pulseNproj", "csv", sep = ".")
            },
            
            content = function(file) {
              sep <- ","
              pulseNprojDL <<- NprojStochPulseFixSet$Nrange
              pulseNprojDL[,2] <- round(NprojStochPulseFixSet$Nrange[,2], 0)
              pulseNprojDL[,3] <- round(NprojStochPulseFixSet$Nrange[,3], 0)
              pulseNprojDL[,4] <- round(NprojStochPulseFixSet$Nrange[,4], 0)
              write.table(pulseNprojDL, file, sep=sep, row.names = F, col.names = T)
            }
          )
          
          reactiveVal({
            output$PrExtPulse <- renderText({
              paste("i. Pr(Ext) = ", round(NprojStochPulseFixSet$PrQext, 4))
            })
          }) # end observe
          
          reactiveVal({
            output$minPopPulse <- renderText({
              paste("ii. min N = ", round(NprojStochPulseFixSet$minMnN, 0), "(range: ", round(NprojStochPulseFixSet$minLoN, 0), " to ", round(NprojStochPulseFixSet$minUpN, 0),")")
            })
          }) # end observe
          
          reactiveVal({
            output$rMnPulse <- renderText({
              paste("iii. mean r = ", round(NprojStochPulseFixSet$mnRmn, 2), "(range: ", round(NprojStochPulseFixSet$mnRlo, 2), " to ", round(NprojStochPulseFixSet$mnRup, 2),")")
            })
          }) # end observe
          
        } # end if
        
        
      }) # end observe
      
    } # end tab 6 if
    
  }) # end tab Events
  
  observeEvent(input$projectPress, {
    
    if(input$tabs == "tab7"){
      
      observe({
        ###################
        ## project by years
        ###################
        if (input$percOrFixPress == "percPress" & input$fullOrSubsetPress == "fullInt" & input$yrsOrGen == "years") {
          
          output$projectPressPlot <- renderPlot({
            
            NprojStochPressPrecRand <<- projStochPress(x=Dmat, initN=input$Nstart/2, projYrs=input$yrsFutureProj, iter=as.numeric(input$iter),
                                                       S_SD=DemRSDdat[,1], F_SD=DemRSDdat[,2], Qthresh=input$Qthresh/2,
                                                       DFparams = c(input$DFa, input$DFb, input$DFc), 
                                                       CatParams = c(as.numeric(input$catMag), as.numeric(input$catMagSD)),
                                                       PercOff=as.numeric(input$percentPress))
            
            Ctheme = theme(
              axis.title.x = element_text(size = 16),
              axis.text.x = element_text(size = 14),
              axis.title.y = element_text(size = 16),
              axis.text.y = element_text(size = 14))
            
            ggplot(NprojStochPressPrecRand$Nrange, aes(x=yrvec, y=Nmd)) +
              geom_ribbon(aes(ymin = Nlo, ymax = Nup), fill = "grey70", alpha=0.4) +
              geom_path(linetype=2) +
              labs(x="years into future", y="N (♀ only)") +
              Ctheme
          })
          
          output$downloadPressN <- downloadHandler(
            filename = function() {
              paste("pressNproj", "csv", sep = ".")
            },
            
            content = function(file) {
              sep <- ","
              pressNprojDL <<- NprojStochPressPrecRand$Nrange
              pressNprojDL[,2] <- round(NprojStochPressPrecRand$Nrange[,2], 0)
              pressNprojDL[,3] <- round(NprojStochPressPrecRand$Nrange[,3], 0)
              pressNprojDL[,4] <- round(NprojStochPressPrecRand$Nrange[,4], 0)
              write.table(pressNprojDL, file, sep=sep, row.names = F, col.names = T)
            }
          )
          
          reactiveVal({
            output$PrExtPress <- renderText({
              paste("i. Pr(Ext) = ", round(NprojStochPressPrecRand$PrQext, 4))
            })
          }) # end observe
          
          reactiveVal({
            output$minPopPress <- renderText({
              paste("ii. min N = ", round(NprojStochPressPrecRand$minMnN, 0), "(range: ", round(NprojStochPressPrecRand$minLoN, 0), " to ", round(NprojStochPressPrecRand$minUpN, 0),")")
            })
          }) # end observe
          
          reactiveVal({
            output$rMnPress <- renderText({
              paste("iii. mean r = ", round(NprojStochPressPrecRand$mnRmn, 2), "(range: ", round(NprojStochPressPrecRand$mnRlo, 2), " to ", round(NprojStochPressPrecRand$mnRup, 2),")")
            })
          }) # end observe
          
        } # end if
        
        if (input$percOrFixPress == "percPress" & input$fullOrSubsetPress == "subsetInt" & input$yrsOrGen == "years") {
          
          output$projectPressPlot <- renderPlot({
            
            NprojStochPressPrecSet <<- projStochPress(x=Dmat, initN=input$Nstart/2, projYrs=input$yrsFutureProj, iter=as.numeric(input$iter),
                                                      S_SD=DemRSDdat[,1], F_SD=DemRSDdat[,2], Qthresh=input$Qthresh/2,
                                                      DFparams = c(input$DFa, input$DFb, input$DFc), 
                                                      CatParams = c(as.numeric(input$catMag), as.numeric(input$catMagSD)),
                                                      PercOff=as.numeric(input$percentPress), IntOff=as.numeric(input$setSubsetPress))
            
            Ctheme = theme(
              axis.title.x = element_text(size = 16),
              axis.text.x = element_text(size = 14),
              axis.title.y = element_text(size = 16),
              axis.text.y = element_text(size = 14))
            
            ggplot(NprojStochPressPrecSet$Nrange, aes(x=yrvec, y=Nmd)) +
              geom_ribbon(aes(ymin = Nlo, ymax = Nup), fill = "grey70", alpha=0.4) +
              geom_path(linetype=2) +
              geom_vline(xintercept=c(round(as.numeric(input$yrsFutureProj)*as.numeric(input$setSubsetPress[1]/100), 0), round(as.numeric(input$yrsFutureProj)*as.numeric(input$setSubsetPress[2]/100), 0)),
                         linetype=3, color="black", size=0.5) +
              labs(x="years into future", y="N (♀ only)") +
              Ctheme
          })
          
          output$downloadPressN <- downloadHandler(
            filename = function() {
              paste("pressNproj", "csv", sep = ".")
            },
            
            content = function(file) {
              sep <- ","
              pressNprojDL <<- NprojStochPressPrecSet$Nrange
              pressNprojDL[,2] <- round(NprojStochPressPrecSet$Nrange[,2], 0)
              pressNprojDL[,3] <- round(NprojStochPressPrecSet$Nrange[,3], 0)
              pressNprojDL[,4] <- round(NprojStochPressPrecSet$Nrange[,4], 0)
              write.table(pressNprojDL, file, sep=sep, row.names = F, col.names = T)
            }
          )
          
          reactiveVal({
            output$PrExtPress <- renderText({
              paste("i. Pr(Ext) = ", round(NprojStochPressPrecSet$PrQext, 4))
            })
          }) # end observe
          
          reactiveVal({
            output$minPopPress <- renderText({
              paste("ii. min N = ", round(NprojStochPressPrecSet$minMnN, 0), "(range: ", round(NprojStochPressPrecSet$minLoN, 0), " to ", round(NprojStochPressPrecSet$minUpN, 0),")")
            })
          }) # end observe
          
          reactiveVal({
            output$rMnPress <- renderText({
              paste("iii. mean r = ", round(NprojStochPressPrecSet$mnRmn, 2), "(range: ", round(NprojStochPressPrecSet$mnRlo, 2), " to ", round(NprojStochPressPrecSet$mnRup, 2),")")
            })
          }) # end observe
          
        } # end if
        
        
        if (input$percOrFixPress == "fixPress" & input$fullOrSubsetPress == "fullInt" & input$yrsOrGen == "years") {
          
          output$projectPressPlot <- renderPlot({
            
            NprojStochPressFixRand <<- projStochPress(x=Dmat, initN=input$Nstart/2, projYrs=input$yrsFutureProj, iter=as.numeric(input$iter),
                                                      S_SD=DemRSDdat[,1], F_SD=DemRSDdat[,2], Qthresh=input$Qthresh/2,
                                                      DFparams = c(input$DFa, input$DFb, input$DFc), 
                                                      CatParams = c(as.numeric(input$catMag), as.numeric(input$catMagSD)),
                                                      FixOff=as.numeric(input$fixedPress))
            
            Ctheme = theme(
              axis.title.x = element_text(size = 16),
              axis.text.x = element_text(size = 14),
              axis.title.y = element_text(size = 16),
              axis.text.y = element_text(size = 14))
            
            ggplot(NprojStochPressFixRand$Nrange, aes(x=yrvec, y=Nmd)) +
              geom_ribbon(aes(ymin = Nlo, ymax = Nup), fill = "grey70", alpha=0.4) +
              geom_path(linetype=2) +
              labs(x="years into future", y="N (♀ only)") +
              Ctheme
          })
          
          output$downloadPressN <- downloadHandler(
            filename = function() {
              paste("pressNproj", "csv", sep = ".")
            },
            
            content = function(file) {
              sep <- ","
              pressNprojDL <<- NprojStochPressFixRand$Nrange
              pressNprojDL[,2] <- round(NprojStochPressFixRand$Nrange[,2], 0)
              pressNprojDL[,3] <- round(NprojStochPressFixRand$Nrange[,3], 0)
              pressNprojDL[,4] <- round(NprojStochPressFixRand$Nrange[,4], 0)
              write.table(pressNprojDL, file, sep=sep, row.names = F, col.names = T)
            }
          )
          
          reactiveVal({
            output$PrExtPress <- renderText({
              paste("i. Pr(Ext) = ", round(NprojStochPressFixRand$PrQext, 4))
            })
          }) # end observe
          
          reactiveVal({
            output$minPopPress <- renderText({
              paste("ii. min N = ", round(NprojStochPressFixRand$minMnN, 0), "(range: ", round(NprojStochPressFixRand$minLoN, 0), " to ", round(NprojStochPressFixRand$minUpN, 0),")")
            })
          }) # end observe
          
          reactiveVal({
            output$rMnPress <- renderText({
              paste("iii. mean r = ", round(NprojStochPressFixRand$mnRmn, 2), "(range: ", round(NprojStochPressFixRand$mnRlo, 2), " to ", round(NprojStochPressFixRand$mnRup, 2),")")
            })
          }) # end observe
          
        } # end if
        
        if (input$percOrFixPress == "fixPress" & input$fullOrSubsetPress == "subsetInt" & input$yrsOrGen == "years") {
          
          output$projectPressPlot <- renderPlot({
            
            NprojStochPressFixSet <<- projStochPress(x=Dmat, initN=input$Nstart/2, projYrs=input$yrsFutureProj, iter=as.numeric(input$iter),
                                                     S_SD=DemRSDdat[,1], F_SD=DemRSDdat[,2], Qthresh=input$Qthresh/2,
                                                     DFparams = c(input$DFa, input$DFb, input$DFc), 
                                                     CatParams = c(as.numeric(input$catMag), as.numeric(input$catMagSD)),
                                                     FixOff=as.numeric(input$fixedPress), IntOff=as.numeric(input$setSubsetPress))
            
            Ctheme = theme(
              axis.title.x = element_text(size = 16),
              axis.text.x = element_text(size = 14),
              axis.title.y = element_text(size = 16),
              axis.text.y = element_text(size = 14))
            
            ggplot(NprojStochPressFixSet$Nrange, aes(x=yrvec, y=Nmd)) +
              geom_ribbon(aes(ymin = Nlo, ymax = Nup), fill = "grey70", alpha=0.4) +
              geom_path(linetype=2) +
              geom_vline(xintercept=c(round(as.numeric(input$yrsFutureProj)*as.numeric(input$setSubsetPress[1]/100), 0), round(as.numeric(input$yrsFutureProj)*as.numeric(input$setSubsetPress[2]/100), 0)),
                         linetype=3, color="black", size=0.5) +
              labs(x="years into future", y="N (♀ only)") +
              Ctheme
          })
          
          output$downloadPressN <- downloadHandler(
            filename = function() {
              paste("pressNproj", "csv", sep = ".")
            },
            
            content = function(file) {
              sep <- ","
              pressNprojDL <<- NprojStochPressFixSet$Nrange
              pressNprojDL[,2] <- round(NprojStochPressFixSet$Nrange[,2], 0)
              pressNprojDL[,3] <- round(NprojStochPressFixSet$Nrange[,3], 0)
              pressNprojDL[,4] <- round(NprojStochPressFixSet$Nrange[,4], 0)
              write.table(pressNprojDL, file, sep=sep, row.names = F, col.names = T)
            }
          )
          
          reactiveVal({
            output$PrExtPress <- renderText({
              paste("i. Pr(Ext) = ", round(NprojStochPressFixSet$PrQext, 4))
            })
          }) # end observe
          
          reactiveVal({
            output$minPopPress <- renderText({
              paste("ii. min N = ", round(NprojStochPressFixSet$minMnN, 0), "(range: ", round(NprojStochPressFixSet$minLoN, 0), " to ", round(NprojStochPressFixSet$minUpN, 0),")")
            })
          }) # end observe
          
          reactiveVal({
            output$rMnPress <- renderText({
              paste("iii. mean r = ", round(NprojStochPressFixSet$mnRmn, 2), "(range: ", round(NprojStochPressFixSet$mnRlo, 2), " to ", round(NprojStochPressFixSet$mnRup, 2),")")
            })
          }) # end observe
          
        } # end if
        
        #########################
        ## project by generations
        #########################
        
        if (input$percOrFixPress == "percPress" & input$fullOrSubsetPress == "fullInt" & input$yrsOrGen == "gens") {
          
          output$projectPressPlot <- renderPlot({
            
            NprojStochPressPrecRand <<- projStochPress(x=Dmat, initN=input$Nstart/2, projGens=input$gensFutureProj, iter=as.numeric(input$iter),
                                                       S_SD=DemRSDdat[,1], F_SD=DemRSDdat[,2], Qthresh=input$Qthresh/2,
                                                       DFparams = c(input$DFa, input$DFb, input$DFc), 
                                                       CatParams = c(as.numeric(input$catMag), as.numeric(input$catMagSD)),
                                                       PercOff=as.numeric(input$percentPress))
            
            Ctheme = theme(
              axis.title.x = element_text(size = 16),
              axis.text.x = element_text(size = 14),
              axis.title.y = element_text(size = 16),
              axis.text.y = element_text(size = 14))
            
            ggplot(NprojStochPressPrecRand$Nrange, aes(x=yrvec, y=Nmd)) +
              geom_ribbon(aes(ymin = Nlo, ymax = Nup), fill = "grey70", alpha=0.4) +
              geom_path(linetype=2) +
              labs(x="years into future", y="N (♀ only)") +
              Ctheme
          })
          
          output$downloadPressN <- downloadHandler(
            filename = function() {
              paste("pressNproj", "csv", sep = ".")
            },
            
            content = function(file) {
              sep <- ","
              pressNprojDL <<- NprojStochPressPrecRand$Nrange
              pressNprojDL[,2] <- round(NprojStochPressPrecRand$Nrange[,2], 0)
              pressNprojDL[,3] <- round(NprojStochPressPrecRand$Nrange[,3], 0)
              pressNprojDL[,4] <- round(NprojStochPressPrecRand$Nrange[,4], 0)
              write.table(pressNprojDL, file, sep=sep, row.names = F, col.names = T)
            }
          )
          
          reactiveVal({
            output$PrExtPress <- renderText({
              paste("i. Pr(Ext) = ", round(NprojStochPressPrecRand$PrQext, 4))
            })
          }) # end observe
          
          reactiveVal({
            output$minPopPress <- renderText({
              paste("ii. min N = ", round(NprojStochPressPrecRand$minMnN, 0), "(range: ", round(NprojStochPressPrecRand$minLoN, 0), " to ", round(NprojStochPressPrecRand$minUpN, 0),")")
            })
          }) # end observe
          
          reactiveVal({
            output$rMnPress <- renderText({
              paste("iii. mean r = ", round(NprojStochPressPrecRand$mnRmn, 2), "(range: ", round(NprojStochPressPrecRand$mnRlo, 2), " to ", round(NprojStochPressPrecRand$mnRup, 2),")")
            })
          }) # end observe
          
        } # end if
        
        if (input$percOrFixPress == "percPress" & input$fullOrSubsetPress == "subsetInt" & input$yrsOrGen == "gens") {
          
          output$projectPressPlot <- renderPlot({
            
            NprojStochPressPrecSet <<- projStochPress(x=Dmat, initN=input$Nstart/2, projGens=input$gensFutureProj, iter=as.numeric(input$iter),
                                                      S_SD=DemRSDdat[,1], F_SD=DemRSDdat[,2], Qthresh=input$Qthresh/2,
                                                      DFparams = c(input$DFa, input$DFb, input$DFc), 
                                                      CatParams = c(as.numeric(input$catMag), as.numeric(input$catMagSD)),
                                                      PercOff=as.numeric(input$percentPress), IntOff=as.numeric(input$setSubsetPress))
            
            Ctheme = theme(
              axis.title.x = element_text(size = 16),
              axis.text.x = element_text(size = 14),
              axis.title.y = element_text(size = 16),
              axis.text.y = element_text(size = 14))
            
            ggplot(NprojStochPressPrecSet$Nrange, aes(x=yrvec, y=Nmd)) +
              geom_ribbon(aes(ymin = Nlo, ymax = Nup), fill = "grey70", alpha=0.4) +
              geom_path(linetype=2) +
              geom_vline(xintercept=c(round(as.numeric(input$gensFutureProj)*G*as.numeric(input$setSubsetPress[1]/100), 0), round(as.numeric(input$gensFutureProj)*G*as.numeric(input$setSubsetPress[2]/100), 0)),
                         linetype=3, color="black", size=0.5) +
              labs(x="years into future", y="N (♀ only)") +
              Ctheme
          })
          
          output$downloadPressN <- downloadHandler(
            filename = function() {
              paste("pressNproj", "csv", sep = ".")
            },
            
            content = function(file) {
              sep <- ","
              pressNprojDL <<- NprojStochPressPrecSet$Nrange
              pressNprojDL[,2] <- round(NprojStochPressPrecSet$Nrange[,2], 0)
              pressNprojDL[,3] <- round(NprojStochPressPrecSet$Nrange[,3], 0)
              pressNprojDL[,4] <- round(NprojStochPressPrecSet$Nrange[,4], 0)
              write.table(pressNprojDL, file, sep=sep, row.names = F, col.names = T)
            }
          )
          
          reactiveVal({
            output$PrExtPress <- renderText({
              paste("i. Pr(Ext) = ", round(NprojStochPressPrecSet$PrQext, 4))
            })
          }) # end observe
          
          reactiveVal({
            output$minPopPress <- renderText({
              paste("ii. min N = ", round(NprojStochPressPrecSet$minMnN, 0), "(range: ", round(NprojStochPressPrecSet$minLoN, 0), " to ", round(NprojStochPressPrecSet$minUpN, 0),")")
            })
          }) # end observe
          
          reactiveVal({
            output$rMnPress <- renderText({
              paste("iii. mean r = ", round(NprojStochPressPrecSet$mnRmn, 2), "(range: ", round(NprojStochPressPrecSet$mnRlo, 2), " to ", round(NprojStochPressPrecSet$mnRup, 2),")")
            })
          }) # end observe
          
        } # end if
        
        
        if (input$percOrFixPress == "fixPress" & input$fullOrSubsetPress == "fullInt" & input$yrsOrGen == "gens") {
          
          output$projectPressPlot <- renderPlot({
            
            NprojStochPressFixRand <<- projStochPress(x=Dmat, initN=input$Nstart/2, projGens=input$gensFutureProj, iter=as.numeric(input$iter),
                                                      S_SD=DemRSDdat[,1], F_SD=DemRSDdat[,2], Qthresh=input$Qthresh/2,
                                                      DFparams = c(input$DFa, input$DFb, input$DFc), 
                                                      CatParams = c(as.numeric(input$catMag), as.numeric(input$catMagSD)),
                                                      FixOff=as.numeric(input$fixedPress))
            
            Ctheme = theme(
              axis.title.x = element_text(size = 16),
              axis.text.x = element_text(size = 14),
              axis.title.y = element_text(size = 16),
              axis.text.y = element_text(size = 14))
            
            ggplot(NprojStochPressFixRand$Nrange, aes(x=yrvec, y=Nmd)) +
              geom_ribbon(aes(ymin = Nlo, ymax = Nup), fill = "grey70", alpha=0.4) +
              geom_path(linetype=2) +
              labs(x="years into future", y="N (♀ only)") +
              Ctheme
          })
          
          output$downloadPressN <- downloadHandler(
            filename = function() {
              paste("pressNproj", "csv", sep = ".")
            },
            
            content = function(file) {
              sep <- ","
              pressNprojDL <<- NprojStochPressFixRand$Nrange
              pressNprojDL[,2] <- round(NprojStochPressFixRand$Nrange[,2], 0)
              pressNprojDL[,3] <- round(NprojStochPressFixRand$Nrange[,3], 0)
              pressNprojDL[,4] <- round(NprojStochPressFixRand$Nrange[,4], 0)
              write.table(pressNprojDL, file, sep=sep, row.names = F, col.names = T)
            }
          )
          
          reactiveVal({
            output$PrExtPress <- renderText({
              paste("i. Pr(Ext) = ", round(NprojStochPressFixRand$PrQext, 4))
            })
          }) # end observe
          
          reactiveVal({
            output$minPopPress <- renderText({
              paste("ii. min N = ", round(NprojStochPressFixRand$minMnN, 0), "(range: ", round(NprojStochPressFixRand$minLoN, 0), " to ", round(NprojStochPressFixRand$minUpN, 0),")")
            })
          }) # end observe
          
          reactiveVal({
            output$rMnPress <- renderText({
              paste("iii. mean r = ", round(NprojStochPressFixRand$mnRmn, 2), "(range: ", round(NprojStochPressFixRand$mnRlo, 2), " to ", round(NprojStochPressFixRand$mnRup, 2),")")
            })
          }) # end observe
          
        } # end if
        
        if (input$percOrFixPress == "fixPress" & input$fullOrSubsetPress == "subsetInt" & input$yrsOrGen == "gens") {
          
          output$projectPressPlot <- renderPlot({
            
            NprojStochPressFixSet <<- projStochPress(x=Dmat, initN=input$Nstart/2, projGens=input$gensFutureProj, iter=as.numeric(input$iter),
                                                     S_SD=DemRSDdat[,1], F_SD=DemRSDdat[,2], Qthresh=input$Qthresh/2,
                                                     DFparams = c(input$DFa, input$DFb, input$DFc), 
                                                     CatParams = c(as.numeric(input$catMag), as.numeric(input$catMagSD)),
                                                     FixOff=as.numeric(input$fixedPress), IntOff=input$setSubsetPress)
            
            Ctheme = theme(
              axis.title.x = element_text(size = 16),
              axis.text.x = element_text(size = 14),
              axis.title.y = element_text(size = 16),
              axis.text.y = element_text(size = 14))
            
            ggplot(NprojStochPressFixSet$Nrange, aes(x=yrvec, y=Nmd)) +
              geom_ribbon(aes(ymin = Nlo, ymax = Nup), fill = "grey70", alpha=0.4) +
              geom_path(linetype=2) +
              geom_vline(xintercept=c(round(as.numeric(input$gensFutureProj)*G*as.numeric(input$setSubsetPress[1]/100), 0), round(as.numeric(input$gensFutureProj)*G*as.numeric(input$setSubsetPress[2]/100), 0)),
                         linetype=3, color="black", size=0.5) +
              labs(x="years into future", y="N (♀ only)") +
              Ctheme
          })
          
          output$downloadPressN <- downloadHandler(
            filename = function() {
              paste("pressNproj", "csv", sep = ".")
            },
            
            content = function(file) {
              sep <- ","
              pressNprojDL <<- NprojStochPressFixSet$Nrange
              pressNprojDL[,2] <- round(NprojStochPressFixSet$Nrange[,2], 0)
              pressNprojDL[,3] <- round(NprojStochPressFixSet$Nrange[,3], 0)
              pressNprojDL[,4] <- round(NprojStochPressFixSet$Nrange[,4], 0)
              write.table(pressNprojDL, file, sep=sep, row.names = F, col.names = T)
            }
          )
          
          reactiveVal({
            output$PrExtPress <- renderText({
              paste("i. Pr(Ext) = ", round(NprojStochPressFixSet$PrQext, 4))
            })
          }) # end observe
          
          reactiveVal({
            output$minPopPress <- renderText({
              paste("ii. min N = ", round(NprojStochPressFixSet$minMnN, 0), "(range: ", round(NprojStochPressFixSet$minLoN, 0), " to ", round(NprojStochPressFixSet$minUpN, 0),")")
            })
          }) # end observe
          
          reactiveVal({
            output$rMnPress <- renderText({
              paste("iii. mean r = ", round(NprojStochPressFixSet$mnRmn, 2), "(range: ", round(NprojStochPressFixSet$mnRlo, 2), " to ", round(NprojStochPressFixSet$mnRup, 2),")")
            })
          }) # end observe
          
        } # end if
        
      }) # end observe
      
    } # end tab 7 if
    
  }) # end tab Events
  
  observeEvent(input$calcMVP, {
    
    if(input$tabs == "tab8"){
      
      observe({
        ############################
        ## project by generations ##
        ############################
        if (input$genOrYrs == "gen") {
          
            output$MVPPlot <- renderPlot({
              
              withProgress(message = "MVP calculation in progress",
                           detail = "(this can take a while ...)",
                           value = 0, {
                           MVPStochGen <<- StochMVP(x=Dmat, initN=input$Nhigh/2, projGens=input$genProj, maxPrQe=(1-input$persistPr), Nint=input$Nstep, Nmin=input$Nlo,
                                       iter=as.numeric(input$iterMVP),S_SD=DemRSDdat[,1], F_SD=DemRSDdat[,2], Qthresh=input$Qthresh/2,
                                                         DFparams = c(input$DFa, input$DFb, input$DFc), 
                                                         CatParams = c(as.numeric(input$catMag), as.numeric(input$catMagSD)))
                           
                           MVPseries <<- MVPStochGen$MVPout[order(MVPStochGen$MVPout$Ninit.vec),]
                           MVPasympOut <<- LogPowFunc(PrPersist=(1-MVPseries[,2]), N=MVPseries[,1], a=0.95, b=MVPStochGen$MVPst,
                                                      c=-17, alpha=0.05, sim=10000)
                           incProgress(1/as.numeric(input$iterMVP))
                           
                           }) # end Progress

              Ctheme = theme(
                axis.title.x = element_text(size = 16),
                axis.text.x = element_text(size = 14),
                axis.title.y = element_text(size = 16),
                axis.text.y = element_text(size = 14))
              
              if(is.na(MVPasympOut$PerPredDat$Nvec[1])==T) {
                ggplot() +
                  geom_path(data=MVPseries, aes(x=Ninit.vec, y=1-Qext)) +
                  geom_vline(xintercept=MVPStochGen$MVP, linetype=2, color="red", size=0.5) +
                  geom_vline(xintercept=MVPStochGen$MVPst, linetype=2, color="red", size=0.5) +
                  geom_hline(yintercept=(input$persistPr), linetype=2, color="red", size=0.5) +			  
                  labs(x="N initial (♀ only)", y="Pr(persistence)") +
                  Ctheme
                
              } else {
                ggplot() +
                  geom_path(data=MVPseries, aes(x=Ninit.vec, y=1-Qext)) +
                  geom_line(data=MVPasympOut$PerPredDat, aes(x=Nvec, y=PerPredMed), col="blue", linetype=2) +
                  geom_ribbon(data=MVPasympOut$PerPredDat, aes(x=Nvec, ymin=PerPredLo, ymax=PerPredUp), alpha=0.1) +
                  annotate("rect", xmin=MVPasympOut$MPVasLo, xmax=MVPasympOut$MVPasUp, ymin=0, ymax=1, fill="red", alpha=0.05) +
                  geom_vline(xintercept=MVPStochGen$MVP, linetype=2, color="red", size=0.5) +
                  geom_vline(xintercept=MVPStochGen$MVPst, linetype=2, color="red", size=0.5) +
                  geom_hline(yintercept=(input$persistPr), linetype=2, color="red", size=0.5) +			  
                  labs(x="N initial (♀ only)", y="Pr(persistence)") +
                  Ctheme
              }
            
              })
       
          reactiveVal({
            output$MVPest <- renderText({
              if (is.na(MVPStochGen$MVP) == T || MVPStochGen$MVP == input$Nhigh/2 || length(MVPStochGen$MVP) == 0) {
                "i. persistence probability not achieved within specified parameter range; lower persistence probability or increase upper initial N"
              } else {
                paste("i. MVP = ", round(MVPStochGen$MVP, 0), " ♀")
              } # end if/else
            })
          }) # end observe
          
          reactiveVal({
            output$MVPstepChange <- renderText({
              if (is.na(MVPStochGen$MVPst) == F) {
                paste("ii. MVP breakpoint = ", round(MVPStochGen$MVPst, 0), " ♀")
              } # end if
            })
          }) # end observe
          
          reactiveVal({
            output$MVPasymptote <- renderText({
              if (is.na(MVPasympOut$MVPasMed) == F & is.na(MVPasympOut$PerPredDat$Nvec[1])==F) {
                paste("iii. MVP asymptotic = ", MVPasympOut$MVPasMed, " ♀ (CI: ", min(c(MVPasympOut$MPVasLo,MVPasympOut$MVPasUp)), " — ",
                      max(c(MVPasympOut$MPVasLo,MVPasympOut$MVPasUp)),
                      ") ➪➪ fit: Pr(Persist) = ", round(MVPasympOut$aLP, 3), " / (1 + (N / ", round(MVPasympOut$bLP, 0), ")^",
                      round(MVPasympOut$cLP, 1), sep="")
              } else {
                "iii. Cannot estimate aymptotic MVP function with set parameters"
                } # end if/else
            })
          }) # end observe
          
          
        } # end if
        
        ######################
        ## project by years ##
        ######################
        
        if (input$genOrYrs == "yr") {
          
          output$MVPPlot <- renderPlot({
            withProgress(message = "MVP calculation in progress",
                         detail = "(this can take a while ...)",
                         value = 0, {
                           
                         MVPStochYr <<- StochMVP(x=Dmat, initN=input$Nhigh/2, projYrs=input$yrProj, maxPrQe=(1-input$persistPr), Nint=input$Nstep, Nmin=input$Nlo,
                                     iter=as.numeric(input$iterMVP),S_SD=DemRSDdat[,1], F_SD=DemRSDdat[,2], Qthresh=input$Qthresh/2,
                                     DFparams = c(input$DFa, input$DFb, input$DFc), 
                                     CatParams = c(as.numeric(input$catMag), as.numeric(input$catMagSD)))
                         
                         MVPseries <<- MVPStochYr$MVPout[order(MVPStochYr$MVPout$Ninit.vec),]
                         MVPasympOut <<- LogPowFunc(PrPersist=(1-MVPseries[,2]), N=MVPseries[,1], a=0.95, b=MVPStochYr$MVPst,
                                                    c=-17, alpha=0.05, sim=10000)
                         incProgress(1/as.numeric(input$iterMVP))
                         }) # end Progress
            
            Ctheme = theme(
              axis.title.x = element_text(size = 16),
              axis.text.x = element_text(size = 14),
              axis.title.y = element_text(size = 16),
              axis.text.y = element_text(size = 14))
            
            if(is.na(MVPasympOut$PerPredDat$Nvec[1])==T) {
              ggplot() +
                geom_path(data=MVPseries, aes(x=Ninit.vec, y=1-Qext)) +
                geom_vline(xintercept=MVPStochYr$MVP, linetype=2, color="red", size=0.5) +
                geom_vline(xintercept=MVPStochYr$MVPst, linetype=2, color="red", size=0.5) +
                geom_hline(yintercept=(input$persistPr), linetype=2, color="red", size=0.5) +			  
                labs(x="N initial (♀ only)", y="Pr(persistence)") +
                Ctheme
              
            } else {
              ggplot() +
                geom_path(data=MVPseries, aes(x=Ninit.vec, y=1-Qext)) +
                geom_line(data=MVPasympOut$PerPredDat, aes(x=Nvec, y=PerPredMed), col="blue", linetype=2) +
                geom_ribbon(data=MVPasympOut$PerPredDat, aes(x=Nvec, ymin=PerPredLo, ymax=PerPredUp), alpha=0.1) +
                annotate("rect", xmin=MVPasympOut$MPVasLo, xmax=MVPasympOut$MVPasUp, ymin=0, ymax=1, fill="red", alpha=0.05) +
                geom_vline(xintercept=MVPStochYr$MVP, linetype=2, color="red", size=0.5) +
                geom_vline(xintercept=MVPStochYr$MVPst, linetype=2, color="red", size=0.5) +
                geom_hline(yintercept=(input$persistPr), linetype=2, color="red", size=0.5) +			  
                labs(x="N initial (♀ only)", y="Pr(persistence)") +
                Ctheme
            }
          })
          
          reactiveVal({
            output$MVPest <- renderText({
              if (is.na(MVPStochYr$MVP) == T || MVPStochYr$MVP == input$Nhigh/2 || length(MVPStochYr$MVP)==0) {
                "i. persistence probability not achieved within specified parameter range; lower persistence probability or increase upper initial N"
              } else {
                 paste("i. MVP = ", round(MVPStochYr$MVP, 0), "total ♀")
              } # end if/else
            })
          }) # end observe
          
          reactiveVal({
            output$MVPstepChange <- renderText({
              if (is.na(MVPStochYr$MVPst) == F) {
                paste("ii. MVP breakpoint = ", round(MVPStochYr$MVPst, 0), " total ♀")
              } # end if
            })
          }) # end observe
          
          reactiveVal({
            output$MVPasymptote <- renderText({
              if (is.na(MVPasympOut$MVPasMed) == F & is.na(MVPasympOut$PerPredDat$Nvec[1])==F) {
                paste("iii. MVP asymptotic = ", MVPasympOut$MVPasMed, " ♀ (CI: ", min(c(MVPasympOut$MPVasLo,MVPasympOut$MVPasUp)), " — ",
                      max(c(MVPasympOut$MPVasLo,MVPasympOut$MVPasUp)),
                      ") ➪➪ fit: Pr(Persist) = ", round(MVPasympOut$aLP, 3), " / (1 + (N / ", round(MVPasympOut$bLP, 0), ")^",
                      round(MVPasympOut$cLP, 1), sep="")
              } else {
                "iii. Cannot estimate aymptotic MVP function with set parameters"
              } # end if/else
            })
          }) # end observe
          
        } # end if
        
        
      }) # end observe
      
    } # end tab 8 if
    
  }) # end tab Events
  
  
  session$onSessionEnded(stopApp)
  
} # end server

shinyApp(ui, server)
