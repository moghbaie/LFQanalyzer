
library(shiny)
library(shinyFiles)
library(bslib)
library(rintrojs)
library(PTXQC)
library(foreach)
library(DT)
library(plotly)

source("src/functions.R",local=TRUE)
source("src/class.R",local=TRUE)

#thematic::thematic_shiny(font = "auto")


ui <- fluidPage(
  
  introjsUI(),
  sidebarLayout(
    sidebarPanel(
      
      selectizeInput("software", "software",
                     choices = c("Maxquant", "Proteom Discoverer"), 
                     multiple = FALSE),
      
      htmlOutput("intensity"), 
      shinyFilesButton("file", "File select", "Please select a file", multiple = TRUE, viewtype = "detail",
                       data.step = 2,
                       data.intro = "This is the sidebar. Look how intro elements can nest"), 
      htmlOutput("import"), 
      htmlOutput("quality"), 
      htmlOutput("select"), 
      htmlOutput("add"), 
      htmlOutput("del"), 
      htmlOutput("submit")
    ),
    mainPanel(
      
      navbarPage("", id = "tabs",
                 tabPanel("Selection", 
                          # verbatimTextOutput("filepaths"),
                          htmlOutput("list"),
                          htmlOutput("software")), 
                 tabPanel("Quality report",htmlOutput("qualityreport")),
                 tabPanel("Volcano Plot", 
                          htmlOutput("comparison"),
                          fluidRow( div(style = "height:60%;",plotlyOutput("volcano")),
                                    div(style = "height:30%;",plotlyOutput("histogram")))
                            
                         )
      )
      
    )
  ), title = 'Options groups for select(ize) input')


server <- function(input, output, session) {
  
  volumes <- c(Home = fs::path_home(), "R Installation" = R.home(), getVolumes()())
  
  # note we use the syntax input[['foo']] instead of input$foo, because we have
  # to construct the id as a character string, then use it to access the value;
  # same thing applies to the output object below
  
  metadata <- reactiveValues()
  output[["a_out"]] <- renderPrint({
    
    res <- foreach(j = 1:loop(), .combine="rbind") %do% {
      
      if(length(input[[paste0('b', j)]])>0 &length(input[[paste0('a', j)]])>0){
        rbind(
          cbind(input[[paste0('b', j)]], rep("case", length(input[[paste0('b', j)]])), rep(input[[paste0('c', j)]],length(input[[paste0('b', j)]]))),
          cbind(input[[paste0('a', j)]], rep("control", length(input[[paste0('a', j)]])), rep(input[[paste0('c', j)]],length(input[[paste0('a', j)]])))
        )
      }
      else if(length(input[[paste0('b', j)]])>0 &length(input[[paste0('a', j)]])==0){
        cbind(input[[paste0('b', j)]], rep("case", length(input[[paste0('b', j)]])), rep(input[[paste0('c', j)]],length(input[[paste0('b', j)]])))
        
      } 
      else if(length(input[[paste0('b', j)]])==0 &length(input[[paste0('a', j)]])>0){
        cbind(input[[paste0('a', j)]], rep("case", length(input[[paste0('a', j)]])), rep(input[[paste0('c', j)]],length(input[[paste0('a', j)]])))
        
      } else{
        matrix(NA, ncol=3, nrow=0)
      }
      
    }
    
    res <- data.frame(res)
    colnames(res) <- c("Condition","Type","Comparison")
    metadata[["targetfile"]] <- res
    
    
    return(res)   
  })
  
  
  
  loop <- reactiveVal(1) 
  lapply(1:2, function(i) {
    output[[paste0('a', i)]] <- renderUI({
      strong(paste(unname (unlist(input[[paste0('a', i)]])), collapse=','))
    })
  })
  
  listen0 <- reactive({
    input[["import"]]
  })
  
  listen1 <- reactive({
    input[["add"]]
  })
  
  listen2 <- reactive({
    input[['del']]
  })
  
  listen3 <- reactive({
    input[['quality']]
  })
  
  
  listen4 <- reactive({
    input[['comparison']]
  })
  
  observeEvent(listen1(),{ nloop <- loop() +1 
  print(nloop)
  loop(nloop)
  })
  
  
  observeEvent(listen2(),{
    
    if(loop()>1){
      nloop <- loop() -1 
    }
    print(nloop)
    loop(nloop)
  })
  
  output$select <- renderUI({
    if(length(input$import)>0){
      if((!is.integer(input$file))& input$import==1){
        lapply(1:loop(), function(i) {
          
          
          div(style="display:inline-block",
              selectizeInput(paste0('b', i), paste0('Select  Case ', i),
                             choices =  getcolumns(), multiple = TRUE,
                             selected = lapply(reactiveValuesToList(input), unclass)[[paste0('b', i)]]),
              selectizeInput(paste0('a', i), paste0('Select  Control ', i),
                             choices =  getcolumns(), multiple = TRUE, 
                             selected = lapply(reactiveValuesToList(input), unclass)[[paste0('a', i)]]),
              textInput(paste0('c', i), label = paste0('Comparison ', i), 
                            value= lapply(reactiveValuesToList(input), unclass)[[paste0('c', i)]]), sep="\n")
        })
      } 
    }
    
  })
  
  output$intensity <- renderUI({
    if(input$software=="Maxquant"){
      selectizeInput("intensity", "Intensity",
                     choices = c("LFQ", "iBAQ"), 
                     multiple = FALSE)
    }
  })
  
  
  
  output$import <- renderUI({
    if(!is.integer(input$file)){
      actionButton("import", "import")
    }
    
  })
  
  
  output$quality <- renderUI({
    if(input$software == "Maxquant" & length(input$import)>0 ){
      if(input$import==1){
        actionButton("quality", "Quality control")
      }
    }
    
  })
  
  
  output$add <- renderUI({
    
    if(length(input$import)>0){
      if((!is.integer(input$file))& input$import==1){
        actionButton("add", "add more comparison")
      }
    }
    
  })
  
  output$del <- renderUI({
    if(loop()>1){
      actionButton("del", "del")
    }
    
  })
  
  
  output$list <- renderUI ({
    if(!is.integer(input$file)){
      verbatimTextOutput('a_out')
    }
  })
  
  
  output$submit <- renderUI({
    
    if(length(input$import)>0){
      if((!is.integer(input$file))& input$import==1){
        actionButton("submit", "Submit")
      }
    }
    
  })
  
  shinyFileChoose(input, "file", roots = volumes, session = session)
  
  observe({
    cat("\ninput$file value:\n\n")
    
    print(file.path(volumes[unlist(input$file[2])], 
                    paste0(c(unlist(input$file[[1]])[-1]),collapse= "//")
    ))
    
  })
  
  ## print to browser
  output$filepaths <- renderPrint({
    if (is.integer(input$file)) {
      cat("No files have been selected (shinyFileChoose)")
    } else {
      file.path(volumes[unlist(input$file[[2]])], 
                paste0(c(unlist(input$file[[1]])[-1]),collapse= "//")
      )
    }
  })
  
  
  proteinGroups <- reactive({
    if(!is.na(input$file[[1]])){
      file.dir <- file.path(volumes[unlist(input$file[[2]])], 
                            paste0(c(unlist(input$file[[1]])[-1]),collapse= "//")
      )
      df <- read.delim(file.dir)
    }
  })
  
  getcolumns <- reactiveVal() 
  input.dir <-  reactiveVal()
  
  observeEvent(listen0(),{ 
    df <- proteinGroups()
    getcolumns(gsub(paste0(input$intensity,".intensity."),"",colnames(df)[grepl(paste0(input$intensity,".intensity."), colnames(df))]))
    
    dir <- dirname(file.path(volumes[unlist(input$file[[2]])], 
                             paste0(c(unlist(input$file[[1]])[-1]),collapse= "//")))
    dir <- gsub("//","/",dir)
    input.dir(dir)
    
  })
  
  
  observeEvent(listen3(),{ 
    
    if(length(input$quality)>0){
      
      html.file <- list.files(input.dir() , pattern =  ".html",full.names = TRUE )
      
      if(length(html.file)==0 & !grepl('report',html.file)){
        runQC(input.dir())
      }
      
      output$qualityreport  <- renderUI({
        includeHTML(list.files(input.dir() , pattern =  ".html",full.names = TRUE ))
      })
    }
  })
  
  
  observeEvent(input$software, {
    
    if(input$software=="Proteom Discoverer"& length(input$import==1)){
      hideTab(inputId = "tabs", target = "Quality report")
    }else{
      showTab(inputId = "tabs", target = "Quality report")
    }
  })
  
  
  observeEvent(input$submit,{
    
    ## Save app reactiveValues in metadata
    metadata[["input.dir"]] <- input.dir()
    metadata[["variables"]] <- getcolumns()
    metadata[["intensity.type"]] <- input$intensity
    
    ## Create putput folder and it's subfolders
    dir.create(file.path(input.dir(),"output"), showWarnings = FALSE)
    dir.create(file.path(input.dir(),"output","backup"), showWarnings = FALSE)
    dir.create(file.path(input.dir(),"output","img"), showWarnings = FALSE)
    dir.create(file.path(input.dir(),"output","table"), showWarnings = FALSE)
    print("output folder was created")
    ## Saving metadata
    if(exists(file.path(input.dir(),"output","backup","targetfile.xlsx"))){
      file.remove(file.path(input.dir(),"output","backup","targetfile.xlsx"))
    }
    
    print(metadata[["targetfile"]])
    xlsx::write.xlsx(metadata[["targetfile"]], file.path(input.dir(),"output", "backup","targetfile.xlsx"), row.names= FALSE)
    
    MSInput <- TemplateProtein$new()$
      importInput(input.dir= metadata[["input.dir"]])$
      removeContaminant(intensity.type = metadata[["intensity.type"]], variables= metadata[["variables"]])$
      transformData()$
      anovaAnalysis()$
      saveTable()$
      visualize() 
    
    save(MSInput, file=file.path(input.dir(),"output","backup","MSInput.RData"))
    
    print("MS object was saved in backup folder ")
  })
  
  
  observeEvent(listen4(),{ 
    
    if(length(input$comparison)>0){
      
      output$volcano <- renderPlotly({
        
        load(file.path(input.dir(),"output","backup","MSInput.RData"))
        if(length(MSInput$df_significant)>0){
          
          df <- MSInput$df_significant[[input$comparison]]
          print(input$comparison)
          y = c(ceiling(10* min(-log10(df$p.adjust.value)))/10, ceiling(10* max(-log10(df$p.adjust.value)))/10)
          x = c(ceiling(10* min(df$log2foldChange, na.rm=T))/10, ceiling(10* max(df$log2foldChange, na.rm=T))/10)
          p <- plot_ly(data = df) %>%
            add_markers( x = ~log2foldChange, y = ~ -log10(p.adjust.value), text = ~ Gene.name, color = ~ Significant, colors = setNames(c("grey","red"), c("No","Yes"))) %>%
            add_lines(x = 1, y = y, line=list( width = 1), name = "log2foldChange = 1") %>%
            add_lines(x= x , y = -log10(0.05), line = list(width = 1), name = "adjusted p.value = 0.05")
          return(p)
        }
      })
    }
  })
  
  
  observeEvent(listen4(),{ 
    
    if(length(input$comparison)>0){
      
      output$histogram <- renderPlotly({
        
        load(file.path(input.dir(),"output","backup","MSInput.RData"))
        if(length(MSInput$df_significant)>0){
          
          df <- MSInput$df_significant[[input$comparison]]
          #fig <- plot_ly(data = df, alpha = 0.6)
          #fig <- fig %>% add_histogram(x = ~p.value, name ="p.value" , color="grey")
          #fig <- fig %>% add_histogram(x = ~p.adjust.value, name ="adjusted p.value", color="red")
          
          
          density1 <- density(df$p.value)
          density2 <- density(df$p.adjust.value)
          
          fig <- plot_ly(x = ~density1$x, y = ~density1$y, type = 'scatter', mode = 'lines', name = 'p.value', fill = 'tozeroy',
                         fillcolor = 'grey',
                         line = list(width = 0.5))%>%
            add_trace(x = ~density2$x, y = ~density2$y, name = 'Adjusted p.value', fill = 'tozeroy',
                                   fillcolor = 'red')%>% 
            layout(xaxis = list(title = ''),
                                yaxis = list(title = 'Density'))
          
          
          return(fig)
        }
      })
    }
  })
  
  
  output$comparison <- renderUI({
    
    selectInput("comparison", "comparison", choices = unique(metadata[["targetfile"]][["Comparison"]]))
  })
  
  
  observeEvent(input$comparison, {
    if(length(input$comparison)>0){
      if(input$submit==1){
        showTab(inputId = "plotlyOutput", target = "volcano")
      }else{
        hideTab(inputId = "plotlyOutput", target = "volcano")
      }
    }
  })
  
  
  observeEvent(input$comparison, {
    if(length(input$comparison)>0){
      if(input$submit==1){
        showTab(inputId = "plotlyOutput", target = "histogram")
      }else{
        hideTab(inputId = "plotlyOutput", target = "histogram")
      }
    }
  })
  
  
}




shinyApp(ui = ui, server = server)