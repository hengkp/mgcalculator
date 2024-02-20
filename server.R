#
# This is the server logic of a Shiny web application. You can run the
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

options(warn=-1)

library(shiny)
library(shinyjs)
library(datasets)
library(dplyr) # for Data table manipulation
library(minpack.lm) # for fitting non-linear equation
library(MESS) # for AUC calculation
library(reshape) # for melt() function
library(ggplot2) # for plotting bar graph

`%then%` <- shiny:::`%OR%`

# Define server logic required to draw a graph
shinyServer(function(input, output, session) {
    
    output$manualTime <- renderUI({
        numinputs <- numericInput("manualTime","Time (hr)",value = 72)
    })
    
    output$manualSlope <- renderUI({
        numinputs <- numericInput("manualSlope","Slope (%)",value = 2)
    })
    
    output$manualRWD <- renderUI({
        numinputs <- numericInput("manualRWD","RWD (%)",value = 100)
    })
    
    observeEvent(eventExpr = input$selectTime, handlerExpr = {
        if (input$selectTime == 1) {
            updateNumericInput(session, "manualTime", value = 72)
            shinyjs::show("manualTime")
            shinyjs::hide("manualSlope")
            shinyjs::hide("manualRWD")
        } else if (input$selectTime == 2) {
            updateNumericInput(session, "manualSlope", value = 2)
            shinyjs::hide("manualTime")
            shinyjs::show("manualSlope")
            shinyjs::hide("manualRWD")
        } else if (input$selectTime == 3) {
            updateNumericInput(session, "manualRWD", value = 100)
            shinyjs::hide("manualTime")
            shinyjs::hide("manualSlope")
            shinyjs::show("manualRWD")
        }
    })
    
    # Import from user
    data <- reactive({
        validate(
            need(input$dataFile != "", "Please load a CSV file") 
        )
        inFile <- input$dataFile
        if (is.null(inFile)) {
            return(NULL)
        }
        TableDATA <- read.csv(inFile$datapath, header = TRUE)
        headertable <- colnames(TableDATA)
        headerformat <- c("Batch","Cell_Name","Density","Replicate","Drug_Name","Dose",
                          "Mitomycin_C","Time","Relative_Wound_Density","Relative_Caspase")
        # Check input file
        validate(
            need(all(headertable %in% headerformat),"Unrecognized data set")
        )
        CellName <- levels(TableDATA[,"Cell_Name"])
        Density <- unique(TableDATA[,"Density"], incomparables = FALSE)
        DrugName <- levels(TableDATA[,"Drug_Name"])
        # Update lists
        updateSelectInput(session, inputId = "selectCell", label = "Cell Line",
                          choices = CellName, selected = CellName[1])
        updateSelectInput(session, inputId = "selectDensity", label = "Cell Density",
                          choices = Density, selected = Density[1])
        updateSelectInput(session, inputId = "selectDrug", label = "Drug",
                          choices = DrugName, selected = DrugName[1])
        return(TableDATA)
    })
    
    # Import from Example
    data0 <- eventReactive(input$loadexample,{
        myfile <- file.path("www", "ExampleMigrationDATA.csv") 
        TableDATA <- read.csv(myfile, header = TRUE)
        CellName <- levels(TableDATA[,"Cell_Name"])
        Density <- unique(TableDATA[,"Density"], incomparables = FALSE)
        DrugName <- levels(TableDATA[,"Drug_Name"])
        updateSelectInput(session, inputId = "selectCell", label = "Cell Line",
                          choices = CellName, selected = CellName[1])
        updateSelectInput(session, inputId = "selectDensity", label = "Cell Density",
                          choices = Density, selected = Density[1])
        updateSelectInput(session, inputId = "selectDrug", label = "Drug",
                          choices = DrugName, selected = DrugName[1])
        return(TableDATA)
    })
    
    # Check whether Load from example or user
    checkexample <- reactiveValues(a = 1)
    observeEvent(input$loadexample,{
        checkexample$a <- 0 })
    observeEvent(input$dataFile,{
        checkexample$a <- 1 })
    
    output$contents <- renderTable({ 
        if (checkexample$a == 0) { #is.null(input$loadexample) && 
            return(data0())
        } else if (checkexample$a == 1) {
            return(data())
        }
    })
    
    plot.RWDT <- reactive({
        if (checkexample$a == 0) {
            TableDATA<-data0()
        } else if (checkexample$a == 1) {
            TableDATA<-data()
        }
        batch <- input$selectBatch
        cellname <- input$selectCell
        density <- input$selectDensity
        drug <- input$selectDrug
        # Select Control group ----
        Control <- TableDATA %>%
            filter( Cell_Name == cellname , Density == density , Dose == 0 , Mitomycin_C == "No" , 
                    Time %in% seq(0, 72) ) %>%
            group_by(Time) %>%
            summarise(avg = mean(Relative_Wound_Density))
        if(is.data.frame(Control) && nrow(Control)==0) next # skip if no match condition found
        # Fit control group ----
        x <- as.matrix(Control$Time)
        y <- as.matrix(Control$avg)
        dat <- data.frame(x,y)
        yfit <- nlsLM( y ~ c - exp(a + b*x) ,dat, start=list(c = 0, a = 0.4, b = -0.1)) # One-phase asscociated function
        coef <- coef(yfit) # [c , a , b]
        RSS <- sum(residuals(yfit)^2)  # Residual sum of squares
        TSS <- sum((y - mean(y))^2)  # Total sum of squares
        rsquare <- 1 - (RSS/TSS)  # R-squared measure
        if (input$selectTime == 1) {
            halftime <- input$manualTime
        } else if (input$selectTime == 2) {
            z <- input$manualSlope/100
            halftime <- 2 * round(min(x[diff(y)/diff(x)<=z])/2) # Minimum Time of 2% Change (0.02)
            if (is.infinite(halftime) || is.null(halftime)) { 
                halftime <- max(x)
            }
        } else if (input$selectTime == 3) {
            z <- input$manualRWD/100
            halftime <- z * 2 * round(as.numeric((log(coef[1])-log(2)-coef[2])/coef[3])/2) # Half Time in even number
            if (is.infinite(halftime) || is.null(halftime)) { 
                halftime <- max(x)
            }
        }
        return(list(x = x, y = y, yf = predict(yfit), halftime = halftime, rsquare = rsquare))
    })
    
    plot.WCDC <- reactive({
        if (checkexample$a == 0) {
            TableDATA<-data0()
        } else if (checkexample$a == 1) {
            TableDATA<-data()
        }
        cellname <- input$selectCell
        density <- input$selectDensity
        drug <- input$selectDrug
        dat <- plot.RWDT()
        halftime <- dat$halftime
        # Calculate AUC from 0 to (half time) ----
        AUC <- cbind(TableDATA %>%
                         filter( Cell_Name == cellname , Density == density , Drug_Name == drug ,
                                 Mitomycin_C == "No" , Time %in% seq(0, halftime, by=2) ) %>%
                         group_by(Dose,Time) %>%
                         summarize(RWD = mean(Relative_Wound_Density)) %>%
                         summarize(RWD = auc(Time,RWD, type = "spline")) , 
                     TableDATA %>%
                         filter( Cell_Name == cellname , Density == density , Drug_Name == drug ,
                                 Mitomycin_C == "Yes" , Time %in% seq(0, halftime, by=2) ) %>%
                         group_by(Dose,Time) %>%
                         summarize(RWD = mean(Relative_Wound_Density)) %>%
                         summarize(RWD_MMC = auc(Time,RWD, type = "spline")) %>%
                         select(RWD_MMC) , 
                     TableDATA %>%
                         filter( Cell_Name == cellname , Density == density , Drug_Name == drug ,
                                 Mitomycin_C == "No" , Time %in% seq(0, halftime, by=2) ) %>%
                         group_by(Dose,Time) %>%
                         summarize(RCS = mean(Relative_Caspase)) %>%
                         summarize(RCS = auc(Time,RCS, type = "spline")) %>%
                         select(RCS) , 
                     TableDATA %>%
                         filter( Cell_Name == cellname , Density == density , Drug_Name == drug ,
                                 Mitomycin_C == "Yes" , Time %in% seq(0, halftime, by=2) ) %>%
                         group_by(Dose,Time) %>%
                         summarize(RCS = mean(Relative_Caspase)) %>%
                         summarize(RCS_MMC = auc(Time,RCS, type = "spline")) %>%
                         select(RCS_MMC))
        # Calculate WC and DC from AUC ----
        Corre <- data.frame( Dose = AUC$Dose ,
                             WC = (AUC$RWD[1]-AUC$RWD_MMC)/AUC$RWD[1] ,
                             DC = (AUC$RCS_MMC-AUC$RCS_MMC[1])/AUC$RCS_MMC[end(AUC$RCS_MMC)[1]] )
        Corre$WC[Corre$WC<0] = 0
        Corre$WC[Corre$WC>1] = 1
        Corre$DC[Corre$DC<0] = 0
        Corre$DC[Corre$DC>1] = 1
        # Change WC and DC to 4-decimal ----
        Corre$WC = round(Corre$WC, 4)
        Corre$DC = round(Corre$DC, 4)
        Corre$WC[1] = 0
        Corre$DC[1] = 0
        x <- as.vector(Corre$DC)
        y <- as.vector(Corre$WC)
        dat <- data.frame(x,y)
        xnew <- seq(from = 0, to = 1, by = 0.01)
        if (sum(x) == 0) {
            rsquare <- 0
        } else if (sum(y) == 0) {
            rsquare <- 0
        } else {
            file <- try(nlsLM( y ~ a*(1-exp(-b*x)) , dat , start=list(a = 0.1, b = 0.1) , lower = c(0,0) , upper = c(110,1000)))
            if (class(file) == "try-error") {
                rsquare <- 0
            } else {
                yfit <- nlsLM( y ~ a*(1-exp(-b*x)) , dat , start=list(a = 0.1, b = 0.1) , lower = c(0,0) , upper = c(110,1000))
                coef <- as.vector(coef(yfit)) # [c , a , b]
                RSS <- sum(residuals(yfit)^2)  # Residual sum of squares
                TSS <- sum((y - mean(y))^2)  # Total sum of squares
                rsquare <- 1 - (RSS/TSS)  # R-squared measure
                y1 = coef[1]*(1-exp(-coef[2]*x))
                y2 = coef[1]*(1-exp(-coef[2]*(x+1E-6)))
                Slope = (y2-y1)/1E-6
                ynew = coef[1]*(1-exp(-coef[2]*xnew))
            }
            remove(file)
        }
        rsquare[is.nan(rsquare)] = 0
        if (rsquare < 0.3) {
            yfit <- lm( y ~ x ,dat) # Linear Algebra function
            coef <- as.vector(coef(yfit)) # [a, b]
            #coef[is.na(coef)] = 1
            RSS <- sum(residuals(yfit)^2)  # Residual sum of squares
            TSS <- sum((y - mean(y))^2)  # Total sum of squares
            rsquare <- 1 - (RSS/TSS)  # R-squared measure
            if (is.na(coef[2])) {
                Slope = rep(1, length(x))*tan(pi/2)
            } else {
                y1 = coef[1]+(coef[2]*x)
                y2 = coef[1]+(coef[2]*(x+1E-6))
                Slope = (y2-y1)/1E-6
            }
            ynew = coef[1]+(coef[2]*xnew)
        }
        return(list(Dose = Corre$Dose, WC = Corre$WC, DC = Corre$DC, 
                    xnew = xnew, ynew = ynew, Slope = Slope, rsquare = rsquare))
    })
    
    plot.AS <- reactive({
        dat <- plot.WCDC()
        x <- as.vector(dat$DC)
        y <- as.vector(dat$WC)
        Slope <- dat$Slope
        Corre <- data.frame( Dose = dat$Dose )
        # Calculate Migrate , Death , and Drug Independent----
        # Corre$Migrate = c(0,atan((diff(predict(yfit))/diff(x)))/(pi/2))*y
        Corre$Migrate = atan(Slope)/(pi/2)*y
        Corre$Migrate[Corre$Migrate<0] = 0
        Corre$Migrate[Corre$Migrate>1] = 1
        Corre$Migrate[is.na(Corre$Migrate)] = 0
        Corre$Migrate[1] = 0
        Corre$Death = y - Corre$Migrate
        Corre$Death[Corre$Death<0] = 0
        Corre$Death[Corre$Death>1] = 1
        Corre$Death[is.na(Corre$Death)] = 0
        Corre$Death[1] = 0
        Corre$DrugIndependent = 1 - Corre$Migrate - Corre$Death
        # Change Migrate, Death, and DI to 4-decimal ----
        Corre$Migrate = round(Corre$Migrate, 4) 
        Corre$Death = round(Corre$Death, 4)
        Corre$DrugIndependent = round(Corre$DrugIndependent, 4)
        return(list(Migrate = Corre$Migrate,
                    Death = Corre$Death,
                    DrugIndependent = Corre$DrugIndependent))
    })
    
    output$graphRWDT <- renderPlot({
        dat <- plot.RWDT()
        x <- dat$x
        y <- dat$y
        yf <- dat$yf
        halftime <- dat$halftime
        rsquare <- dat$rsquare
        cellname <- input$selectCell
        density <- input$selectDensity
        drug <- input$selectDrug
        # Plot Fitted Control ----
        p <- plot(x,y,col="black", pch = 16, cex = 1.5, ylim=c(0, 100), xlab= "Time(hr)", ylab="RWD",
                  main=sprintf("%s, %s, %s", cellname, drug, density))
        lines(x,yf, col="blue")
        lines(c(halftime,halftime),c(0,100), col="red",lty = 5)
        output$textRWD.halftime <- renderText(sprintf("Half time is %.0f",halftime/2))
        output$textRWD.rsquare <- renderText(sprintf("R-square is %.2f",rsquare))
        
    }, height = function() {
        session$clientData$output_graphRWDT_width * 0.7
    })
    
    output$graphWCDC <- renderPlot({
        Corre <- plot.WCDC()
        x <- as.vector(Corre$DC)
        y <- as.vector(Corre$WC)
        xnew <- Corre$xnew
        ynew <- Corre$ynew
        rsquare <- Corre$rsquare
        cellname <- input$selectCell
        density <- input$selectDensity
        drug <- input$selectDrug
        # Plot Fitted DC Versus WC ----
        p <- plot(x,y,col = rgb(1, 1-(1:8)/8, 0), pch = 16, cex = 2,
                  xlim=c(0,1), ylim=c(0, 1), xlab= "Death change (DC)", ylab="Wound closure change (WC)",
                  main=sprintf("%s, %s, %s", cellname, drug, density))
        if (sum(x) == 0) {
            lines(c(0,0),c(0,1), col="blue")
        } else if (sum(y) == 0) {
            lines(c(0,1),c(0,0), col="blue")
        } else {
            lines(xnew,ynew, col="blue")
        }
        output$textWCDC.rsquare <- renderText(sprintf("R-square is %.2f",rsquare))
    }, height = function() {
        session$clientData$output_graphWCDC_width * 0.7
    })
    
    output$graphAS <- renderPlot({
        cellname <- input$selectCell
        density <- input$selectDensity
        drug <- input$selectDrug
        datw <- plot.WCDC()
        dats <- plot.AS()
        Corre <- data.frame(Dose = datw$Dose,
                            WC = datw$WC,
                            DC = datw$DC,
                            Migration = dats$Migrate,
                            Death = dats$Death,
                            DrugIndependent = dats$DrugIndependent)
        # Plot Correlation between Migrate , Death , and Drug Independent using bar graph ----
        data <- melt(Corre[,c(1,4:6)],1)
        data <- data[order(data$variable,decreasing = TRUE),]
        data$variable <- factor(data$variable, levels = rev(levels(data$variable)))
        x <- c(seq(1,nrow(data)/3,by = 1),seq(1,nrow(data)/3,by = 1),seq(1,nrow(data)/3,by = 1))
        
        ggplot(data, aes(x=x, y=value, fill=variable, order=-as.numeric(variable))) +
            theme(text = element_text(size=14)) +
            geom_bar( colour="black", stat="identity", position="fill") +
            scale_fill_manual("Variable", values = c("DrugIndependent" = "gray", "Death" = "firebrick3", "Migration" = "royalblue4")) +
            scale_x_continuous(breaks = 1:nrow(Corre),
                               labels = unique(data$Dose) ) +
            labs( title = sprintf("%s, %s", cellname, density) ) +
            xlab( sprintf("Dose of %s",drug) ) +
            ylab("Affected score")
    }, height = function() {
        session$clientData$output_graphAS_width * 0.7
    })
    
    output$outputs <- renderTable({
        cellname <- input$selectCell
        density <- input$selectDensity
        drug <- input$selectDrug
        datw <- plot.WCDC()
        dats <- plot.AS()
        Corre <- data.frame(Dose = datw$Dose,
                            WC = datw$WC,
                            DC = datw$DC,
                            Migration = dats$Migrate,
                            Death = dats$Death,
                            DrugIndependent = dats$DrugIndependent)
        
        # Save Corre DATA to Excel ----
        n = nrow(Corre)
        outputs <- matrix(ncol=9, nrow=n)
        outputs[1:n,1] <- cellname
        outputs[1:n,2] <- density
        outputs[1:n,3] <- drug
        outputs[1:n,4:9] <- as.matrix(Corre)
        colnames(outputs) <- c("CellName","Density","DrugName",colnames(Corre))
        outputs <- outputs[complete.cases(outputs), ]
        return(outputs)
    })
})