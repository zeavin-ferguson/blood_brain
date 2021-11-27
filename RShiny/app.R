library(shiny)
library(ggplot2)
library(Hmisc)
library(corrplot)
library(shinythemes)

all_genes<-read.csv("all_genes.csv")
ui <- fluidPage(theme = shinytheme("superhero"),
                titlePanel("Gene Expression Correlation with Alcohol Traits"),
                sidebarLayout(
                  sidebarPanel(selectInput(inputId = "x", 
                                           label = "Genes:", 
                                           choices = all_genes$x,
                                           selected = "Ccl5", 
                                           multiple = F,
                                           selectize = T),
                               selectInput(inputId = "y",
                                           label = "Alcohol Trait:",
                                           choices = c("IncreaseFromBaseline", 
                                                       "LastCIE_g_kg_2hr"),
                                           selected = "LastCIE_g_kg_2hr")),
                  mainPanel(
                    tabsetPanel(type = "tabs",
                                tabPanel("Blood (F)", plotOutput(outputId = "scatterplot1")),
                                tabPanel("Blood (M)", plotOutput(outputId = "scatterplot2")),
                                tabPanel("Hyp (F)", plotOutput(outputId = "scatterplot3")),
                                tabPanel("Hyp (M)", plotOutput(outputId = "scatterplot4")),
                                tabPanel("Amyg (F)", plotOutput(outputId = "scatterplot5")),
                                tabPanel("Amyg (M)", plotOutput(outputId = "scatterplot6")),
                                tabPanel("PFC (F)", plotOutput(outputId = "scatterplot7")),
                                tabPanel("PFC (M)", plotOutput(outputId = "scatterplot8")),
                                tabPanel("About", textOutput("text1"))
                    )
                    
                  )
                )
)




server <- function(input, output) {

  output$scatterplot1 <- renderPlot({
    dataset<-read.csv("blood_f.csv")
    p = ggplot(data = dataset) +
      aes_string(x = paste0("`",input$x,"`"), y = input$y) + #i tried putting quotes around input x (gene name) so that the Riken gene would work....throw errors due to starting with a digit...this did not work
      geom_point(aes(color=Group, shape=Group))+
      geom_smooth(method="lm")+
      labs(subtitle = paste0("corr = ", round(cor.test(dataset[, names(dataset) == input$x], dataset[, names(dataset) == input$y], method="pearson")$estimate, digits = 2), ", p = ", round(cor.test(dataset[, names(dataset) == input$x], dataset[, names(dataset) == input$y], method="pearson")$p.value, digits = 4)))
    plot(p)})
  
  output$scatterplot2 <- renderPlot({
    dataset<-read.csv("blood_m.csv")
    p = ggplot(data = dataset) +
      aes_string(x = paste0("`",input$x,"`"), y = input$y) + #i tried putting quotes around input x (gene name) so that the Riken gene would work....throw errors due to starting with a digit...this did not work
      geom_point(aes(color=Group, shape=Group))+
      geom_smooth(method="lm")+
      labs(subtitle = paste0("corr = ", round(cor.test(dataset[, names(dataset) == input$x], dataset[, names(dataset) == input$y], method="pearson")$estimate, digits = 2), ", p = ", round(cor.test(dataset[, names(dataset) == input$x], dataset[, names(dataset) == input$y], method="pearson")$p.value, digits = 4)))
    plot(p)})

  
  output$scatterplot3 <- renderPlot({
    dataset<-read.csv("hyp_f.csv")
    p = ggplot(data = dataset) +
      aes_string(x = paste0("`",input$x,"`"), y = input$y) + #i tried putting quotes around input x (gene name) so that the Riken gene would work....throw errors due to starting with a digit...this did not work
      geom_point(aes(color=Group, shape=Group))+
      geom_smooth(method="lm")+
      labs(subtitle = paste0("corr = ", round(cor.test(dataset[, names(dataset) == input$x], dataset[, names(dataset) == input$y], method="pearson")$estimate, digits = 2), ", p = ", round(cor.test(dataset[, names(dataset) == input$x], dataset[, names(dataset) == input$y], method="pearson")$p.value, digits = 4)))
    plot(p)})
  
  output$scatterplot4 <- renderPlot({
    dataset<-read.csv("hyp_m.csv")
    p = ggplot(data = dataset) +
      aes_string(x = paste0("`",input$x,"`"), y = input$y) + #i tried putting quotes around input x (gene name) so that the Riken gene would work....throw errors due to starting with a digit...this did not work
      geom_point(aes(color=Group, shape=Group))+
      geom_smooth(method="lm")+
      labs(subtitle = paste0("corr = ", round(cor.test(dataset[, names(dataset) == input$x], dataset[, names(dataset) == input$y], method="pearson")$estimate, digits = 2), ", p = ", round(cor.test(dataset[, names(dataset) == input$x], dataset[, names(dataset) == input$y], method="pearson")$p.value, digits = 4)))
    plot(p)})
  
  output$scatterplot5 <- renderPlot({
    dataset<-read.csv("amy_f.csv")
    p = ggplot(data = dataset) +
      aes_string(x = paste0("`",input$x,"`"), y = input$y) + #i tried putting quotes around input x (gene name) so that the Riken gene would work....throw errors due to starting with a digit...this did not work
      geom_point(aes(color=Group, shape=Group))+
      geom_smooth(method="lm")+
      labs(subtitle = paste0("corr = ", round(cor.test(dataset[, names(dataset) == input$x], dataset[, names(dataset) == input$y], method="pearson")$estimate, digits = 2), ", p = ", round(cor.test(dataset[, names(dataset) == input$x], dataset[, names(dataset) == input$y], method="pearson")$p.value, digits = 4)))
    plot(p)})
  
  output$scatterplot6 <- renderPlot({
    dataset<-read.csv("amy_m.csv")
    p = ggplot(data = dataset) +
      aes_string(x = paste0("`",input$x,"`"), y = input$y) + #i tried putting quotes around input x (gene name) so that the Riken gene would work....throw errors due to starting with a digit...this did not work
      geom_point(aes(color=Group, shape=Group))+
      geom_smooth(method="lm")+
      labs(subtitle = paste0("corr = ", round(cor.test(dataset[, names(dataset) == input$x], dataset[, names(dataset) == input$y], method="pearson")$estimate, digits = 2), ", p = ", round(cor.test(dataset[, names(dataset) == input$x], dataset[, names(dataset) == input$y], method="pearson")$p.value, digits = 4)))
    plot(p)})
  
  output$scatterplot7 <- renderPlot({
    dataset<-read.csv("pfc_f.csv")
    p = ggplot(data = dataset) +
      aes_string(x = paste0("`",input$x,"`"), y = input$y) + #i tried putting quotes around input x (gene name) so that the Riken gene would work....throw errors due to starting with a digit...this did not work
      geom_point(aes(color=Group, shape=Group))+
      geom_smooth(method="lm")+
      labs(subtitle = paste0("corr = ", round(cor.test(dataset[, names(dataset) == input$x], dataset[, names(dataset) == input$y], method="pearson")$estimate, digits = 2), ", p = ", round(cor.test(dataset[, names(dataset) == input$x], dataset[, names(dataset) == input$y], method="pearson")$p.value, digits = 4)))
    plot(p)})
  
  output$scatterplot8 <- renderPlot({
    dataset<-read.csv("pfc_m.csv")
    p = ggplot(data = dataset) +
      aes_string(x = paste0("`",input$x,"`"), y = input$y) + #i tried putting quotes around input x (gene name) so that the Riken gene would work....throw errors due to starting with a digit...this did not work
      geom_point(aes(color=Group, shape=Group))+
      geom_smooth(method="lm")+
      labs(subtitle = paste0("corr = ", round(cor.test(dataset[, names(dataset) == input$x], dataset[, names(dataset) == input$y], method="pearson")$estimate, digits = 2), ", p = ", round(cor.test(dataset[, names(dataset) == input$x], dataset[, names(dataset) == input$y], method="pearson")$p.value, digits = 4)))
    plot(p)})
  
  output$text1 <- renderText({
    "We used the chronic intermittent ethanol (CIE) procedure to induce alcohol dependence in male and female C57BL/6J mice. One group was exposed to intermittent ethanol vapor and the other to control air in identical chambers. Mice underwent 4 vapor cycles which involved 4 days of vapor (16 hours vapor on, 8 hours vapor off) followed by 3 days of rest in the homecage. All mice were given 5 days of access to water and ethanol to measure voluntary alcohol intake in between the vapor sessions. We measured gene expression levels in whole blood, amygdala, prefrontal cortex, and hypothalamus one week following a final (5th) vapor exposure. Details can be found in the publication (publication link). E-mail: lferguson@utexas.edu"
  })
  
}

shinyApp(ui, server)
