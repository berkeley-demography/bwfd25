## this is the user-interface
## a slider with time "t"

library(shiny)

shinyUI(
    fluidPage(
        ## Application title
        titlePanel("Peak Population"),
        h4("Explore the timing of peak births and peak population in an age-structured population with falling fertility"),
        p("This app shows the dynamics of populations undergoing declining fertility rates. It is useful for investigating the case of exponentially declining fertility with constant mortality, what we call the 'Coale Scenario'. It can also be used to see the impact of increasing longevity and of fertility that levels-out below replacement"),
        sidebarLayout(
            sidebarPanel(
        ## layout:
        ## advance time
        ## k
        ## k2
        ## rho
        ## show NRR slope
        ## show predictions
        ## advance time
                p("Advance time with 'play' button or with slider to project the population"),
                sliderInput("myt",
                            "time",
                            min = -50,
                            max = +50,
                            step = 5,
                            value = -50,
                            animate = T),
                h4("Options:"),
                radioButtons("k1",
                             label = "k: rate of NRR change",
                             choices = list("-2%" = -.02,
                                            "-1.5%" = -.015,
                                            "-1%" = -.01),
                             selected = -.015),
        ## k2
        checkboxInput("curvature",
                      label = "Extra-curvature",
                      value = FALSE),
        conditionalPanel(
            condition = "input.curvature == true",
            sliderInput("k2",
                        "Extra-curvature (k2)",
                        min = 0,
                        max = .0003,
                        step = .00005,
                        value = 0)),
        ## rho
        radioButtons("rho",
                     label = "rho: rate of e0 change",
                     choices = list("0%" = 0,
                                    "0.2% (about 1.5 years per decade)" = .002,
                                    "0.3% (about 2.5 years per decade)" = .003),
                     selected = 0),

        ## show NRR slope
        checkboxInput("nrr_slope",
                      label = "Plot generational 'slope'",
                      value = FALSE),

        ## show predictions
        checkboxInput("show_prediction",
                      label = "Show predicted peaks",
                      value = FALSE),
        tags$hr(),
        h4("Exercises:"),
        p("1. Where will births peak -- before or after NRR reaches replacement? Extra-credit: does this change with 'k'?"),
        p("2. To undertand why, check 'plot generational slope' box and notice that when the slope is zero, the peak is about half-way between t = -30 and t = 0."),
        p("3. When does population peak? More extra-credit: does this change with 'k'"),
        p("4. How much does longevity increase delay peak population?"),
        p("5. How much does leveling of NRR below replacement delay poplation?"),
        p("6. How accurate are the analytical approximations?"),
        p("t_B = t_0 - mu_0/2    (Coale)"),
        p("t_N = t_B + A_0       (Coale)"),
        p("t_N^+ = t_N + (-rho/k1)*mu_0  (Goldstein and Cassidy)"),
      p("t_N2 = t_N +  (k2/k1)*sigma_0^2   (Goldstein and Cassidy)")),
        mainPanel(
            plotOutput("peakPlot", height = "800px"),
            ##                h6("Note: insert note here"),
            width = 8
        )
        )
    )
)

##         ## p("To peak when fertility falls below replacement in the specific case of exponentially falling fertility (the 'Coale Scenario')."),
##         ## h4("See how additional changes in longevity and the time path of fertility decline affect the timing of birth and population peaks"),
##         ## p("The app computes a population projection using single-years of time and age with time-varying demographic rates. The NRR panel shows the NRR used for the projection. The 'Births' panel shows the number of births individuals aged '0' in each year of the projection. The 'Population' panel shows the total projected population size."),
##         ## h4("The main outcomes of the 'Coale Scenario' are:"),
##         ## tags$ul(
##         ##          tags$li("Births peak about half-a-generation -- about 15 years -- before fertility reaches replacement"),
##         ##          tags$li("Population peaks about the mean age of the stationary population -- about 40 years -- after births peak")),
##         ## h4("Our extensions suggest that expected increases in longevity and a leveling of fertility decline delay population decline by another 10 to 20 years."),
##         ## Sidebar with a slider input for the number of bins
##         sidebarLayout(
##             sidebarPanel(
##                 h4("Instructions"),
##                 h5("(1) See births peak before replacement fertility"),
##                 p("(Move time slider to the right)"),
##                 sliderInput("myt",
##                             "time (0 when NRR = 1)",
##                             min = -50,
##                             max = +50,
##                             step = 5,
##                             value = -50,
##                             animate = T),
##                 h5("(2) See how birth peak timing does not depend on rate of fertility decline"),
##                 p("(Choose a different rate 'k' of NRR change)"),
##                 radioButtons("k1",
##                              label = "k: rate of NRR change",
##                              choices = list("-2%" = -.02,
##                                             "-1.5%" = -.015,
##                                             "-1%" = -.01),
##                              selected = -.015),
##                 h5("(3) Visualize the link between NRR(t) and B(t)"),
##                 p("(Check the box below, and move the time slider. Check k = -2% for clearer view.)"),
##                 checkboxInput("nrr_slope",
##                               label = "Plot generational 'slope'",
##                               value = FALSE),
##                 h5("(4) See how peak population is about 40 years after peak births"),
##                 ## p("(Check the box below, and move the time slider)"),
##                 ## checkboxInput("show_population",
##                 ##               label = "Plot population",
##                 ##               value = FALSE),
##                 p("(Choosing a different rate of NRR change has only a small effect on the timing of peak population)"),
##                 h5("(5) See how increasing longevity delays peak population"),

##                 p("(Move the time slider to t = 50 and check a rate of e(0) increase below"),
##                 radioButtons("rho",
##                              label = "rho: rate of e0 change",
##                              choices = list("0%" = 0,
##                                             "0.2% (about 1.5 years per decade)" = .002,
##                                             "0.3% (about 2.5 years per decade)" = .003),
##                              selected = 0),
##                 h5("(6) See how flattening of the NRR below replacement delays peak population"),
##                 p("(Check the box below and move the curvature slider"),
##                 checkboxInput("curvature",
##                               label = "Extra-curvature",
##                               value = FALSE),
##                 conditionalPanel(
##                     condition = "input.curvature == true",
##                     sliderInput("k2",
##                                 "Extra-curvature (k2)",
##                                 min = 0,
##                                 max = .0003,
##                                 step = .00005,
##                                 value = 0)),
##                 h5("(7) Compare simulation (black dashed lines) with analytical predictions (orange dashed lines)"),
##                 p("(Try various combinations of k1 (k), k2, and rho)."),
##                 checkboxInput("show_prediction",
##                               label = "Show predicted peaks",
##                               value = FALSE),
##                 h5("Notes: Population is normalized to be 1000 at time 0 across scenarios. Simulation starts at year -100, with stable population for NRR(-100)"),
##             width = 4),
##             ## Show a plot of the generated distribution
 ##             mainPanel(
 ##                 plotOutput("peakPlot", height = "800px"),
 ##                 ##                h6("Note: insert note here"),
 ##                 width = 8
 ##                 )
 ##         )
 ##     )
 ## )



