require(shiny)
require(shinydashboard)
require(shinymanager)

source('body.R')
source('serv.R')
source('tcga.R')

USE_PASSWORD = FALSE

cred = data.table( user     = c('MPAC'),
                   password = c('HNSC'),
                   stringsAsFactors = FALSE )

header = dashboardHeader(disable=TRUE)

sidebar = sidebarMenu(
    menuItem('TCGA-HNSCC results', startExpanded = TRUE,
             radioButtons( inputId  = 'tcga_hpv',
                           label    = 'Select an HPV type:',
                           choices  = c('HPV+', 'HPV-'),
                           inline   = TRUE,
                           selected = 'HPV+' ),
             
             radioButtons( inputId  = 'tcga_set',
                           label    = 'Select a dataset:',
                           choices  = c('exploratory', 'validation'),
                           inline   = FALSE,
                           selected = 'exploratory' ),
             
             menuSubItem('Pathway',  tabName='tcga_pathway'),
             menuSubItem('Protein',  tabName='tcga_protein')),

    menuItem('Docs', startExpanded = TRUE,
             menuSubItem('Acknowledgement',    tabName='acknowledgement'),
             menuSubItem('Future development', tabName='futuredev'      ),
             menuSubItem('Manuscript figures', tabName='manufig'        ),
             menuSubItem('References',         tabName='references'     ))
) %>% dashboardSidebar(width=200)

body = tabItems(
    tabItem(tabName = 'tcga_pathway',
            h4('Pathway features in TCGA patient samples', align='center'), 
            body_tcga_pathway),
    tabItem(tabName = 'tcga_protein',
            h4('Protein features in TCGA patient samples', align='center'), 
            body_tcga_protein),

    tabItem(tabName = 'acknowledgement',
            h4('Acknowledgement',    align='center'), body_acknowledgement),
    tabItem(tabName = 'futuredev',
            h4('Future Development', align='center'), body_futuredev      ),
    tabItem(tabName = 'manufig',
            h4('Manuscript figures', align='center'), body_manufig        ),
    tabItem(tabName = 'references',
            h4('References',         align='center'), body_references     )
) %>% dashboardBody()


if (USE_PASSWORD == TRUE) {
    ui = secure_app(ui=dashboardPage(header, sidebar, body))
} else {
    ui = dashboardPage(header, sidebar, body)
}

server <- function(input, output, session) {
    if ( USE_PASSWORD == TRUE ) {
        res_auth = secure_server(check_credentials=check_credentials(cred))
        output$res_auth = renderPrint({reactiveValuesToList(res_auth)})
    }

    ## Panel A ##
    servTcgaOvr(input, output, session)

    ## Panel B ##
    servTcgaIpl(input, output, session)

    ## Panel C ##
    servTcgaProtTabs(input, output, session)

    ## Panel D ##
    servTcgaProtNei(input, output, session)
}

shinyApp(ui=ui, server=server)
