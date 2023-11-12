body_tcga_pathway = fluidRow(
    box(title = '(A) Pathway over-representation in TCGA patient samples',
        width=12, collapsible=TRUE,
        selectizeInput( inputId  = 'tcga_ovr_gonames',
                        label    = 'Enter pathway name(s):',
                        choice   = NULL,
                        multiple = TRUE ),
        plotOutput('tcga_ovr')),

    box(title = paste0('(B) Inferred pathway levels (IPLs) of proteins in the ',
                       'selected pathway below'),
        width=12, collapsible=TRUE, height=1600,
        selectizeInput( inputId  = 'tcga_ipl_goname',
                        label    = 'Enter a pathway name:',
                        choice   = NULL,
                        multiple = FALSE ),
        plotOutput('tcga_ipl'))
)

body_tcga_protein = fluidRow(
    box(title = "(C) Features of protein groups",
        width=12, collapsible=TRUE, 
        selectizeInput( inputId  = 'tcga_prots',
                        label    = 'Enter protein name(s):',
                        choice   = NULL,
                        multiple = TRUE ),

        tabBox(title=h4('Click tab name for figures'), width=12, side='left',
            tabPanel('Protein inferred pathway levels (IPLs)',   
                     plotOutput('tcga_prot_ipl')),
            tabPanel('Clinical relevance',       plotOutput('tcga_prot_clin')),
            tabPanel('Immune cell infiltration', plotOutput('tcga_prot_infil')))
    ),

    box(title = "(D) A protein's pathway neighbors and states",
        width=12, collapsible=TRUE, 
        selectizeInput( inputId  = 'tcga_nei_stt_prot',
                        label    = 'Enter a protein name:',
                        choice   = NULL,
                        multiple = FALSE ),
        plotOutput('tcga_nei_stt')
    )
)

serveTcga <- function(input, output, session) {
    updateSelectizeInput( 
        session, 
        inputId  = 'tcga_ovr_gonames',
        choices  = sort(rownames(inpl$tcgal$ovrmat)),
        selected = c(
            'formation of primary germ layer',
            'anatomical structure formation involved in morphogenesis',
            'positive regulation of interleukin-2 biosynthetic process',
            'immune effector process',
            'leukocyte activation',
            'regulation of immune effector process'
        ),
        server = TRUE 
    )
    output$tcga_ovr = renderPlot(hmTcgaOvr(input))

    updateSelectizeInput( 
        session, 
        inputId  = 'tcga_ipl_goname',
        choices  = sort(rownames(inpl$tcgal$ovrmat)),
        selected = 'regulation of immune effector process',
        server = TRUE 
    )
    output$tcga_ipl = renderPlot(hmTcgaIpl(input), height=1500)


    updateSelectizeInput( 
        session, 
        inputId  = 'tcga_prots',
        choices  = sort(rownames(inpl$tcgal$omic_iplmat)),
        selected = c('CD28', 'CD86', 'LCP2', 'IL12RB1', 'TYK2', 'CD247', 
                     'FASLG'),
        server = TRUE 
    )
    output$tcga_prot_ipl   = renderPlot(hmTcgaProtIpl(input))
    output$tcga_prot_clin  = renderPlot(srTcgaProtIplOsPfi(input))
    output$tcga_prot_infil = renderPlot(bxTcgaProtIplInfil(input))


    updateSelectizeInput( 
        session, 
        inputId  = 'tcga_nei_stt_prot',
        choices  = sort(rownames(inpl$tcgal$omic_iplmat)),
        selected = 'CD28',
        server   = TRUE 
    )
    output$tcga_nei_stt = renderPlot(hmTcgaNeiStt(input), height=1000)
}

body_acknowledgement = fluidRow(
    box(title = "", width=12, collapsible=TRUE, 
        htmltools::includeMarkdown('ack.md') )
)

body_futuredev = fluidRow(
    box(title = "", width=12, collapsible=TRUE, 
        htmltools::includeMarkdown('todo.md') )
)

body_manufig = fluidRow(
    box(title = "", width=12, collapsible=TRUE, 
        htmltools::includeMarkdown('fig.md') )
)

body_references = fluidRow(
    box(title = "", width=12, collapsible=TRUE, 
        htmltools::includeMarkdown('ref.md') )
)
