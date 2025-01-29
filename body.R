body_tcga_pathway = fluidRow(
    box(title = '(i) Pathway over-representation in TCGA patient samples',
        width=12, collapsible=TRUE,
        selectizeInput( inputId  = 'tcga_ovr_gonames',
                        label    = 'Enter pathway name(s):',
                        choice   = NULL,
                        multiple = TRUE ),
        uiOutput('tcga_ovr')
    ),

    box(title =paste0('(ii) Inferred pathway levels (IPLs) of proteins in the ',
                       'selected pathway below'),
        width=12, collapsible=TRUE,
        selectizeInput( inputId  = 'tcga_ipl_goname',
                        label    = 'Enter a pathway name:',
                        choice   = NULL,
                        multiple = FALSE ),
        uiOutput('tcga_ipl')
    )
)

body_tcga_protein = fluidRow(
    box(title = "(iii) Features of protein groups",
        width=12, collapsible=TRUE,
        selectizeInput( inputId  = 'tcga_prots',
                        label    = 'Enter protein name(s):',
                        choice   = NULL,
                        multiple = TRUE ),

        tabBox(title=h4('Click tab name for figures'), width=12, side='left',
            tabPanel('Protein inferred pathway levels (IPLs)',
                   # plotOutput('tcga_prot_ipl')),
                     uiOutput('tcga_prot_ipl')),
            tabPanel('Clinical relevance',       plotOutput('tcga_prot_clin')),
            tabPanel('Immune cell infiltration', plotOutput('tcga_prot_infil')))
    ),

    box(title = "(iv) A protein's pathway neighbors and states",
        width=12, collapsible=TRUE,
        selectizeInput( inputId  = 'tcga_nei_stt_prot',
                        label    = 'Enter a protein name:',
                        choice   = NULL,
                        multiple = FALSE ),
      # plotOutput('tcga_nei_stt')
        uiOutput('tcga_nei_stt')
    )
)

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
