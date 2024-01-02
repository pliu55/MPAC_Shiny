## Panel A ##
servTcgaOvr <- function(input, output, session) {
    updateSelectizeInput( 
        session, 
        inputId  = 'tcga_ovr_gonames',
        choices  = sort(rownames(inpl$tcgal$ovrmat)),
        selected = c(
          # 'formation of primary germ layer',
          # 'anatomical structure formation involved in morphogenesis',
            'positive regulation of interleukin-2 biosynthetic process',
            'immune effector process',
            'leukocyte activation',
            'regulation of immune effector process'
        ),
        server = TRUE 
    )

    getHmTcgaOvrHeight = reactive({
        req(input$tcga_ovr_gonames)
        return(40 * length(input$tcga_ovr_gonames) + 100)
    })

    output$tcga_ovr_hm = hmTcgaOvr(input) |> renderPlot()
    output$tcga_ovr = plotOutput('tcga_ovr_hm', height=getHmTcgaOvrHeight()) |> 
                      renderUI()
}

## Panel B ##
servTcgaIpl <- function(input, output, session) {
    updateSelectizeInput( 
        session, 
        inputId  = 'tcga_ipl_goname',
        choices  = sort(rownames(inpl$tcgal$ovrmat)),
      # selected = 'regulation of immune effector process',
        selected = 'positive regulation of interleukin-2 biosynthetic process',
        server = TRUE 
    )

    getHmTcgaIplHeight = reactive({
        req(input$tcga_ipl_goname)
        n_gns = intersect(inpl$gol$gmtl[[input$tcga_ipl_goname]],
                          rownames(inpl$tcgal$omic_iplmat)) |> length()
        return(15 * n_gns + 100)
    })

    output$tcga_ipl_hm = hmTcgaIpl(input) |> renderPlot()
    output$tcga_ipl = plotOutput('tcga_ipl_hm', height=getHmTcgaIplHeight()) |> 
                      renderUI()
}

## Panel C ##
servTcgaProtTabs <- function(input, output, session) {
    updateSelectizeInput( 
        session, 
        inputId  = 'tcga_prots',
        choices  = sort(rownames(inpl$tcgal$omic_iplmat)),
        selected = c('CD28', 'CD86', 'LCP2', 'IL12RB1', 'TYK2', 'CD247', 
                     'FASLG'),
        server = TRUE 
    )

    getHmTcgaProtIplHeight = reactive({
        req(input$tcga_prots)
        return(20 * length(input$tcga_prots) + 100)
    })
    output$tcga_prot_ipl_hm = hmTcgaProtIpl(input) |> renderPlot()
    output$tcga_prot_ipl = plotOutput('tcga_prot_ipl_hm', 
                                      height=getHmTcgaProtIplHeight()) |> 
                           renderUI()

    output$tcga_prot_clin  = renderPlot(srTcgaProtIplOsPfi(input))
    output$tcga_prot_infil = renderPlot(bxTcgaProtIplInfil(input))
}

## Panel D ##
servTcgaProtNei <- function(input, output, session) {
    updateSelectizeInput( 
        session, 
        inputId  = 'tcga_nei_stt_prot',
        choices  = sort(rownames(inpl$tcgal$omic_iplmat)),
        selected = 'CD86',
        server   = TRUE 
    )

    getHmTcgaNeiSttHeight = reactive({
        req(input$tcga_nei_stt_prot)

        n_up = nrow(inpl$pthl$edgedt[ to   == input$tcga_nei_stt_prot])
        n_dn = nrow(inpl$pthl$edgedt[ from == input$tcga_nei_stt_prot])
        return(15 * (n_up + n_dn + 3) + 100)
    })
    output$tcga_nei_stt_hm = hmTcgaNeiStt(input) |> renderPlot()
    output$tcga_nei_stt = plotOutput('tcga_nei_stt_hm', 
                                     height=getHmTcgaNeiSttHeight()) |> 
                          renderUI()
}
