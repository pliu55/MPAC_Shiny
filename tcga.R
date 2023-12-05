#
#  pliu 2022-10-17
#
#  - cannot use ggrastr because it requires textshaping package, which cannot
#    be installed on data-viz.it.wisc.edu
#
#

require(data.table)
require(magrittr)
require(ggplot2)
require(patchwork)
suppressMessages(require(ComplexHeatmap))

source('param.R')
inpl = readRDS('inpl.rds')

main <- function() {
  # hmTcgaOvr(input=list(tcga_hpv='HPV+', tcga_set='exploratory',
  #     tcga_ovr_gonames=c(
  #         'immune effector process',
  #         'leukocyte activation',
  #         'immune response-activating signal transduction',
  #         'immune response-activating cell surface receptor signaling pathway'
  #     )))
  # hmTcgaIpl(input=list(tcga_hpv='HPV+', tcga_set='exploratory',
  #     tcga_ipl_goname='immune response-activating signal transduction'))
  # hmTcgaProtIpl(input=list(tcga_hpv='HPV+', tcga_set='exploratory',
  #     tcga_prots=c('CD28', 'CD86', 'LCP2', 'IL12RB1', 'TYK2', 'CD247',
  #                  'FASLG')))
  # srTcgaProtIplOsPfi(input=list(tcga_hpv='HPV+', tcga_set='exploratory',
  #     tcga_prots=c('CD28', 'CD86', 'LCP2', 'IL12RB1', 'TYK2', 'CD247',
  #                  'FASLG')))
  # bxTcgaProtIplInfil(input=list(tcga_hpv='HPV+', tcga_set='exploratory',
  #     tcga_prots=c('CD28', 'CD86', 'LCP2', 'IL12RB1', 'TYK2', 'CD247',
  #                  'FASLG')))
  # hmTcgaNeiStt(input=list(tcga_hpv='HPV+', tcga_set='exploratory',
  #              tcga_nei_stt_prot_name='CD28'))
}

hmTcgaNeiStt <- function(input) {
    shiny::req(input$tcga_hpv)
    shiny::req(input$tcga_set)
    shiny::req(input$tcga_nei_stt_prot)

    prot = input$tcga_nei_stt_prot
    cldt = inpl$tcgal$patdt[(hpv == input$tcga_hpv) & (set == input$tcga_set)]
    pats = cldt$pat

    prot_iplmat = inpl$tcgal$ipl_sttmat[prot, pats, drop=FALSE]
    prot_cnmat  = inpl$tcgal$cn_sttmat[prot,  pats, drop=FALSE]
    prot_rnamat = inpl$tcgal$rna_sttmat[prot, pats, drop=FALSE]
    rownames(prot_cnmat)  = 'CNA'
    rownames(prot_rnamat) = 'RNA-seq'

    grpdt = data.table( rbind(
        data.table(ent = rownames(prot_cnmat),  loc = 'omic' ),
        data.table(ent = rownames(prot_rnamat), loc = 'omic' ),
        data.table(ent = prot,                  loc = prot   )
    )) 

    up_edgedt = inpl$pthl$edgedt[ to == prot, .(from, title_lab) ] %>% 
                setnames(1, 'ent')
    dn_edgedt = inpl$pthl$edgedt[ from == prot, .(to, title_lab) ] %>%
                setnames(1, 'ent')

    up_ents = dn_ents = up_iplmat = dn_iplmat = NULL
    edgedt = data.table()
    if ( nrow(up_edgedt) > 0 ) {
        up_ents = sort(up_edgedt$ent)
        edgedt = rbind(edgedt, up_edgedt)
        up_iplmat = inpl$tcgal$ipl_sttmat[up_ents, pats, drop=FALSE]
        grpdt = data.table(ent = up_ents, loc = 'upstream') %>% rbind(grpdt, .)
    }
    if ( nrow(dn_edgedt) > 0 ) {
        dn_ents = sort(dn_edgedt$ent)
        edgedt = rbind(edgedt, dn_edgedt)
        dn_iplmat = inpl$tcgal$ipl_sttmat[dn_ents, pats, drop=FALSE]
        grpdt = data.table(ent = dn_ents, loc = 'downstream') %>%rbind(grpdt, .)
    }

    grpdt = merge(grpdt, edgedt, by='ent', all.x=TRUE) %>%
            .[, grp := ifelse(is.na(title_lab), loc, 
                              paste0(loc, ":\n", title_lab))] %>%
            .[, loc := factor(loc, levels=c('upstream', 'omic', prot, 
                                            'downstream'))] %>%
            .[order(loc, grp, ent)]

    grpdt[, grp := factor(grp, levels=unique(grpdt$grp))]

    pltmat = rbind(up_iplmat, prot_cnmat, prot_rnamat, prot_iplmat, 
                   dn_iplmat) %>% .[ grpdt$ent, ]
    rownames(pltmat) = rownames(pltmat) %>% stringr::str_wrap(width=60)

    Heatmap(pltmat,
        col = structure(c(MAGENTA, 'gray85', CYAN), names=c(1, 0, -1)),
        rect_gp            = gpar(col='white'),
        border             = 'black',
        column_split       = cldt$icl_lab,
        row_split          = grpdt$grp,
        cluster_column_slices = FALSE,
        cluster_row_slices = FALSE,
        row_gap            = unit(0.2, 'cm'),
        row_title_rot      = 0,
        row_title_gp       = gpar(fontsize=FONT_SIZE-2),
        row_title_side     = 'right',
        cluster_rows       = TRUE,
        cluster_columns    = TRUE,
        clustering_distance_rows    = 'manhattan',
        clustering_distance_columns = 'manhattan',
        show_row_dend      = FALSE,
        show_column_dend   = FALSE,
        show_row_names     = TRUE,
        row_names_side     = 'left',
        show_column_names  = FALSE,
        row_names_gp       = gpar(fontsize=FONT_SIZE-2),
      # column_names_gp    = gpar(fontsize=FONT_SIZE-3),
      # column_names_side  = 'top',
        column_names_rot   = 30,
        height = unit(0.8*nrow(pltmat), 'cm'),
      # width  = unit(0.2*ncol(pltmat), 'cm'),
        width  = unit(20, 'cm'),
        heatmap_legend_param = list(
            title       = "omic or\npathway\nstate",
            title_gp    = gpar(fontsize=FONT_SIZE),
            labels_gp   = gpar(fontsize=FONT_SIZE),
            grid_height = grid::unit(3.5, 'mm'),
            grid_width  = grid::unit(3.5, 'mm'),
            at          = 1:-1,
            labels      = c('activated', 'normal', 'repressed'),
            border      = 'black',
            title_position = 'lefttop',
            direction   = 'horizontal'
        )
    ) %>% draw( heatmap_legend_side = 'bottom',
                background   = 'transparent',
              # row_names_gp = gpar(fontsize=FONT_SIZE-2),
                padding      = unit(c(0, 7, 0, 0), 'cm') )
}


bxTcgaProtIplInfil <- function(input) {
    shiny::req(input$tcga_hpv)
    shiny::req(input$tcga_set)
    shiny::req(input$tcga_prots)

    cldt = inpl$tcgal$patdt[(hpv == input$tcga_hpv) & (set == input$tcga_set)]
    iplmat = inpl$tcgal$omic_iplmat[input$tcga_prots, cldt$pat, drop=FALSE]
    act_pats = getActPats(iplmat)

    n_gns = length(input$tcga_prots)
    prt_lab = ifelse(n_gns > 1, paste0('all ', n_gns, ' proteins'), 
                     input$tcga_prots[1])
    act_lab = paste0(prt_lab, ' with IPLs > 0 (n=', length(act_pats), ')') %>%
              stringr::str_wrap(width=40)
    other_lab = paste0('others (n=', ncol(iplmat) - length(act_pats), ')')

    frcdt = inpl$tcgal$infilmat[, cldt$pat, drop=FALSE] %>% 
            as.data.table(keep.rownames='ct') %>%
            melt(id='ct', variable='pat', value='frc') %>%
            .[, grp := ifelse(pat %in% act_pats, act_lab, other_lab)]

    title = paste0(ncol(iplmat), ' patient samples in ', input$tcga_hpv,
                   ' ', input$tcga_set, ' set')

    p = ggplot(frcdt, aes(x=ct, y=frc, fill=grp)) +
        geom_boxplot( outlier.shape=NA, width=0.5, alpha=0.7) +
        theme_classic() + SQ_TM +
        theme( aspect.ratio = 1/4, 
               axis.title.x = element_blank(),
               axis.text.x = element_text(angle=40, hjust=1),
               plot.margin = margin(l=0, r=1.0, b=0, t=0, unit='cm'),
               legend.position = 'top',
               legend.title = element_blank(),
               legend.direction = 'vertical' ) +
        scale_fill_manual(values = c(MAGENTA, 'gray90')) +
        scale_y_continuous(labels=scales::percent) +
        ylab('cell fraction')

    print(p)
}

getActPats <- function(iplmat) {
    stt_sums = t(iplmat) %>% sign() %>% rowSums()
    act_pats = names(stt_sums[stt_sums == nrow(iplmat)])
    return(act_pats)
}

srTcgaProtIplOsPfi <- function(input) {
    shiny::req(input$tcga_hpv)
    shiny::req(input$tcga_set)
    shiny::req(input$tcga_prots)

    cldt = inpl$tcgal$patdt[(hpv == input$tcga_hpv) & (set == input$tcga_set)]
    iplmat = inpl$tcgal$omic_iplmat[input$tcga_prots, cldt$pat, drop=FALSE]
    os_frm  = survival::Surv(OS_days, OS)~group
    pfi_frm = survival::Surv(PFI_days, PFI)~group

    os_title  = 'overall survival'
    os_ylab   = 'overall survival probability'

    pfi_title = 'progression-free survival'
    pfi_ylab  = 'progression-free survival probability'

    os_p  = srByActIpls(input$tcga_prots, os_title, iplmat, cldt, os_frm, 
                        os_ylab)

    pfi_p = srByActIpls(input$tcga_prots, pfi_title, iplmat, cldt, pfi_frm, 
                        pfi_ylab)

    wrap_plots(os_p, plot_spacer(), pfi_p, widths=c(1, 0.1, 1), ncol=3) %>%
    return()
}

srByActIpls <- function(gns, title, iplmat, cdrdt, frm, ylab) {
    act_pats = getActPats(iplmat)
    n_gns = length(gns)
    prt_lab = ifelse(n_gns > 1, paste0('all ', n_gns, ' proteins'), gns[1])
    act_lab = paste0(prt_lab, ' with IPLs > 0 (n=', length(act_pats), ')') %>%
              stringr::str_wrap(width=40)
    other_lab = paste0('others (n=', ncol(iplmat) - length(act_pats), ')')

    cdrdt[, group := ifelse(pat %in% act_pats, act_lab, other_lab) %>%
                     factor(levels=c(act_lab, other_lab))]

    p = survminer::ggsurvplot(
        fit = do.call(survival::survfit, args=list(formula=frm, data=cdrdt)),
        data = cdrdt,
        palette = c(MAGENTA, GREY),
        pval = TRUE,
        pval.size = 5,
        pval.coord = c(0.1, 0.05),
        xlab = 'days',
        ylab = ylab,
        legend = 'top',
        legend.title = '',
        title = title,
        ggtheme = theme_classic() + SQ_TM ) + 
        guides(color=guide_legend(nrow=2))

    return(p$plot)
}

hmTcgaProtIpl <- function(input) {
    shiny::req(input$tcga_hpv)
    shiny::req(input$tcga_set)
    shiny::req(input$tcga_prots)

    cldt = inpl$tcgal$patdt[(hpv == input$tcga_hpv) & (set == input$tcga_set)]
    pltmat = inpl$tcgal$omic_iplmat[input$tcga_prots, cldt$pat, drop=FALSE]

    col_title = paste0(ncol(pltmat), ' ', 'patient samples in ',
                       input$tcga_hpv, ' ', input$tcga_set, ' set ')
    hmIplWithSampCl(pltmat, cldt, col_title)
}


hmTcgaIpl <- function(input) {
    shiny::req(input$tcga_hpv)
    shiny::req(input$tcga_set)
    shiny::req(input$tcga_ipl_goname)

    cldt = inpl$tcgal$patdt[(hpv == input$tcga_hpv) & (set == input$tcga_set)]
    gnnames = intersect(inpl$gol$gmtl[[input$tcga_ipl_goname]],
                        rownames(inpl$tcgal$omic_iplmat))
    pltmat = inpl$tcgal$omic_iplmat[gnnames, cldt$pat, drop=FALSE]

    col_title = paste0(ncol(pltmat), ' ', 'patient samples in ',
                       input$tcga_hpv, ' ', input$tcga_set, ' set ')
    hmIplWithSampCl(pltmat, cldt, col_title)
}

hmIplWithSampCl <- function(pltmat, cldt, col_title) {
    Heatmap( pltmat,
        col = STT_CLRS,
        na_col = 'black',
        border            = 'black',
        rect_gp           = gpar(col='white'),
        column_split      = cldt$icl_lab,
        column_gap        = unit(0.2, 'cm'),
        show_row_names    = TRUE,
        show_column_names = FALSE,
        column_names_side = 'top',
        row_names_gp      = gpar(fontsize=10-floor(length(nrow(pltmat))/50)),
        column_title_gp   = gpar(fontsize=FONT_SIZE),
        row_title_rot     = 0,
        row_title_gp      = gpar(fontsize=FONT_SIZE),
        show_row_dend     = FALSE,
        show_column_dend  = FALSE,
        cluster_column_slices = FALSE,
        cluster_row_slices    = FALSE,
        cluster_columns   = TRUE,
        cluster_rows      = FALSE,
        use_raster        = FALSE,
      # heatmap_height    = unit(length(gnnames)/10, 'cm'), 
        heatmap_legend_param = list(
            title       = 'IPL',
            title_gp    = gpar(fontsize=FONT_SIZE),
            labels_gp   = gpar(fontsize=FONT_SIZE),
            grid_height = grid::unit(3, 'mm'),
            grid_width  = grid::unit(3, 'mm'),
            border      = 'black',
            title_position = 'lefttop',
            direction   = 'horizontal'
        )
    ) %>% draw(
        heatmap_legend_side = 'bottom',
        column_title        = col_title,
        column_title_gp     = gpar(fontsize=FONT_SIZE),
        background          = 'transparent'
    )
}

hmTcgaOvr <- function(input) {
    shiny::req(input$tcga_hpv)
    shiny::req(input$tcga_set)
    shiny::req(input$tcga_ovr_gonames)

    cldt = inpl$tcgal$patdt[(hpv == input$tcga_hpv) & (set == input$tcga_set)]
    pltmat = inpl$tcgal$ovrmat[input$tcga_ovr_gonames, cldt$pat, drop=FALSE]
    rownames(pltmat) = rownames(pltmat) %>% stringr::str_wrap(width=40)

    set.seed(123456)
    p = Heatmap( pltmat,
        col = OVR_CLRS,
        na_col = 'black',
        border            = 'black',
        column_split      = cldt$icl_lab,
        column_title_gp   = gpar(fontsize=FONT_SIZE),
        cluster_column_slices = FALSE,
        show_row_names    = TRUE,
        row_names_gp      = gpar(fontsize=(FONT_SIZE-1)),
        show_column_names = FALSE,
        column_names_side = 'top',
        column_names_gp   = gpar(fontsize=FONT_SIZE),
        show_row_dend     = FALSE,
        show_column_dend  = FALSE,
        cluster_columns   = TRUE,
        cluster_rows      = TRUE,
        use_raster        = ifelse( (nrow(pltmat) > 20) | (ncol(pltmat) > 20),
                                    TRUE, FALSE),
        raster_quality    = 5,
        heatmap_legend_param = list(
            title       = expression(paste('log'[10], '(adjusted P)')),
            title_gp    = gpar(fontsize=FONT_SIZE),
            labels_gp   = gpar(fontsize=FONT_SIZE),
            legend_width = grid::unit(0.8, 'inch'),
            title_position = 'lefttop',
            direction   = 'horizontal'
        )
    ) %>% draw(
        background      = 'transparent',
        heatmap_legend_side = 'bottom',
        column_title    = paste0(ncol(pltmat), ' ', 'patient samples in ',
                                 input$tcga_hpv, ' ', input$tcga_set, ' set '),
        column_title_gp = gpar(fontsize=FONT_SIZE),
        padding = unit(c(2, 2, 2, 35), 'mm')
    )

    return(p)
}

main()
