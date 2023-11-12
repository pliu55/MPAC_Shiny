require(ggplot2)

SQ_TM = theme( aspect.ratio  = 1,
               plot.title    = element_text(size=15, hjust=0.5),
               strip.text    = element_text(size=15),
               legend.title  = element_text(size=15),
               legend.text   = element_text(size=14),
               axis.title    = element_text(size=15),
               axis.text     = element_text(size=14),
               legend.margin = margin(l=-0.3, unit='cm') )

FONT_SIZE = 15

DT_OPT=list(dom='Bfrtip', buttons=c('excel'),
            columnDefs=list(list(className = 'dt-center', targets = '_all')))

## no difference in using colorRampPalette or not
EXPR_COLORS = rev(RColorBrewer::brewer.pal(11, 'RdYlBu'))

## color-blind friendly
## from http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/#a-colorblind-friendly-palette
BLACK   = '#000000'
GREY    = '#999999'
ORANGE  = '#E69F00'
CYAN    = '#56B4E9'
GREEN   = '#009E73'
YELLOW  = '#F0E442'
BLUE    = '#0072B2'
RED     = '#D55E00'
MAGENTA = '#CC79A7'
CBF_COLORS = c(RED, CYAN, GREEN, YELLOW, MAGENTA, BLUE, ORANGE, GREY, BLACK)
OVR_CLRS = circlize::colorRamp2(seq(-4, 0, 0.1), rev(viridis::cividis(41)))
STT_CLRS = circlize::colorRamp2(c(-0.5, 0, 0.5), c(CYAN, 'gray95', MAGENTA))

PALETTE2COLORS = list(
    'Set1'       = RColorBrewer::brewer.pal(n=9, name='Set1'),
    'Set2'       = RColorBrewer::brewer.pal(n=8, name='Set2'),
    'Tableau 10' = ggthemes::tableau_color_pal(palette='Tableau 10')(10),
    'Tableau 20' = ggthemes::tableau_color_pal(palette='Tableau 20')(20),

    ## https://github.com/scverse/scanpy/blob/master/scanpy/plotting/palettes.py
    ## vega_10
    'Scanpy10' = c( '#1f77b4', '#ff7f0e', '#279e68', '#d62728', '#aa40fc',
                    '#8c564b', '#e377c2', '#7f7f7f', '#b5bd61', '#17becf' ),
    ## default_20
    'Scanpy20' = c( '#1f77b4', '#ff7f0e', '#279e68', '#d62728', '#aa40fc',
                    '#8c564b', '#e377c2', '#b5bd61', '#17becf', '#aec7e8',
                    '#ffbb78', '#98df8a', '#ff9896', '#c5b0d5', '#c49c94',
                    '#f7b6d2', '#dbdb8d', '#9edae5', '#ad494a', '#8c6d31' )
)
