

haploplot.vp <- function(y = unit(5/6, 'npc'), ...){
  grid::viewport(height = unit(1, 'snpc'), width = unit(1,'snpc'), angle = 45, y = y,
                 just = 'center', ...)
}

mapping.vp <- function(x, y = unit(5/6, 'npc'), just = 'center', xscale, ...){

  grid::viewport(height = unit(0.05, 'npc'), y = y, just = just,
                 width = unit(sqrt(2), 'snpc'), xscale = xscale, yscale = c(0,1), ...)
}

center.vp <- function(xscale = range(x, na.rm = T), y = unit(1/6, 'npc'), ...){

  grid::viewport(name = 'vp_center',height = unit(0.15, 'npc'), y = y, just = 'center',
                 width = unit(sqrt(2), 'snpc'), xscale = xscale, ...)
}

genes2.vp <- function(xscale = range(x, na.rm = T), yscale = range(y, na.rm = T)){

  grid::viewport(name = 'vp_genes', height = unit(0.75,'npc'),  width = unit(1, 'npc'),
                 xscale = xscale, yscale = yscale+c(0,0.5))
}
