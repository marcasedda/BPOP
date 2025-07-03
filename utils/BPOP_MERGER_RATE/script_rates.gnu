
set log xy

p   sprintf("SET%s/Merger_rates_SET%s_none_massive_no.txt", Nset, Nset) u 1:2 smooth bezier w l lw 4 t 'IB'     #, '' u 1:2 w lp notitle
rep sprintf("SET%s/Merger_rates_SET%s_young_massive_no.txt", Nset, Nset) u 1:2 smooth bezier w l lw 4 t 'YC'    #, '' u 1:2 w lp notitle
rep sprintf("SET%s/Merger_rates_SET%s_globular_massive_no.txt", Nset, Nset) u 1:2 smooth bezier w l lw 4 t 'GC' #, '' u 1:2 w lp notitle
rep sprintf("SET%s/Merger_rates_SET%s_nuclear_massive_no.txt", Nset, Nset) u 1:2 smooth bezier w l lw 4 t 'NC'  #, '' u 1:2 w lp notitle
