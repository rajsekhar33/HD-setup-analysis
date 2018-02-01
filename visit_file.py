import sys
sys.path.append("/home/pkg/pub/gnu/dellvs/visit2_10_0.linux-x86_64/2.10.0/linux-x86_64/lib/site-packages")
import visit
visit.Launch()
 
visit.OpenDatabase("/home/pkg/pub/gnu/dellvs/visit2_10_0.linux-x86_64/data/noise.silo")
visit.AddPlot("Pseudocolor", "hardyglobal")
visit.DrawPlots()
visit.SaveWindow()
