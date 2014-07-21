
from thermocouples_reference import thermocouples

fts = "{:<10s} {:<16s} {:<37s} {: >6.1f} {: >6.1f}"
# Table
print "SORTED BY MIN"
tlist = sorted([(t.minT_C,t.maxT_C,k,t) for k,t in thermocouples.items()])
for min,max,k,t in tlist:
    print fts.format(
           k,   t.type,t.composition,  min,    max)

# Table
print "SORTED BY MAX"
tlist = sorted([(t.maxT_C,t.minT_C,k,t) for k,t in thermocouples.items()])
for max,min,k,t in tlist:
    print fts.format(
           k,   t.type,t.composition,  min,    max)


