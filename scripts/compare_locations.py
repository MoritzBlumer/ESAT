import pybedtools

original_insertions = pybedtools.BedTool("10Mbp_insertions.sorted.bed")
discovered_insertions = pybedtools.BedTool("test.vcf")


intersection = original_insertions.window(discovered_insertions, w=100)
print(intersection)
#header
intersection.saveas("Intersections.vcf")

