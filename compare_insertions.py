import pybedtools
import argparse



#def find_intersection(original_insertion, detected_insertion):
     #original_insertions.window(detected_insertion, w=100)

#print(intersection)
#header



if __name__ == '__main__':

    parser = argparse.ArgumentParser(description="Arguments for de program")
    parser.add_argument("original_ins", metavar="original_insertions", type=argparse.FileType('r'), nargs="?", help="Original Insertions. Output .bed file from sequence_simulation.py")
    parser.add_argument("detected_ins", metavar="detected_insertions", type=argparse.FileType('r'), nargs='?', help="Detected Insertions. Output .vcf from MOBSTER")
    parser.add_argument("title", metavar="title", type=str, nargs='?',
                        help="title as ID")



    args = parser.parse_args()

original_insertions = pybedtools.BedTool(args.original_ins)
detected_insertions = pybedtools.BedTool(args.detected_ins)
#header = args.header



intersection = original_insertions.window(detected_insertions, w=100)
intersection.saveas(args.title+".comparison.vcf")
