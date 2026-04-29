import sys
import subprocess   
import HTSeq, time, argparse, re
from Bio.Emboss.Applications import NeedleCommandline
from Bio import AlignIO
from pybedtools import BedTool
import pandas as pd

parser = argparse.ArgumentParser(
    description = 'This script takes in an annotation file of human Alu and\
        MIR elements and the corresponding human reference genome sequence\
    to produce the index for SINEs_find.py. The index file is the annotation \
        input file with added columns reporting the lenght, start and end of \
            the alignment of the element sequence to its consensus sequence \
                on the reference genome',
    epilog = 'Written by Davide Carnevali davide.carnevali@crg.eu')
parser.add_argument("annotation", help="Annotation file in GTF format \
                    (either gzipped or not). Should refer to the same version \
                        of the human reference genome file")
parser.add_argument("genome", help="Human reference genome sequence. Should refer \
                    to the same version of the annotation file")
parser.add_argument("output", help="output filename")
args = parser.parse_args()

genome = BedTool(args.genome)
B2_Mm1a = 'GGGGCTGGTGAGATGGCTCAGTGGGTAAGAGCACCCGACTGCTCTTCCGAAGGTCCGGAGTTCAAATCCCAGCAACCACATGGTGGCTCACAACCATCCGTAACGAGATCTGACTCCCTCTTCTGGAGTGTCTGAAGACAGCTACAGTGTACTTACATATAATAAATAAATAAATCTTTAAAAAAAAAAAAAA'
B2_Mm1t = 'GGGGCTGGTGAGATGGCTCAGCGGGTAAGAGCACCCGACTGCTCTTCCGAAGGTCCGGAGTTCAAATCCCAGCAACCACATGGTGGCTCACAACCATCCGTAACGAGATCTGACGCCCTCTTCTGGTGTGTCTGAAGACAGCTACAGTGTACTTACATATAATAAATAAATAAATCTTTAAAAAAAAAAAAAA'
B2_Mm2 = 'GGGGCTGGAGAGATGGCTCAGCGGTTAAGAGCACTGACTGCTCTTCCAGAGGTCCTGAGTTCAATTCCCAGCAACCACATGGTGGCTCACAACCATCTGTAATGGGATCTGATGCCCTCTTCTGGTGTGTCTGAAGACAGCTACAGTGTACTCACATACATAAATAAATAAATAAATCTTTAAAAAAAAAAAAAA'
B3 = 'GGGGCTGGAGAGATGGCTCAGCGGTTAAGAGCACTGGCTGCTCTTCCAGAGGACCCGGGTTCGATTCCCAGCACCCACATGGCGGCTCACAACCGTCTGTAACTCCAGTTCCAGGGGATCCGACGCCCTCTTCTGGCCTCCGCGGGCACCAGGCACGCACGTGGTGCACAGACATACATGCAGGCAAAACACCCATACACATAAAATAAAAATAAA'
B3A = 'GGGGCTGGAGAGATGGCTCAGCGGTTAAGAGCACTTGCTGCTCTTGCAGAGGACCCGAGTTCGGTTCCCAGCACCCACGTCGGGCGGCTCACAACCGCCTGTAACTCCAGCTCCAGGGGATCCGACGCCCTCTTCTGGCCTCCGCGGGCACCCGCACNCACACGCGCACACACACACACACAATAAAATAAAAATAAA'
B1_Mm = 'AGCCGGGCGTGGTGGCGCACGCCTTTAATCCCAGCACTCGGGAGGCAGAGGCAGGCGGATTTCTGAGTTCGAGGCCAGCCTGGTCTACAAAGTGAGTTCCAGGACAGCCAGGGCTATACAGAGAAACCCTGTCTCGAAAAAACAAAA'
B1_Mus1 = 'AGCCGGGCGGTGGTGGCGCACGCCTTTAATCCCAGCACTTGGGAGGCAGAGGCAGGCGGATTTCTGAGTTCGAGGCCAGCCTGGTCTACAGAGTGAGTTCCAGGACAGCCAGGGCTACACAGAGAAACCCTGTCTCGAAAAAACAAAA'
B1_Mus2 = 'AGCCGGGCGTGGTGGCGCACGCCTTTAATCCCAGCACTCGGGAGGCAGAGGCAGGCGGATTTCTGAGTTCGAGGCCAGCCTGGTCTACAGAGTGAGTTCCAGGACAGCCAGGGCTACACAGAGAAACCCTGTCTCGAAAAAACAAAA'
MIR = 'ACAGTATAGCATAGTGGTTAAGAGCACGGGCTCTGGAGCCAGACTGCCTGGGTTCGAATCCCGGCTCTGCCACTTACTAGCTGTGTGACCTTGGGCAAGTTACTTAACCTCTCTGTGCCTCAGTTTCCTCATCTGTAAAATGGGGATAATAATAGTACCTACCTCATAGGGTTGTTGTGAGGATTAAATGAGTTAATACATGTAAAGCGCTTAGAACAGTGCCTGGCACATAGTAAGCGCTCAATAAATGTTAGCTATTATT'
MIRb = 'CAGAGGGGCAGCGTGGTGCAGTGGAAAGAGCACGGGCTTTGGAGTCAGACAGACCTGGGTTCGAATCCCGGCTCTGCCACTTACTAGCTGTGTGACCTTGGGCAAGTTACTTAACCTCTTGAGCCTCAGTTTCCTCATCTGTAAAATGGGGATAATAATACCTACCTCGCAGGGTTGTTGTGAGGATTAAATGAGATAATGCATGTAAAGCGCTTAGCACAGTGCCTGGCACATAGTAAGCGCTCAATAAATGGTAGCTCTATTATT'


start_time = time.time()
annotation = HTSeq.GFF_Reader(args.annotation)
alu_list = []
char = re.compile('-*')
char2 = re.compile('[-NATGCatgc]*')
    
def needle(chrom, start, end, name, score, strand):
    item=BedTool([(chrom, start, end, name, score, strand)])
    item = item.sequence(fi=genome, s=True)
    temp = open(item.seqfn).read().split('\n')[1]
    
    if name == "B2_Mm1a":
        sine_length = 193
        needle_cline = NeedleCommandline(asequence="asis:"+B2_Mm1a, bsequence="asis:"+temp, gapopen=10, gapextend=0.5, outfile='tmpaln.txt')
        
    elif name == "B2_Mm1t":
        sine_length = 193
        needle_cline = NeedleCommandline(asequence="asis:"+B2_Mm1t, bsequence="asis:"+temp, gapopen=10, gapextend=0.5, outfile='tmpaln.txt')
                    
    elif name == "B2_Mm2":
        sine_length = 195
        needle_cline = NeedleCommandline(asequence="asis:"+B2_Mm2, bsequence="asis:"+temp, gapopen=10, gapextend=0.5, outfile='tmpaln.txt')

    elif name == "B3":
        sine_length = 216
        needle_cline = NeedleCommandline(asequence="asis:"+B3, bsequence="asis:"+temp, gapopen=10, gapextend=0.5, outfile='tmpaln.txt')

    elif name == "B3A":
        sine_length = 198
        needle_cline = NeedleCommandline(asequence="asis:"+B3A, bsequence="asis:"+temp, gapopen=10, gapextend=0.5, outfile='tmpaln.txt')
    
    elif name == "B1_Mm":
        sine_length = 147
        needle_cline = NeedleCommandline(asequence="asis:"+B1_Mm, bsequence="asis:"+temp, gapopen=10, gapextend=0.5, outfile='tmpaln.txt')
    
    elif name == "B1_Mus1":
        sine_length = 148
        needle_cline = NeedleCommandline(asequence="asis:"+B1_Mus1, bsequence="asis:"+temp, gapopen=10, gapextend=0.5, outfile='tmpaln.txt')
    
    elif name == "B1_Mus2":
        sine_length = 147
        needle_cline = NeedleCommandline(asequence="asis:"+B1_Mus2, bsequence="asis:"+temp, gapopen=10, gapextend=0.5, outfile='tmpaln.txt')
        
    elif name == "MIR":
        sine_length = 262
        needle_cline = NeedleCommandline(asequence="asis:"+MIR, bsequence="asis:"+temp, gapopen=10, gapextend=0.5, outfile='tmpaln.txt')    
        
    elif name == "MIRb":
        sine_length = 267
        needle_cline = NeedleCommandline(asequence="asis:"+MIRb, bsequence="asis:"+temp, gapopen=10, gapextend=0.5, outfile='tmpaln.txt')    
        
    else:
        raise ValueError(f"Unknown SINE name: {name}")
    
    # Run the alignment
    child = subprocess.Popen(str(needle_cline), stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=(sys.platform!="win32"))
    child.wait()
    
    # Read the alignment
    align = AlignIO.read('tmpaln.txt', "emboss")
    aln_start = char.search(str(align[1, :].seq)).end()
    aln_end = char2.search(str(align[1, :].seq)).end()
        
    return (aln_start, aln_end, sine_length)

# To mantain GTF notation, which is 1 based, we add +1 to element.iv.start
for element in annotation:
    aln_start, aln_end, sine_length = needle(element.iv.chrom, element.iv.start, element.iv.end, element.attr['gene_id'], int(element.score), element.iv.strand)
    if element.iv.strand == "+":
        alu_list.append([element.iv.chrom, element.source, "exon", element.iv.start + 1, element.iv.end, "0", element.iv.strand, ".", "gene_id " + "\""+element.attr['transcript_id']+"\"; " + "transcript_id " + "\""+element.attr['transcript_id']+"\";",element.iv.start - aln_start, (element.iv.start - aln_start) + aln_end])
    else:
        alu_list.append([element.iv.chrom, element.source, "exon", element.iv.start + 1, element.iv.end, "0", element.iv.strand, ".", "gene_id " + "\""+element.attr['transcript_id']+"\"; " + "transcript_id " + "\""+element.attr['transcript_id']+"\";",(element.iv.end + aln_start) - aln_end, element.iv.end + aln_start])

final_list = pd.DataFrame(alu_list)
final_list.to_csv(args.output,sep="\t", header=False,index=False,doublequote=False,quotechar='\'',escapechar='')

print("Finished!")
print("Time elapsed {}".format(time.time() - start_time))
