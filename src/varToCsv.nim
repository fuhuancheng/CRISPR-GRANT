# Convert varscan2 output to csv
import parseopt
import strutils
import os

const help = """

Usage: varToCsv -i=varscan2_output -o=csv_file

-h, --help      print this help and exit.
-i, --input     input file, varscan2 output file.
-o, --output    output csv file name.

Output csv columns: Chrom,Position,Freq,P-value

"""

var
  varscan: string
  output: string

# command arguments
for kind, key, val in parseopt.getopt():
  case kind
  of cmdArgument:
    continue
  of cmdShortOption, cmdLongOption:
    case key
    of "h", "help":
      echo help
      quit()
    # of "v", "version":
    #   echo version
    #   quit()
    of "i", "input":
      varscan = val
    of "o", "output":
      output = val
  else: assert(false) #not needed

# test if input file not exits, quit.
if not existsFile(varscan):
  echo(help, "\n")
  echo("Cannot open input file.")
  quit(errorcode = 1)

let outputCsv = open(output, mode = fmWrite)

outputCsv.write("Chrom,Position,Freq,P-value\n")

# TODO: rewrite a more elegant way to parse and handle varscan2 output format, like the header and the last blank line
# if there is a header, skip
var lineStart = 1
# Chrom,Position,Ref,Cons,Reads1,Reads2,VarFreq,Strands1,Strands2,Qual1,Qual2,Pvalue,MapQual1,MapQual2,Reads1Plus,Reads1Minus,Reads2Plus,Reads2Minus,VarAllele
let file = open(varscan).readAll().split("\n")[lineStart..^1]
for l in file:
  if l != "":
    let line = l.split("\t")
    outputCsv.write(line[0], ",")
    outputCsv.write(line[1], ",")
    outputCsv.write(line[6], ",")
    outputCsv.write(line[11], "\n")

outputCsv.close()
