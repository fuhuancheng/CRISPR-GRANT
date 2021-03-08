import ggplotnim, chroma
import strutils
import ginger
import parseopt
import os
import tables

const help = """
usage: snp_plot -i=var.snp -o=plot.pdf

-h, --help      print this help and exit.
-i, --input     input file, var frequency csv file, containing 3 columns: Chrom, Position, Frequency.
-o, --output    output plot file name, the extension indicate the plot format (.pdf, .png).
--width         set output figure width.
--height        set output figure height.
"""

var
  var_snp: string
  output: string
  width = 500.0
  height = 300.0

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
      var_snp = val
    of "o", "output":
      output = val
    of "width":
      width = parseFloat(val)
    of "height":
      height = parseFloat(val)
  else: assert(false) #not needed

# test if input file not exits, quit.
if not existsFile(var_snp):
  echo(help, "\n")
  echo("Cannot open input file.")
  quit(errorcode = 1)

# test if output exits or has an extension name
if output == "":
  output = var_snp & ".pdf"

if splitFile(output).ext == "":
  output &= ".pdf"

# type
#   SnpCount = object
#     AT, AG, AC, TA, TG, TC, GA, GT, GC, CA, CT, CG: int

proc parseSNP(filename: string, header = true): Table[string, int] =
  result = {"A>T": 0, "A>G": 0, "A>C": 0, "T>A": 0, "T>G": 0, "T>C": 0, "G>A": 0, "G>T": 0, "G>C": 0, "C>A": 0, "C>T": 0, "C>G": 0}.toTable
  var lineStart = 0
  # if there is a header, skip
  if header:
    lineStart = 1
  let file = open(filename).readAll().split("\n")[lineStart..^1]
  for l in file:
    if l != "":
      let
        line = l.split("\t")
        refBase = line[2]
        varBase = line[18].strip()
        varCount = parseInt(line[16])
        baseMut = refBase & ">" & varBase
      result[baseMut] += varCount

var
  baseMutationSeqSub: seq[string]
  baseMutationSeqNum: seq[int]
  baseMutationSeqType: seq[string]

let baseMutation = parseSNP(var_snp, header = true)

for k in baseMutation.keys:
  baseMutationSeqSub.add(k)
  baseMutationSeqNum.add(baseMutation[k])
  baseMutationSeqType.add(k[0..0])

let  baseMutationDf = seqsToDf({"substitution": baseMutationSeqSub, "counts": baseMutationSeqNum, "type": baseMutationSeqType})

# theme to remove background
func theme_classic(): Theme =
  result = Theme(canvasColor: some(parseHex("FFFFFF")), 
                # hideTicks: some(true), 
                # hideTickLabels: some(true), 
                # hideLabels: some(true), 
                plotBackgroundColor: some(parseHex("FFFFFF"))
                )

# set min width and height
width = max(width, 400.0)
height = max(height, 250.0)

ggplot(baseMutationDf, aes("substitution", "counts", fill = "type")) + 
geom_bar(stat = "identity", position = "identity") +
scale_y_continuous() +
theme_classic() +
xlab("Substitution", rotate = -45.0) +
ylab("Number of substitutions") +
legendPosition(1.0, 1.0) +
ggsave(output, width = width, height = height)
