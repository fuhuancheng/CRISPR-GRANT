# Visualization of reads counts
import ggplotnim, chroma
import ginger
import parseopt
import os
import strutils
import tables

const help = """

Usage: reads_count_plot -i=count-file.txt -o=plot.pdf

-h, --help      print this help and exit.
-i, --input     input text file, reads counting output.
-o, --output    output plot file name, the extension indicates the plot format (.pdf, .png).
--width         output figure width.
--height        output figure height.

"""

var
  countFile: string
  output: string
  width: float
  height: float

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
      countFile = val
    of "o", "output":
      output = val
    of "width":
      width = parseFloat(val)
    of "height":
      height = parseFloat(val)
  else: assert(false) #not needed

# if input file not exits, quit.
if not existsFile(countFile):
  echo(help, "\n")
  echo("Cannot open input file.")
  quit(errorcode = 1)

# test if output exits or has an extension name
if output == "":
  output = countFile & ".pdf"

if splitFile(output).ext == "":
  output &= ".pdf"

# custom colors of reads types
const typeCol = {"Total reads": parseHtmlHex("#EA4335"), #red
           "Not-mapped reads": parseHtmlHex("#ff751b"), #orange
           "Mapped reads": parseHtmlHex("#4285F4"), #blue
           "Unmodified reads": parseHtmlHex("#6269d8"), #slate blue
           "Modified reads": parseHtmlHex("#6964ff")}.toTable() #light blue

type
  CountTab = object
    ids: seq[int]
    readsType: seq[string]
    readsNumber: seq[int]

proc parseCountTable(filename: string): CountTab =
  let counts = open(filename).readAll().strip().split("\n")
  var i = 0
  for l in counts:
    let line = l.split(":")
    result.ids.add(i)
    result.readsType.add(line[0].strip())
    result.readsNumber.add(parseInt(line[1].strip()))
    i += 1

# TODO: make CountTab, parseCountTable() use openarray data for easy parsing
let
  readsCounts = parseCountTable(countFile)
  countDf = seqsToDf({"ids": readsCounts.ids,
                      "type": readsCounts.readsType,
                      "number": readsCounts.readsNumber})

# theme to remove background
func theme_classic(): Theme =
  result = Theme(canvasColor: some(parseHex("FFFFFF")), 
                # hideTicks: some(true), 
                # hideTickLabels: some(true), 
                # hideLabels: some(true), 
                plotBackgroundColor: some(parseHex("FFFFFF"))
                )

# TODO: make the plot width and height more scalable
if width == 0:
  width = 600
if height == 0:
  height = 300

proc xLabel(x: int): string =
  result = readsCounts.readsType[x]

ggplot(countDf, aes("ids", "number")) + 
geom_bar(aes(fill = "type"), stat = "identity", position = "identity") + 
scale_fill_manual(typeCol) + 
# theme_void() + 
theme_classic() + 
xlab("Reads Types") + 
ylab("Reads Number") + 
scale_x_discrete(labels = proc(x: Value): string =
                    xLabel(parseInt($x))) +
scale_y_continuous(labels = proc(x: float): string =
                    x.formatFloat(ffDecimal, 0)) + 
legendPosition(0.0, 1.0) + # remove legend
ggsave(output, width = width, height = height)
