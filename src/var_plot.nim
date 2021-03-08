# visualization of output from varscan2
import ggplotnim, chroma
import strutils
import ginger
import parseopt
import os

const help = """

usage: var_plot -i=var-frequency.csv -o=plot.pdf

-h, --help      print this help and exit.
-i, --input     input file, var frequency csv file, containing 3 columns: Chrom, Position, Frequency.
-o, --output    output plot file name, the extension indicate the plot format (.pdf, .png).
--width         set output figure width.
--height        set output figure height.
--only          seqeuence name, only output sequences aligned to this reference.

"""

var
  varscan: string
  output: string
  width: float
  height: float
  refName: string

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
    of "width":
      width = parseFloat(val)
    of "height":
      height = parseFloat(val)
    of "only":
      refName = val
  else: assert(false) #not needed

# test if input file not exits, quit.
if not existsFile(varscan):
  echo(help, "\n")
  echo("Cannot open input file.")
  quit(errorcode = 1)

# test if output exits or has an extension name
if output == "":
  output = varscan & ".pdf"

if splitFile(output).ext == "":
  output &= ".pdf"

type
  Varcount = object
    chrome: seq[string]
    position: seq[int]
    freq: seq[float]

# TODO: a better way to parse varscan2 format
# In varscan2 result: when Var is N, Freq will be '-'
proc parsePercent(percent: string): float =
  var pct = replace(percent, "%", "")
  pct = replace(pct, "-", "0")
  result = parseFloat(pct) / 100.0

# TODO: rewrite a more elegant way to parse and handle varscan2 output format, like the header and the last blank line
proc parseVarCsv(filename: string, header = true, only = ""): Varcount =
  var lineStart = 0
  # if there is a header, skip
  if header:
    lineStart = 1
  let file = open(filename).readAll().split("\n")[lineStart..^1]
  for l in file:
    if l != "":
      let line = l.split(",")
      let freq = parsePercent(line[2])
      if only == "" or only == line[0]:
        result.chrome.add(line[0])
        result.position.add(parseInt(line[1]))
        result.freq.add(freq)

# round a float number to near multiples of 10
proc roundTen(x: float): float =
  result = (floor(x / 10) + 1) * 10

let
  varscanPos = parseVarCsv(varscan, header = true, only = refName)
  varscanDf = seqsToDf({"chrome": varscanPos.chrome,
                        "position": varscanPos.position,
                        "frequency": varscanPos.freq})
  # add axis to plot, for by default there is none
  yMax = 100.0
  axis_y = seqsToDf({"x": [0.0, 0.0], "y": [0.0, yMax]})
  xMax = roundTen(toFloat(varscanPos.position.max))
  axis_x = seqsToDf({"x": [0.0, xMax], "y": [0.0, 0.0]})

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
  width = xMax * 2
if height == 0:
  height = 200

# let width larger than 200 in case reference length too small
# let height at least 150
width = max(width, 200)
height = max(height, 150)

ggplot(varscanDf, aes("position", f{`frequency`*100.0})) + 
geom_line() + 
xlim(0, xMax) + 
ylim(0.0, yMax) + 
geom_line(data = axis_y, aes = aes(x = "x", y = "y")) + 
geom_line(data = axis_x, aes = aes(x = "x", y = "y")) + 
# theme_void() + 
theme_classic() + 
xlab("Position (bp)") + 
ylab("Frequency (%)") + 
scale_y_continuous() + 
ggsave(output, width = width, height = height)
