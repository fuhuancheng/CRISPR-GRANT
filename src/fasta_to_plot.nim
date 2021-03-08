import ggplotnim, seqmath, tables, chroma
import os, strutils
import ginger
import parseopt
import os
import sequtils

const help = """
usage: fasta_to_plot -i:fasta_file -o:output_plot_name

-i        input file, aligned fasta format.
-o        output plot file name, the extension indicate the plot format (.pdf, .png).
--width   set output figure width.
--height  set output figure height.
--diff    make additional plot with only highlighting different bases.
-h        print this help and exit.
"""

var
  input: string
  output: string
  width: float
  height: float
  diffAign = false

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
      input = val
    of "o", "output":
      output = val
    of "width":
      width = parseFloat(val)
    of "height":
      height = parseFloat(val)
    of "diff":
      diffAign = true
  else: assert(false) #not needed

# test if input file exits
if not existsFile(input):
  echo(help, "\n")
  echo("Cannot open input file.")
  quit(errorcode = 1)

# test if output exits or has an extension name
if output == "":
  output = input & ".pdf"

if splitFile(output).ext == "":
  output &= ".pdf"

# convert hex to color table
proc hexToBaseCol(a, t, g, c: string, del = "#EEEEEE"): Table[string, Color] =
  result = {"A": parseHtmlHex(a),
           "T": parseHtmlHex(t),
           "G": parseHtmlHex(g),
           "C": parseHtmlHex(c),
           "a": parseHtmlHex(a),
           "t": parseHtmlHex(t),
           "g": parseHtmlHex(g),
           "c": parseHtmlHex(c),
           "-": parseHtmlHex(del)}.toTable()

# custom colors of nucleitide
const
  baseCol = hexToBaseCol(a = "#FF6E6E", 
                        t = "#64C8AA", 
                        g = "#8096FA", 
                        c = "#FADC96")
  # color palette same as CRISPResso
  crispressoCol = hexToBaseCol(a = "#91E185", 
                              t = "#C3ADDF", 
                              g = "#FAF19D", 
                              c = "#FFC38A")

# TODO: add another alignment color palette: if different to ref, colored background, else white background.
const alignCol = {"same": parseHtmlHex("#FFFFFF"), 
                  "diff": parseHtmlHex("#14BCBC"), 
                  "del": parseHtmlHex("#D6D8D8"), 
                  "A": parseHtmlHex("#000000"), 
                  "T": parseHtmlHex("#000000"), 
                  "G": parseHtmlHex("#000000"), 
                  "C": parseHtmlHex("#000000"), 
                  "-": parseHtmlHex("#FF3333")}.toTable()

type
  SEQUENCE = object
    name: string
    sequence: string
  FASTA = object
  # used for ggplotnim
    position: seq[int]
    nameIdx: seq[int]
    names: seq[string]
    bases: seq[string]
  AlignNotion = object
    x: seq[int]
    nameIdx: seq[int]
    text: seq[string]
  SameNotion = object
    pos: seq[int]
    nameIdx: seq[int]
    baseDiff: seq[string]

proc parseFasta(str: string): seq[SEQUENCE] =
  let fasta = str.split(">")[1..^1]
  for f in fasta:
    let fs = f.split("\n")
    var tmp: SEQUENCE
    tmp.name = fs[0]
    tmp.sequence = join(fs[1..^1], "").toUpper
    result.add(tmp)

proc fastaToSeq(fasta: seq[SEQUENCE]): FASTA =
  for j in 0..<fasta.len:
    var f = fasta[j]
    for i in 0..<f.sequence.len:
      result.position.add(i)
      result.nameIdx.add(j)
      result.names.add(f.name)
      result.bases.add(f.sequence[i..i])

proc fastaNote(fasta: seq[SEQUENCE]): AlignNotion =
  for f in 0..<fasta.len:
    result.x.add(0)
    result.nameIdx.add(f)
    result.text.add(fasta[f].name)

# test if bases same as in reference
proc sameBase(base: seq[char]): seq[string] =
  if all(base, proc(c: char): bool = return c == base[0]):
    result = sequtils.repeat("same", base.len)
  else:
    for i in base:
      if i == base[0] or i == '-':
        result.add("del")
      else:
        result.add("diff")

proc fastaAlignSame(fasta: seq[SEQUENCE]): SameNotion =
  for i in 0..<fasta[0].sequence.len:
    var bases: seq[char]
    for j in 0..<fasta.len:
      bases.add(fasta[j].sequence[i])
      result.pos.add(i)
      result.nameIdx.add(j)
    for b in sameBase(bases):
      result.baseDiff.add(b)

let
  alignFasta = open(input).readAll()
  savePlotFile = output

  fasta = parseFasta(alignFasta)
  fastaSeq = fastaToSeq(fasta)
  fastaDf = seqsToDf({"position": fastaSeq.position,
                    "names": fastaSeq.names,
                    "nameIdx": fastaSeq.nameIdx,
                    "bases": fastaSeq.bases})

  label = fastaNote(fasta)
  labelDf = seqsToDf({"x": label.x,
                      "num": label.nameIdx,
                      "label": label.text})

  diffLabel = fastaAlignSame(fasta)
  diffLabelDf = seqsToDf({"x": diffLabel.pos,
                          "num": diffLabel.nameIdx,
                          "diff": diffLabel.baseDiff})

# TODO: remove x.ticks y.ticks;
# set plot width and height automatically: at least col.num * 20, base.num * 10
if width == 0:
  width = toFloat(max(fastaSeq.position))*30.0
if height == 0:
  height = toFloat(max(fastaSeq.nameIdx))*30.0

const
  tileOffset = 0.5
  ylabelOffset = 5 # TODO: to set a more resonable value
  seqNameFont = some(font(family = "monospace"))
  # subplotAdd = 100

# TODO: add nucleitide position to x axis: at each 10s' site.

let
  align = ggplot(fastaDf, aes("position", "nameIdx")) +
  geom_tile(aes(x = f{`position` - tileOffset}, fill = "bases")) +
  geom_text(data = fastaDf, aes = aes(x = "position", y = "nameIdx", text = "bases")) +
  geom_text(data = labelDf, aes = aes(x = f{`x` - ylabelOffset}, y = "num", text = "label"), alignKind = taRight, font = seqNameFont) + # add custom seqeuence name
  # scale_x_discrete() +
  scale_x_continuous() +
  scale_y_discrete() +
  theme_void() +
  # theme_opaque() +
  legendPosition(1.0, 1.0) + # remove legend
  scale_fill_manual(baseCol) # use costum fill color

align + ggsave(savePlotFile, width, height)

if diffAign == true:
  # another alignment visualization type: only show up different bases
  let align2 = ggplot(fastaDf, aes("position", "nameIdx")) +
  geom_tile(data = diffLabelDf, aes = aes(x = f{`x` - tileOffset}, y = "num", fill = "diff")) +
  geom_text(data = fastaDf, aes = aes(x = "position", y = "nameIdx", text = "bases", color = "bases")) +
  geom_text(data = labelDf, aes = aes(x = f{`x` - ylabelOffset}, y = "num", text = "label"), alignKind = taRight) + # add custom seqeuence name
  scale_x_continuous() +
  scale_y_discrete() +
  theme_void() +
  legendPosition(1.0, 1.0) + # remove legend
  scale_fill_manual(alignCol) + # use costum fill color
  scale_color_manual(alignCol) # use costum text color

  # TODO: add alternative plot
  align2 + ggsave(savePlotFile & ".pdf", width, height)

## below used subplot to avoid ylabel exceeding
#   alignPlot = ggcreate(align, width = width, height = height)

# var plt = initViewport(wImg = width + subplotAdd, hImg = height)
# plt.layout(cols = 2, rows = 1, colwidths = @[quant(subplotAdd, ukPoint), quant(0.0, ukRelative)])
# plt.embedAt(1, alignPlot.view)
# plt.draw(savePlotFile)
