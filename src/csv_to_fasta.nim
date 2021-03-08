import strutils
import strformat
import parseopt

const
  help = """
  usage: csv_to_fasta -i:csv-file -o:output-fasta

  -i input file, csv format.
  -n number of top sequences in input file for output, default 10.
  -o output fasta file name.
  -r reference file, fasta format. If provided, this sequence will be added to output fasta.
  --only $1 seqeuence name, only output sequences aligned to this reference.
  -h print this help and exit.
  """ % ["\t"]

var
  align_seq:seq[string]
  top_number = 10
  fasta: File
  csv_in: string
  outFile: string
  reference: string
  refName: string

for kind, key, val in getopt():
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
      csv_in = val
    of "o", "output":
      outFile = val
    of "r", "reference":
      reference = val
    of "n":
      top_number = parseInt(val)
    of "only":
      refName = val
  else: assert(false) #not needed

# test if input and output names were provided
if csv_in == "" or outFile == "":
  echo(help)
  echo("Please provide input and output file names.")
  quit()

fasta = open(outFile, mode = fmWrite)

# to align reads and percent in sequence names of output fasta
var
  percentWidth = 7
  readsWidth = 0
  firstReads = true

# add counts of each reference to sub-directory

for l in lines(csv_in):
  if l[0] != '#' and top_number > 0:
    let line = l.split(",")
    let line_seq = line[0].strip()
    let line_ref = line[1]
    var line_reads = line[2]
    if firstReads:
      firstReads = false
      readsWidth = len(line_reads)
    line_reads = align(line_reads, readsWidth)
    let line_percent = align(line[3], percentWidth)
    if refName == "" or refName == line_ref:
      align_seq.add(fmt(">{line_ref} {line_reads} {line_percent}\n{line_seq}"))
      top_number -= 1

# add reference to output, if possible
if reference != "":
  fasta.write(readFile(reference).strip() & "\n\n")

for s in align_seq:
  fasta.write(s, "\n")
fasta.close()