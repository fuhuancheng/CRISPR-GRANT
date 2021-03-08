import tables, os, strutils, algorithm, strformat
import parseopt

const
  help = """
usage: sam_count -i:csv-file -o:output-fasta

-i            input file, SAM format.
-o            output csv file name.
-h            print this help and exit.
--reference   reference file, for counting reads exactly the same as it.
--countTab    counting table of all reads type: total, not-mapped and mapped to reference, unmodified and modified.
"""

var 
  sam: string
  outFile: string
  refFile: string
  reference: seq[string]
  countFile: string

for kind, key, val in getopt():
  case kind
  of cmdArgument:
    continue
  of cmdShortOption, cmdLongOption:
    case key
    of "h", "help":
      echo(help)
      quit()
    # of "v", "version":
    #   echo version
    #   quit()
    of "i", "input":
      sam = val
    of "o", "output":
      outFile = val
    of "reference":
      refFile = val
    of "countTab":
      countFile = val
  else: assert(false) #not needed

# test if input and output files not given
if sam == "" or outFile == "":
  echo(help)
  echo("Input file or output name wasn't given.")
  quit()

# extract reference sequence
proc parseFasta(fa: string): seq[string] =
  let fasta = fa.split(">")[1..^1]
  for f in fasta:
    let fs = f.split("\n")
    var fsSeq = join(fs[1..^1], "").toUpper
    result.add(fsSeq)

if refFile != "":
  let refFasta = open(refFile).readAll()
  reference = parseFasta(refFasta)

var
  csv_out = open(outFile, mode = fmWrite)
  totalReads = 0
  no_reference = 0
  mappedReads = 0
  unmodified = 0
  modified = 0
  output = initTable[(string, string), int]()

for l in lines(sam):
  if l[0] != '@':
    totalReads += 1
    let line = l.split("\t")
    let fastq_seq = line[9]
    let ref_name = line[2]
    if ref_name != "*":
      # counting reads mapped to reference
      mappedReads += 1
      # counting reads unmodified and modified seperately
      if fastq_seq in reference:
        unmodified += 1
      else:
        modified += 1
      # add mapped reads to dict
      try:
        output[(fastq_seq, ref_name)] += 1
      except:
        output[(fastq_seq, ref_name)] = 1
    else:
      no_reference += 1

var output_s: seq[(int, (string, string))]
for (k, v) in pairs(output):
  output_s.add((v, k))

output_s.sort(system.cmp, order = SortOrder.Descending)

# write to file
csv_out.write(&"#Number of reads not mapped to reference: {no_reference}\n")
csv_out.write("#sequence,reference,count,percent\n")
for s in output_s:
  var
    sequence = s[1][0]
    ref_name = s[1][1]
    count = s[0]
    readsRatio:float = (s[0]/mappedReads)*100
  csv_out.write(fmt("{sequence},{ref_name},{count},{readsRatio:2.2f}%\n"))
csv_out.close()

# write all reads counts to file
if countFile != "":
  var allReadsCount = &"""
Total reads: {totalReads}
    Not-mapped reads: {no_reference}
    Mapped reads: {mappedReads}
      Unmodified reads: {unmodified}
      Modified reads: {modified}
"""
  writeFile(countFile, allReadsCount)
