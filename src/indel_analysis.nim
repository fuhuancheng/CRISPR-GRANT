# call external tools for each step in analysis
# fastqc
# flash for merge pair-end reads
# bwa for fast mapping reads to reference
# needle for global alignment
# result visualization

import parseopt
import os
import osproc
import strutils
import strformat
import tables
import times

# version
const version = "v0.5.0"

# help
const help = """
This tool will mapping indel amplicons to reference library, and output basic statistics: number of aligned reads to each reference, alignment plot of top number reads, indel distribution along each reference.

Usage:

indel_analysis -1:fastq-1.gz -2:fastq-2.gz -r:reference -o:output

    -h, --help       get this help and exit.
    -1               fastq file 1.
    -2               fastq file 2.
    -r, --reference  amplicon reference file, fasta format.
    -o, --output     output file name.
    -t, --thread     number of threads.
    -n, --topNumber  number of top reads for alignment plot.
    --plotSNP        plot number of substitutions.
    -q               quality phred value in Phred quality score, default is 30.
    --fastp          additional parameters passed to fastp.
    --flash          additional parameters passed to FLASH.
    --bwa            additional parameters passed to BWA-MEM.
    --varscan        additional parameters passed to VarScan2.
    --javaPath       path to java, default using "java" within bin/jre11.
    -v, --version    print version ($#) and exit.
""" % [version]

# get tool directory
let bin_dir = getAppDir() & "/"

var
  # java path, using jre11-lite within
  java = bin_dir & "jre11/bin/java"
  # analysis type
  analysis_type = "amplicon"
  # use current time for default output file name
  output_dir = getCurrentDir()
  output = ""
  # output = getTime().format("yyyyMMdd-HHmmss")
  # outputTime = output
  # default threads
  threads = "1"
  # default quality score
  phred = "30"
  # default visualization sequences
  top_number = "10"
  # fastq file name
  fastq1:string
  fastq2:string
  # reference file name
  reference:string
  # additional parameters
  fastpArgs: string
  flashArgs: string
  bwaArgs: string
  varscanArgs: string
  # single end or pair end
  pairEnd = true
  # plot SNP
  plotSNP = false

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
    of "v", "version":
      echo version
      quit()
    of "o", "output":
      output = val
      output_dir = parentDir(output) & "/"
      output = splitPath(output).tail
      output_dir &= output & "/"
    of "t", "thread":
      threads = val
    of "r", "reference":
      reference = val
    of "1":
      fastq1 = val
    of "2":
      fastq2 = val
    of "q":
      phred = val
    of "n", "topNumber":
      top_number = val
    of "plotSNP":
      plotSNP = true
    of "fastp":
      fastpArgs = val
    of "flash":
      flashArgs = val
    of "bwa":
      bwaArgs = val
    of "varscan":
      varscanArgs = val
    of "javaPath":
      java = val
  else: assert(false) #not needed

# test if fastq files and reference were provided
if reference == "" or fastq1 == "":
  # echo(help)
  echo("Please provide valid reference and FASTQ files.")
  quit(errorcode = 2)

# whether single end according to the number of input files
if fastq1 != "" and fastq2 == "":
  pairEnd = false

# if output not given, use current time as output name
if output == "":
  output = getTime().format("yyyyMMdd-HHmmss")

# test if there were duplicate ref-names in reference
# to count ref-name, not ref-name and sequence
proc countFasta(file: string): Table[string, (string, int)] =
  let fasta = open(file).readAll().split(">")
  for l in fasta[1 .. ^1]:
    let line = l.split("\n")
    try:
      result[line[0]][1] += 1
    except:
      result[line[0]] = (join(line[1..^1], ""), 1)

let referenceCount = countFasta(reference)
# if referenceCount.len < 1:
#   echo("Reference is empty or format is unvalid.")
#   quit(errorcode = 2)
for k in referenceCount.keys:
  if referenceCount[k][1] > 1:
    echo(&"Sequence name '{k}' is duplicate in provided reference file.")
    quit(errorcode = 2)

# check if folder exists
if dirExists(output_dir):
  echo("Output directory exists.")
  quit(errorcode = 2)
else:
  createDir(output_dir)

# list of output files
type
  OutFileList = object
    qcReport, readsCount, topSeq, topAlignment, varFreq, sam, bam, consensus, snp: string

# proc NumberOutFile(outName: string): outFileList =
#   result.qcReport = "0.QC-report"
#   result.readsCount = &"1.{outName}_readsCount"
#   result.topSeq = &"{outName}_top"
#   result.topAlignment = &"2.{outName}_top_alignment"
#   result.varFreq = &"3.{outName}_varFreq"
#   result.sam = "Mapping"
#   result.bam = &"{outName}_sorted"
#   result.consensus = &"{outName}"
#   result.indel = &"{outName}"

var outputFileNames: OutFileList
outputFileNames.qcReport = "0.QC-report"
outputFileNames.readsCount = "1.readsCount"
outputFileNames.topSeq = "Top_sequences"
outputFileNames.topAlignment = "2.Top_sequences_alignment"
outputFileNames.varFreq = "3.varFreq"
outputFileNames.sam = "Mapping"
outputFileNames.bam = "4.Mapping_sorted"
outputFileNames.consensus = "Consensus"
outputFileNames.snp = "Base_variation"

var
  # external commands: fastp, flash, bwa, samtools, varscan2, mafft
  fastqc = fmt"{bin_dir}fastp/fastp --html {output_dir}{outputFileNames.qcReport}.html --json {output_dir}QC-report.json --qualified_quality_phred {phred} --in1 {fastq1} --in2 {fastq2} --out1 {fastq1}-qc.fq.gz --out2 {fastq2}-qc.fq.gz --thread {threads} {fastpArgs} "
  flash = bin_dir & "flash/flash --allow-outies --max-overlap=150 $1 --output-directory $2 -o $3 -z -t $4 $5-qc.fq.gz $6-qc.fq.gz" % [flashArgs, output_dir, output, threads, fastq1, fastq2]
  reference_index = fmt"{bin_dir}bwa/bwa index $1" #% [reference]
  mapping = fmt"{bin_dir}bwa/bwa mem {bwaArgs} -o {output_dir}{outputFileNames.sam}.sam -t {threads} {reference} {output_dir}{output}.extendedFrags.fastq.gz"
  seq_count = fmt"{bin_dir}sam_count -i=$1{outputFileNames.sam}.sam -o=$1{outputFileNames.sam}.csv --reference=$2 --countTab={output_dir}{outputFileNames.readsCount}.txt" % [output_dir, reference]
  readsPlot = &"{bin_dir}reads_count_plot -i={output_dir}{outputFileNames.readsCount}.txt -o={output_dir}{outputFileNames.readsCount}.pdf"
  csv_to_fasta = &"{bin_dir}csv_to_fasta -i:{output_dir}{outputFileNames.sam}.csv -n:{top_number} -o:$1{outputFileNames.topSeq}.fasta -r:$2 --only:$3" #% [output_dir, reference, ref_name]

  # multiple sequence alignment
  # TODO: write mafft stdout to file, not use ">" pipe
  msa = &"{bin_dir}mafft/mafft.bat --anysymbol --quiet --thread {threads} --maxiterate 1000 --globalpair $1{outputFileNames.topSeq}.fasta > $#{outputFileNames.topAlignment}.fasta" #% [output_dir]
  fastaPlot = &"{bin_dir}/fasta_to_plot -i=$#{outputFileNames.topAlignment}.fasta -o=$#{outputFileNames.topAlignment}.pdf" #% [output_dir, output_dir]

  # use samtools and varscan2 for calculating indel ratio along sequence
  samSort = &"{bin_dir}samtools/samtools sort --threads {threads} -l 4 -O BAM -o {output_dir}{outputFileNames.bam}.bam {output_dir}{outputFileNames.sam}.sam"

  # bam index
  bamIndex = &"{bin_dir}samtools/samtools index -@ {threads} {output_dir}{outputFileNames.bam}.bam"

  mpileup = &"{bin_dir}samtools/samtools mpileup -f {reference} -o {output_dir}{output}.mpileup {output_dir}{outputFileNames.bam}.bam"

  callIndel = &"{java} -jar {bin_dir}varscan2/VarScan2.jar pileup2cns {output_dir}{output}.mpileup --min-avg-qual {phred} --p-value 0.05 {varscanArgs} > {output_dir}{output}.var"
  varToCsv = &"{bin_dir}varToCsv -i={output_dir}{output}.var -o={output_dir}{outputFileNames.varFreq}.csv "
  # plot indel ratio along sequence
  freqPlot = &"{bin_dir}var_plot -i=$#{outputFileNames.varFreq}.csv -o=$#{outputFileNames.varFreq}.pdf --only=$#" #% [output_dir, output_dir, ref_name]

  # call SNP
  callSNP = &"{java} -jar {bin_dir}varscan2/VarScan2.jar pileup2snp {output_dir}{output}.mpileup --min-avg-qual {phred} --min-coverage 1 --min-reads2 1 --p-value 0.05 {varscanArgs} > {output_dir}{outputFileNames.snp}.snp"
  # SNP plot
  snp_plot = &"{bin_dir}snp_plot -i={output_dir}{outputFileNames.snp}.snp -o={output_dir}{outputFileNames.snp}.pdf"

  # analysis log
  analysis_log = open(&("{output_dir}{output}.log"), mode = fmAppend)

# produce a line of "#"
const star_n = repeat("#", 20) & "\n"

# convert "/" in Unix path to "\" in Windows path
proc slashCov(s: string): string =
  result = s.replace("/", "\\")

# call external commands and save log to file
proc extCall(cmd:string, log: File) {.discardable.} =
  var cmd = cmd
  when defined(windows):
    cmd = cmd.slashCov()
  log.write(star_n)
  log.write(cmd & "\n")
  log.write(star_n)
  let (cmdLog, cmdExitcode) = execCmdEx(cmd)
  if cmdExitcode != 0:
    log.write("Command failed:\n" & cmdLog & "\n")
    quit(errorcode = 3)
  else:
    log.write(cmdLog & "\n")

# call external commands through shell
proc extShellCall(cmd: string, log: File) {.discardable.} =
  var cmd = cmd
  when defined(windows):
    cmd = cmd.slashCov()
  log.write(star_n)
  log.write(cmd & "\n")
  log.write(star_n)
  let cmdExitcode = execShellCmd(cmd)
  if cmdExitcode != 0:
    log.write("Command failed.\n")
    quit(errorcode = 3)

# write stdout of command to file, not use ">", which is not compatible on Windows
proc writeExtCmdStdout(cmd: string, log: File) {.discardable.} =
  log.write(star_n)
  log.write(cmd & "\n")
  log.write(star_n)
  let
    cmdSplit = cmd.split('>')
    cmdPre = cmdSplit[0].strip()
    cmdSuf = cmdSplit[1].strip()
  let (cmdOut, cmdExitcode) = execCmdEx(cmdPre)
  if cmdExitcode == 0:
    writeFile(cmdSuf, cmdOut)
  else:
    log.write("Command failed:\n" & cmdOut & "\n")
    quit(errorcode = 3)

# processing single end reads
if pairEnd:
  extCall(fastqc, analysis_log)
  extCall(flash, analysis_log)
else:
  fastqc = fmt"{bin_dir}fastp/fastp --html {output_dir}{outputFileNames.qcReport}.html --json {output_dir}QC-report.json --qualified_quality_phred {phred} --in1 {fastq1} --out1 {output_dir}{output}.extendedFrags.fastq.gz --thread {threads} {fastpArgs} "
  extCall(fastqc, analysis_log)

extCall(reference_index % [reference], analysis_log)
extCall(mapping, analysis_log)
extCall(seq_count, analysis_log)
extCall(readsPlot, analysis_log)
extCall(samSort, analysis_log)
extCall(bamIndex, analysis_log)
extCall(mpileup, analysis_log)
when defined(windows):
  extShellCall(callIndel, analysis_log)
  if plotSNP:
    extShellCall(callSNP, analysis_log)
    extCall(snp_plot, analysis_log)
else:
  extCall(callIndel, analysis_log)
  if plotSNP:
    extCall(callSNP, analysis_log)
    extCall(snp_plot, analysis_log)
extCall(varToCsv, analysis_log)

# add results to seperate folders according to ref-name
if referenceCount.len <= 1:
  let outPrefix = output_dir
  extCall(csv_to_fasta % [outPrefix, reference, ""], analysis_log)
  writeExtCmdStdout(msa % [outPrefix], analysis_log)
  extCall(fastaPlot % [outPrefix, outPrefix], analysis_log)
  extCall(freqPlot % [outPrefix, outPrefix, ""], analysis_log)
else:
  for k in referenceCount.keys:
    let
      refName = k
      refSeq = referenceCount[k][0]
      refDir = output_dir & refName & "/"
      refPrefix = refDir
    createDir(refDir)
    writeFile(refPrefix & refName & ".fa", &">{refName}\n{refSeq}\n")

    extCall(csv_to_fasta % [refPrefix, refPrefix & refName & ".fa", refName], analysis_log)
    writeExtCmdStdout(msa % [refPrefix], analysis_log)
    extCall(fastaPlot % [refPrefix, refPrefix], analysis_log)
    extCall(freqPlot % [output_dir, refPrefix, refName], analysis_log)

# close log file
analysis_log.close()
