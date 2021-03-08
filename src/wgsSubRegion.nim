# mapping, sam to bam, extracting reference and reads
import parseopt
import os
import osproc
import times
import strutils
import strformat

# version
const version = "v0.1.0"

# help
const help = """
Usage: wgsSubRegion -1=fastq-1.gz -2=fastq-2.gz -r=reference --region=region.txt -o=output

-h, --help       get this help and exit.
-1               fastq file 1.
-2               fastq file 2.
-r, --reference  amplicon reference file, fasta format.
-o, --output     output file name.
-t, --thread     number of threads.
-n, --topNumber  number of top reads for alignment plot.
-q               quality phred value, default is 30.
--region         text file include specific region for analysis (format: chrome[:start[-end]], one site per line). For example, chr1:42-2020.
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
  # specific region to extract
  region, regionName, regionRanges: string
  # fastq file name
  fastq1, fastq1Name:string
  fastq2, fastq2Name:string
  # reference file name
  reference:string
  # additional parameters
  fastpArgs: string
  flashArgs: string
  bwaArgs: string
  varscanArgs: string
  # single end or pair end
  pairEnd = true

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
      fastq1Name = splitFile(fastq1).name
    of "2":
      fastq2 = val
      fastq2Name = splitFile(fastq2).name
    of "q":
      phred = val
    of "n", "topNumber":
      top_number = val
    of "region":
      region = val
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

createDir(output_dir)

# if region file provided, process
if region != "":
  regionName = splitFile(region).name
  regionRanges = open(region).readAll().strip().split('\n').join(" ")

var
  # external commands: fastp, bwa, samtools
  fastqc = fmt"{bin_dir}fastp/fastp --html {output_dir}0.QC-report.html --json {output_dir}QC-report.json --qualified_quality_phred {phred} --in1 {fastq1} --in2 {fastq2} --out1 $1{fastq1Name}-qc.fq.gz --out2 $1{fastq2Name}-qc.fq.gz --thread {threads} {fastpArgs} " % [output_dir]
  reference_index = fmt"{bin_dir}bwa/bwa index {reference}" #% [reference]
  mapping = fmt"{bin_dir}bwa/bwa mem -t {threads} {bwaArgs} -o {output_dir}{output}.sam {reference} $1{fastq1Name}-qc.fq.gz $1{fastq2Name}-qc.fq.gz" % [output_dir]
  samSort = &"{bin_dir}samtools/samtools sort --threads {threads} -l 4 -O BAM -o {output_dir}{output}_sorted.bam {output_dir}{output}.sam"

  # call mutation on whole genome
  mpileup = &"{bin_dir}samtools/samtools mpileup -f {reference} -o {output_dir}{output}.mpileup {output_dir}{output}_sorted.bam"
  callIndel = &"{java} -jar {bin_dir}varscan2/VarScan2.jar pileup2indel {output_dir}{output}.mpileup --min-avg-qual {phred} --p-value 0.05 --variants {varscanArgs} > {output_dir}{output}.var"
  varToCsv = &"{bin_dir}varToCsv -i={output_dir}{output}.var -o={output_dir}{output}_var.csv "
  # plot indel ratio along sequence
  # freqPlot = &"{bin_dir}var_plot -i={output_dir}{output}_var.csv -o={output_dir}{output}_varFreq.pdf" #% [output_dir & output, ref_name]

  # extract reference in specific region
  refIndex = fmt"{bin_dir}samtools/samtools faidx {reference}"
  extractRef = fmt"{bin_dir}samtools/samtools faidx -o {output_dir}reference-{regionName}.fa $1 $2" % [reference, regionRanges]
  # extract bam mapping to specific region
  bamIndex = fmt"{bin_dir}samtools/samtools index -@ {threads} $1$2_sorted.bam " % [output_dir, output]
  extractBam = fmt"{bin_dir}samtools/samtools view -o $1$2_sorted-sub.bam -O BAM $1$2_sorted.bam $3" % [output_dir, output, regionRanges]
  # convert bam to fastq
  bamToFastq = fmt"{bin_dir}samtools/samtools fastq -o $1$2_mapped-sub.fq $1$2_sorted-sub.bam" % [output_dir, output]

  # call indel_analysis
  indelAnalysis = fmt"{bin_dir}indel_analysis -1={output_dir}{output}_mapped-sub.fq -r={output_dir}reference-{regionName}.fa -o={output_dir}{output}-{regionName} -t={threads} -n={top_number} -q={phred} --javaPath={java}"

  # analysis log
  analysis_log = open(&("{output_dir}{output}.log"), mode = fmWrite)

if fastpArgs != "":
  indelAnalysis &= fmt" --fastp={fastpArgs} "
if flashArgs != "":
  indel_analysis &= fmt" --flash={flashArgs} "
if bwaArgs != "":
  indel_analysis &= fmt" --bwa={bwaArgs} "
if varscanArgs != "":
  indel_analysis &= fmt" --varscan={varscanArgs} "

# produce a line of "#"
const star_n = repeat("#", 20) & "\n"

# call external commands and save log to file
proc extCall(cmd: string, log: File) {.discardable.} =
  let (cmdLog, cmdExitcode) = execCmdEx(cmd)
  log.write(star_n)
  log.write(cmd & "\n")
  log.write(star_n)

  if cmdExitcode != 0:
    log.write("Command failed:\n" & cmdLog & "\n")
    quit(errorcode = 3)
  log.write(cmdLog & "\n")

extCall(reference_index, analysis_log)

# processing single end reads
if pairEnd:
  extCall(fastqc, analysis_log)
  extCall(mapping, analysis_log)
else:
  fastqc = fmt"{bin_dir}fastp/fastp --html {output_dir}0.QC-report.html --json {output_dir}QC-report.json --qualified_quality_phred {phred} --in1 {fastq1} --out1 $1{fastq1Name}-qc.fq.gz --thread {threads} {fastpArgs} " % [output_dir]
  mapping = fmt"{bin_dir}bwa/bwa mem -t {threads} {bwaArgs} -o {output_dir}{output}.sam {reference} $1{fastq1Name}-qc.fq.gz" % [output_dir]
  extCall(fastqc, analysis_log)
  extCall(mapping, analysis_log)

extCall(samSort, analysis_log)
# call mutation on whole genome
extCall(mpileup, analysis_log)
extCall(callIndel, analysis_log)
extCall(varToCsv, analysis_log)
# extCall(freqPlot, analysis_log)

if region != "":
  extCall(refIndex, analysis_log)
  extCall(extractRef, analysis_log)
  extCall(bamIndex, analysis_log)
  extCall(extractBam, analysis_log)
  extCall(bamToFastq, analysis_log)

  extCall(indel_analysis, analysis_log)

# close log file
analysis_log.close()
