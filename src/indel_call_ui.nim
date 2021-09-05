import ui, os, strutils, osproc, strformat
import cpuinfo
from ui/rawui import Form, newForm, formSetPadded, formAppend

# add icon to windows exe
when defined(windows):
  {.passC: "-include src/windows/icon.h".}
  {.link: "bin/windows/ico.o".}

# Form control to align label and entries
type
  Form = ref object of Widget
    impl: ptr rawui.Form

proc newForm(): Form =
  result = new(Form)
  result.impl = rawui.newForm()

proc padded(f: Form; padded: int) =
  formSetPadded(f.impl, cint(padded))

proc append[SomeWidget: Widget](f: Form; label: string; c: SomeWidget; stretchy: int) =
  formAppend(f.impl, label = cstring(label), c = c.impl, stretchy = cint(stretchy))

# get tool directory
let bin_dir = getAppDir() & "/bin/"

const
  windowTitle = "CRISPR indel analysis"
  winWidth = 500
  winHeight = 500
  MAX_THREAD = 99
  MAX_ALIGNMENT = 200 # for high accuracy mafft
  PhredScore = [0, 128] # Phred score

init()

var mainwin: Window
# add menu
var menu = newMenu("File")
menu.addQuitItem(proc(): bool {.closure.} =
    mainwin.destroy()
    return true)

# title: string; width, height: int; hasMenubar: bool
mainwin = newWindow(windowTitle, winWidth, winHeight, false)
# Margined returns whether the Window has margins around its child. 
mainwin.margined = true
# OnClosing registers f to be run when the user clicks the Window's close button.
mainwin.onClosing = (proc ():bool = return true)

# padded = false
# Padded returns whether there is space between each control of the Box. 
let box = newVerticalBox(true)

# use tabs for different setting
var tab = newTab()
tab.add("Basic inputs", box)
# SetChild sets the Window's child to child. If child is nil, the Window will not have a child.
mainwin.setChild(tab)

var fastq1 = newHorizontalBox(true)
# widget, strechy = bool
# box.add(fastq1, false)
# fastq1.add(newLabel("FASTQ1 file: "), false)
var fastq1_name = newEntry("")
fastq1.add(fastq1_name, true)
fastq1.add newButton("Browse", proc() = 
  var file_open = ui.openFile(mainwin)
  fastq1_name.text = file_open
  )

# slect fastq file 2
var fastq2 = newHorizontalBox(true)
# box.add(fastq2, false)
# fastq2.add(newLabel("FASTQ2 file: "), false)
var fastq2_name = newEntry("")
fastq2.add(fastq2_name, true)
fastq2.add newButton("Browse", proc() =
  var file_open = ui.openFile(mainwin)
  fastq2_name.text=file_open)

var ref_box = newHorizontalBox(true)
# box.add(ref_box)
# ref_box.add(newLabel("Reference: "), false)
# var ref_seq = newNonWrappingMultilineEntry()
var ref_entry = newEntry("")
ref_box.add(ref_entry, true)
ref_box.add newButton("Browse", proc() =
  var file_open = ui.openFile(mainwin)
  ref_entry.text=file_open)

# slect output directory
var output = newHorizontalBox(true)
# box.add(output, false)
# output.add(newLabel("Output directory: "), false)
var output_dir = newEntry("")
output.add(output_dir, true)
output.add newButton("Browse", proc() =
  var file_open = ui.saveFile(mainwin)
  output_dir.text=file_open)

# specific region in WGS
var
  region = newHorizontalBox(true)
  regionFile = newEntry("")
region.add(regionFile, true)
region.add newButton("Browse", proc() =
  var file_open = ui.openFile(mainwin)
  regionFile.text = file_open)

# a form to contain all inputs
var inputForm = newForm()
box.add(inputForm)
inputForm.padded(1)
inputForm.append("FASTQ1:", fastq1, 0)
inputForm.append("FASTQ2:", fastq2, 0)
inputForm.append("Reference:", ref_box, 0)
inputForm.append("Output folder:", output, 0)
inputForm.append("Region file:", region, 0)
inputForm.append("", newLabel(""), 0)

# analysis type and if plot substitution
box.add(newHorizontalSeparator())

var analysis_type_box = newVerticalBox(true)
box.add(analysis_type_box)

var type_label = newHorizontalBox(true)
type_label.add(newLabel(""), true)
type_label.add(newLabel("Analysis type"))
type_label.add(newLabel(""), true)
analysis_type_box.add(type_label)

var
  radio_check_box_horizon = newHorizontalBox(true)
  analysis_type_radio_box = newVerticalBox(true)
  analysis_type_radio_buttons = newRadioButtons()
  substitution_check_box = newVerticalBox(true)
  substitution_check = newCheckbox("Substitution")

analysis_type_box.add(radio_check_box_horizon)
radio_check_box_horizon.add(newLabel(""), true)
radio_check_box_horizon.add(analysis_type_radio_box)
radio_check_box_horizon.add(newLabel(""), true)
radio_check_box_horizon.add(substitution_check_box)
radio_check_box_horizon.add(newLabel(""), true)

analysis_type_radio_buttons.add("Amplicon(s)")
analysis_type_radio_buttons.add("Whole Genome Sequencing")
analysis_type_radio_box.add(analysis_type_radio_buttons)

substitution_check_box.add(substitution_check)
substitution_check_box.add(newLabel(""))

# add separator
box.add(newHorizontalSeparator())

# paramters
# var param = newGroup("")
# box.add(param)
var param_list = newVerticalBox(true)
box.add(param_list)

var opt_label = newHorizontalBox(true)
opt_label.add(newLabel(""), true)
opt_label.add(newLabel("Options"))
opt_label.add(newLabel(""), true)
param_list.add(opt_label)

var qc_trim = newHorizontalBox(true)
param_list.add(qc_trim)

# var qc = newHorizontalBox(true)
# qc_trim.add(qc)

qc_trim.add(newLabel("Phred quality score"))
qc_trim.add(newLabel(""), true)
var qc_score = newSpinbox(min = PhredScore[0], max = PhredScore[1])
qc_score.value = 30
qc_trim.add(qc_score, false)

# if disabe default trim provided by fastp
var trim = newHorizontalBox(true)
param_list.add(trim)
trim.add(newLabel("Disable adapter trimming (enabled by default)"))
trim.add(newLabel(""), true)
var trim_disable = newCheckbox("Disable")
trim.add(trim_disable, true)

# threads
# detect total CPU cores
var totalCPU = cpuinfo.countProcessors()
var thread_box = newHorizontalBox(true)
param_list.add(thread_box)
thread_box.add(newLabel(&"CPU threads ({totalCPU} cores found)"))
thread_box.add(newLabel(""), true)
var thread = newSpinbox(min = 1, max = MAX_THREAD)
# default cpu thread used is max cores minus 1?
# thread.value = cpuinfo.countProcessors()-1
# change default cpu thread to 1?
thread.value = 1
thread_box.add(thread, false)

# number of reads for global alignment
var number_box = newHorizontalBox(true)
param_list.add(number_box)
number_box.add(newLabel("Number of reads for alignment plot"))
number_box.add(newLabel(""), true)
var top_reads = newSpinbox(min = 1, max = MAX_ALIGNMENT)
top_reads.value=10
number_box.add(top_reads, true)

# additional parameters passed to external tools
var addBox = newVerticalBox(true)
tab.add("Additional Options", addBox)
# add separator
# box.add(newHorizontalSeparator())

# fastp
var fastpBox = newVerticalBox(true)
addBox.add(fastpBox)
var fastp_label = newHorizontalBox(true)
fastp_label.add(newLabel(""), true)
fastp_label.add(newLabel("Additional QC options"))
fastp_label.add(newLabel(""), true)
fastpBox.add(fastp_label)

var fastpParam = newHorizontalBox(true)
fastpBox.add(fastpParam)
fastpParam.add(newLabel("fastp options:"))
var fastpArgs = newEntry("")
fastpParam.add(fastpArgs, true)

var fastp_example = newHorizontalBox(true)
fastp_example.add(newLabel(""), true)
fastp_example.add(newLabel("Example: --thread 2"))
fastp_example.add(newLabel(""), true)
fastp_example.add(newLabel(""), true)
fastpBox.add(fastp_example)

# add separator
addBox.add(newHorizontalSeparator())
# flash
var flashBox = newVerticalBox(true)
addBox.add(flashBox)
var flash_label = newHorizontalBox(true)
flash_label.add(newLabel(""), true)
flash_label.add(newLabel("Additional merge options"))
flash_label.add(newLabel(""), true)
flashBox.add(flash_label)

var flashParam = newHorizontalBox(true)
flashBox.add(flashParam)
flashParam.add(newLabel("FLASH options:"))
var flashArgs = newEntry("")
flashParam.add(flashArgs, true)

var flash_example = newHorizontalBox(true)
flash_example.add(newLabel(""), true)
flash_example.add(newLabel("Example: --max-overlap=150"))
flash_example.add(newLabel(""), true)
flash_example.add(newLabel(""), true)
flashBox.add(flash_example)

# add separator
addBox.add(newHorizontalSeparator())
# bwa
var bwaBox = newVerticalBox(true)
addBox.add(bwaBox)

var bwa_label = newHorizontalBox(true)
bwa_label.add(newLabel(""), true)
bwa_label.add(newLabel("Additional mapping options"))
bwa_label.add(newLabel(""), true)
bwaBox.add(bwa_label)

var bwaParam = newHorizontalBox(true)
bwaBox.add(bwaParam)
bwaParam.add(newLabel("BWA-MEM options:"))
var bwaArgs = newEntry("")
bwaParam.add(bwaArgs, true)

var bwa_example = newHorizontalBox(true)
bwa_example.add(newLabel(""), true)
bwa_example.add(newLabel("Example: -t 2"))
bwa_example.add(newLabel(""), true)
bwa_example.add(newLabel(""), true)
bwaBox.add(bwa_example)

# add separator
addBox.add(newHorizontalSeparator())
# varscan
var varscanBox = newVerticalBox(true)
addBox.add(varscanBox)

var varscan_label = newHorizontalBox(true)
varscan_label.add(newLabel(""), true)
varscan_label.add(newLabel("Additional indel calling options"))
varscan_label.add(newLabel(""), true)
varscanBox.add(varscan_label)

var varscanParam = newHorizontalBox(true)
varscanBox.add(varscanParam)
varscanParam.add(newLabel("VarScan2 options:"))
var varscanArgs = newEntry("")
varscanParam.add(varscanArgs, true)

var varscan_example = newHorizontalBox(true)
varscan_example.add(newLabel(""), true)
varscan_example.add(newLabel("Example: --min-avg-qual 30"))
varscan_example.add(newLabel(""), true)
# varscan_example.add(newLabel(""), true)
varscanBox.add(varscan_example)

# confirmation window
proc conv_done() = 
  var complete_info = &"Analysis result saved to {output_dir.text}"
  msgBox(mainwin, "Analysis complete", complete_info)

# failure message
proc failMsg(msg: string) =
  msgBox(mainwin, "Error occurs!\n", msg)

# call external indel_analysis tool
proc call_indel_analysis() =
  var indel_analysis: string
  # analysis type
  var analysis_type = "amplicons"
  if regionFile.text != "" or analysis_type_radio_buttons.radioButtonsSelected == 1:
    analysis_type = "WGS"

  if analysis_type == "amplicons":
    indel_analysis = fmt"{bin_dir}indel_analysis "
    if substitution_check.checked:
      indel_analysis &= " --plotSNP "
  else:
    indel_analysis = fmt"{bin_dir}wgsSubRegion"
    if regionFile.text != "":
      indel_analysis &= &" --region={regionFile.text} "
  indel_analysis &= fmt" -1={fastq1_name.text} -r={ref_entry.text} -o={output_dir.text} -q={qc_score.value} -t={thread.value} -n={top_reads.value} "
  if fastq2_name.text != "":
    indel_analysis &= &" -2={fastq2_name.text} "
  
  # if disable adapter trimming is checked, add the option to fastpArgs
  var fastpOpt = ""
  if trim_disable.checked:
    fastpOpt &= " --disable_adapter_trimming "
  if fastpArgs.text != "":
    fastpOpt &= fastpArgs.text
  if fastpOpt != "":
    indel_analysis &= &" --fastp=\"{fastpOpt}\" "
  
  if flashArgs.text != "":
    indel_analysis &= &" --flash=\"{flashArgs.text}\" "
  if bwaArgs.text != "":
    indel_analysis &= &" --bwa=\"{bwaArgs.text}\" "
  if varscanArgs.text != "":
    indel_analysis &= &" --varscan=\"{varscanArgs.text}\" "
  # msgBox(mainwin, "Convert", indel_analysis)
  # echo(indel_analysis)
  # use poDaemon in execCmdEx options to hide console on Windows
  let (call_indel_analysis_log, call_indel_analysis_exitcode) = execCmdEx(indel_analysis, options = {poDaemon})
  if call_indel_analysis_exitcode == 0:
    conv_done()
  else:
    failMsg(call_indel_analysis_log)

# convert
var conv_button = newButton("Begin analysis", proc() = call_indel_analysis())

# box contain analysis button to align the button center
var cbBox = newHorizontalBox(true)
cbBox.add(newLabel(""), true)
cbBox.add(conv_button, false)
cbBox.add(newLabel(""), true)

box.add(cbBox, false)

show(mainwin)
mainLoop()

