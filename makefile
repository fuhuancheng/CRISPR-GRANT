APPNAME = CRISPR-GRANT
version = v1.1.0

buildDir = bin/
binDir = ${buildDir}bin/
WINDOWSDIR = ${buildDir}windows
macBUILDDIR = ${buildDir}macos/$(APPNAME)-$(version).app

debug: indel_analysis indel_call_ui sam_count csv_to_fasta fasta_to_plot var_plot varToCsv reads_count_plot wgsSubRegion snp_plot

indel_call_ui: src/indel_call_ui.nim
	nim c --app:gui --outdir:${buildDir} src/indel_call_ui.nim

indel_analysis: src/indel_analysis.nim
	nim c --outdir:${binDir} src/indel_analysis.nim

sam_count: src/sam_count.nim
	nim c --outdir:${binDir} --opt:speed src/sam_count.nim

csv_to_fasta: src/csv_to_fasta.nim
	nim c --outdir:${binDir} src/csv_to_fasta.nim

fasta_to_plot: src/fasta_to_plot.nim
	nim c --outdir:${binDir} src/fasta_to_plot.nim

var_plot: src/var_plot.nim
	nim c --outdir:${binDir} src/var_plot.nim

varToCsv: src/varToCsv.nim
	nim c --outdir:${binDir} src/varToCsv.nim

reads_count_plot: src/reads_count_plot.nim
	nim c --outdir:${binDir} src/reads_count_plot.nim

wgsSubRegion: src/wgsSubRegion.nim
	nim c --outdir:${binDir} src/wgsSubRegion.nim

snp_plot: src/snp_plot.nim
	nim c --outdir:${binDir} src/snp_plot.nim

unix-gui: src/indel_call_ui.nim
	nim c -d:release --opt:speed --app:gui --outdir:${buildDir} src/indel_call_ui.nim

release: src/indel_analysis.nim src/sam_count.nim src/csv_to_fasta.nim src/fasta_to_plot.nim src/varToCsv.nim src/reads_count_plot.nim src/wgsSubRegion.nim src/snp_plot.nim
	nim c -d:release --opt:speed --outdir:${binDir} src/indel_analysis.nim
	nim c -d:release --opt:speed --outdir:${binDir} src/sam_count.nim
	nim c -d:release --opt:speed --outdir:${binDir} src/csv_to_fasta.nim
	nim c -d:release --opt:speed --outdir:${binDir} src/fasta_to_plot.nim
	nim c -d:release --opt:speed --outdir:${binDir} src/var_plot.nim
	nim c -d:release --opt:speed --outdir:${binDir} src/varToCsv.nim
	nim c -d:release --opt:speed --outdir:${binDir} src/reads_count_plot.nim
	nim c -d:release --opt:speed --outdir:${binDir} src/wgsSubRegion.nim
	nim c -d:release --opt:speed --outdir:${binDir} src/snp_plot.nim

windows-ico: src/windows/app.ico src/windows/app.rc
	windres src/windows/app.rc $(WINDOWSDIR)/ico.o

windows-gui: windows-ico src/indel_call_ui.nim
	nim cpp -d:release -d:mingw --app:gui --passL:-static --outdir:${buildDir} src/indel_call_ui.nim

windows: windows-gui release

unix: unix-gui release

mac: src/macos/Info.plist src/macos/App.icns
	# rm -rf bin/build
	install_name_tool -add_rpath @executable_path/ ${binDir}fasta_to_plot
	install_name_tool -add_rpath @executable_path/ ${binDir}var_plot
	install_name_tool -add_rpath @executable_path/ ${binDir}reads_count_plot
	install_name_tool -add_rpath @executable_path/ ${binDir}snp_plot
	mkdir -p $(macBUILDDIR)/Contents/Resources
	mkdir -p $(macBUILDDIR)/Contents/MacOS
	cp src/macos/Info.plist $(macBUILDDIR)/Contents/
	cp src/macos/App.icns $(macBUILDDIR)/Contents/Resources/
	cp ${buildDir}indel_call_ui $(macBUILDDIR)/Contents/MacOS/
	cp -r ${binDir} $(macBUILDDIR)/Contents/MacOS/bin
	cp -r ${buildDir}example-data ${buildDir}macos/
	cp doc/user-guide.pdf ${buildDir}macos/
	cd ${buildDir}macos && zip -r -9 $(APPNAME)-$(version)_x64_macOS.zip $(APPNAME)-$(version).app example-data/* user-guide.pdf

linux:
	cp doc/user-guide.pdf ${buildDir}
	cd ${buildDir} &&	zip -r $(APPNAME)-$(version)_x64_linux.zip indel_call_ui bin/* example-data/* user-guide.pdf

win:
	cd ${buildDir} &&	7z a $(APPNAME)-$(version)_x64_windows.zip indel_call_ui.exe bin/* example-data/* user-guide.pdf
