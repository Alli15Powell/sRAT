import sys, subprocess, os
from PyQt5 import QtWidgets, QtCore

class MainWin(QtWidgets.QWidget):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("sRAT – small RNA Alignment Tool")
        self.resize(560, 220)

        self.ref_path = None
        self.reads_path = None
        self.outdir = os.path.abspath("out")

        # widgets
        self.ref_btn = QtWidgets.QPushButton("Select Genome FASTA…")
        self.ref_label = QtWidgets.QLabel("No file selected")
        self.reads_btn = QtWidgets.QPushButton("Select Reads (FASTQ/FASTA)…")
        self.reads_label = QtWidgets.QLabel("No file selected")
        self.out_btn = QtWidgets.QPushButton("Select Output Folder…")
        self.out_label = QtWidgets.QLabel(self.outdir)
        self.run_btn = QtWidgets.QPushButton("Run sRAT")

        # layout
        grid = QtWidgets.QGridLayout(self)
        grid.addWidget(self.ref_btn,     0, 0)
        grid.addWidget(self.ref_label,   0, 1)
        grid.addWidget(self.reads_btn,   1, 0)
        grid.addWidget(self.reads_label, 1, 1)
        grid.addWidget(self.out_btn,     2, 0)
        grid.addWidget(self.out_label,   2, 1)
        grid.addWidget(self.run_btn,     3, 0, 1, 2)

        # signals
        self.ref_btn.clicked.connect(self.choose_ref)
        self.reads_btn.clicked.connect(self.choose_reads)
        self.out_btn.clicked.connect(self.choose_out)
        self.run_btn.clicked.connect(self.run_job)

    def choose_ref(self):
        path, _ = QtWidgets.QFileDialog.getOpenFileName(self, "Choose Genome FASTA", "", "FASTA (*.fa *.fasta);;All Files (*)")
        if path:
            self.ref_path = path
            self.ref_label.setText(path)

    def choose_reads(self):
        path, _ = QtWidgets.QFileDialog.getOpenFileName(self, "Choose Reads (FASTQ/FASTA)", "", "FASTQ/FASTA (*.fastq *.fq *.fa *.fasta);;All Files (*)")
        if path:
            self.reads_path = path
            self.reads_label.setText(path)

    def choose_out(self):
        path = QtWidgets.QFileDialog.getExistingDirectory(self, "Choose Output Folder", self.outdir)
        if path:
            self.outdir = path
            self.out_label.setText(path)

    def run_job(self):
        if not self.ref_path or not self.reads_path:
            QtWidgets.QMessageBox.warning(self, "Missing files", "Please select both Genome FASTA and Reads.")
            return
        cmd = [sys.executable, "srat_cli.py", "--genome", self.ref_path, "--reads", self.reads_path, "--out", self.outdir]
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
        out, _ = proc.communicate()
        QtWidgets.QMessageBox.information(self, "sRAT finished", out or "Done.")
        if QtWidgets.QMessageBox.question(self, "Open output folder?", "Open the output folder?") == QtWidgets.QMessageBox.Yes:
            QtWidgets.QDesktopServices.openUrl(QtCore.QUrl.fromLocalFile(self.outdir))

if __name__ == "__main__":
    app = QtWidgets.QApplication(sys.argv)
    win = MainWin()
    win.show()
    sys.exit(app.exec_())
