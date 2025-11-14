import os, sys, threading
from multiprocessing import Process, Queue
from PyQt6.QtWidgets import (
    QApplication,
    QMainWindow,
    QWidget,
    QVBoxLayout,
    QPushButton,
    QFileDialog,
    QLabel,
    QTextEdit,
    QLineEdit,
    QHBoxLayout,
    QMessageBox,
)
from PyQt6.QtCore import QProcess, Qt
from PyQt6.QtWebEngineWidgets import QWebEngineView
import qdarkstyle

import py3Dmol
import pipeline

root_dir = os.path.dirname(os.path.abspath(__file__))
work_dir = os.path.join(root_dir, "work_dir")
log_dir = os.path.join(work_dir, "log")

os.makedirs(work_dir, exist_ok=True)
os.makedirs(log_dir, exist_ok=True)

class ProteinViewer(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Protein 3D Structure Prediction Pipeline")
        self.setGeometry(100, 100, 800, 600)

        central_widget = QWidget()
        self.setCentralWidget(central_widget)
        layout = QVBoxLayout(central_widget)

        self.file_label = QLabel("No sequence file selected")
        layout.addWidget(self.file_label)

        file_layout = QHBoxLayout()
        self.file_button = QPushButton("Select Sequence File")
        self.file_button.clicked.connect(self.pick_file)
        file_layout.addWidget(self.file_button)

        layout.addLayout(file_layout)

        layout.addWidget(QLabel("Or enter sequence directly:"))
        self.seq_input = QTextEdit()
        self.seq_input.setPlaceholderText(">seq1\nMENSDSENQKAVKLLIAGHGEVQGVSVKSAMKGVKGGVKAQML")
        layout.addWidget(self.seq_input)

        self.db_label = QLabel("No database selected")
        layout.addWidget(self.db_label)

        db_layout = QHBoxLayout()
        self.db_button = QPushButton("Select Database Prefix")
        self.db_button.clicked.connect(self.pick_database)
        db_layout.addWidget(self.db_button)
        layout.addLayout(db_layout)

        self.run_button = QPushButton("Run MSA Pipeline")
        self.run_button.clicked.connect(self.run_pipeline_process)
        layout.addWidget(self.run_button)

        self.stop_button = QPushButton("Stop Pipeline")
        self.stop_button.clicked.connect(self.stop_pipeline)
        self.stop_button.setEnabled(False)
        layout.addWidget(self.stop_button)

        self.status_output = QTextEdit()
        self.status_output.setReadOnly(True)
        layout.addWidget(self.status_output)

        self.seq_file = None
        self.db_prefix = None
        self.pipeline_proc = None
        try:
            self.log_file = open(os.path.join(log_dir, "pipeline.log"), 'w')
        except OSError:
            self.log_file = None
            self.log("Could not write to file:", fname)

        # Button to load PDB
        # self.load_button = QPushButton("Load PDB")
        # self.load_button.clicked.connect(self.load_pdb)
        # layout.addWidget(self.load_button)

        # WebEngineView for 3Dmol
        # self.web_view = QWebEngineView()
        # layout.addWidget(self.web_view)

    def pick_file(self):
        file_path, _ = QFileDialog.getOpenFileName(
            self, "Open FASTA File", "", "(*.fasta *.fa *.fna *.aln);;All Files (*)"
        )
        if file_path:
            self.seq_file = file_path
            self.file_label.setText(f"Sequence file: {file_path}")
        else:
            self.seq_file = None
            self.file_label.setText("No sequence file selected")

    def pick_database(self):
        # Ask user to pick any file in the HHblits database
        db_prefix, _ = QFileDialog.getOpenFileName(
            self,
            "Select any HHblits DB file (e.g. pdb100_a3m.ffdata)",
            "",
            "HHsuite files (*.ffdata *.ffindex)"
        )

        if not db_prefix:
            return

        # Determine prefix by stripping known HHblits suffixes
        suffixes = [
            "_a3m.ffdata", "_a3m.ffindex",
            "_cs219.ffdata", "_cs219.ffindex",
            "_hhm.ffdata", "_hhm.ffindex"
        ]
        prefix = None
        for suf in suffixes:
            if db_prefix.endswith(suf):
                prefix = db_prefix[: -len(suf)]
                break

        if not prefix:
            QMessageBox.warning(self, "Error", "Not a valid HHblits database file!")
            return

        # Check required files exist
        required_files = [prefix + s for s in suffixes]
        missing = [str(f) for f in required_files if not os.path.exists(f)]

        if missing:
            QMessageBox.warning(
                self,
                "Missing files",
                "The following required HHblits DB files are missing:\n" + "\n".join(missing)
            )
            return

        self.db_prefix = prefix
        self.db_label.setText(f"HHblits database prefix: {prefix}")

    def load_pdb(self):
        pdb_prefix, _ = QFileDialog.getOpenFileName(self, "Open PDB File", "", "PDB Files (*.pdb)")
        if pdb_prefix:
            with open(pdb_prefix, 'r') as f:
                pdb_str = f.read()

            # Create py3Dmol HTML
            view = py3Dmol.view(js='https://3dmol.org/build/3Dmol.js')
            view.addModel(pdb_str, 'pdb')

            view.setStyle({'cartoon': {'color':'spectrum'}})
            view.zoomTo()

            html = view._make_html()
            self.web_view.setHtml(html)

    def log(self, msg):
        if self.log_file:
            self.log_file.write(msg)
            self.log_file.write('\n')
            self.log_file.flush()
        self.status_output.append(msg)
        QApplication.processEvents()  # update GUI
    
    def run_pipeline_process(self):
        if self.pipeline_proc and self.pipeline_proc.poll() is None:
            return
        if self.seq_file:
            query_file = self.seq_file
        elif self.seq_input.toPlainText().strip():
            query_file = os.path.join(work_dir, "seq.fasta")
            with open(query_file, "w") as f:
                f.write(self.seq_input.toPlainText())
        else:
            self.log("No sequence input provided.")
            return

        if not self.db_prefix:
            self.log("No database selected!")
            return

        self.pipeline_proc = QProcess(self)
        self.pipeline_proc.setProcessChannelMode(QProcess.ProcessChannelMode.MergedChannels)

        self.pipeline_proc.readyReadStandardOutput.connect(self.handle_stdout)
        self.pipeline_proc.finished.connect(self.on_pipeline_finished)

        # Start the process (replace with your Python command or shell script)
        cmd = ["python3", "-u", "-c",
               f"import pipeline; pipeline.run_pipeline('{root_dir}', '{work_dir}', '{log_dir}', '{query_file}', '{self.db_prefix}', top_hits=200)"]

        self.log(f"Starting pipeline: {' '.join(cmd)}")
        self.pipeline_proc.start(cmd[0], cmd[1:])

        self.run_button.setEnabled(False)
        self.stop_button.setEnabled(True)

    def handle_stdout(self):
        data = self.pipeline_proc.readAllStandardOutput().data().decode()
        for line in data.splitlines():
            self.log(line)

    def stop_pipeline(self):
        if self.pipeline_proc and self.pipeline_proc.state() != QProcess.ProcessState.NotRunning:
            self.pipeline_proc.kill()  # hard kill
            self.pipeline_proc = None
            self.log("Pipeline process killed by user.")
            self.run_button.setEnabled(True)
            self.stop_button.setEnabled(False)
    
    def on_pipeline_finished(self, exitCode, exitStatus):
        self.log(f"Pipeline finished with code {exitCode}")
        self.run_button.setEnabled(True)
        self.stop_button.setEnabled(False)

    def closeEvent(self, event):
        self.stop_pipeline()
        event.accept()

if __name__ == "__main__":
    app = QApplication(sys.argv)
    app.setApplicationName("ProteinViewer")
    app.setApplicationDisplayName("Protein 3D Viewer")
    app.setApplicationVersion("1.0")
    app.setStyleSheet(qdarkstyle.load_stylesheet_pyqt6())
    window = ProteinViewer()
    window.show()
    sys.exit(app.exec())

