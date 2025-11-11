import os, sys, threading
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
    QHBoxLayout
)

from PyQt6.QtWebEngineWidgets import QWebEngineView
import qdarkstyle

from Bio import SeqIO
import py3Dmol

import pipeline.database as db
import pipeline.search as search
import pipeline.align as align
import pipeline.profold as profold

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
        self.db_button = QPushButton("Select BLAST Database")
        self.db_button.clicked.connect(self.pick_database)
        db_layout.addWidget(self.db_button)
        layout.addLayout(db_layout)


        self.run_button = QPushButton("Run MSA Pipeline")
        self.run_button.clicked.connect(self.run_pipeline_threaded)
        layout.addWidget(self.run_button)

        self.stop_button = QPushButton("Stop Pipeline")
        # self.stop_button.clicked.connect(self.stop_pipeline)
        layout.addWidget(self.stop_button)

        self.status_output = QTextEdit()
        self.status_output.setReadOnly(True)
        layout.addWidget(self.status_output)

        self.seq_file = None
        self.db_file = None
        self.pipeline_thread = None

        # Button to load PDB
        # self.load_button = QPushButton("Load PDB")
        # self.load_button.clicked.connect(self.load_pdb)
        # layout.addWidget(self.load_button)

        # WebEngineView for 3Dmol
        # self.web_view = QWebEngineView()
        # layout.addWidget(self.web_view)

    def pick_file(self):
        file_path, _ = QFileDialog.getOpenFileName(
            self, "Open FASTA File", "", "FASTA or GZIP Files (*.fasta *.fa *.fna *.gz);;All Files (*)"
        )
        if file_path:
            self.seq_file = file_path
            self.file_label.setText(f"Selected: {file_path}")

    def pick_database(self):
        db_file, _ = QFileDialog.getOpenFileName(
            self, "Open BLAST Database File", "", 
            "FASTA or gzipped FASTA (*.fasta *.fa *.fna *.gz);;All Files (*)"
        )
        if db_file:
            self.db_file = db_file
            self.db_label.setText(f"Selected: {db_file}")

    def load_pdb(self):
        pdb_file, _ = QFileDialog.getOpenFileName(self, "Open PDB File", "", "PDB Files (*.pdb)")
        if pdb_file:
            with open(pdb_file, 'r') as f:
                pdb_str = f.read()

            # Create py3Dmol HTML
            view = py3Dmol.view(js='https://3dmol.org/build/3Dmol.js')
            view.addModel(pdb_str, 'pdb')

            view.setStyle({'cartoon': {'color':'spectrum'}})
            view.zoomTo()

            html = view._make_html()
            self.web_view.setHtml(html)
    def log(self, msg):
        self.status_output.append(msg)
        QApplication.processEvents()  # update GUI

    def run_pipeline_threaded(self):
        if self.pipeline_thread and self.pipeline_thread.is_alive():
            self.log("Pipeline already running, please stop it before start a new one.")
            return
        self.pipeline_thread = threading.Thread(target=self.run_pipeline, daemon=True)
        self.pipeline_thread.start()

    def run_pipeline(self):
        # --- Prepare query file ---
        if self.seq_file:
            query_file = self.seq_file
        elif self.seq_input.toPlainText().strip():
            query_file = "temp_query.fasta"
            with open(query_file, "w") as f:
                f.write(self.seq_input.toPlainText())
        else:
            self.log("No sequence input provided.")
            return

        if not self.db_file:
            self.log("No database selected!")
            return

        self.log(f"Preparing BLAST database: {self.db_file} ...")
        try:
            blast_db = db.prepare_blast_db(self.db_file)
            self.log(f"Database ready: {blast_db}")
        except Exception as e:
            self.log(f"Database error: {e}")
            return

        # --- Run BLAST search ---
        self.log("Running BLAST search...")
        try:
            hits = search.blast_search(query_file, blast_db, top_hits=10)
            self.log(f"Top hits: {hits}")
        except Exception as e:
            self.log(f"BLAST error: {e}")
            return

        # --- Extract sequences ---
        self.log("Extracting hit sequences...")

        # Use the original FASTA, decompressed if needed
        if self.db_file.endswith(".gz"):
            import gzip, shutil
            fasta_for_db = os.path.splitext(self.db_file)[0]  # remove .gz
            if not os.path.exists(fasta_for_db):
                with gzip.open(self.db_file, "rt") as f_in, open(fasta_for_db, "w") as f_out:
                    shutil.copyfileobj(f_in, f_out)
        else:
            fasta_for_db = self.db_file

        # Parse sequences from the actual FASTA file
        records = [r for r in SeqIO.parse(fasta_for_db, "fasta") if r.id in hits]

        if not records:
            self.log("No hits found.")
            return

        hits_file = "hits.fasta"
        SeqIO.write(records, hits_file, "fasta")
        self.log(f"Extracted {len(records)} sequences -> {hits_file}")

        # --- Run MAFFT ---
        self.log("Running MAFFT alignment...")
        alignment_fasta = "aligned.fasta"
        try:
            alignment = align.run_mafft(hits_file, alignment_fasta)
            self.log(f"Alignment complete: {len(alignment)} sequences, length {alignment.get_alignment_length()}")
            self.log(f"MSA exported to '{alignment_fasta}'")
        except Exception as e:
            self.log(f"MAFFT error: {e}")
            return

        aln_file, aln_content = align.write_aln_file(alignment_fasta, "hits.aln")
        self.log(f"Write {len(aln_content)} sequences to {aln_file}")

        profold.run(aln_file, output_dir="predictions", log_func=self.log)

if __name__ == "__main__":
    app = QApplication(sys.argv)
    app.setApplicationName("ProteinViewer")
    app.setApplicationDisplayName("Protein 3D Viewer")
    app.setApplicationVersion("1.0")
    app.setStyleSheet(qdarkstyle.load_stylesheet_pyqt6())
    window = ProteinViewer()
    window.show()
    sys.exit(app.exec())

