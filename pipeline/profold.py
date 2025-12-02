import os
import subprocess

def run_profold(root_dir: str, fasta_file: str, n_worker: int, n_struct: int,
                n_iter: int, output_dir: str = None):
    """
    Run ProFOLD and stream output to both GUI and optionally a log file.
    """
    profold_script = os.path.join(root_dir, "run_ProFOLD.sh")
    if not output_dir: output_dir = os.path.join(root_dir, "predictions")
    os.makedirs(output_dir, exist_ok=True)
    print(f"Running ProFOLD: {profold_script} {fasta_file} {output_dir}")

    process = subprocess.Popen(
        [profold_script, fasta_file, output_dir, n_worker, n_struct, n_iter],
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        universal_newlines=True,
        bufsize=1
    )

    # Stream output line by line
    for line in process.stdout:
        line = line.rstrip()
        print(line)

    process.stdout.close()
    return_code = process.wait()

    if return_code != 0:
        print(f"ProFOLD failed. Return code: {return_code}")
        raise RuntimeError(f"ProFOLD failed.")

    print(f"ProFOLD completed successfully.")

    return output_dir

