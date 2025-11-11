import os
import subprocess

def run(aln_file: str, output_dir="predictions", log_func=None):
    """
    Run ProFOLD and stream output to both GUI and optionally a log file.
    """
    profold_script = "./run_ProFOLD.sh"
    os.makedirs(output_dir, exist_ok=True)
    if log_func:
        log_func(f"Running ProFOLD: {profold_script} {aln_file} {output_dir}")

    process = subprocess.Popen(
        [profold_script, aln_file, output_dir],
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        universal_newlines=True,
        bufsize=1
    )

    # Stream output line by line
    if log_func:
        for line in process.stdout:
            line = line.rstrip()
            log_func(line)

    process.stdout.close()
    return_code = process.wait()

    if return_code != 0:
        if log_func:
            log_func(f"ProFOLD failed. Return code: {return_code}")
        raise RuntimeError(f"ProFOLD failed.")

    if log_func:
        log_func(f"ProFOLD completed successfully.")

    return output_dir

