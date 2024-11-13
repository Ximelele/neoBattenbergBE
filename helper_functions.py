import subprocess
import shlex

def safe_run(command: str) -> None:
    args = shlex.split(command)
    subprocess.run(args, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
