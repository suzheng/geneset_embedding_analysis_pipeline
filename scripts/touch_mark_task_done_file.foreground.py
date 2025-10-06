import time
import os
from pathlib import Path

def main():
    # Sleep for 900 seconds (15 minutes)
    time.sleep(900)

    # Get the directory of the current script
    script_dir = Path(__file__).resolve().parent

    # File path for MARK_TASK_DONE
    done_file = script_dir / "MARK_TASK_DONE"

    # "touch" the file (create if not exists, update timestamp if exists)
    done_file.touch()

    print(f"Created or updated: {done_file}")

if __name__ == "__main__":
    main()
