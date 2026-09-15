#!/usr/bin/env python3
"""Regression check for retrying a full-sized FASTQ partial with a bad MD5."""
import csv
import hashlib
import os
import subprocess
import sys
import tempfile
from pathlib import Path


SCRIPT = Path(__file__).parents[1] / "server_download_ena_fastq_20260914.py"


def main() -> None:
    good = b"GOOD"
    bad = b"FAIL"
    with tempfile.TemporaryDirectory() as temporary:
        root = Path(temporary) / "server"
        manifest = root / "metadata" / "manifest.tsv"
        manifest.parent.mkdir(parents=True)
        fields = [
            "TaskID", "GroupIDs", "GSE", "SelectedAccession", "SRR", "ENAFastqURL",
            "ExpectedMD5", "ExpectedBytes", "ResolutionState",
        ]
        with manifest.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
            writer.writeheader()
            writer.writerow({
                "TaskID": "test", "GroupIDs": "KLA31_TEST", "GSE": "GSE_TEST",
                "SelectedAccession": "GSM_TEST", "SRR": "SRR_TEST",
                "ENAFastqURL": "https://example.invalid/SRR_TEST.fastq.gz",
                "ExpectedMD5": hashlib.md5(good).hexdigest(), "ExpectedBytes": str(len(good)),
                "ResolutionState": "resolved",
            })
        partial = root / "raw" / "ena_fastq" / "GSE_TEST" / "SRR_TEST_SRR_TEST.fastq.gz.part"
        partial.parent.mkdir(parents=True)
        partial.write_bytes(bad)
        bin_dir = root / "bin"
        bin_dir.mkdir()
        fake_wget = bin_dir / "wget"
        fake_wget.write_text(
            "#!/usr/bin/env python3\n"
            "import pathlib, sys\n"
            "out = pathlib.Path(sys.argv[sys.argv.index('-O') + 1])\n"
            "if out.exists(): raise SystemExit(17)\n"
            "out.write_bytes(b'GOOD')\n",
            encoding="utf-8",
        )
        fake_wget.chmod(0o755)
        env = os.environ | {"PATH": f"{bin_dir}:{os.environ['PATH']}"}
        subprocess.run([sys.executable, str(SCRIPT), str(root), str(manifest)], check=True, env=env)
        destination = partial.with_suffix("")
        assert destination.read_bytes() == good
        assert list(partial.parent.glob("*.part.checksum-mismatch-*"))


if __name__ == "__main__":
    main()
