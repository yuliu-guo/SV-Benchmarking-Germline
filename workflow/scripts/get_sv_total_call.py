#!/usr/bin/env python3

import argparse
import json
from pathlib import Path
from statistics import mean

SAMPLES = ["NA12878", "NA12879", "NA12881", "NA12882"]


def read_sv_call_count(summary_json: Path) -> int:
    """Return TP-base + FP from a summary.json file."""
    with summary_json.open("r") as f:
        data = json.load(f)

    try:
        tp_base = data["TP-base"]
        fp = data["FP"]
    except KeyError as e:
        raise KeyError(f"Missing key {e} in {summary_json}")

    return tp_base + fp


def collect_sv_counts(root_folder: Path):
    """
    Expected structure:

    root_folder/
        Dragen/
            NA12878/summary.json
            NA12879/summary.json
            NA12881/summary.json
            NA12882/summary.json
        Manta/
            NA12878/summary.json
            ...
    """
    results = {}

    for tool_dir in sorted(root_folder.iterdir()):
        if not tool_dir.is_dir():
            continue

        tool_name = tool_dir.name
        sample_counts = {}

        for sample in SAMPLES:
            summary_json = tool_dir / sample / "summary.json"

            if summary_json.exists():
                sample_counts[sample] = read_sv_call_count(summary_json)
            else:
                sample_counts[sample] = None

        results[tool_name] = sample_counts

    return results


def write_tsv(results, output_tsv: Path):
    with output_tsv.open("w") as out:
        header = ["calling_tool"] + SAMPLES + ["Average"]
        out.write("\t".join(header) + "\n")

        for tool_name, sample_counts in results.items():
            values = [sample_counts[sample] for sample in SAMPLES]

            numeric_values = [v for v in values if v is not None]
            avg = round(mean(numeric_values)) if numeric_values else ""

            row = [tool_name]
            row.extend("" if v is None else str(v) for v in values)
            row.append(str(avg))

            out.write("\t".join(row) + "\n")


def main():
    parser = argparse.ArgumentParser(
        description="Calculate total SV calls per calling tool and sample from summary.json files."
    )
    parser.add_argument(
        "root_folder",
        help="Root folder containing one subfolder per calling tool."
    )
    parser.add_argument(
        "-o",
        "--output",
        default="sv_call_counts.tsv",
        help="Output TSV file name. Default: sv_call_counts.tsv"
    )

    args = parser.parse_args()

    root_folder = Path(args.root_folder)
    output_tsv = Path(args.output)

    if not root_folder.exists():
        raise FileNotFoundError(f"Root folder does not exist: {root_folder}")

    results = collect_sv_counts(root_folder)
    write_tsv(results, output_tsv)

    print(f"Wrote TSV to: {output_tsv}")


if __name__ == "__main__":
    main()
