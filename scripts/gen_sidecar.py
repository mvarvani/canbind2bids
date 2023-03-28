"""Handle conversion of CANBIND documentation to BIDS"""

import json
from typing import Any, Optional

import pandas as pd


def gen_description(row: pd.Series):
    """Generate a column description.

    The generated description will come from question text and special
    instructions as applicable.
    """
    question_text = (
        f"Question text: {row['Question text']}"
        if pd.notna(row["Question text"])
        else None
    )
    special_instructions = (
        f"Special instructions: {row['Special Instructions']}"
        if pd.notna(row["Special Instructions"])
        else None
    )
    return ", ".join(
        [
            desc_part
            for desc_part in [question_text, special_instructions]
            if desc_part is not None
        ]
    )


def handle_levels(row: pd.Series):
    """Generate a Levels object for a row if applicable."""
    if row["Variable Type"] not in ["single-select", "radio"]:
        return None

    for delimiter in [";", ","]:
        try:
            return {
                key.strip(): val.strip()
                for key, val in [
                    item.strip().split("=")
                    for item in row["Value Codes"].split(delimiter)
                ]
            }
        except ValueError:
            # The delimiter was probably wrong
            pass
        except AttributeError as err:
            # We didn't get the kind of value we expected
            raise AttributeError(f"Could not parse a cell from row:\n{row}") from err
    raise ValueError(f"Could not parse row:\n{row}")


def row_to_dict(row: pd.Series):
    """Generate a BIDS TSV sidecar-formatted dict from a row."""
    return (
        row["Variable Name"]
        if row["Variable Name"] != "SUBJLABEL"
        else "participant_id",
        dict(
            [("LongName", row["Item Name"]), ("Description", gen_description(row))]
            + ([("Levels", levels)] if (levels := handle_levels(row)) else [])
        ),
    )


def main(doc_path: str, out_path: str) -> None:
    """Read a CANBIND excel documentation file and reformat it to BIDS JSON."""
    table = pd.read_excel(doc_path)
    with open(out_path, "w", encoding="utf-8") as out_file:
        json.dump(dict(list(table.apply(row_to_dict, axis=1).values)), out_file)


if __name__ == "__main__":
    main(snakemake.input.doc_excel, snakemake.output.json)
