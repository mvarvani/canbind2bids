import pandas as pd


def reformat_table(table: pd.DataFrame):
    """Reformat a table to conform to the BIDS standard."""
    return table.assign(
        participant_id=table.SUBJLABEL.str.split("_").str[1:3].str.join("")
    ).drop(columns=["SUBJLABEL"])


def gen_description(row: pd.Series) -> str:
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
    if row["Variable Type"] not in ["single-select", "radio"]:
        return None
    return dict([item.split(" = ") for item in row["Value Codes"].split("; ")])


def row_to_dict(row: pd.Series):
    return (
        row["Variable Name"],
        dict(
            [("LongName", row["Item Name"]), ("Description", gen_description(row))]
            + ([("Levels", levels)] if (levels := handle_levels(row)) else [])
        ),
    )


def main(doc_path, out_path):
    table = pd.read_excel(doc_path)


if __name__ == "__main__":
    main(snakemake.input.control_csv, snakemake.input.patient_csv, snakemake.output.tsv)
