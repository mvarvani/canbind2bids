import pandas as pd


def reformat_table(table: pd.DataFrame):
    """Reformat a table to conform to the BIDS standard."""
    return table.assign(
        participant_id=table.SUBJLABEL.str.split("_").str[1:3].str.join("")
    ).drop(columns=["SUBJLABEL"])


def write_combined_table(combined_table: pd.DataFrame, out_path):
    combined_table.to_csv(out_path, sep="\t", na_rep="n/a", index=False)


def main(csv_path, out_path):
    write_combined_table(
        reformat_table(pd.read_csv(csv_path)),
        out_path,
    )


if __name__ == "__main__":
    main(snakemake.input.csv, snakemake.output.tsv)
