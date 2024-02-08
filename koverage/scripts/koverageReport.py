import pandas as pd
import datapane as dp
from plotly.subplots import make_subplots
import plotly.graph_objects as go
from datetime import datetime


def create_title(sample_names, ref_fasta):
    """Generate the title for the html report incorporating information from the run

    Args:
        sample_names (list): list of sample names (str)
        ref_fasta (str): reference fasta file path
    """
    current_date = f"{datetime.now():%d-%m-%Y}"
    sample_list_str = ", ".join(str(item) for item in sample_names)
    wonk = "\n".join(
        [
            "# Koverage Report",
            "### Report Date: " + current_date,
            "### Reference Sequence: " + ref_fasta,
            "### Reads Processed: " + sample_list_str + " ###",
        ]
    )
    title = dp.Group(
        dp.Text(
            "![](https://raw.githubusercontent.com/beardymcjohnface/Koverage/main/koverage.png)"
        ),
        dp.Text(wonk),
        columns=2,
    )
    return title


def qualgraph(sample_name, df, ref_fa):
    """Create graph for report

    Args:
        sample_name (str): name of sample
        df (dataframe): pandas dataframe from which to plot
    """
    df = df[(df["Sample"] == sample_name)]
    df = df.sort_values("Count", ascending=False)
    fig = make_subplots(specs=[[{"secondary_y": True}]])
    fig.add_trace(
        go.Bar(x=df["Contig"], y=df["Count"], name="Count"), secondary_y=False
    )
    fig.add_trace(
        go.Scatter(x=df["Contig"], y=df["Mean"], name="Mean Depth"), secondary_y=True
    )
    fig.update_xaxes(title_text=ref_fa + " Contig Number")
    fig.update_yaxes(title_text="Count", secondary_y=False)
    fig.update_yaxes(title_text="Mean", secondary_y=True)
    fig.update_layout(title=sample_name)
    return dp.Group(dp.Plot(fig), dp.DataTable(df), label=sample_name)


def create_layout(sample_names, graphs):
    """Decide whether to use tabbed layout or not (i.e. if more than one sample)

    Args:
        sample_names (list): list of sample names (str)

    Returns:
        sample_coverage (datapane.Block): dp.block object for laying out samples
    """
    if len(sample_names) > 1:
        sample_coverage = dp.Blocks(
            dp.Text("## Sample Coverage"), dp.Select(blocks=[*graphs])
        )
    else:
        sample_coverage = dp.Blocks(dp.Text("## Sample Coverage"), blocks=[*graphs])
    return sample_coverage


def generate_buttons(all_df):
    """Generate update buttons for use in figure pane generation

    Args:
        headers (list): List of columns in dataframe

    Returns:
        buttons (list): List of dict() objects for button creation in plotly.graph_object.Figure()
    """
    buttons = []
    head = list(all_df)
    for i in head:
        if i != "Contig":
            plonk = dict(label=i, method="update", args=[{"y": [all_df[i]]}])
            buttons.append(plonk)
    return buttons


def generate_figure(all_df, ref_fa, buttons):
    """Generate the figure panel

    Args:
        all_df (dataframe): pandas dataframe of all_coverage.tsv
        ref_fa (str): reference fasta file
        buttons (list): list of button objects dict()

    Retruns:
        fig (go.Figure()): Figure object generated with plotly.graph_object.Figure()
    """
    fig = go.Figure()
    fig.add_trace(go.Bar(x=all_df["Contig"], y=all_df["Count"], name="Count"))
    fig.update_xaxes(title_text=ref_fa + " Contig Number")
    fig.update_layout(
        title_text="All Coverage",
        autosize=True,
        updatemenus=[
            dict(
                type="buttons",
                bgcolor="mediumspringgreen",
                bordercolor="black",
                xanchor="left",
                yanchor="top",
                direction="left",
                pad={"r": 10, "t": 10},
                x=0.11,
                y=1.4,
                showactive=True,
                buttons=buttons,
            )
        ],
    )
    return fig


def main(**kwargs):
    # if kwargs["pyspy"]:
    #     subprocess.Popen(
    #         [
    #             "py-spy",
    #             "record",
    #             "-s",
    #             "-o",
    #             kwargs["pyspy_svg"],
    #             "--pid",
    #             str(os.getpid()),
    #         ]
    #     )
    # Read in data frames
    sample_df = pd.read_csv(kwargs["sample_cov"], sep="\t")
    all_df = pd.read_csv(kwargs["all_cov"], sep="\t")

    # shorten contig IDs ### TODO update to work with Phables output
    # sample_df["Contig"] = sample_df["Contig"].str.extract(
    #     r"([A-Za-z0-9]+_[A-Za-z0-9]+)"
    # )
    # all_df["Contig"] = all_df["Contig"].str.extract(r"([A-Za-z0-9]+_[A-Za-z0-9]+)")

    # Filter for top n contigs
    all_df = all_df.sort_values("Count", ascending=False).head(kwargs["max_ctg"])
    sample_df = sample_df[sample_df["Contig"].isin(all_df["Contig"])]

    # generate graphs
    graphs = []
    for sample in kwargs["sample_names"]:
        graphs.append(qualgraph(sample, sample_df, kwargs["ref_fasta"]))

    # layout
    report_title = create_title(kwargs["sample_names"], kwargs["ref_fasta"])
    sample_layout = create_layout(kwargs["sample_names"], graphs)
    panels = [
        dp.Group(
            dp.Text(kwargs["sample_cov_desc"]), sample_layout, label="Sample Coverage"
        )
    ]
    buttons = generate_buttons(all_df)
    figure = generate_figure(all_df, kwargs["ref_fasta"], buttons)
    all_cov_pane = dp.Group(
        dp.Text(str(kwargs["all_cov_desc"])),
        dp.Plot(figure),
        dp.DataTable(all_df),
        label="All Coverage",
    )
    panels.append(all_cov_pane)
    report = dp.Blocks(report_title, dp.Select(type=dp.SelectType.TABS, blocks=panels))
    dp.save_report(report, kwargs["out_file"])


if __name__ == "__main__":
    main(
        log_file=snakemake.log.err,
        sample_cov=snakemake.input.smpl,
        all_cov=snakemake.input.all,
        out_file=snakemake.output.html,
        sample_cov_desc=snakemake.params.sample_cov_desc,
        all_cov_desc=snakemake.params.all_cov_desc,
        sample_names=snakemake.params.sample_names,
        ref_fasta=snakemake.params.ref_fasta,
        max_ctg=snakemake.params.max_ctg,
        # pyspy=snakemake.params.pyspy,
        # pyspy_svg=snakemake.log.pyspy,
    )
