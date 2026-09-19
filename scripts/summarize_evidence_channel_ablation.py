"""Read-only scientific summaries for validated channel-ablation artefacts."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd


def _policy(channels: list[str]) -> str:
    direct = any(not value.endswith("_transferred") for value in channels)
    transferred = any(value.endswith("_transferred") for value in channels)
    return "direct+transferred" if direct and transferred else ("direct-only" if direct else "transferred-only")


def run(output: Path) -> None:
    summary = json.loads((output / "evidence_channel_ablation_summary.json").read_text(encoding="utf-8"))
    conditions = pd.DataFrame(summary["conditions"]).rename(columns={"id": "Condition", "family": "Type"})
    conditions["Channels_Retained"] = conditions["included_channels"].map(lambda values: ";".join(values))
    conditions["Channels_Removed"] = conditions["excluded_channels"].map(lambda values: ";".join(values))
    conditions["Direct_Transferred_Policy"] = conditions["included_channels"].map(_policy)
    conditions = conditions.drop(columns=["included_channels", "excluded_channels"])

    bias = pd.read_csv(output / "evidence_channel_ablation_bias.csv")
    d10 = bias.loc[(bias["Rank_Domain"] == "absolute") & (bias["Response_Scope"] == "meaningful") & (bias["Stratum"] == "D10")]
    degree = d10.loc[d10["Audit"] == "Degree", ["Condition", "Tissue", "Target", "Median_Promotion"]].rename(columns={"Median_Promotion": "Degree_D10_Promotion"})
    pagerank = d10.loc[d10["Audit"] == "PageRank", ["Condition", "Tissue", "Target", "Median_Promotion"]].rename(columns={"Median_Promotion": "PageRank_D10_Promotion"})
    target = degree.merge(pagerank, on=["Condition", "Tissue", "Target"], validate="one_to_one")

    matched = pd.read_csv(output / "evidence_channel_matched_controls.csv")
    matched = matched.loc[(matched["Audit"] == "Coverage_Matched_Control") & (matched["Response_Scope"] == "meaningful_response")]
    matched = matched.groupby(["Condition", "Tissue", "Target"], as_index=False)["Promotion_Difference_High_Minus_Low"].median().rename(columns={"Promotion_Difference_High_Minus_Low": "Annotation_Matched_Difference"})
    annotation = pd.read_csv(output / "evidence_channel_ablation_annotation.csv")
    annotation = annotation.loc[(annotation["Metric"] == "Annotation_Coverage") & (annotation["Response_Scope"] == "meaningful_response"), ["Condition", "Tissue", "Target", "Spearman_Rho"]].rename(columns={"Spearman_Rho": "Promotion_vs_Annotation_Spearman"})
    target = target.merge(matched, on=["Condition", "Tissue", "Target"], validate="one_to_one").merge(annotation, on=["Condition", "Tissue", "Target"], validate="one_to_one")
    c1 = target.loc[target["Condition"] == "C1", ["Tissue", "Target", "Degree_D10_Promotion", "PageRank_D10_Promotion", "Annotation_Matched_Difference"]].rename(columns={
        "Degree_D10_Promotion": "C1_Degree_D10_Promotion", "PageRank_D10_Promotion": "C1_PageRank_D10_Promotion", "Annotation_Matched_Difference": "C1_Annotation_Matched_Difference",
    })
    target = target.merge(c1, on=["Tissue", "Target"], validate="many_to_one")
    target["Degree_Bias_Reduction_vs_C1"] = target["C1_Degree_D10_Promotion"] - target["Degree_D10_Promotion"]
    target["PageRank_Bias_Reduction_vs_C1"] = target["C1_PageRank_D10_Promotion"] - target["PageRank_D10_Promotion"]
    target["Annotation_Bias_Reduction_vs_C1"] = target["C1_Annotation_Matched_Difference"] - target["Annotation_Matched_Difference"]

    aggregate = target.groupby("Condition", as_index=False).agg(
        Degree_D10_Aggregate_Mean=("Degree_D10_Promotion", "mean"), Degree_D10_Aggregate_Median=("Degree_D10_Promotion", "median"), Degree_D10_Positive_Targets=("Degree_D10_Promotion", lambda s: int((s > 0).sum())),
        PageRank_D10_Aggregate_Mean=("PageRank_D10_Promotion", "mean"), PageRank_D10_Aggregate_Median=("PageRank_D10_Promotion", "median"), PageRank_D10_Positive_Targets=("PageRank_D10_Promotion", lambda s: int((s > 0).sum())),
        Annotation_Matched_Aggregate_Median=("Annotation_Matched_Difference", "median"), Annotation_Matched_Positive_Targets=("Annotation_Matched_Difference", lambda s: int((s > 0).sum())),
        Promotion_vs_Annotation_Spearman_Median=("Promotion_vs_Annotation_Spearman", "median"),
        Degree_Bias_Reduction_Mean=("Degree_Bias_Reduction_vs_C1", "mean"), Degree_Bias_Improved_Targets=("Degree_Bias_Reduction_vs_C1", lambda s: int((s > 0).sum())),
        PageRank_Bias_Reduction_Mean=("PageRank_Bias_Reduction_vs_C1", "mean"), PageRank_Bias_Improved_Targets=("PageRank_Bias_Reduction_vs_C1", lambda s: int((s > 0).sum())),
        Annotation_Bias_Reduction_Median=("Annotation_Bias_Reduction_vs_C1", "median"), Annotation_Bias_Improved_Targets=("Annotation_Bias_Reduction_vs_C1", lambda s: int((s > 0).sum())),
    )

    benchmark = pd.read_csv(output / "evidence_channel_ablation_benchmarks.csv")
    benchmark = benchmark.groupby(["Condition", "Reference"], as_index=False).median(numeric_only=True)
    benchmark = benchmark.pivot(index="Condition", columns="Reference")
    benchmark.columns = [f"{metric}_{reference}" for metric, reference in benchmark.columns]
    benchmark = benchmark.reset_index()
    health = pd.read_csv(output / "evidence_channel_ablation_network_health.csv").groupby("Condition", as_index=False).median(numeric_only=True)
    health = health.drop(columns=[column for column in ("Nodes", "Edges") if column in health])
    master = conditions.merge(aggregate, on="Condition", validate="one_to_one").merge(benchmark, on="Condition", validate="one_to_one").merge(health, on="Condition", validate="one_to_one")
    master.to_csv(output / "evidence_channel_ablation_master_condition_table.csv", index=False, encoding="utf-8-sig")
    target.to_csv(output / "evidence_channel_ablation_target_bias_reduction.csv", index=False, encoding="utf-8-sig")

    availability_rows = []
    for axis in ("degree", "pagerank"):
        values = pd.read_csv(output / f"evidence_channel_availability_by_{axis}.csv")
        d1 = values.loc[values["Stratum"] == "D1"].set_index(["Tissue", "Channel"])
        d10 = values.loc[values["Stratum"] == "D10"].set_index(["Tissue", "Channel"])
        joined = d10.join(d1, lsuffix="_D10", rsuffix="_D1").reset_index()
        joined["Axis"] = axis
        joined["Availability_D10_minus_D1"] = joined["Mean_Node_Nonzero_Fraction_D10"] - joined["Mean_Node_Nonzero_Fraction_D1"]
        joined["MedianScore_D10_minus_D1"] = joined["Median_Node_Channel_Score_D10"] - joined["Median_Node_Channel_Score_D1"]
        availability_rows.append(joined)
    pd.concat(availability_rows, ignore_index=True).to_csv(output / "evidence_channel_ablation_availability_d10_vs_d1.csv", index=False, encoding="utf-8-sig")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    run(parser.parse_args().output)
