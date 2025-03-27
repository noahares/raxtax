import matplotlib.pyplot as plt
import pandas as pd
import argparse
import seaborn as sns
from matplotlib.ticker import FuncFormatter

sns.set_context("paper")
sns.set_theme(style="whitegrid")
plt.rcParams.update(
        {
            "pgf.texsystem": "pdflatex",
            "font.family": "sans-serif",
            "font.size": 10,
            "text.usetex": True,
            "pgf.rcfonts": False,
            }
        )

def compare_speedup(df1, df2, out_path, weak):
    program_colors = ['tomato', 'tab:cyan', 'black']
    labels = ['raxtax', 'raxtax TP']
    fig, ax1 = plt.subplots(figsize=(9, 6))
    for i, df in enumerate([df1, df2]):

        baseline = df[df["Threads"] == 1].set_index("Rep")["RuntimeSeconds"]

        if weak:
            df["Speedup"] = df.apply(lambda row: baseline[row["Rep"]] / row["RuntimeSeconds"], axis=1)
        else:
            df["Speedup"] = df.apply(lambda row: (baseline[row["Rep"]] / row["RuntimeSeconds"]) / row['Threads'], axis=1)

        program_stats = df.groupby("Threads")["Speedup"].agg(["mean", "std"]).reset_index()
        print(program_stats)
        sns.lineplot(data=program_stats, x="Threads", y="mean", color=program_colors[i], label=labels[i], ax=ax1, marker='o')

        plt.fill_between(program_stats["Threads"],
                        program_stats["mean"] - program_stats["std"],
                        program_stats["mean"] + program_stats["std"],
                        alpha=0.2)
    plt.axhline(y=1, color=program_colors[2], linestyle='--', label='ideal')
    ax1.set_xlabel("Threads")
    if weak:
        ax1.set_ylabel('rel. Speedup')
    else:
        ax1.set_ylabel('Efficiency')

    ax1.set_xscale("log", base=2)
    ax1.set_ylim(0.0, 1.1)
    ax1.set_xticks(df1['Threads'].unique())
    ax1.xaxis.set_major_formatter(FuncFormatter(lambda x, _: f'{int(x)}'))
    ax1.legend(loc='upper center', bbox_to_anchor=(0.5, 1.15), ncol=2, frameon=False, handletextpad=1.5)
    plt.savefig(out_path)
    plt.close()

def plot_speedup(df, out_path, weak):
    fig, ax1 = plt.subplots(figsize=(9, 6))

    program_colors = ['tomato', 'tab:cyan']
    baseline = df[df["Threads"] == 1].set_index("Rep")["RuntimeSeconds"]

    if weak:
        df["Speedup"] = df.apply(lambda row: baseline[row["Rep"]] / row["RuntimeSeconds"], axis=1)
    else:
        df["Speedup"] = df.apply(lambda row: (baseline[row["Rep"]] / row["RuntimeSeconds"]) / row['Threads'], axis=1)

    program_stats = df.groupby("Threads")["Speedup"].agg(["mean", "std"]).reset_index()
    print(program_stats)
    sns.lineplot(data=program_stats, x="Threads", y="mean", color=program_colors[0], label="raxtax", ax=ax1, marker='o')

    plt.fill_between(program_stats["Threads"],
                    program_stats["mean"] - program_stats["std"],
                    program_stats["mean"] + program_stats["std"],
                    alpha=0.2)
    plt.axhline(y=1, color=program_colors[1], linestyle='--', label='ideal')
    ax1.set_xlabel("Threads")
    if weak:
        ax1.set_ylabel('rel. Speedup')
    else:
        ax1.set_ylabel('Efficiency')


    ax1.set_xscale("log", base=2)
    ax1.set_ylim(0.0, 1.1)
    ax1.set_xticks(df['Threads'].unique())
    ax1.xaxis.set_major_formatter(FuncFormatter(lambda x, _: f'{int(x)}'))
    ax1.legend(loc='upper center', bbox_to_anchor=(0.5, 1.15), ncol=2, frameon=False, handletextpad=1.5)
    plt.savefig(out_path)
    plt.close()

def plot_results(stats, field, y_field, out_path):
    fig, ax1 = plt.subplots(figsize=(9, 6))
    program_colors = ['tomato', 'tab:cyan']

    programs = stats['Tool'].unique()

    for i, program in enumerate(programs):
        program_stats = stats[stats['Tool'] == program]
        if "Bytes" in y_field:
            program_stats[y_field] = program_stats[y_field] / (1024**3)
        sns.lineplot(data=program_stats, x=field, y=y_field, errorbar='sd', color=program_colors[i % len(program_colors)], label=program, ax=ax1, marker='o')

    ax1.set_xlabel("Dataset Size")
    if y_field == "RuntimeSeconds":
        ax1.set_ylabel("Time (s)")
    elif y_field == "MaxMemoryBytes":
        ax1.set_ylabel('Max Memory (GiB)')
    ax1.set_xticks(stats[field].unique())
    ax1.tick_params(axis='x', labelrotation=45)
    ax1.xaxis.set_major_formatter(FuncFormatter(lambda x, _: f'{int(x):,}'))
    ax1.yaxis.set_major_formatter(FuncFormatter(lambda x, _: f'{int(x):,}'))
    ax1.legend(loc='upper center', bbox_to_anchor=(0.5, 1.15), ncol=2, frameon=False, handletextpad=1.5)
    plt.tight_layout()
    plt.savefig(out_path)
    plt.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Plot runtime and memory results")
    parser.add_argument("-i", dest="input_csv", type=str, help="CSV file")
    parser.add_argument("-j", dest="input_csv_2", type=str, help="CSV file")
    parser.add_argument("-f", dest="field", type=str, help="Group Field")
    parser.add_argument("-y", dest="y_field", type=str, help="Group Field")
    parser.add_argument("-o", dest="out_path", type=str, help="Output path")
    parser.add_argument("-w", dest="weak", type=str, help="Weak Speedup")
    args = parser.parse_args()
    stats = pd.read_csv(args.input_csv)
    if args.field == "Threads":
        if args.input_csv_2:
            stats2 = pd.read_csv(args.input_csv_2)
            compare_speedup(stats, stats2, args.out_path, args.weak)
        else:
            plot_speedup(stats, args.out_path, args.weak)
    else:
        plot_results(stats, args.field, args.y_field, args.out_path)
