#######Use of matplotlib, seaborn and scatter plots####
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Load data from Tuesday
df = pd.read_csv('metada.csv')

# Clean resolution column
df['resolution'] = pd.to_numeric(df['resolution'], errors='coerce')
df['year'] = df['year'].astype(str)

# Set seaborn style
sns.set_style("whitegrid")
sns.set_palette("muted")

# ── Figure 1: Resolution histogram ──────────────────────────────
plt.figure(figsize=(8, 5))
sns.histplot(df['resolution'].dropna(), bins=30, kde=True, color='steelblue')
plt.xlabel('Resolution (Å)', fontsize=12)
plt.ylabel('Count', fontsize=12)
plt.title('Resolution Distribution of PDB Structures', fontsize=14)
plt.tight_layout()
plt.savefig('resolution_histogram.png', dpi=150)
plt.show()
print("Saved resolution histogram")

# ── Figure 2: Depositions per year (bar chart) ──────────────────
year_counts = df.groupby('year')['pdb_id'].count().reset_index()
year_counts.columns = ['year', 'count']

plt.figure(figsize=(12, 5))
sns.barplot(data=year_counts, x='year', y='count', color='steelblue')
plt.xticks(rotation=45, ha='right')
plt.xlabel('Year', fontsize=12)
plt.ylabel('Number of Structures', fontsize=12)
plt.title('PDB Structures Deposited Per Year', fontsize=14)
plt.tight_layout()
plt.savefig('depositions_per_year.png', dpi=150)
plt.show()
print("Saved deposition year chart")

# ── Figure 3: Experimental method pie chart ─────────────────────
method_counts = df['method'].value_counts()

plt.figure(figsize=(7, 7))
plt.pie(method_counts.values,
        labels=method_counts.index,
        autopct='%1.1f%%',
        startangle=90,
        colors=sns.color_palette("muted", len(method_counts)))
plt.title('PDB Structures by Experimental Method', fontsize=14)
plt.tight_layout()
plt.savefig('method_pie_chart.png', dpi=150)
plt.show()
print("Saved method pie chart")

# ── Figure 4: Resolution vs Year scatter ────────────────────────
plt.figure(figsize=(10, 5))
sns.boxplot(data=df.dropna(subset=['resolution']),
            x='year', y='resolution', color='lightblue')
plt.xticks(rotation=45, ha='right')
plt.xlabel('Year', fontsize=12)
plt.ylabel('Resolution (Å)', fontsize=12)
plt.title('Resolution Improvement Over Years', fontsize=14)
plt.tight_layout()
plt.savefig('resolution_vs_year.png', dpi=150)
plt.show()
print("Saved resolution vs year boxplot")
