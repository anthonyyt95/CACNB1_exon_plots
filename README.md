# Exon Expression Visualization and Analysis

This MATLAB script suite provides tools to visualize and analyze exon-level expression data across conditions (e.g., brain vs. immune cells). The primary functions include bar plots with error bars, exon structure diagrams, and transcript-level comparisons.

---

## Features

### 1. **Exon Expression Bar Plots**
- Generates side-by-side bar plots for each exon comparing two conditions.
- Applies `log10(x + 1)` transformation to raw counts.
- Includes standard deviation error bars.
- Optional: overlays individual replicate scatter points.

### 2. **Exon Structure Diagram**
- Creates a linear schematic of exon start and end positions.
- Useful for illustrating transcript structure visually.

### 3. **Transcript Exon Extraction**
- Aggregates unique exons across multiple transcripts.
- Tracks exon usage across isoforms.
- Outputs exon ID, genomic coordinates, and frequency.

### 4. **ENSM Exon ID Mapping**
- Extracts transcript-specific exon IDs and matches them to corresponding ENSM exon identifiers.
- Produces matrix format for cross-transcript comparison.

### 5. **Data Aggregation by Reference**
- Averages replicate exon expression values using a reference table.
- Merges data across replicates for clean bar and error plots.

---

## Inputs

- `.csv` files with exon count data (first column: exon number, remaining columns: sample expression).
- Reference tables for exon structures or transcript assignments.
- Transcript ID list (optional).

---

## Outputs

- Bar plots of exon expression (`log10+1`) with error bars and optional replicates.
- Schematic diagrams of exon structure.
- Aggregated exon frequency and mapping tables.

---

## Requirements

- MATLAB with access to:
  - `importdata`
  - `scatter`, `bar`, `plot`, `fill` functions
- CSV-formatted input data

---

## Notes

- Customize `cond1`, `cond2`, and `bColor` variables to fit your dataset and visual style.
- Modify Y-axis tick labels to represent original (unlogged) expression values if needed.
- Designed for flexible reuse with multiple datasets in a structured directory.

