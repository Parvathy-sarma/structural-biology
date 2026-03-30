##SciPy Stats & Structural Analysis

###Do B-factors differ between high and low resolution structures?

**Biological question:** Proteins that are more ordered and rigid tend to 
diffract better in X-ray crystallography, giving higher resolution structures. 
I tested whether this is statistically confirmed by comparing Cα B-factors 
between high-resolution (<1.5 Å) and low-resolution (>3.0 Å) PDB structures.

**Data:**
- High-resolution group: 22 structures (resolution < 1.5 Å)
- Low-resolution group: 45 structures (resolution > 3.0 Å)

**Results:**

| Group | Mean B-factor (Ų) | Std | Min | Max |
|---|---|---|---|---|
| High-res (<1.5 Å) | 21.05 | 4.09 | 11.11 | 29.50 |
| Low-res (>3.0 Å) | 54.01 | 33.92 | 15.00 | 148.30 |

**Statistical tests:**

| Test | Statistic | p-value | Result |
|---|---|---|---|
| t-test | -4.47 | 3.18e-05 | Significant |
| Mann-Whitney U | 166.0 | 1.16e-05 | Significant |
| Shapiro-Wilk (high-res) | — | 0.259 | Normal |
| Shapiro-Wilk (low-res) | — | 0.001 | NOT normal |

**Interpretation:**

Low-resolution structures have significantly higher B-factors 
(mean 54 Ų) compared to high-resolution structures (mean 21 Ų). 
This difference is highly significant (p < 0.0001).

Since the low-resolution group failed the normality test 
(Shapiro-Wilk p = 0.001), the Mann-Whitney U test is the more 
appropriate result here — and it confirms the same conclusion.

The wide spread in the low-resolution group (std = 33.9, max = 148.3 Ų) 
reflects genuine biological diversity — some low-resolution structures 
are flexible proteins that are simply hard to crystallise well, 
while others are large complexes where parts of the structure are 
inherently disordered.

**Conclusion:** B-factors are a reliable proxy for structural order, 
and high-resolution crystal structures do indeed correspond to 
more rigid, well-ordered proteins. This has direct implications 
for ML models — B-factors can be used as a feature for predicting 
crystallisability and structural quality.


