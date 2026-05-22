# Q: Should RIPPLE use `log(distance)` instead of linear distance in the Poisson GLM?

**Short answer:** No, RIPPLE stays with linear distance. The current model is `glm(counts ~ distance + offset(log(total_counts)), family = poisson)`, which under the canonical log link gives `rate ∝ exp(β · distance)` — an exponential decay/growth. A user suggested switching to `log(distance + 1)`, which would make the rate a power-law in distance. We considered it and decided against it. This note records the reasoning so we don't have to re-derive it every time the question comes up.

---

## Where the question comes from

The intuition is right that **pure diffusion in free space** produces a `1/r` (power-law) concentration profile from a point source. Anyone who has worked through Fick's second law for a continuous point source will reach for `log(distance)` as the "physical" choice for a paracrine-signaling model.

## Why the intuition doesn't translate to tissue

Very little real paracrine signaling is pure free-space diffusion. Three near-universal modifiers pull the profile toward **exponential**:

| Regime                                              | Profile                              | Functional form           |
|-----------------------------------------------------|---------------------------------------|---------------------------|
| Pure diffusion, no degradation (free 3D)             | `c(r) ~ 1/r`                          | Power-law                 |
| Diffusion + first-order degradation (steady state)   | `c(r) ~ exp(-r/λ)`, `λ = sqrt(D/k)`   | **Exponential**           |
| Restricted diffusion + receptor binding (tissue)     | near-exponential                      | **~Exponential**          |
| Source-sink steady state                             | exponential to bi-exponential         | **Exponential**           |

Every measured morphogen gradient I am aware of fits the exponential form: Bicoid (Driever & Nüsslein-Volhard 1988; Gregor et al. 2007), Hedgehog (Briscoe & Therond 2013), Wg/Wnt (Stamataki et al. 2005), Spätzle/Toll. Characteristic length scales 100–300 µm. Reviews: Lander 2007 (*Cell*), Wartlick et al. 2009 (*CSH Perspect. Biol.*), Müller et al. 2013 (*Development*).

For **chemokines/cytokines in lymphoid tissue** — RIPPLE's actual target — the published spatial-proteomics literature shows exponential or near-exponential profiles with length scales 20–150 µm (CXCL13 ~50 µm in B-cell follicles, CCL19/CCL21 ~100 µm in T zone, IL-2 ~50 µm from Treg sinks). Reviews: Oyler-Yaniv et al. 2017 (*Immunity*), Krummel et al. 2016 (*Nat. Rev. Immunol.*).

## The deeper problem: RIPPLE doesn't measure cytokine concentration

Even if the cytokine profile were power-law, **RIPPLE's dependent variable is gene expression in receiver cells, not cytokine concentration**. The receptor-binding → signaling → transcription pipeline introduces multiple non-linearities between the diffusion physics and what RIPPLE sees:

1. **Receptor density** varies by cell type and activation state (non-linear filter).
2. **Signaling thresholds** are switch-like below threshold, saturating above (sigmoidal filter).
3. **Transcription factor dynamics** add delay, half-life, and feedback (temporal smoothing).
4. **mRNA half-life** differs per transcript (further distortion of steady-state).

So the observed RNA-vs-distance curve is several non-linear transforms removed from the diffusion physics. **Neither linear nor log-distance is mechanistically correct** for the gene-expression readout. Both are empirical approximations.

## What matters for RIPPLE's job

RIPPLE's design intent is to **detect** whether expression depends on distance, reproducibly across replicates — not to fit the correct mechanistic decay law. For detection, the two models differ in:

- **Sensitivity to decay shape.** Linear fits exponentials well, power-laws less well. Log-distance fits power-laws well, exponentials less well. For empirical gene-expression curves (which are neither), both lose some power.
- **Robustness to short-distance outliers.** Log-distance puts disproportionate weight on the 0–10 µm bin because that's where the log-axis stretches out. Three noisy cells very close to query can drive the slope. Linear is more robust here.
- **The "per µm" interpretation.** Linear gives `β = log-rate change per µm` — directly interpretable in physical units, easy to communicate. Log-distance gives a unitless decay exponent (`β = −0.4` means "doubling distance halves the rate to the 0.5 power"), which is harder to explain and harder to compare across studies with different distance scales.
- **Reproducibility with the existing pipeline.** All cached results in `inst/extdata/*_cached/`, the benchmarks in `bench_*_results.rds`, and the paper figures use linear-distance gradient scores. Switching defaults invalidates all of them and requires re-validating null/power calibration under the new model.

## Decision: stay with linear

Two reasons:

1. **The empirical case for log-distance is weak.** Tissue paracrine signaling is well-approximated by exponentials in the published literature, and the gene-expression readout is too far removed from the diffusion physics for either functional form to be "correct". Linear-on-log-rate at least aligns with the dominant morphogen-with-degradation framework.
2. **The cost is high for unclear benefit.** Switching the default would invalidate all cached results, all paper figure numbers, and require re-validating the null and power benchmarks. The benchmarks currently show clean FDR control (1 FP / 7500 in null calibration) under the linear model. Without a concrete dataset where linear demonstrably mis-fits, that's a lot of churn for little gain.

## If you do want log-distance for a specific analysis

The change is a one-liner. The Poisson GLM lives at `R/glm.R:60-108` in `fit_poisson()`:

```r
glm(counts ~ distances + offset(log_total), family = poisson)
```

To fit with log-distance, transform the input vector before calling:

```r
fit_poisson(counts, distances = log1p(distances), total_counts = ...)
```

(Use `log1p(x) = log(x + 1)` to handle distance = 0 at the query cells; the `+ 1` µm pseudocount is a standard Wald-stable choice for count-data modelling.)

The returned `beta` is now a power-law decay exponent — interpret accordingly. The pipeline's downstream meta-analysis, FDR correction, and sign-consistency gate all still work; only the gradient-score *interpretation* changes.

## When to revisit

Open this question again if any of:

- A specific dataset shows clear evidence the linear model mis-fits (e.g. systematic non-zero residuals as a function of distance, or known long-range gradients where linear-fit's `β` saturates as `max_distance_um` increases).
- The paper reviewers ask for a sensitivity analysis under an alternative model.
- A user reports that switching to `log1p` ad-hoc gives qualitatively different top hits on their data — and the difference is biologically meaningful.

If any of these land, the next step would be to expose `distance_transform = c("linear", "log")` as an opt-in argument on `run_ripple()` (and `fit_poisson()`), re-run the benchmark suite under `"log"` to confirm FDR control, and add a vignette section walking through the comparison. Estimated cost: ~half a day of implementation + ~30 min of benchmark re-run + ~2 hours of vignette.

---

*This note was added 2026-05-04 as part of the contamination → broad rename PR. It records the reasoning behind a deliberate non-change.*
