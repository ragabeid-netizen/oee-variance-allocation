# Variance allocation for Overall Equipment Effectiveness

Data and code for the manuscript

> **Why Overall Equipment Effectiveness Tracks Its Performance Term: Variance Allocation and
> Operating-Range Asymmetry in a Six-Year Plant Record**
> R. E. R. Shehata, Y. Shaban, N. Ayoub, M. M. Maghawry, S. A. Abdelwahab

Everything reported in the paper — every table, every figure and every number in the text —
is reproduced by one script from the two data files in this repository.

```
python code/reproduce.py
```

Requires Python 3.9+, `pandas` and `numpy`. Runtime is about two minutes, most of it the
25,000 bootstrap replicates.

---

## What the study asks

OEE is the product of an Availability, a Performance and a Quality rate, and is read across
manufacturing as a three-way diagnostic: a plant reporting a low figure is expected to be able
to say whether the shortfall came from stoppages, from running below rate, or from defective
output. That reading holds only if the three rates actually differentiate.

Decomposing the variance of `ln OEE` over 2,189 machine-shifts gives Performance 100.99% of it,
Availability −1.10% and Quality 0.11%. The composite tracks its Performance term at r = 0.9995.

The obvious objection is that Availability only looks inert because the plant records it once a
year while Performance is recorded every shift. That objection is tested rather than assumed:
Availability is rebuilt from the same maintenance register at monthly and then at daily
resolution, a 365-fold refinement that quadruples its dispersion, and its allocated share moves
only to +0.35%. What operates instead is an asymmetry of operating range.

---

## Contents

| Path | What it is |
|---|---|
| `data/shift_records_2020_2025.csv` | 2,190 machine-shift records, January 2020 – December 2025. Five packaging machines, two plants, three product families. |
| `data/maintenance_daily_register.csv` | The daily maintenance register the plant derives Availability from, 2020–2025. |
| `data/openalex_search_2026-08-15.csv` | The twelve bibliographic queries reported in Section 2.3, with the record counts each returned and the date they were run. |
| `code/reproduce.py` | Reproduces every reported quantity. |
| `figures/` | The four figures, 400 dpi. |

Product and plant names are replaced by neutral labels (`Product A`–`C`, `Plant 1`–`2`). Nothing
else is altered: no record is removed, no outlier is trimmed and no value is imputed.

Two of the 2,190 records carry a missing field — one lacks a product label, one lacks a waste
entry and therefore a Quality rate. Both are excluded from the decomposition, giving n = 2,189.
Retaining them with the missing fields imputed at the machine mean changes no reported result.

---

## Data provenance, and one open item

The three rates do not come from the same system, and this matters enough to state here as well
as in the paper.

**Performance** and **Quality** are computed for each machine and each shift from that shift's
own counts.

**Availability** is not. The plant derives a single annual figure from the maintenance record of
one unit and applies it unchanged to every machine and every shift of that year. The Availability
column therefore holds six distinct values across six years — at shift resolution it is a year
effect, not a shift-varying quantity.

> **Open item.** Those six annual figures are being re-verified against the plant's original
> paper registers, and this repository will be updated when that is complete. Readers should
> treat the Availability *levels* as provisional. The result reported in the paper is a statement
> about *dispersion* and does not depend on them: shifting all six values down by twenty
> percentage points moves the allocated share of Availability from −1.10% to −1.36%, and the six
> figures would have to run from roughly 71% to 100% before Availability was allocated even 5% of
> the variance. `code/reproduce.py` recomputes both of those checks.

---

## Robustness, in brief

- **23 strata** formed by six stratification variables: the Performance share stays between
  98.73% and 101.29%.
- **Five resampling schemes**, 5,000 replicates each — i.i.d., three cluster bootstraps
  (machine-month, machine-year, calendar-year) and a moving-block bootstrap. The shifts are
  serially dependent (lag-1 ρ = +0.65, effective n ≈ 460), which is why the ordinary bootstrap is
  not reported alone. The most conservative scheme still places Performance above 98.9%.
- **Covariance sweep**: holding the coefficients of variation at their observed values and
  sweeping the correlation between the log rates across −0.9 to +0.9 never lifts the magnitude of
  the Availability share above 2.54%.

---

## Licence

Data and code are released under [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/).
Please cite the manuscript above.

## Contact

Ragab Eid Ragab Shehata — ragabeid2029@gmail.com
Production Technology Department, Faculty of Technology and Education, Helwan University,
Cairo, Egypt
