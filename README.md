# SimCOVID-19

Simulation of the SARS-CoV-2 pandemic virus outbreak which cause the corona virus disease 2019 (COVID-19) by using a modified S. E. I. R. model.
**This repository is currently in a preview state! It should be used with caution!** Currently only German data has been processed, other countries will follow. Further questions can answered by <soeren.thiering@hs-anhalt.de>

## Data

Thanks for sharing the data by Johns Hopkins CSSE (<https://systems.jhu.edu/research/public-health/ncov/>).

## Model

**parameter:**

```R
k       = 0.3009135
x_max   = 83019213
te      = 5.1
ti      = 12
th      = 8
thi     = 12
nh_max  = 28031
kh      = 0.05262713
kd_1    = 0.002
kd_2    = 0.05262713
```

```R
if (u[4] > nh_max) kd = kd_1 * nh_max / u[4] + kd_2 * (u[4] - nh_max) / u[4]
else kd = kd_1
```

**Suspected:**

```R
du[1] = - k * u[1] / x_max * u[2]
```

**Exposed / Incubating:**

```R
du[2] = + k * u[1] / x_max * u[2] - k * h(p, t - te)[1] / x_max * h(p, t - te)[2] 
```

**Infected:**

```R
du[3] = + k * h(p, t - te)[1] / x_max * h(p, t - te)[2] - k * h(p, t - te - ti)[1] / x_max * h(p, t - te - ti)[2] * (1-kh)  - k * h(p, t - te - th)[1] / x_max * h(p, t - te - th)[2] * kh
```

**Hostspitalisation:**

```R
du[4] = + k * h(p, t - te - th)[1] / x_max * h(p, t - te - th)[2] * kh - k * h(p, t - te - th - thi)[1] / x_max * h(p, t - te - th - thi)[2] * kh
```

**Recovered:**

```R
du[5] = + k * h(p, t - te - ti)[1] / x_max * h(p, t - te - ti)[2] * (1-kh) + k * h(p, t - te - th - thi)[1] / x_max * h(p, t - te - th - thi)[2] * (kh - kd)
```

**Deaths:**

```
du[6] = + k * h(p, t - te - th - thi)[1] / x_max * h(p, t - te - th - thi)[2] * kd
```

## Results

### Current Situation

![Situation](Situation-1.png)

The growths rate is determined by logarithmic fitting of the confirmed cases in the last 10 days:

![Situation - logarithmic scale](Situation-2.png)

Basic Reproductive Number R0 has been determined from growths rate and median incubation time:

```R
R0 = exp(k * te) -1
```

![Situation - Basic Reproductive Number](Situation-3.png)

### Model vs Situation

![Model vs Situation](Model_vs_Situation-1.png)

### Forecast (100 Days)

![Forecast as Fraction](Forecast-2.png)

![Forecast ARDS](Forecast-ARDS-2.png)

## References

* <https://www.rki.de/DE/Content/InfAZ/N/Neuartiges_Coronavirus/Steckbrief.html>
