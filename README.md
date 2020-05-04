# SimCOVID-19

Simulation of the SARS-CoV-2 pandemic virus outbreak which cause the corona virus disease 2019 (COVID-19) by using a modified S. E. I. R. D. model.
**This repository is currently in a preview state! It should be used with caution!** Currently only German data has been processed, other countries will follow. Further questions can answered by <soeren.thiering@hs-anhalt.de>

## Data

Thanks for sharing the data by Johns Hopkins CSSE (<https://systems.jhu.edu/research/public-health/ncov/>).

## Model

**Parameter:**

| Name                                                                                                            |         Value | Description                                | Source |
|-----------------------------------------------------------------------------------------------------------------|--------------:|--------------------------------------------|--------|
| ![$k(t)$](https://render.githubusercontent.com/render/math?math=%24k(t)%24)                                     | ~ 0.15 to 0.3 | growth rate                                |        |
| ![$n_{max}$](https://render.githubusercontent.com/render/math?math=%24n_%7Bmax%7D%24)                           |      83019213 | max inahbitants                            |        |
| ![$t_{e}$](https://render.githubusercontent.com/render/math?math=%24t_%7Be%7D%24)                               |           5.1 | exposed / incubation time                  |        |
| ![$t_{i}$](https://render.githubusercontent.com/render/math?math=%24t_%7Bi%7D%24)                               |            12 | time of a moderate  infection              |        |
| ![$t_{h}$](https://render.githubusercontent.com/render/math?math=%24t_%7Bh%7D%24)                               |             8 | time until hostspitalisation               |        |
| ![$t_{h,i}$](https://render.githubusercontent.com/render/math?math=%24t_%7Bh%2Ci%7D%24)                         |            14 | time for hostspitalisation with ARDS       |        |
| ![$t_{h,d}$](https://render.githubusercontent.com/render/math?math=%24t_%7Bh%2Cd%7D%24)                         |             2 | time until death                           |        |
| ![$n_{h,max}$](https://render.githubusercontent.com/render/math?math=%24n_%7Bh%2Cmax%7D%24)                     |         28031 | max intensive beds for treating ARDS       |        |
| ![$\varphi_{h}$](https://render.githubusercontent.com/render/math?math=%24%5Cvarphi_%7Bh%7D%24)                 |          0.14 | propotion hostspitalisation with ARDS      |        |
| ![$\varphi_{d}(u_h(t))$](https://render.githubusercontent.com/render/math?math=%24%5Cvarphi_%7Bd%7D(u_h(t))%24) |        ~ 0.05 | propotion deaths                           |        |
| ![$\varphi_{d,min}$](https://render.githubusercontent.com/render/math?math=%24%5Cvarphi_%7Bd%2Cmin%7D%24)       |    0.02282158 | propotion deaths with hostspitalisation    |        |
| ![$\varphi_{d,max}$](https://render.githubusercontent.com/render/math?math=%24%5Cvarphi_%7Bd%2Cmax%7D%24)       |    0.09687916 | propotion deaths without hostspitalisation |        |

**Growth Rate:**

![formula](https://render.githubusercontent.com/render/math?math=k(t)%20%3D%20%20%20%5Cbegin%7Bcases%7D%20%20%20%20%20k_%7Bi%2B1%7D%20%20%20%20%20%26%20%5Cquad%20%5Ctext%7B%20if%20%7D%20t%20%3E%20t_%7Bk%2Cn%2B1%7D%20%5C%5C%20%20%20%20%20k_%7Bi%2B1%7D%20%20%20%20%20%26%20%5Cquad%20%5Ctext%7B%20if%20%7D%20t%20%3E%20t_%7Bk%2Cn%7D%20%5Cwedge%20t%20%5Cleq%20t_%7Bk%2Cn%2B1%7D%20%5C%5C%20%20%20%20%20%5Ccdots%20%20%20%20%20%20%26%20%5Cquad%20%5C%5C%20%20%20%20%20k_%7B1%7D%20%20%20%20%20%20%20%26%20%5Cquad%20%5Ctext%7B%20if%20%7D%20t%20%3E%20t_%7Bk%2C1%7D%20%5Cwedge%20t%20%5Cleq%20t_%7Bk%2C2%7D%20%20%20%5Cend%7Bcases%7D)

**Deaths:**

![formula](https://render.githubusercontent.com/render/math?math=%5Cvarphi_%7Bd%7D(u_h(t))%20%20%5Cbegin%7Bcases%7D%20%20%5Cvarphi_%7Bd%2Cmin%7D%20%5Ccfrac%20%7B1-u_h(t)%7D%20%7Bn_%7Bh%2Cmax%7D%7D%20%2B%20%5Cvarphi_%7Bd%2Cmax%7D%20%5Ccfrac%20%7Bu_h(t)%7D%20%7Bn_%7Bh%2Cmax%7D%7D%20%26%20%5Cquad%20%5Ctext%7B%20if%20%7D%20%5Ccfrac%20%7Bu_h(t)%7D%20%7Bn_%7Bh%2Cmax%7D%7D%20%3C%201%20%5C%5C%20%20%5Cvarphi_%7Bd%2Cmax%7D%20%26%20%5Cquad%20%5Ctext%7B%20if%20%7D%20%5Ccfrac%20%7Bu_h(t)%7D%20%7Bn_%7Bh%2Cmax%7D%7D%20%5Cgeq%201%20%5C%5C%20%20%5Cend%7Bcases%7D)

**Susceptibles:**

![formula](https://render.githubusercontent.com/render/math?math=%5Ccfrac%7Bdu_7(t)%7D%20%7Bdt%7D%20%3D%20-%20k(t)%20%5Ccfrac%7Bu_7(t)%7D%20%7Bn_%7Bmax%7D%7D%20u_7(t))

**Exposed / Incubating:**

![formula](https://render.githubusercontent.com/render/math?math=%5Ccfrac%7Bdu_6(t)%7D%20%7Bdt%7D%20%20%3D%20%20%20%20%20%2B%20k(t)%20%5Ccfrac%7Bu_7(t)%7D%20%7Bn_%7Bmax%7D%7D%20u_6(t)%20%20%20%20%20-%20k(t-t_e)%20%5Ccfrac%7Bu_7(t-t_e)%7D%20%7Bn_%7Bmax%7D%7D%20u_6(t-t_e))

**Infected:**

![formula](https://render.githubusercontent.com/render/math?math=%5Ccfrac%7Bdu_5(t)%7D%20%7Bdt%7D%20%20%3D%20%20%20%20%20%2B%20k(t-t_e)%20%5Ccfrac%7Bu_7(t-t_e)%7D%20%7Bn_%7Bmax%7D%7D%20u_6(t-t_e)%20%20%20%20%20-%20k(t-t_e-t_i)%20%5Ccfrac%7Bu_7(t-t_e-t_i)%7D%20%7Bn_%7Bmax%7D%7D%20u_6(t-t_e-t_i)%20(1-%5Cvarphi_%7Bh%7D)%20%20%20%20%20-%20k(t-t_e-t_h)%20%5Ccfrac%7Bu_7(t-t_e-t_h)%7D%20%7Bn_%7Bmax%7D%7D%20u_6(t-t_e-t_h)%20(%5Cvarphi_%7Bh%7D))

**Hostspitalisation:**

![formula](https://render.githubusercontent.com/render/math?math=%5Ccfrac%7Bdu_4(t)%7D%20%7Bdt%7D%20%20%3D%20%20%20%20%20%2B%20k(t-t_e-t_h)%20%5Ccfrac%7Bu_7(t-t_e-t_h)%7D%20%7Bn_%7Bmax%7D%7D%20u_6(t-t_e-t_h)%20(%5Cvarphi_%7Bh%7D)%20%20%20%20%20-%20k(t-t_e-t_h-t_%7Bh%2Ci%7D)%20%5Ccfrac%7Bu_7(t-t_e-t_h-t_%7Bh%2Ci%7D)%7D%20%7Bn_%7Bmax%7D%7D%20u_6(t-t_e-t_h-t_%7Bh%2Ci%7D)%20(%5Cvarphi_%7Bh%7D-%5Cvarphi_%7Bd%7D(u_4(t-t_%7Bh%2Ci%7D%2Bt%7Bh%2Cd%7D)))%20%20%20%20%20-%20k(t-t_e-t_h-t_%7Bh%2Cd%7D)%20%5Ccfrac%7Bu_7(t-t_e-t_h-t_%7Bh%2Cd%7D)%7D%20%7Bn_%7Bmax%7D%7D%20u_6(t-t_e-t_h-t_%7Bh%2Cd%7D)%20%5Cvarphi_%7Bd%7D(u_4(t)))

**Recovered:**

![formula](https://render.githubusercontent.com/render/math?math=%5Ccfrac%7Bdu_3(t)%7D%20%7Bdt%7D%20%20%3D%20%20%20%20%20%2B%20k(t-t_e-t_i)%20%5Ccfrac%7Bu_7(t-t_e-t_i)%7D%20%7Bn_%7Bmax%7D%7D%20u_6(t-t_e-t_i)%20(1-%5Cvarphi_%7Bh%7D)%20%20%20%20%20%2B%20k(t-t_e-t_h-t_%7Bh%2Ci%7D)%20%5Ccfrac%7Bu_7(t-t_e-t_h-t_%7Bh%2Ci%7D)%7D%20%7Bn_%7Bmax%7D%7D%20u_6(t-t_e-t_h-t_%7Bh%2Ci%7D)%20(%5Cvarphi_%7Bh%7D-%5Cvarphi_%7Bd%7D(u_4(t-t_%7Bh%2Ci%7D%2Bt%7Bh%2Cd%7D))))

**Deaths:**

![formula](https://render.githubusercontent.com/render/math?math=%5Ccfrac%7Bdu_2(t)%7D%20%7Bdt%7D%20%20%3D%20%20%20%20%20%2B%20k(t-t_e-t_h-t_%7Bh%2Cd%7D)%20%5Ccfrac%7Bu_7(t-t_e-t_h-t_%7Bh%2Cd%7D)%7D%20%7Bn_%7Bmax%7D%7D%20u_6(t-t_e-t_h-t_%7Bh%2Cd%7D)%20%5Cvarphi_%7Bd%7D(u_4(t)))

**Confirmed:**

![formula](https://render.githubusercontent.com/render/math?math=%5Ccfrac%7Bdu_1(t)%7D%20%7Bdt%7D%20%20%3D%20%20%20%20%20%2B%20%5Ccfrac%7Bdu_2(t)%7D%20%7Bdt%7D%20%20%20%20%20%2B%20%5Ccfrac%7Bdu_3(t)%7D%20%7Bdt%7D%20%20%20%20%20%2B%20%5Ccfrac%7Bdu_4(t)%7D%20%7Bdt%7D%20%20%20%20%20%2B%20%5Ccfrac%7Bdu_5(t)%7D%20%7Bdt%7D)

## Results

### Current Situation

![Situation](Situation-1.png)

The growths rate is determined by logarithmic fitting of the confirmed cases:

![Situation - logarithmic scale](Situation-2.png)

Basic Reproductive Number R0 has been determined from growths rate and median incubation time:

![formula](https://render.githubusercontent.com/render/math?math=R_0%20%3D%20e%5E%7Bk%20%20t_e%7D)

![Situation - Basic Reproductive Number](Situation-3.png)

![Death Casses](Situation-5.png)

Hidden cases are calculated by the amount of death casses per day.

![Hidden Casses](Situation-6.png)

### Model vs Situation

![Model vs Situation](Model_vs_Situation-1.png)

The predicted begin of a new phase is linked to a real event!

|   Date   |         Event         |        Phase        |       Estimate       |
|----------|-----------------------|---------------------|----------------------|
| 03/01/20 | begin of simulation   | uncontrolled growth |                    0 |
| 03/12/20 | who declares pandemic | social distancing   | 11.713200 (03/12/20) |
| 03/22/20 | contact restriction   | shutdown            | 21.946100 (03/22/20) |

### Forecast (250 Days)

![Forecast as Fraction](Forecast-1.png)

![Forecast as Fraction](Forecast-2.png)

![Forecast ARDS](Forecast-ARDS-2.png)

### Scenario 1 - no social distancing and no lock down

![Forecast as Fraction](Forecast-2-scenario-1.png)

![Forecast ARDS](Forecast-ARDS-2-scenario-1.png)

The day of health system crash is close to the timepiont the shutdown if an $R_0$ of 1 after shutdown is assumed.

![Forecast as Fraction](Forecast-ARDS-3.png)

![Forecast as Fraction](Forecast-ARDS-4.png)

### Scenario 2 - exit lock down after day 50 without social distancing (worest case)

![Forecast as Fraction](Forecast-2-scenario-2.png)

![Forecast ARDS](Forecast-ARDS-2-scenario-2.png)

### Scenario 3 - exit lock down after day 50 with social distancing (better case)

![Forecast as Fraction](Forecast-2-scenario-3.png)

![Forecast ARDS](Forecast-ARDS-2-scenario-3.png)

## References

* <https://www.rki.de/DE/Content/InfAZ/N/Neuartiges_Coronavirus/Steckbrief.html>
* ...
