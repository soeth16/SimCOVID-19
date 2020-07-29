# SimCOVID-19

Simulation of the SARS-CoV-2 pandemic virus outbreak which cause the corona virus disease 2019 (COVID-19) by using a modified S. E. I. R. D. model.
**This repository is currently in a preview state! It should be used with caution!** Currently only German data has been processed, other countries will follow. Further questions can answered by <soeren.thiering@hs-anhalt.de>

## Data

Thanks for sharing the data by Johns Hopkins CSSE (<https://systems.jhu.edu/research/public-health/ncov/>).

## Model

![Model](flowchart.svg)

**Parameter:**

| Name                                                                                                            |         Value | Description                                | Source |
|-----------------------------------------------------------------------------------------------------------------|--------------:|--------------------------------------------|--------|
| ![$k(t)$](https://render.githubusercontent.com/render/math?math=%24k(t)%24)                                     | ~ 0.15 to 0.4 | growth rate                                |        |
| ![$n_{max}$](https://render.githubusercontent.com/render/math?math=%24n_%7Bmax%7D%24)                           |      83019213 | max inahbitants                            |        |
| ![$t_{e}$](https://render.githubusercontent.com/render/math?math=%24t_%7Be%7D%24)                               |           5.1 | exposed / incubation time                  |        |
| ![$t_{i}$](https://render.githubusercontent.com/render/math?math=%24t_%7Bi%7D%24)                               |            13 | time of a moderate  infection              |        |
| ![$t_{h}$](https://render.githubusercontent.com/render/math?math=%24t_%7Bh%7D%24)                               |             8 | time until hostspitalisation               |        |
| ![$t_{h,i}$](https://render.githubusercontent.com/render/math?math=%24t_%7Bh%2Ci%7D%24)                         |            10 | time for hostspitalisation with ARDS       |        |
| ![$t_{h,ii}$](https://render.githubusercontent.com/render/math?math=%24t_%7Bh%2Cii%7D%24)                       |            10 | time after ARDS                            |        |
| ![$t_{h,d}$](https://render.githubusercontent.com/render/math?math=%24t_%7Bh%2Cd%7D%24)                         |           2.5 | time until death while ARDS                |        |
| ![$n_{h,max}$](https://render.githubusercontent.com/render/math?math=%24n_%7Bh%2Cmax%7D%24)                     |         28031 | max intensive beds for treating ARDS       |        |
| ![$\varphi_{h}$](https://render.githubusercontent.com/render/math?math=%24%5Cvarphi_%7Bh%7D%24)                 |          0.14 | propotion hostspitalisation with ARDS      |        |
| ![$\varphi_{d}(u_h(t))$](https://render.githubusercontent.com/render/math?math=%24%5Cvarphi_%7Bd%7D(u_h(t))%24) |        ~ 0.05 | propotion deaths                           |        |
| ![$\varphi_{d,min}$](https://render.githubusercontent.com/render/math?math=%24%5Cvarphi_%7Bd%2Cmin%7D%24)       |    0.02282158 | propotion deaths with hostspitalisation    |        |
| ![$\varphi_{d,max}$](https://render.githubusercontent.com/render/math?math=%24%5Cvarphi_%7Bd%2Cmax%7D%24)       |    0.09687916 | propotion deaths without hostspitalisation |        |

**Growth Rate:**

![Growth Rate](https://render.githubusercontent.com/render/math?math=k(t)%20%3D%5Cbegin%7Bcases%7Dk_%7Bn%2B1%7D%26%20%5Cquad%20%5Ctext%7B%20if%20%7D%20t%20%3E%20t_%7Bk%2Cn%2B1%7D%20%5C%5Ck_%7Bn%7D%26%20%5Cquad%20%5Ctext%7B%20if%20%7D%20t%20%3E%20t_%7Bk%2Cn%7D%20%5Cwedge%20t%20%5Cleq%20t_%7Bk%2Cn%2B1%7D%20%5C%5C%5Ccdots%26%20%5Cquad%20%5C%5Ck_%7B1%7D%26%20%5Cquad%20%5Ctext%7B%20if%20%7D%20t%20%3E%20t_%7Bk%2C1%7D%20%5Cwedge%20t%20%5Cleq%20t_%7Bk%2C2%7D%20%5C%5C%20%5Cend%7Bcases%7D%20%5C%5C)

**Deaths:**

![Deaths](https://render.githubusercontent.com/render/math?math=%5Cvarphi_%7Bd%7D(H(t))%20%3D%20%20%20%5Cbegin%7Bcases%7D%20%20%20%20%20%5Cvarphi_%7Bd%2Cmin%7D%20%5Ccfrac%20%7B1-H(t)%7D%20%7Bn_%7Bh%2Cmax%7D%7D%20%2B%20%5Cvarphi_%7Bd%2Cmax%7D%20%5Ccfrac%20%7BH(t)%7D%20%7Bn_%7Bh%2Cmax%7D%7D%20%20%26%20%5Cquad%20%5Ctext%7B%20if%20%7D%20%5Ccfrac%20%7BH(t)%7D%20%7Bn_%7Bh%2Cmax%7D%7D%20%3C%201%20%5C%5C%20%20%20%20%20%5Cvarphi_%7Bd%2Cmax%7D%20%26%20%5Cquad%20%5Ctext%7B%20if%20%7D%20%5Ccfrac%20%7BH(t)%7D%20%7Bn_%7Bh%2Cmax%7D%7D%20%5Cgeq%201%20%5C%5C%20%20%20%5Cend%7Bcases%7D)

**Susceptibles:**

![Susceptibles](https://render.githubusercontent.com/render/math?math=S%27(t)%20%3D%20-%20k(t)%20%5Ccfrac%7BS(t)%7D%20%7Bn_%7Bmax%7D%7D%20S(t))

**Exposed / Incubating:**

![Exposed](https://render.githubusercontent.com/render/math?math=E%27(t)%20%3D%20-%20S%27(t)%20%2B%20S%27(t-t_e))

**Infected:**

![Infected](https://render.githubusercontent.com/render/math?math=I%27(t)%20%3D%20-%20S%27(t%20-%20t_%7Be%7D)%20%2B%20S%27(t%20-%20t_%7Be%7D%20-%20t_%7Bi%7D)%20(1-%5Cvarphi_%7Bh%7D)%20%2B%20S%27(t%20-%20t_%7Be%7D%20-%20t_%7Bh%7D)%20%5Cvarphi_%7Bh%7D%20-%20S%27(t%20-%20t_%7Be%7D%20-%20t_%7Bh%7D%20-%20t_%7Bh%2Ci%7D)%20(%5Cvarphi_%7Bh%7D%20-%20%5Cvarphi_%7Bd%7D(H(t%20-%20t_%7Bh%2Ci%7D%20%2B%20t_%7Bh%2Cd%7D)))%20%2B%20S%27(t%20-%20t_%7Be%7D%20-%20t_%7Bh%7D%20-%20t_%7Bh%2Ci%7D%20-%20t_%7Bh%2Cii%7D)(%5Cvarphi_%7Bh%7D%20-%20%5Cvarphi_%7Bd%7D(H(t%20-%20t_%7Bh%2Ci%7D%20%2B%20t_%7Bh%2Cd%7D%20-%20t_%7Bh%2Cii%7D))))

**Hostspitalisation:**

![Hostspitalisation](https://render.githubusercontent.com/render/math?math=H%27(t)%20%3D%20-%20S%27(%20t%20-%20t_%7Be%7D%20-%20t_%7Bh%7D)%20%20%5Cvarphi_%7Bh%7D%20%20%2B%20S%27(%20t%20-%20t_%7Be%7D%20-%20t_%7Bh%7D%20-%20t_%7Bh%2Ci%7D)%20(%5Cvarphi_%7Bh%7D%20-%20%5Cvarphi_%7Bd%7D(H(%20t%20-%20t_%7Bh%2Ci%7D%20%2B%20t_%7Bh%2Cd%7D)))%20%2B%20S%27(%20t%20-%20t_%7Be%7D%20-%20t_%7Bh%7D%20-%20t_%7Bh%2Cd%7D)%20%5Cvarphi_%7Bd%7D(H(t)))

**Recovered:**

![Recovered](https://render.githubusercontent.com/render/math?math=R%27(t)%20%3D%20-%20S%27(%20t%20-%20t_%7Be%7D%20-%20t_%7Bi%7D)%20%20(1-%5Cvarphi_%7Bh%7D)%20-%20S%27(%20t%20-%20t_%7Be%7D%20-%20t_%7Bh%7D%20-%20t_%7Bh%2Ci%7D%20-%20t_%7Bh%2Cii%7D)%20(%5Cvarphi_%7Bh%7D%20-%20%5Cvarphi_%7Bd%7D(H(%20t%20-%20t_%7Bh%2Ci%7D%20%2B%20t_%7Bh%2Cd%7D%20-%20t_%7Bh%2Cii%7D))))

**Deaths:**

![Deaths](https://render.githubusercontent.com/render/math?math=D%27(t)%20%3D%20-%20S%27(%20t%20-%20t_%7Be%7D%20-%20t_%7Bh%7D%20-%20t_%7Bh%2Cd%7D)%20%20%5Cvarphi_%7Bd%7D(H(t)))

**Confirmed:**

![Confirmed](https://render.githubusercontent.com/render/math?math=C%27(t)%20%3D%20%2B%20D%27(t)%20%2B%20R%27(t)%20%2B%20H%27(t)%20%2B%20I%27(t))

## Results

### Current Situation

![Situation](Situation-1.png)

The growths rate is determined by logarithmic fitting of the confirmed cases:

![Situation - logarithmic scale](Situation-2.png)

Basic Reproductive Number R0 has been determined by the nowcast method:

![Situation - Basic Reproductive Number](Situation-3.png)

![Death Casses](Situation-5.png)

Hidden cases are calculated by the amount of death casses per day.

![Hidden Casses](Situation-6.png)

### Model vs Situation

![Model vs Situation](Model_vs_Situation-1.png)

The predicted begin of a new phase is linked to a real event!

|   date   |                   event                |        phase        |     estimated date     | estimated growth rate |
|----------|----------------------------------------|---------------------|------------------------|-----------------------|
| 03/01/20 | begin of simulation                    | uncontrolled growth |                      0 |             0.3003057 |
| 03/12/20 | who declares pandemic                  | social distancing   |  11.6507903 (03/12/20) |             0.2818724 |
| 03/22/20 | contact restriction                    | shutdown            |  21.8597928 (03/22/20) |             0.1863497 |
| 04/22/20 | mask requirement in all federal states |                     |  53.7949697 (04/23/20) |             0.1765100 |
| 04/30/20 | end of shutdown                        | re open             |  61.2062883 (05/01/20) |             0.1691815 |
| 06/18/20 | endemic outbreak in Gütersloh          | begin               | 101.7158569 (06/10/20) |             0.1992506 |
|          |                                        | end                 | 110.9930793 (06/19/20) |             0.1647673 |
| 07/29/20 | 2nd phase                              | begin               | 131.2207343 (07/10/20) |             0.1743329 |

### Forecast (250 Days)

![Forecast as Fraction](Forecast-1.png)

![Forecast as Fraction](Forecast-2.png)

![Forecast ARDS](Forecast-ARDS-2.png)

### Scenario 1 - no social distancing and no lock down

![Forecast as Fraction](Forecast-2-scenario-1.png)

![Forecast ARDS](Forecast-ARDS-2-scenario-1.png)

The day of health system crash is close to the time point the shutdown if an ![$R_0$](https://render.githubusercontent.com/render/math?math=%24R_0%24) of 1 after shutdown is assumed.

![Forecast as Fraction](Forecast-ARDS-3.png)

![Forecast as Fraction](Forecast-ARDS-4.png)

### Scenario 2 - End of lock down after day 50 without social distancing (worst case)

![Forecast as Fraction](Forecast-2-scenario-2.png)

![Forecast ARDS](Forecast-ARDS-2-scenario-2.png)

### Scenario 3 - End of lock down after day 50 with social distancing (better case)

![Forecast as Fraction](Forecast-2-scenario-3.png)

![Forecast ARDS](Forecast-ARDS-2-scenario-3.png)

## References

* <https://www.rki.de/DE/Content/InfAZ/N/Neuartiges_Coronavirus/Steckbrief.html>
* <https://www.rki.de/DE/Content/InfAZ/N/Neuartiges_Coronavirus/Projekte_RKI/R-Wert-Erlaeuterung.pdf?__blob=publicationFile>
* <https://www.recoverytrial.net/files/recovery_dexamethasone_statement_160620_v2final.pdf>
* ...
