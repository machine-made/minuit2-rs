iminuit: 2.33.0.
Bundled Minuit2: ROOT 6.37.01 development snapshot, commit `56be091a18cc07bc12242d65447d63611b4ebda4`.
ROOT reference tag: `v6-36-08`.
Rust base commit: `e50460b`. The caller commits these harness files later.

# Matched accuracy

The CSV contains 90 rows in a fixed order. Values and gradients are counted inside each user function.
Each fit runs Migrad once, then Hesse. Explicit Hesse calls are the difference between the two counter readings.
iminuit uses one Migrad iteration, with wrapper retries and Simplex disabled. Both sides use the same call budgets.
Value calls determine call ratios. Gradient calls are shown separately. Failed fits retain their measurements and a message.
True distance uses the objective evaluated at the returned parameters with 40-digit arithmetic, minus its reference minimum.
Covariance error is the maximum absolute entry error divided by the geometric mean of the corresponding reference diagonal entries.
Reference covariance is twice the error definition times the inverse reference Hessian on free parameters.
The NLL error definition is 0.5. Other cases use 1.0. Fixed parameters are removed before covariance comparisons.
A ratio with two zero errors is one. A positive error divided by zero is infinite. Empty medians say not reached.
Verdicts include rows where both true distances are at most 1e-6. Each median uses those paired rows.
A group says difference if its call median leaves [0.9, 1.1], its covariance median leaves [0.5, 2], or bounds success differs.
Otherwise it says no difference. These measurements describe differences; they do not establish a general speed advantage.

Gaussian data have seeded noise. Poisson bins contain expected counts, with totals of 10000 and 1000000 events.
The seed is 20261113. All reference minima and Hessians use mpmath at 40 decimal digits.
The one-sided bound distance scale is five units. Double limits also span five units.
Bounds success means every free parameter is within 0.01 of its constrained reference minimum, regardless of the validity flag.

Version sources: [iminuit source tree](https://github.com/scikit-hep/iminuit/tree/v2.33.0/extern) and [bundled ROOT header](https://github.com/root-project/root/blob/56be091a18cc07bc12242d65447d63611b4ebda4/core/foundation/inc/ROOT/RVersion.hxx).

## Hessian

Only non-quadratic cases enter this group's verdict: Rosenbrock, the Gaussian fit, and both Poisson fits.
Strategy 3 rows are iminuit only. Rust strategy 3 is not implemented.
The pinned iminuit setter rejects strategy 3. These rows use its bundled MnStrategy constructor and check the resulting level.

Each paired cell gives Rust / iminuit. Full parameters, covariance matrices, and counter readings are in the CSV.

| case | strategy | tolerance | status | value calls | gradient calls | Hesse calls | EDM | true distance | covariance error | success |
|---|---|---|---|---|---|---|---|---|---|---|
| dense2 | 0 | 0.1 | ok / ok | 22 / 23 | 0 / 0 | 5 / 6 | 1.7513277221717818e-26 / 3.667496490707142e-27 | 1.7513282571514718e-26 / 3.667495173723496e-27 | 6.716100729279464e-14 / 3.964678522351168e-14 |  |
| dense2 | 1 | 0.1 | ok / ok | 29 / 39 | 0 / 0 | 5 / 10 | 4.531133717370879e-27 / 8.29233772581379e-28 | 4.531136014458611e-27 / 8.292337729525389e-28 | 7.26653671576597e-14 / 1.9189366279121915e-13 |  |
| dense2 | 2 | 0.1 | ok / ok | 35 / 42 | 0 / 0 | 5 / 10 | 1.407454594746429e-19 / 6.317327676963028e-20 | 1.4074545947321644e-19 / 6.31732767698661e-20 | 1.0871130706975044e-14 / 7.738487123110742e-14 |  |
| dense2 | 3 | 0.1 | not implemented / ok | not reached / 93 | not reached / 0 | not reached / 37 | not reached / 8.292333621183789e-28 | not reached / 8.292337729525389e-28 | not reached / 4.051898132228871e-16 |  |
| dense5 | 0 | 0.1 | ok / ok | 77 / 78 | 0 / 0 | 20 / 21 | 1.7365063300542297e-06 / 1.7365063726930433e-06 | 1.7365063300542305e-06 / 1.7365063726930428e-06 | 8.375888481060834e-14 / 4.9368968104543165e-14 |  |
| dense5 | 1 | 0.1 | ok / ok | 145 / 127 | 0 / 0 | 40 / 31 | 1.7365063674575203e-06 / 1.7365063576959599e-06 | 1.7365063674574528e-06 / 1.7365063576942376e-06 | 7.093459850155435e-14 / 4.810524828812577e-13 |  |
| dense5 | 2 | 0.1 | ok / ok | 112 / 125 | 0 / 0 | 20 / 31 | 3.0068213432632512e-18 / 2.452793893800051e-18 | 3.00682134318678e-18 / 2.4527938938019156e-18 | 5.781716291186023e-14 / 1.878577760651111e-13 |  |
| dense5 | 3 | 0.1 | not implemented / ok | not reached / 309 | not reached / 0 | not reached / 121 | not reached / 1.7365063590483826e-06 | not reached / 1.7365063590482852e-06 | not reached / 1.2325237018922224e-15 |  |
| dense10 | 0 | 0.1 | ok / ok | 216 / 217 | 0 / 0 | 65 / 66 | 8.066089323670796e-05 / 8.06608936688021e-05 | 8.066089323669296e-05 / 8.066089366881502e-05 | 2.228763460228227e-13 / 1.0308191466720878e-13 |  |
| dense10 | 1 | 0.1 | ok / ok | 361 / 323 | 0 / 0 | 105 / 86 | 8.06608932366993e-05 / 8.066089366871302e-05 | 8.066089323669296e-05 / 8.066089366881502e-05 | 3.447513923478433e-13 / 4.046065335597054e-12 |  |
| dense10 | 2 | 0.1 | ok / ok | 297 / 320 | 0 / 0 | 65 / 86 | 2.9400676918224244e-16 / 8.330319260156563e-17 | 2.9400676961033753e-16 / 8.330319260164856e-17 | 8.08485676759831e-14 / 3.94278694496564e-13 |  |
| dense10 | 3 | 0.1 | not implemented / ok | not reached / 833 | not reached / 0 | not reached / 341 | not reached / 8.066089366881418e-05 | not reached / 8.066089366881502e-05 | not reached / 3.883856400557528e-15 |  |
| dense20 | 0 | 0.1 | ok / ok | 602 / 603 | 0 / 0 | 230 / 231 | 2.00390225794735e-05 / 2.0039022498263553e-05 | 2.003902257947239e-05 / 2.0039022498263668e-05 | 9.960257228485692e-14 / 1.0350108548618187e-13 |  |
| dense20 | 1 | 0.1 | ok / ok | 992 / 914 | 0 / 0 | 310 / 271 | 2.0039022579474032e-05 / 2.00390224982444e-05 | 2.003902257947239e-05 / 2.0039022498263668e-05 | 1.3473718450418407e-13 / 1.0078273252829128e-12 |  |
| dense20 | 2 | 0.1 | ok / ok | 892 / 935 | 0 / 0 | 230 / 271 | 4.474557411599843e-16 / 3.536365092044338e-16 | 4.474557408873442e-16 / 3.5363650920448984e-16 | 1.1072055362801121e-13 / 5.648582817026248e-13 |  |
| dense20 | 3 | 0.1 | not implemented / ok | not reached / 2534 | not reached / 0 | not reached / 1081 | not reached / 2.0039022498263912e-05 | not reached / 2.0039022498263668e-05 | not reached / 4.780909489023107e-15 |  |
| rosenbrock | 0 | 0.1 | ok / ok | 124 / 125 | 0 / 0 | 5 / 6 | 1.3190530708951375e-05 / 1.3190530741822684e-05 | 1.3231417635634005e-05 / 1.3231417668691429e-05 | 0.010286249132200154 / 0.010286249150155977 |  |
| rosenbrock | 1 | 0.1 | ok / ok | 156 / 150 | 0 / 0 | 13 / 10 | 3.451792444266796e-05 / 3.4254032651655184e-05 | 3.459740883846543e-05 / 3.459747545738712e-05 | 0.013912781408871875 / 0.02153431006603148 |  |
| rosenbrock | 2 | 0.1 | ok / ok | 181 / 172 | 0 / 0 | 13 / 10 | 2.9803328483419318e-06 / 2.958617302328516e-06 | 2.9758239687521107e-06 / 2.9758226451137953e-06 | 0.0016990556998048706 / 0.009445868056594878 |  |
| rosenbrock | 3 | 0.1 | not implemented / ok | not reached / 210 | not reached / 0 | not reached / 37 | not reached / 3.418967947583445e-05 | not reached / 3.460139011873646e-05 | not reached / 0.023474774516818866 |  |
| quadratic2 | 0 | 0.1 | ok / ok | 23 / 24 | 0 / 0 | 5 / 6 | 3.30648779426829e-05 / 3.3064878495403884e-05 | 3.306487793905337e-05 / 3.3064878489688526e-05 | 2.1929555005004563e-10 / 2.0889667676451542e-10 |  |
| quadratic2 | 1 | 0.1 | ok / ok | 44 / 38 | 0 / 0 | 13 / 10 | 3.306487793494943e-05 / 3.3064878609280166e-05 | 3.306487793905337e-05 / 3.3064878489688526e-05 | 8.527275030725434e-11 / 5.627136897142916e-09 |  |
| quadratic2 | 2 | 0.1 | ok / ok | 37 / 44 | 0 / 0 | 5 / 10 | 3.573662517193599e-19 / 3.779548062334506e-18 | 3.5740806649741747e-19 / 3.782315205600814e-18 | 4.943933414071304e-10 / 5.823627321074467e-09 |  |
| quadratic2 | 3 | 0.1 | not implemented / ok | not reached / 94 | not reached / 0 | not reached / 37 | not reached / 3.306487786113649e-05 | not reached / 3.3064877860996814e-05 | not reached / 3.2122624240716587e-13 |  |
| quadratic2_fixed | 0 | 0.1 | ok / ok | 10 / 11 | 0 / 0 | 2 / 3 | 9.79678043496136e-22 / 4.535465274137128e-18 | 9.796929768705238e-22 / 4.53547927040122e-18 | 5.840594674566546e-11 / 7.317835226672287e-11 |  |
| quadratic2_fixed | 1 | 0.1 | ok / ok | 12 / 18 | 0 / 0 | 2 / 5 | 9.79678043496136e-22 / 4.534893118657347e-18 | 9.796929768705238e-22 / 4.53547927040122e-18 | 5.840594674566546e-11 / 3.8525560519531155e-10 |  |
| quadratic2_fixed | 2 | 0.1 | ok / ok | 18 / 23 | 0 / 0 | 2 / 5 | 9.99964194926884e-22 / 4.533462885500546e-18 | 1.0007101230030648e-21 / 4.532944982629981e-18 | 1.0555734064610078e-10 / 1.6855787876579598e-09 |  |
| quadratic2_fixed | 3 | 0.1 | not implemented / ok | not reached / 42 | not reached / 0 | not reached / 17 | not reached / 4.535612797099679e-18 | not reached / 4.53547927040122e-18 | not reached / 2.731148640577885e-14 |  |
| gaussian200 | 0 | 0.1 | ok / ok | 62 / 63 | 0 / 0 | 9 / 10 | 1.073358590552285e-05 / 1.0733585687786578e-05 | 1.0733384090237475e-05 / 1.0733383971284701e-05 | 3.656073654435664e-05 / 3.656648570725632e-05 |  |
| gaussian200 | 1 | 0.1 | ok / ok | 95 / 85 | 0 / 0 | 21 / 16 | 1.0733585840284747e-05 / 1.0733514867770894e-05 | 1.0733384090237475e-05 / 1.0733383971284701e-05 | 3.6568560979010436e-05 / 2.879602064744411e-05 |  |
| gaussian200 | 2 | 0.1 | ok / ok | 112 / 97 | 0 / 0 | 21 / 16 | 1.4965049796474662e-07 / 1.4965109267078875e-07 | 1.4965137746912974e-07 / 1.496512763053865e-07 | 1.4365926577313213e-05 / 5.392405805399705e-06 |  |
| gaussian200 | 3 | 0.1 | not implemented / ok | not reached / 175 | not reached / 0 | not reached / 61 | not reached / 1.0733493256965485e-05 | not reached / 1.0733383971284701e-05 | not reached / 2.6532240922020072e-05 |  |
| poisson10000 | 0 | 0.1 | ok / ok | 124 / 125 | 0 / 0 | 20 / 21 | 4.093436876032087e-06 / 4.093435832389763e-06 | 4.093398290471228e-06 / 4.093397247091046e-06 | 0.0001101092267172601 / 0.00011010892400135987 |  |
| poisson10000 | 1 | 0.1 | ok / ok | 184 / 166 | 0 / 0 | 40 / 31 | 4.093436873077739e-06 / 4.093438735212173e-06 | 4.093398290471228e-06 / 4.093397247091046e-06 | 0.00011010727787042258 / 0.00010605401180551369 |  |
| poisson10000 | 2 | 0.1 | ok / ok | 209 / 182 | 0 / 0 | 40 / 31 | 3.1272719957251886e-07 / 3.1272724633218794e-07 | 3.1272750681749784e-07 / 3.127277280081266e-07 | 1.9482162214167433e-05 / 2.229123247492859e-05 |  |
| poisson10000 | 3 | 0.1 | not implemented / ok | not reached / 348 | not reached / 0 | not reached / 121 | not reached / 4.093438993160723e-06 | not reached / 4.093397236918288e-06 | not reached / 0.00010502898906350494 |  |
| poisson1000000 | 0 | 0.1 | ok / ok | 133 / 134 | 0 / 0 | 20 / 21 | 1.891596878681953e-07 / 1.8915964059224766e-07 | 1.891598039160504e-07 / 1.891597566632023e-07 | 1.5393199548439396e-06 / 1.5393499094524938e-06 |  |
| poisson1000000 | 1 | 0.1 | ok / ok | 195 / 177 | 0 / 0 | 40 / 31 | 1.8915967325885614e-07 / 1.8915967949341896e-07 | 1.8915978923504715e-07 / 1.8915974455612928e-07 | 1.5392756376290368e-06 / 1.1910812004522249e-06 |  |
| poisson1000000 | 2 | 0.1 | ok / ok | 203 / 176 | 0 / 0 | 40 / 31 | 3.070654108013947e-05 / 3.070655874903748e-05 | 3.070658135879605e-05 / 3.070659481302611e-05 | 2.3434584563007317e-05 / 2.357689921318864e-05 |  |
| poisson1000000 | 3 | 0.1 | not implemented / ok | not reached / 357 | not reached / 0 | not reached / 121 | not reached / 1.891596800443235e-07 | not reached / 1.8915974455612928e-07 | not reached / 1.076367682781607e-06 |  |

difference; median Rust/iminuit calls = 1.2903225806451613; median Rust/iminuit covariance errors = 1.1461576504973254; paired rows = 4.

## EDM

Each paired cell gives Rust / iminuit. Full parameters, covariance matrices, and counter readings are in the CSV.

| case | strategy | tolerance | status | value calls | gradient calls | Hesse calls | EDM | true distance | covariance error | success |
|---|---|---|---|---|---|---|---|---|---|---|
| rosenbrock | 1 | 1.0 | ok / ok | 130 / 127 | 0 / 0 | 13 / 10 | 0.0016467319157285176 / 0.0016467384678486574 | 0.0021443287733624726 / 0.002144329477598134 | 0.38260045087028266 / 0.38260048569330696 |  |
| rosenbrock | 1 | 0.1 | ok / ok | 143 / 140 | 0 / 0 | 13 / 10 | 3.4517099511865e-05 / 3.4518507168690824e-05 | 3.459740883846543e-05 / 3.459747545738712e-05 | 0.013936476288027012 / 0.013936475194053222 |  |
| rosenbrock | 1 | 0.01 | ok / ok | 149 / 146 | 0 / 0 | 13 / 10 | 1.6868090457168244e-08 / 1.6879299363496475e-08 | 1.685947641451356e-08 / 1.685962537937131e-08 | 0.007453395688580944 / 0.0074533820926370176 |  |
| rosenbrock | 1 | 0.001 | ok / ok | 149 / 146 | 0 / 0 | 13 / 10 | 1.6868090457168244e-08 / 1.6879299363496475e-08 | 1.685947641451356e-08 / 1.685962537937131e-08 | 0.007453395688580944 / 0.0074533820926370176 |  |
| rosenbrock | 1 | 0.0001 | ok / ok | 149 / 146 | 0 / 0 | 13 / 10 | 1.6868090457168244e-08 / 1.6879299363496475e-08 | 1.685947641451356e-08 / 1.685962537937131e-08 | 0.007453395688580944 / 0.0074533820926370176 |  |
| dense2 | 1 | 1.0 | ok / ok | 24 / 29 | 0 / 0 | 5 / 10 | 4.531133717370879e-27 / 8.29233652768636e-28 | 4.531136014458611e-27 / 8.292337729525389e-28 | 7.26653671576597e-14 / 4.334693380514508e-14 |  |
| dense2 | 1 | 0.1 | ok / ok | 24 / 29 | 0 / 0 | 5 / 10 | 4.531133717370879e-27 / 8.29233652768636e-28 | 4.531136014458611e-27 / 8.292337729525389e-28 | 7.26653671576597e-14 / 4.334693380514508e-14 |  |
| dense2 | 1 | 0.01 | ok / ok | 24 / 29 | 0 / 0 | 5 / 10 | 4.531133717370879e-27 / 8.29233652768636e-28 | 4.531136014458611e-27 / 8.292337729525389e-28 | 7.26653671576597e-14 / 4.334693380514508e-14 |  |
| dense2 | 1 | 0.001 | ok / ok | 24 / 29 | 0 / 0 | 5 / 10 | 4.531133717370879e-27 / 8.29233652768636e-28 | 4.531136014458611e-27 / 8.292337729525389e-28 | 7.26653671576597e-14 / 4.334693380514508e-14 |  |
| dense2 | 1 | 0.0001 | ok / ok | 24 / 29 | 0 / 0 | 5 / 10 | 4.531133717370879e-27 / 8.29233652768636e-28 | 4.531136014458611e-27 / 8.292337729525389e-28 | 7.26653671576597e-14 / 4.334693380514508e-14 |  |
| dense10 | 1 | 1.0 | ok / ok | 256 / 237 | 0 / 0 | 105 / 86 | 8.066089323670677e-05 / 8.066089366878947e-05 | 8.066089323669296e-05 / 8.066089366881502e-05 | 2.228763460228227e-13 / 1.0308191466720878e-13 |  |
| dense10 | 1 | 0.1 | ok / ok | 256 / 237 | 0 / 0 | 105 / 86 | 8.066089323670677e-05 / 8.066089366878947e-05 | 8.066089323669296e-05 / 8.066089366881502e-05 | 2.228763460228227e-13 / 1.0308191466720878e-13 |  |
| dense10 | 1 | 0.01 | ok / ok | 278 / 259 | 0 / 0 | 105 / 86 | 5.346706092037349e-06 / 5.346706118935724e-06 | 5.346706092037085e-06 / 5.346706118935966e-06 | 9.000977163155978e-14 / 7.858117239577392e-14 |  |
| dense10 | 1 | 0.001 | ok / ok | 300 / 281 | 0 / 0 | 105 / 86 | 4.2976421098451103e-07 / 4.297642127919101e-07 | 4.2976421098445836e-07 / 4.297642127918784e-07 | 1.2790103148961412e-13 / 6.653170773420468e-14 |  |
| dense10 | 1 | 0.0001 | ok / ok | 322 / 303 | 0 / 0 | 105 / 86 | 6.232796404457407e-09 / 6.2327965265840024e-09 | 6.232796404456802e-09 / 6.232796526583347e-09 | 1.162386950414538e-13 / 9.457874211165546e-14 |  |
| quadratic2_fixed | 1 | 1.0 | ok / ok | 10 / 13 | 0 / 0 | 2 / 5 | 9.79678043496136e-22 / 4.535465274137128e-18 | 9.796929768705238e-22 / 4.53547927040122e-18 | 5.840594674566546e-11 / 7.317835226672287e-11 |  |
| quadratic2_fixed | 1 | 0.1 | ok / ok | 10 / 13 | 0 / 0 | 2 / 5 | 9.79678043496136e-22 / 4.535465274137128e-18 | 9.796929768705238e-22 / 4.53547927040122e-18 | 5.840594674566546e-11 / 7.317835226672287e-11 |  |
| quadratic2_fixed | 1 | 0.01 | ok / ok | 10 / 13 | 0 / 0 | 2 / 5 | 9.79678043496136e-22 / 4.535465274137128e-18 | 9.796929768705238e-22 / 4.53547927040122e-18 | 5.840594674566546e-11 / 7.317835226672287e-11 |  |
| quadratic2_fixed | 1 | 0.001 | ok / ok | 10 / 13 | 0 / 0 | 2 / 5 | 9.79678043496136e-22 / 4.535465274137128e-18 | 9.796929768705238e-22 / 4.53547927040122e-18 | 5.840594674566546e-11 / 7.317835226672287e-11 |  |
| quadratic2_fixed | 1 | 0.0001 | ok / ok | 10 / 13 | 0 / 0 | 2 / 5 | 9.79678043496136e-22 / 4.535465274137128e-18 | 9.796929768705238e-22 / 4.53547927040122e-18 | 5.840594674566546e-11 / 7.317835226672287e-11 |  |
| poisson10000 | 1 | 1.0 | ok / ok | 144 / 135 | 0 / 0 | 40 / 31 | 4.093436873802133e-06 / 4.093436522704517e-06 | 4.093398290471228e-06 / 4.093397247091046e-06 | 0.0001101092267172601 / 0.00011010892400135987 |  |
| poisson10000 | 1 | 0.1 | ok / ok | 144 / 135 | 0 / 0 | 40 / 31 | 4.093436873802133e-06 / 4.093436522704517e-06 | 4.093398290471228e-06 / 4.093397247091046e-06 | 0.0001101092267172601 / 0.00011010892400135987 |  |
| poisson10000 | 1 | 0.01 | ok / ok | 144 / 135 | 0 / 0 | 40 / 31 | 4.093436873802133e-06 / 4.093436522704517e-06 | 4.093398290471228e-06 / 4.093397247091046e-06 | 0.0001101092267172601 / 0.00011010892400135987 |  |
| poisson10000 | 1 | 0.001 | ok / ok | 155 / 146 | 0 / 0 | 40 / 31 | 1.8926572654140052e-09 / 1.892673835271905e-09 | 1.892679292321105e-09 / 1.892679372990712e-09 | 9.553113484867428e-06 / 9.553077919736225e-06 |  |
| poisson10000 | 1 | 0.0001 | ok / ok | 155 / 146 | 0 / 0 | 40 / 31 | 1.8926572654140052e-09 / 1.892673835271905e-09 | 1.892679292321105e-09 / 1.892679372990712e-09 | 9.553113484867428e-06 / 9.553077919736225e-06 |  |

Matched-accuracy calls use the first tolerance in the listed order that reaches a true distance of at most 1e-6.

| case | Rust calls | Rust tolerance | iminuit calls | iminuit tolerance |
|---|---|---|---|---|
| rosenbrock | 149 | 0.01 | 146 | 0.01 |
| dense2 | 24 | 1.0 | 29 | 1.0 |
| dense10 | 300 | 0.001 | 281 | 0.001 |
| quadratic2_fixed | 10 | 1.0 | 13 | 1.0 |
| poisson10000 | 155 | 0.001 | 146 | 0.001 |

The verdict uses all paired tolerance rows, following the same row rule as the other groups.

difference; median Rust/iminuit calls = 0.8275862068965517; median Rust/iminuit covariance errors = 1.000003722897636; paired rows = 17.

## gradient

Each paired cell gives Rust / iminuit. Full parameters, covariance matrices, and counter readings are in the CSV.

| case | strategy | tolerance | status | value calls | gradient calls | Hesse calls | EDM | true distance | covariance error | success |
|---|---|---|---|---|---|---|---|---|---|---|
| quadratic2 | 1 | 0.0001 | ok / ok | 17 / 29 | 3 / 3 | 5 / 10 | 1.302320849383499e-26 / 8.842714103621333e-26 | 1.311140362111822e-26 / 6.985048507126285e-32 | 2.3266570423441596e-13 / 2.578002535082128e-10 |  |
| rosenbrock | 1 | 0.0001 | ok / ok | 60 / 74 | 13 / 17 | 13 / 10 | 8.401780790099119e-09 / 2.7301804845750188e-08 | 7.742229623451493e-09 / 2.721136126830823e-08 | 0.6620472038093936 / 0.0382382791292768 |  |
| quadratic2_no_g2_migrad | ROOT only | | ok | 2 | | | 0 | 0 | | |

The gradient rows use tolerance 0.0001 so both libraries can enter the matched-accuracy verdict.
The ROOT-only Hessian row is copied from the existing built runner, using its quadratic2_no_g2_migrad workload.
Its historical count is the runner's reported count. This row is excluded from the paired medians.

difference; median Rust/iminuit calls = 0.6985088536812675; median Rust/iminuit covariance errors = 8.657315772012263; paired rows = 2.

## bounds

Each paired cell gives Rust / iminuit. Full parameters, covariance matrices, and counter readings are in the CSV.

| case | strategy | tolerance | status | value calls | gradient calls | Hesse calls | EDM | true distance | covariance error | success |
|---|---|---|---|---|---|---|---|---|---|---|
| lower_0.0 | 1 | 0.1 | ok / ok | 23 / 37 | 0 / 0 | 6 / 5 | 5.6334968265256994e-08 / 2.1322293390822565e-08 | 5.632776899242308e-08 / 2.1322576910301374e-08 | 3.951853089867452e-05 / 2.4311552464628683e-05 | True / True |
| lower_1e-09 | 1 | 0.1 | ok / ok | 23 / 36 | 0 / 0 | 6 / 5 | 5.6334968265256994e-08 / 2.0982491110246334e-08 | 5.632776899242308e-08 / 2.0982764175167417e-08 | 3.951853089867452e-05 / 2.411696067960989e-05 | True / True |
| lower_1e-06 | 1 | 0.1 | ok / ok | 37 / 36 | 0 / 0 | 6 / 5 | 1.2731006261841464e-08 / 1.2726869282147737e-08 | 1.2726763433624815e-08 / 1.272692709306154e-08 | 1.877951530993016e-05 / 1.8779636436039304e-05 | True / True |
| lower_0.001 | 1 | 0.1 | ok / ok | 36 / 35 | 0 / 0 | 6 / 5 | 6.343457610672827e-11 / 6.313191801845887e-11 | 6.311857566119587e-11 / 6.311930011538789e-11 | 1.3057827601414829e-06 / 1.3057903992530484e-06 | True / True |
| lower_0.1 | 1 | 0.1 | ok / ok | 21 / 20 | 0 / 0 | 6 / 5 | 2.493285026692007e-12 / 2.43340724686149e-12 | 2.430936230886946e-12 / 2.4309363098205367e-12 | 2.41698460468065e-07 / 2.416977725738789e-07 | True / True |
| upper_0.0 | 1 | 0.1 | ok / ok | 24 / 38 | 0 / 0 | 6 / 5 | 5.2903352382987145e-08 / 1.2209161212502362e-06 | 5.290299685268299e-08 / 1.2209049588922654e-06 | 1.9160497322090464e-06 / 9.214633733156674e-06 | True / True |
| upper_1e-09 | 1 | 0.1 | ok / ok | 24 / 37 | 0 / 0 | 6 / 5 | 5.2903352382987145e-08 / 1.2183124211427053e-06 | 5.290299685268299e-08 / 1.2183012946595997e-06 | 1.9160497322090464e-06 / 9.204798417039584e-06 | True / True |
| upper_1e-06 | 1 | 0.1 | ok / ok | 38 / 37 | 0 / 0 | 6 / 5 | 1.1404651091802815e-06 / 1.1404667497224987e-06 | 1.1404570762899714e-06 / 1.1404566780743503e-06 | 8.905685824345255e-06 / 8.905683086313232e-06 | True / True |
| upper_0.001 | 1 | 0.1 | ok / ok | 39 / 38 | 0 / 0 | 6 / 5 | 8.162039186911339e-07 / 8.162046255119097e-07 | 8.16199566341872e-07 / 8.16198548634731e-07 | 7.533295085071501e-06 / 7.5332900033586725e-06 | True / True |
| upper_0.1 | 1 | 0.1 | ok / ok | 28 / 27 | 0 / 0 | 6 / 5 | 7.129505416667621e-07 / 7.129559772998012e-07 | 7.129472009697846e-07 / 7.129510248705801e-07 | 7.040459487184236e-06 / 7.0404765728504515e-06 | True / True |
| both_0.0 | 1 | 0.1 | ok / ok | 39 / 38 | 0 / 0 | 6 / 5 | 4.332679321507613e-10 / 4.3148743257053567e-10 | 4.3141235660947733e-10 / 4.3141599479770496e-10 | 7.777070577308365e-06 / 7.777103573025634e-06 | True / True |
| both_1e-09 | 1 | 0.1 | ok / ok | 39 / 38 | 0 / 0 | 6 / 5 | 4.332679321507613e-10 / 4.3148743257053567e-10 | 4.3141235660947733e-10 / 4.3141599479770496e-10 | 7.777070577308365e-06 / 7.777103573025634e-06 | True / True |
| both_1e-06 | 1 | 0.1 | ok / ok | 39 / 38 | 0 / 0 | 6 / 5 | 4.332679321507613e-10 / 4.2635973612642615e-10 | 4.3141235660947733e-10 / 4.26288685950662e-10 | 7.777070577308365e-06 / 7.730682867190808e-06 | True / True |
| both_0.001 | 1 | 0.1 | ok / ok | 36 / 35 | 0 / 0 | 6 / 5 | 4.295525961805282e-06 / 4.295496538693033e-06 | 4.298671273431668e-06 / 4.298818558078236e-06 | 0.000774604285964875 / 0.0007746175039264891 | True / True |
| both_0.1 | 1 | 0.1 | ok / ok | 21 / 20 | 0 / 0 | 6 / 5 | 4.4206596184168756e-08 / 4.422456372310871e-08 | 4.422191586299312e-08 / 4.42219154939722e-08 | 7.890133743804206e-05 / 7.890133662602494e-05 | True / True |
| at_limit_0.0 | 1 | 0.1 | ok / ok | 18 / 10 | 0 / 0 | 2 / 5 | 0.0 / 0.0 | 0.0 / 0.0 | 1.0 / 1.0 | True / True |
| at_limit_1e-09 | 1 | 0.1 | ok / ok | 18 / 13 | 0 / 0 | 2 / 5 | 0.0 / 6.462348731172794e-25 | 0.0 / 0.0 | 1.0 / 1.0 | True / True |
| at_limit_1e-06 | 1 | 0.1 | ok / ok | 10 / 13 | 0 / 0 | 2 / 5 | 2.0938009888999857e-24 / 6.46234872816353e-25 | 0.0 / 0.0 | 1.0 / 1.0 | True / True |
| at_limit_0.001 | 1 | 0.1 | ok / ok | 10 / 13 | 0 / 0 | 2 / 5 | 2.2840525472452567e-22 / 1.6808569151532745e-21 | 0.0 / 0.0 | 1.0 / 1.0 | True / True |
| at_limit_0.1 | 1 | 0.1 | ok / ok | 10 / 13 | 0 / 0 | 2 / 5 | 3.697539055277227e-20 / 1.8635411653519656e-19 | 0.0 / 2.168404344971009e-19 | 1.0 / 1.0 | True / True |
| quadratic2_y-1.0 | 1 | 0.1 | ok / ok | 77 / 74 | 0 / 0 | 13 / 10 | 6.532518807247125e-11 / 6.61034194364606e-11 | 6.614520426434177e-11 / 6.613663475088805e-11 | 0.15000893514512476 / 8.93506752003307e-06 | True / True |
| quadratic2_y-1.2 | 1 | 0.1 | ok / ok | 63 / 59 | 0 / 0 | 13 / 10 | 1.448496047250657e-05 / 1.4485377126665103e-05 | 1.447398061858056e-05 / 1.4473940925386564e-05 | 0.1500925884630277 / 0.0010635663736521229 | True / True |

quadratic2_y-1.0: Rust fval -0.6995600250752753, success True; iminuit fval -0.6995600250752838, success True.
quadratic2_y-1.2: Rust fval -0.6995455511608019, success True; iminuit fval -0.6995455512004952, success True.
Both libraries succeed from both starts. The earlier failure is absent at this Rust base commit.

no difference; median Rust/iminuit calls = 1.0263157894736843; median Rust/iminuit covariance errors = 1.0; paired rows = 17.

## Reference counts

Rust measured 31 Migrad calls and 13 additional Hesse calls.
Earlier measurements on cca0b9d were 31 Migrad calls and 13 additional Hesse calls.
ROOT measured 28 Migrad calls and 38 total calls after Hesse.
Earlier ROOT measurements were 28 Migrad calls and 38 total calls after Hesse.
The existing runner lacks this unbounded workload. A small counter probe links its existing ROOT libraries.

## Failed fits

No paired fit failed.
