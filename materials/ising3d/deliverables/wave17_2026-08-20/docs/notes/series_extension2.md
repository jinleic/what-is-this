# Series extension, resource wall, and strict holdout analysis

[COMPUTATION] The resource profile below was printed before either coefficient computation. The finite-lattice computation then extended the low-temperature reduced free energy through `x^32`; the same run reached `v^22` on the high-temperature side.

## Profile first

[LEMMA] The HT box bound is `sum(side_i-1) <= order/2`; the LT bound is `sum(side_i) <= floor((order+6)/4)`. Sorting a box puts its largest side in the transfer direction, so the transfer cross-section is the product of its two smallest sides. The memory column is the CRT engine's conservative three-array estimate `3*2^cross_section*(full_degree+1)*8` bytes.

| target | canonical boxes | new budget-frontier canonical boxes | wall box | N | cross-section | peak memory | guard <=22 |
|---|---:|---|---|---:|---:|---:|---|
| HT `v^22` | 83 | `[[1, 1, 12], [1, 2, 11], [1, 3, 10], [1, 4, 9], [1, 5, 8], [1, 6, 7], [2, 2, 10], [2, 3, 9], [2, 4, 8], [2, 5, 7], [2, 6, 6], [3, 3, 8], [3, 4, 7], [3, 5, 6], [4, 4, 6], [4, 5, 5]]` | `[4, 5, 5]` | 100 | 20 | 5939134464 (5.531 GiB) | yes |
| HT `v^24` | 102 | `[[1, 1, 13], [1, 2, 12], [1, 3, 11], [1, 4, 10], [1, 5, 9], [1, 6, 8], [1, 7, 7], [2, 2, 11], [2, 3, 10], [2, 4, 9], [2, 5, 8], [2, 6, 7], [3, 3, 9], [3, 4, 8], [3, 5, 7], [3, 6, 6], [4, 4, 7], [4, 5, 6], [5, 5, 5]]` | `[5, 5, 5]` | 125 | 25 | 242397216768 (225.750 GiB) | no |
| LT `x^30` | 23 | `[[1, 1, 7], [1, 2, 6], [1, 3, 5], [1, 4, 4], [2, 2, 5], [2, 3, 4], [3, 3, 3]]` | `[3, 3, 3]` | 27 | 9 | 1339392 (1.277 MiB) | yes |
| LT `x^32` | 23 | `[[1, 1, 7], [1, 2, 6], [1, 3, 5], [1, 4, 4], [2, 2, 5], [2, 3, 4], [3, 3, 3]]` | `[3, 3, 3]` | 27 | 9 | 1339392 (1.277 MiB) | yes |

[COMPUTATION] `HT v^22` fits the current guard at cross-section 20. `HT v^24` does not: the required canonical `5x5x5` box has cross-section 25, beyond the hard guard 22, and its hypothetical three-array peak is 242,397,216,768 bytes (225.750 GiB). Thus `v^24` is the quantified resource wall.

[COMPUTATION] The frozen `x^28` FLM artifact actually has canonical wall `2x3x3`, cross-section 6. The 20-spin item stored in the frozen HT artifact is explicitly a separate `4x5x2` engine probe, not an LT-required box. This corrects the resource-wall inventory without changing any coefficient.

## Exact extensions and witnesses

[COMPUTATION] The locally derived LT coefficients are `c_30=346966/5` and `c_32=-427509/2` in `phi=3K+sum c_n x^n`. They were obtained from the local plus-boundary finite-lattice calculation before the external table was instantiated. Every coefficient through frozen order `x^28` agrees exactly.

[EXTERNAL] Guttmann and Enting, *Series studies of the Potts model. I. The simple cubic Ising model*, J. Phys. A 26 (1993) 807-821, DOI `10.1088/0305-4470/26/4/010`, Table 2 (arXiv PDF page 21/24), print `lambda_15=65337` and `lambda_16=-200934` for `Lambda_0(u)=sum lambda_n u^n`, `u=x^2`. An exact local formal logarithm of their printed table gives `[346966/5, -427509/2]`, agreeing coefficient-by-coefficient with the local `x^30,x^32` result. This is an external primary-source witness, not the source of the local coefficients.

[LEMMA] The internal `4x4xc` periodic-slab route in `experiments/e26_lt_series_independent.py` is rigorous only below its first transverse wrapping term `x^(4L)`. To witness through `x^32` without wrapping requires `4L>32`, hence `L>=9` in both periodic transverse directions and an 81-spin cross-section. That exceeds the current guard 22. Therefore no independent internal slab witness is available for `x^30,x^32`; the honest witness status is external-primary-source only.

[COMPUTATION] The locally derived HT `v^22` coefficient is `360734896503/22`. When reached, the same Table 2 value `a_22=16809862992` for `Phi(v)` gives the exact formal-log total-free-energy witness `360734896503/22` after adding the known `3/22` from `3 log(cosh K)`. The recorded process peak RSS was 6066749440 bytes against the propagation-array estimate 5939134464 bytes.

## Strict holdout prediction

[COMPUTATION] For the homogeneous Euler-ODE holdout, each model is fitted only to the first `m` coefficients of `G(z)=(phi-log 2)/z`, `z=v^2`; all remaining coefficients are recursively predicted. The fixed selector minimizes `|order-degree|` and then uses lexicographic degrees, without seeing any holdout.

| m | held out | selected (order,degree) | absolute errors | relative errors |
|---:|---:|---|---|---|
| 5 | 6 | `[1, 2]` | `['53357477/142000', '15229658787509/1270332000', '69562111721452531/300645240000', '25438916719193510206061/7684492334400000', '4670996780894755407546919961/147311718050448000000', '35954934955161565577978074039729/719065323733749300000000']` | `['53357477/3414283500', '15229658787509/405452136654000', '69562111721452531/1357340483661592500', '25438916719193510206061/514887647750954200800000', '4670996780894755407546919961/152092974230284433451084000000', '35954934955161565577978074039729/11790543415272283784711304450000000']` |
| 6 | 5 | `[6, 0]` | `['4468383/14', '72236127/16', '134006939/2', '4129826907/4', '360734896503/22']` | `['1', '1', '1', '1', '1']` |
| 7 | 4 | `[1, 3]` | `['2668642054204683591/36603229222000', '50695057118625065399178173/15943313970110962400', '1430606150847025147631879831161958532731/14589419395157160337555581320000', '470285743355821423696448208288539907390978921999403/187315640488046285141954471314911819625000']` | `['889547351401561197/55084906556052149875', '50695057118625065399178173/1068257351325253780784046800', '1430606150847025147631879831161958532731/15062944193906926563937560585840644310000', '470285743355821423696448208288539907390978921999403/3071422190220387867243480585132854435136171785062500']` |
| 8 | 3 | `[2, 2]` | `['81388058555809444387999649/1508530517741320788820', '55287748227499919436253263611342948329712803131/17053668814232961767987193747738250751800', '226802575385410864810044120739338955107895085044036767440516245464693/1850020131969046709978632186402019797165146541367609873505500']` | `['81388058555809444387999649/101076778535299796363416810990', '18429249409166639812084421203780982776570934377/5869058361007172506311316997569299862324718223550', '75600858461803621603348040246446318369298361681345589146838748488231/10111618497489704038238099077321775992129184703016259936295264610625250']` |
| 9 | 2 | `[1, 4]` | `['112537174927757132161289848117/94149691672341237090480', '9614618862385248336014593990726154047369386866083/113735394751413885602133680623078198413600']` | `['112537174927757132161289848117/97205482488547167155482674386340', '9614618862385248336014593990726154047369386866083/1864923902471778974345580687069844111727531581256400']` |
| 10 | 1 | `[10, 0]` | `['360734896503/22']` | `['1']` |

[COMPUTATION] Balanced first-order inhomogeneous differential approximants were likewise fitted to only the first `m` coefficients of `F(z)=phi-log 2`; the table reports genuine recursive predictions of every remaining coefficient.

| m | held out | selected `(deg Q1,deg Q0,deg P)` | absolute errors | relative errors |
|---:|---:|---|---|---|
| 7 | 5 | `[1, 1, 2]` | `['96845007262/50970325', '1835074294578386761/25449541524300', '158759925078478552313125421/83389590166793554125', '31849040983739934835890942816287711/728639059271703713426801250', '656389057208330581005017179688297475354559/700334856422576601167623127403375']` | `['193690014524/32536419104925', '7340297178313547044/459594078410277096525', '28865440923359736784204622/1015889429337863966656552125', '63698081967479869671781885632575422/1504576596235724909720910488595616875', '1312778114416661162010034359376594950709118/22966838359003775914072003438033242159808875']` |
| 8 | 4 | `[1, 2, 2]` | `['13923495197985453/1605217228370', '24375777481221477369281716219/55215478934054244275505', '184143956255562069128193433072879892560571/12661834008234953714737296397439550', '18130861145056345346580317079490805697325735907222458/45627537759301861578568208980289000933024565']` | `['37129320527961208/19325779261853887165', '48751554962442954738563432438/7399257317371592135318697729195', '368287912511124138256386866145759785121142/26145591389588185714510844549290009747985925', '36261722290112690693160634158981611394651471814444916/1496313191935316505643639551766420719873080626366508745']` |
| 9 | 3 | `[2, 2, 2]` | `['179010023987436470443/2006664793391265', '1674813411141134903756270147451246049/383495580289143642820893171450', '1937885118901910710740067281436411513802565209863006/13436536306763737730624029385845380622779725']` | `['358020047974872940886/268907006561430851987835', '3349626822282269807512540294902492098/791885183096842128054861000613387102575', '3875770237803821421480134562872823027605130419726012/440638866725383526284554121282524002919916889276527425']` |
| 10 | 2 | `[2, 2, 3]` | `['21435764283495735519186362/74481208705869790545', '5347206216014081605249896622102368055093286780/348696885446632415267099983796936487813']` | `['85743057133982942076745448/307594499779383709831195194315', '10694412432028163210499793244204736110186573560/11435194080228126477387475432358260081789715983449']` |
| 11 | 1 | `[2, 3, 3]` | `['1144039080326151555755332353760356/975229798709497641127252897']` | `['254230906739144790167851634168968/3553529500041536914496587980323251709']` |

## Finite D-finiteness budget

[COMPUTATION] Every `(order,degree)` pair with parameter budget `(order+1)(degree+1) <= 11` was tested by exact rational rank on all known coefficients of `G`. Writing `F=zG` conjugates `theta` on `F` to `theta+1` on `G`, so this change preserves ODE order, polynomial degree, and the budget. The result is: **no D-finite recurrence of budget <= 11 is consistent with the known coefficients**. Budget 12 is the first vacuous budget because it has more operator coefficients than the 11 data. This finite negative result is not a proof that the full series is non-D-finite.

## Untuned series-only critical estimate

[COMPUTATION] With the model family and pole filters fixed without external input, the extended HT series gives `K_c=0.22217230511543145888577136509626297 +/- 0.0058034119744074035537324910340241222`. The uncertainty is the largest displacement from the longest-series median among every retained approximant at the recorded truncations; it is a method-spread envelope, not a statistical confidence interval.

[EXTERNAL] Final comparison only: the benchmark `K_c=0.221654626` differs from the untuned series-only central estimate by `0.00051767911543145888577136509626297` and was not used in fitting, filtering, or model selection.
