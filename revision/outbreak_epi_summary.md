# Outbreak epidemiologic & transmission-classification summary

Route indicators: ● primary/strong  ◐ contributory  – not evidenced.
theta = common-source fraction (source_decomposition); 1-theta = propagated.
Model: additive renewal-with-source for resolution-eligible series; theta allocates source and propagated incidence inside the recursion. Monthly/fortnightly -> static only.

| Outbreak | Location, country | Period | Suspected | Vehicle / source | Pt-src | Env-water | P2P | theta (mode) | Res. | Model | Independence |
|---|---|---|---|---|:--:|:--:|:--:|---|---|---|---|
| Lewis 2005 | Kathmandu, Nepal | 2002 | ~5,900 | Single municipal water supply | ● | ● | – | 0.80 (A common-source) | daily | Additive renewal-with-source (source-dominant) | independent |
| Aye 2004 | Madaya, Mandalay, Myanmar | 2000 | 49 | Myaung River water + contact | ◐ | ● | ● | 0.50 (C mixed) | weekly | Additive renewal-with-source | independent |
| Muyembe-Tamfum 2009 | Kinshasa, DR Congo | 2004 | peritonitis series | Presumed municipal water | ◐ | ◐ | – | 0.65 (A common-source) | daily | exclude (complications) | independent |
| Ali 2017 | Kikwit, DR Congo | 2017 | ~1,500 | Camp water + overcrowding | ● | ● | ◐ | 0.50 (C mixed) | weekly | Additive renewal-with-source (+spatial) | independent |
| Kabwama 2017 | Kampala, Uganda | Jan-Jun 2015 | 10,230 | Vended/spring water (recurrent) | ◐ | ● | ◐ | 0.60 (C mixed) | daily | Additive renewal-with-source | independent |
| Lutterloh 2012 | Neno/Tsangano (border), Malawi/Mozambique | Mar-Nov 2009 | 303 | Spring water | ◐ | ● | ◐ | 0.55 (C mixed) | fortnightly | Static only | independent |
| Hendriksen 2015 | Lusaka, Zambia | 2010-2012 | ~2,040 | Water/sanitation (endemic) | – | ● | ● | 0.45 (C mixed) | monthly | Static only (+ODE) | independent |
| Davis 2018 | Harare, Zimbabwe | 2016-2017 | ~880 | Boreholes | ◐ | ● | – | 0.72 (A common-source) | daily | Additive renewal-with-source (+ODE) | Harare-recurrent (season 1) |
| N'Cho 2019 | Harare, Zimbabwe | 2017-2018 | ~580 | Boreholes/shallow wells | ◐ | ● | – | 0.80 (A common-source) | daily | Additive renewal-with-source (+ODE) | Harare-recurrent (season 2) |
| Neil 2012 | Kasese, Uganda | 2008-2009 | 709 (+ IP-heavy) | No vehicle implicated | – | ◐ | ● | 0.20 (B propagated) | weekly | Additive renewal-with-source (propagated-dominant) | KASESE-UGANDA (canonical) |
| Walters 2014 (Kasese) | Kasese, Uganda | 2009-2011 | 709 | Piped water (taps/rivers) | – | ● | ◐ | 0.20 (B propagated (=Neil)) | monthly | Static (redundant) | KASESE-UGANDA (=Neil) |
| Walters 2014 (Bundibugyo) | Bundibugyo, Uganda | 2011 | 333 | Gravity-flow-scheme piped water | ● | ● | ◐ | 0.70 (A common-source) | monthly | Static only | KASESE-UGANDA (linked spread) |
| Polonsky 2014 (Dzivaresekwa) | Harare (Dzivaresekwa), Zimbabwe | 2012 | 916 | Borehole water | – | ● | ● | 0.50 (C mixed (propagated-lean)) | weekly | Additive renewal-with-source (propagated-lean) | HARARE 2011-12 (canonical) |
| Polonsky 2014 (Kuwadzana) | Harare (Kuwadzana), Zimbabwe | 2012 | 1,860 | Borehole water | – | ● | ● | 0.50 (C mixed) | weekly | Additive renewal-with-source | HARARE 2011-12 |
| Muti 2014 | Harare (Dzivaresekwa), Zimbabwe | 2011 | 1,078 | Wells | – | ● | ● | 0.50 (C mixed) | daily | Additive renewal-with-source (redundant) | HARARE 2011-12 |
| Imanishi 2014 (Dzivaresekwa) | Harare (Dzivaresekwa), Zimbabwe | 2011-2012 | 4,181 (Harare total) | Boreholes/shallow wells | ◐ | ● | – | 0.50 (A/C water) | monthly | Static (redundant) | HARARE 2011-12 |
| Imanishi 2014 (Kuwadzana) | Harare (Kuwadzana), Zimbabwe | 2011-2012 | 4,181 (Harare total) | Boreholes/shallow wells | ◐ | ● | – | 0.50 (A/C water) | monthly | Static (redundant) | HARARE 2011-12 |
| Qamar 2018 | Hyderabad, Pakistan | 2016-2017 | ~586 | Sewage-contaminated water | – | ● | ● | 0.50 (C mixed (propagated, data theta 0.05)) | weekly | Additive renewal-with-source (propagated-lean) | HYDERABAD-XDR (canonical) |
| Yousafzai 2019 | Hyderabad, Pakistan | 2016-2018 | 101 (confirmed) | Sewage-contaminated water | – | ● | ● | 0.50 (C mixed (propagated, data theta 0.38)) | weekly | Additive renewal-with-source (redundant) | HYDERABAD-XDR |

**Independent epidemics: ~12** (from 19 series). Groups collapse: Kasese-Uganda (Neil+Walters x2), Harare 2011-12 (Polonsky x2 + Muti + Imanishi x2), Hyderabad-XDR (Qamar+Yousafzai); Muyembe excluded (complications).
