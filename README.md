# Research Repository

**Title:** On the estimation of spatial density from mobile network operator data

**Authors (Article):** Fabio Ricciato and Angelo Coluccia

**Authors (Code Repository):** Marco Ramljak, Fabio Ricciatio and Angelo Coluccia

**Abstract:** We tackle the problem of estimating the spatial distribution of mobile phones from Mobile Network Operator (MNO) data, namely CDR or signalling data. The process of transforming MNO data to a density map requires to perform geolocation of radio cells to determine their spatial footprint. Traditional geolocation solutions rely on Voronoi tessellations and approximate cell footprints by mutually disjoint regions. Recently, some pioneering work started to consider more elaborated geolocation methods with partially overlapping (non-disjont) cell footprints coupled with a probabilistic model for phone-to-cell association. Estimating the spatial density in such a probabilistic setup is currently an open research problem and is the focus of the present work. We start by reviewing three different estimation methods proposed in literature and provide novel analytical insights that unveil some key aspects of their mutual relationships and properties. Furthermore, we develop a novel estimation approach for which a closed-form solution can be given. Numerical results based on semi-synthetic data are presented to assess the relative accuracy of each method. Our results indicate that the estimators based on overlapping cells have the potential to improve spatial accuracy over traditional approaches based on Voronoi tessellations.

**Journal:** ... (under review)

**Full article:** <https://arxiv.org/pdf/2009.05410.pdf>

**Supplementary Material and code for reproduction:** <https://r-ramljak.github.io/MNO_mobdensity/>

|           |                                                                                                                                                             |
|-----------|-------------------------------------------------------------------------------------------------------------------------------------------------------------|
| *authors* | Ramljak M., Ricciato F., Coluccia, A.                                                                                                                       |
| *version* | 1.0                                                                                                                                                         |
| *status*  | 2021 - closed                                                                                                                                               |
| *license* | [EUPL](https://joinup.ec.europa.eu/sites/default/files/custom-page/attachment/eupl_v1.2_en.pdf) *(concerning the source code, please cite this repository)* |
