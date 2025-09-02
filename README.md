# flood-risk

![R](https://img.shields.io/badge/R-4.3.0-blue)

# Flood Risk Prediction in Buncombe County, NC

## Authors
- Katie Miller — [katie.miller6572@gmail.com](mailto:katie.miller6572@gmail.com)  
- Maya Hall — [mayalhall@gmail.com](mailto:mayalhall@gmail.com)  
- Truman Anarella — [truman.anarella@gmail.com](mailto:truman.anarella@gmail.com)  

---

## Project Description
This project aimed to leverage different machine learning algorithms to predict flood risk areas in Buncombe County, NC. With more frequent and severe storms and flood events, understanding the potential for flood risk is crucially important. We combined several machine learning techniques to achieve our goal of predicting flood risk areas.  

As a feasibility study, we chose to examine conditions in Buncombe County in September 2024 due to the recent and tragic flood events. Based on existing literature, we selected the following variables to include in our models: **Topographic Wetness Index (TWI)**, **distance to streams (DTS)**, **Normalized Difference Vegetation Index (NDVI)** from *Landsat 9 OLI-2*, **National Landcover Database (NLCD)** data, and **social vulnerability** from the *U.S. Census Bureau*.  

For our flood risk layer, we used the **National Flood Hazard Layer (NFHL)** 100-year floodplain data from *FEMA*, and decided to incorporate a **1 km buffer area** as well. For modeling, we conducted a **principal component analysis (PCA)**, and developed a **random forest classification** and a **support vector machine (SVM) classification**.  

---

## Open-Source Workflow
To ensure that this work is adaptable and accessible for practitioners, we used completely open-source data and software. The workflow created can be replicated for other locations in the United States.  

We hope that our code can be used to inform decision-making around disaster mitigation and resilience planning.  

---

## Workflow
*Insert workflow diagram here*

---

## License
This project is licensed under the **MIT License** — see the [LICENSE](./LICENSE) file for details.
