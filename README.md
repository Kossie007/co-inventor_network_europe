# What Does It Mean to Be a Broker?
## Brokers in the European Co-Inventor Network

### Course: Applied Network Science – Group Project

### Team Members:
-   Gellért Banai
-   Patrik Bruinsma
-   Hunor Kuti
-   Ákos Virág

### Mentors:
-   Sándor Juhász
-   László Lőrincz

### Project Overview
This repository contains the code and analysis for a research project mapping the "brokerage" effect within the European co-inventor network. Innovation is often described as a recombinant process, relying on the transfer of knowledge across disconnected clusters. Actors who span these structural holes—known as brokers—are theorized to be vital for innovation, acting as gatekeepers and translators of knowledge.

Using the **OECD REGPAT database** (Jan 2024 edition) spanning 1978–2023, this project analyzes the structural role of brokers across Europe. We identify "Weak" and "Strong" brokers using **Burt’s Constraint**. The constraint $C_i$ for an inventor $i$ is calculated as:

$$
C_i = \sum_{j \in V_i, i \neq j} \left( p_{ij} + \sum_{q \in V_i, q \neq i,j} p_{iq} p_{qj} \right)^2
$$

Where $p_{ij}$ is the proportional strength of $i$'s tie to $j$. Inventors exhibiting low constraint ($C_i < \mu$) are defined as weak brokers, while those with significantly low constraint ($C_i \leq \mu - \sigma$) are strong brokers.

The study further employs feature-based role discovery to classify brokers into four distinct typologies, revealing that European brokerage is an elite, highly centralized phenomenon.

---
### Repository Structure
```text
.
├── data/
│   └── Confidential               # We are waiting for the green light to publish it.
├── codes/
│   ├── 01_brokerage.r            # under construction - backbone extraction and constraint calculation
│   ├── 02_attacking.r            # under construction - 
│   ├── 03_community.r            # under construction - 
│   ├── 04_core-periphery.r       # under construction - 
│   ├── 05_broker_typology.r      # under construction - 
│   ├── 06_clustering-egos.r      # under construction - 
│   └── 01-06_brokerige_network.r # whole analysis codes
├── figures/
│   ├── dist_coinv_across_countries.png           # interindustry based Lorenz-curve of cumulative brokerage distribution
│   ├── dist_internat_coinv.png                   # international based Lorenz-curve of cumulative brokerage distribution
│   └── coinventor_map.png                        # map of co-inventorship by country
├── documents/
│   └── network_groupwork.pdf     # under review - final research paper 
└── README.md                     # project description and usage
```

---
### Scripts
The analysis is divided into sequential notebooks corresponding to the methodology sections of the paper.

file: `codes/01_brokerage.r`
-   

file: codes/02_attacking.r`
-   

file: `codes/03_community.r`
-   

file: `codes/04_core-periphery.r`
-   

file: `codes/05_broker_typology.r`
-   

file: `codes/06_clustering-egos.r`
-   

file: `codes/01-06_brokerige_network.r`
-   

---
### Data
The project utilizes data derived from the **OECD REGPAT database** (January 2024 edition).

file: `data/confidential`
-   Represents the co-inventor network where nodes are inventors and edges are weighted by shared patents.
-   The analysis is restricted to inventors with at least two co-inventor links (Degree $\ge 2$) to exclude incidental one-time inventors.

### Figures and Reports

file: `figures/dist_coinv_across_countries.png`
-   

file: `figures/dist_internat_coinv.png `
-   

file: `documents/network_groupwork.pdf`**
-   

---
### Licence
MIT License (MIT): see the [License File](https://github.com/sensiolabs/GotenbergBundle/blob/1.x/LICENSE) for more details.

