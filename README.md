# What Does It Mean to Be a Broker?
## Brokers in the European Co-Inventor Network

### Course: Applied Network Science (ADIN150NAMB) – Group Project

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

[cite_start]Using the **OECD REGPAT database** (Jan 2024 edition) spanning 1978–2023, this project analyzes the structural role of brokers across Europe. We identify "Weak" and "Strong" brokers using **Burt’s Constraint**. The constraint $C_i$ for an inventor $i$ is calculated as:

$$
C_i = \sum_{j \in V_i, i \neq j} \left( p_{ij} + \sum_{q \in V_i, q \neq i,j} p_{iq} p_{qj} \right)^2
$$

[cite_start]Where $p_{ij}$ is the proportional strength of $i$'s tie to $j$. Inventors exhibiting low constraint ($C_i < \mu$) are defined as weak brokers, while those with significantly low constraint ($C_i \leq \mu - \sigma$) are strong brokers.

The study further employs feature-based role discovery to classify brokers into four distinct typologies, revealing that European brokerage is an elite, highly centralized phenomenon.

### Repository Structure
```text
.
├── data/
│   ├── regpat_nodes.csv           # anonymized inventor data (active inventors, degree >= 2)
│   ├── regpat_edges.csv           # co-inventor weights and ties
│   └── country_codes.csv          # mapping of ISO codes to country names
├── codes/
│   ├── 01_network_build.ipynb     # backbone extraction and constraint calculation
│   ├── 02_inequality.ipynb        # HHI and Lorenz curve analysis
│   ├── 03_robustness.ipynb        # targeted node removal simulations
│   ├── 04_community.ipynb         # Louvain community & core-periphery detection
│   ├── 05_broker_typology.ipynb   # K-means clustering of broker roles
│   └── 06_main_plots.py           # generation of maps and network visualizations
├── figures/
│   ├── backbone_map.png           # max spanning tree of the network
│   ├── broker_geography.png       # map of broker profiles by country
│   └── robustness_plot.png        # impact of removal on path length
├── documents/
│   └── network_groupwork.pdf      # final research paper
├── requirements.txt               # Python dependencies
└── README.md                      # project description and usage
```

### Scripts
The analysis is divided into sequential notebooks corresponding to the methodology sections of the paper.

**file: `01_network_build.ipynb`**
-   Constructs the cumulative co-inventor network from the OECD REGPAT database.
-   Filters for the "active" inventor population (Degree $\geq 2$) to exclude isolates.
-   Calculates Burt’s Constraint (Eq. 1) to identify Weak ($C < \mu$) and Strong ($C \leq \mu - \sigma$) brokers.

**file: `02_inequality.ipynb`**
-   Calculates the Herfindahl-Hirschman Index (HHI) to measure the concentration of brokerage power across 31 European countries and IPC4 technological fields.
-   Generates Lorenz curves to visualize distributional inequality, finding high geographic concentration (HHI $\approx 0.79$ for strong brokers)

**file: `03_robustness.ipynb`**
-   Simulates reverse percolation-style targeted node removal on the giant connected component (GCC).
-   Compares the impact of removing top-degree brokers ($k=100, 500, \dots, 10000$) versus random nodes on the average shortest path length to assess network fragility.
-   Targeted removal increases path length by $\approx 40\%$, whereas random removal shows minimal impact.

**file: `04_community.ipynb`**
-   Applies the Louvain algorithm to detect community structures, finding a high modularity of $0.994$.
-   Determines the core-periphery placement of brokers using k-core decomposition.
-   Definitions used: Periphery ($k_{core} \leq 2$), Core ($k_{core} \geq 5$), Intermediate (otherwise).

**file: `05_broker_typology.ipynb`**
-   Implements feature-based role discovery using K-means clustering ($k=4$) on strong brokers.
-   Identifies four broker profiles:
    -   **Type 1:** High clustering (embedded in cohesive teams, mean local clustering $\approx 0.95$).
    -   **Type 2:** High degree (super-star hubs, mean degree $\approx 125.5$).
    -   **Type 3:** Low constraint (bridging connectors, constraint $\approx 0.16$).
    -   **Type 4:** Other (peripheral/weakly embedded).

### Data
The project utilizes data derived from the **OECD REGPAT database** (January 2024 edition).

**file: `regpat_nodes.csv` & `regpat_edges.csv`**
-   Represents the co-inventor network where nodes are inventors and edges are weighted by shared patents.
-   The analysis is restricted to inventors with at least two co-inventor links (Degree $\ge 2$) to exclude incidental one-time inventors.

### Figures and Reports

**file: `figures/backbone_map.png`**
-   Visualizes the maximum spanning tree of NUTS2-level collaboration ties, highlighting high-activity regions like Germany and France (Figure 1).

**file: `figures/broker_geography.png`**
-   Maps the clustering of countries based on their broker role composition, revealing distinct groups (e.g., Western/Northern economies vs. Eastern Europe) (Figure 2).

**file: `documents/network_groupwork.pdf`**
-   The final research paper detailing the finding that European brokers form a dense "Elite Club" and that brokerage is technologically specialized rather than a generalist phenomenon.
