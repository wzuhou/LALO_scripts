# LALO Plots and scripts
Date: 18/09/2024

Author: *Zhou Wu* (Institute: The Roslin institute, University of Edinburgh, UK)

The scripts used to create plots in the LALO genome and stress transcriptomic paper.

## Graphic abstract

<p align="center">
  <img src="https://github.com/wzuhou/LALO_scripts/blob/main/Graphic%20abstract.png">
</p>


## Pipeline

```mermaid
%%{init: {'theme':'forest'}}%%

stateDiagram
    #[*] --> Lapland_longspur 
    #genome --> [*]
    OmniC --> Lapland_longspur
    PacBio --> Lapland_longspur
    RNAseq --> Lapland_longspur
    Lapland_longspur --> genome
    genome --> DEG_extreme_events
    RNAseq --> DEG_extreme_events
    DEG_extreme_events --> Extreme_Spring
    DEG_extreme_events --> Snowstorm
Extreme_Spring --> Stress_related_DEGs(e.g.,FKBP5)
Snowstorm --> Stress_related_DEGs(e.g.,FKBP5)

```    
