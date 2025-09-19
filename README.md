# Causal Tool


## Prerequisites
To run the tool, make sure you have installed, at least in the **MAIN** folder of this repository, a **R version $\eq$ 4.4.3** (used to write the procedure).
To visualize and interact with the tool output, make sure you have installed **Cytoscape version 3.10.3**, available at https://cytoscape.org/download.html.

## Input preparation

Pick a _Vitis vinifera_ or a _Homo Sapiens_ gene and find its OneGenE expansion list.

### Vitis vinifera

To check if your Vitis vinifera (Vv) gene has been already expanded, go to [VvOneGenE](http://ibdm.disi.unitn.it/onegene/vv/onegene-vv.php) and under **Gene name(s)** type the _Ordered Locus Name_ (_VIT_XXsYYYYgZZZZZ_) or _gene name_  of your Vv gene. You can also type multiple names (space separated) and get multiple expansion lists as a result.
For example: 
1. Type VIT_04s0008g06000 in the **gene name(s)** box; this name corresponds to the transcription factor VvERF045. 
2. You are then redirected to the output page where you can **check** the expansion list of VIT_04s0008g06000 and press the **download** button. 
3. The expansion list will appear in your **Donwnload** folder as a zip compressed folder, extract it and choose the corresponding **.exp.csv** file (_54651_Vv-VIT_04s0008g06000.exp.csv_). The expansion list is already annotated with additional information about the candidate genes that could be useful for a biologist.
4. Move the espansion list to the **MAIN** folder of this project to provide it as input to the ***causal tool***. 


### Homo sapiens

To check if your human gene has been already expanded, go to [HsOneGenE](http://ibdm.disi.unitn.it/onegene/fantom/onegene-fantom.php) and under **Gene name(s)** type the _gene name_ or _symbol_ of your Hs gene. You can also type multiple names (space separated) and get multiple expansion lists as a result.
For example:
1. Type U2AF1 in the **gene name(s)** box;
2. You are then redirected to the output page where you can **check** the expansion list of the chosen isoform of U2AF1, such as p1@U2AF1 and press the **download** button.
3. The expansion list will appear in your **Donwnload** folder as a zip compressed folder. Unzip the folder to get access to the corresponding **.csv** file (_T105641_p1@<!-- -->U2AF1.csv_)
4. Move the espansion list to the **MAIN** folder of this project to provide it as input to the ***causal tool***.


## Transcriptomic dataset preparation

### Homo sapiens

Right now, the file _fantom_mat.csv_ is a placeholder for the actual FANTOM-full transcriptomic dataset, a gene@home version of the FANTOM5 transcriptomic dataset. This file should be replaced and renamed accordingly in order for the procedure to work. The FANTOM-full transcriptomic dataset can be downoloaded from: [Human OneGenE download page](https://gene.disi.unitn.it/test/download/bc/hgnc_data_mat.csv.gz). The file will need to be extracted, renamed (_fantom_mat.csv_) and placed in the **MAIN** folder.

### Vitis vinifera

Right now, the file _vespucci_mat.csv_ is a placeholder for the actual VESPUCCI transcriptomic dataset. This file should be replaced and renamed accordingly in order for the procedure to work. The VESPUCCI transcriptomic dataset can be downoloaded from: [Vitis OneGenE downolad page](http://ibdm.disi.unitn.it/download-2/). The file will need to be extracted, renamed (_vespucci_mat.csv_) and placed in the **MAIN** folder.


## Input submission

To run the causal, make sure you have a terminal panel open in the **MAIN** folder and type the following command:

  ```
  % Rscript install_packages.R
  ```

In this way, all the necessary R pacakges will be installed, if not present.

Next, you can type the actual command that runs the tool:

  ```
  % Rscript main.R explist.csv organism_type n pc_version
  ```
The arguments you have to provide are:

- **explist.csv**, which corresponds to the annotated expansion list of the gene under investigation. For Vv gene _VIT_04s0008g06000_ it is _54651_Vv-VIT_04s0008g06000.exp.csv_ and for Hs gene _p1@<!-- -->U2AF1_, it is _T105641_p1@<!-- -->U2AF1.csv_;
- **organism_type**, which corresponds to the organism to which the gene under investigation belongs, either _Hs_ or _Vv_;
- **n**, which can be the first **n** genes you select from the expansion list: make sure that **n** is not greater than the expansion list length and that **n** is an integer number followed by the letter _L_. For example, _100L_ means the first 100 genes of the chosen expansion list. This value can be substituted by the relative frequency threshold, according to which you can cut the expansion list, by selecting only the candidate genes with relative frequency >= **n** (0 <**n** $\leq$ 1): make sure that **n** is a real number;
- **pc_version**, which can be selected from _mrf_ (maj.rule fast), _cf_ (conservative fast), _o_ (original), _s_ (stable), _sf_ (stable.fast), _mr_ (maj.rule), _c_ (conservative);

Here are some examples:

#### Vitis vinifera

 ```
  % Rscript main.R 54651_Vv-VIT_04s0008g06000.exp.csv Vv 0.7 o

  ```
  
  ```
  % Rscript main.R 54651_Vv-VIT_04s0008g06000.exp.csv Vv 150L cf

  ```

#### Homo Sapiens

  ```
  % Rscript main.R T105641_p1@U2AF1.csv Hs 0.5 mrf

  ```
  
  ```
  % Rscript main.R T105641_p1@U2AF1.csv Hs 200L s

  ```
  
**Note:** A modified expansion list of p1@U2AF1, used in the biological validation of the paper is made available for testing and reproducibility.

## Output visualisation in Cytoscape

The causal tool has two output files:

 - a list of edges (_\_edges.csv_), which represents the interactions retrieved by pc() between the surviving input gene nodes, divided into _source_ and _target_, and the direction of their interaction, --- if undirected,  --> if directed, or <-> if bidirected. Also, the pearson correlation (_zero_order_r_) computed between the input genes, as the zero-order conditional independence test, is provided, along with its sign (_r_sign_).
 
##### Homo sapiens example _T105641_p1@<!-- -->U2AF1.csv_

| source | interaction | target | zero_order_r | r_sign |
| :---   |     :---  |   :--- |  :---  |:---    | 
| T005391| --- | T126388 | 0.831005760842777| + |
| T017334 | --> |T017335 | 0.893899138354596 | + |

 - a list of nodes (_\_nodes.csv_), which represents the input gene nodes that survived after pc() application and for which an interaction was found in the output graph. Additional information, extracted from human and grapevine annotation files, is added for the biological interpretation of the results.

##### Homo sapiens example _T105641_p1@<!-- -->U2AF1.csv_

| transcript | gene_name | coords | entrezgene_id | hgnc_id | uniprot_id | description | type | rank | Frel | TID |
| :--- | :---  | :--- | :--- |:--- | :--- | :--- | :--- | :--- | :--- | :--- |
| p1@SCARB2 | SCARB2   | chr4:77134970..77134985;-   | e:950   | h:1665 | u:Q14108 | scavenger receptor class B member 2   | gene with protein product | 2 | 1 | T145731|


These two files are contained in the **Vv** folder inside the **MAIN** folder, if a grapevine expansion list was chosen as input to the tool, or they are contained in the **Hs** folder inside the **MAIN** folder, if a human expansion list was chosen as input.


To visualize the pc() ouptut graph on Cytoscape, do the following steps:

1. Open Cytoscape and allow the app to accept incoming network connections;
2. Select the _network_ icon from the main horizontal toolbar, which stands for _Import Network from File System_ (or from File -> Import -> Network from file) and select the _\_edges.csv_ file from the **Vv** folder or **Hs** folder inside **MAIN**;
3. Click the OK button in the _Import Network from Table_ panel and wait for the network to load;
4. Select the _table_ icon from the main horizontal toolbar, which stands for _Import Table from File_ (or from File -> Import -> Table from file) and select the _\_nodes.csv_ file from the **Vv** folder or **Hs** folder inside **MAIN**;
5. In the _Import Columns from Table_ panel, make sure that the **TID** column has the **key** icon (click on **TID** and select it from the menu), then click the OK button and wait for the table to load;
6. The **Vv** and **Hs** folders contain respectively a _Vv_style.xml_ and a _Hs_style.xml_ that can be uploaded in Cytoscape to customize the network appearance (From File -> Import -> Styles from file...). This feature is managed by the _Style_ panel (under _Network_ in the main vertical toolbar), from which you can selected the uploaded style and visualize the network in a more human-friendly and enriched way.

## Reproducibility

The relevant information from the R ```sessionInfo()``` command is presented here to reproduce the output:

```
Output: R version 4.4.3 (2025-02-28)
Matrix products: default
BLAS:   /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/lib/libRblas.0.dylib 
LAPACK: /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.0

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
[1] pcalg_2.7-12 readr_2.1.5 

loaded via a namespace (and not attached):
 [1] bit_4.6.0           compiler_4.4.3      BiocManager_1.30.25
 [4] crayon_1.5.3        tidyselect_1.2.1    Rcpp_1.0.14        
 [7] parallel_4.4.3      fastICA_1.2-7       cluster_2.1.8.1    
[10] R6_2.6.1            igraph_2.1.4        RBGL_1.82.0        
[13] robustbase_0.99-4-1 BiocGenerics_0.52.0 bdsmatrix_1.3-7    
[16] graph_1.84.1        tibble_3.2.1        sfsmisc_1.1-20     
[19] pillar_1.10.2       tzdb_0.5.0          rlang_1.1.6        
[22] bit64_4.6.0-1       cli_3.6.5           magrittr_2.0.3     
[25] vroom_1.6.5         hms_1.1.3           lifecycle_1.0.4    
[28] clue_0.3-66         DEoptimR_1.1-3-1    vctrs_0.6.5        
[31] glue_1.8.0          ggm_2.5.1           corpcor_1.6.10     
[34] abind_1.4-8         stats4_4.4.3        tools_4.4.3        
[37] pkgconfig_2.0.3
```
