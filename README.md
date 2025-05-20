# JBI-Special-Issue-BITS-2022


## Prerequisites

To run the data processing procedure make sure you have installed in the **MAIN** folder the following:

- **R version ≥ 3.5.0**. The R version used to write this data processing procedure is **4.1.1 (2021-08-10) -- "Kick Things"** and it is available at https://cran.r-project.org/src/base/R-4/ for **MacOs**, https://cran.r-project.org/bin/windows/base/old/4.1.1/ for **Windows** and https://cran.r-project.org/doc/manuals/r-patched/R-admin.html for **Unix**.
- **PHP**, available at https://www.php.net/manual/en/install.php

To visualize and interact with the data processing procedure, make sure you have installed **Cytoscape version 3.9.1**, available at https://cytoscape.org/download.html

## Input preparation for the data processing procedure

Pick a _Vitis vinifera_ or a _Homo Sapiens_ gene and find its OneGenE expansion list.

### Vitis vinifera

To check if your Vitis vinifera (Vv) gene has been already expanded, go to [VvOneGenE](http://ibdm.disi.unitn.it/onegene/vv/onegene-vv.php) and under **gene name(s)** type the _Ordered Locus Name_ (_VIT_XXsYYYYgZZZZZ_) or _gene name_  of your Vv gene. You can also type multiple names (space separated) and get multiple expansion lists as a result.
For example: 
1. Type VIT_04s0008g06000 in the **gene name(s)** box; this name corresponds to the transcription factor VvERF045. 
2. You are then redirected to the output page where you can **check** the expansion list of VIT_04s0008g06000 and press the **download** button. 
3. The expansion list will appear in your **Donwnload**n folder as a zip compressed file, just extract it (_54651_Vv-VIT_04s0008g06000.exp.csv_) to use it. The expansion list is already annotated with additional information about the candidate genes that could be useful for a biologist.
4. Move the espansion list to the **MAIN** folder of this project to provide it as input to the ***data processing procedure***. 


### Homo sapiens

To check if your human gene has been already expanded, go to [HsOneGenE](https://gene.disi.unitn.it/test/gene\_history-z.php), choose _Homo Sapiens (Hs)_ in the **Organism** box and leave **Tile size** and **Iterations** blank; the significance level **alpha** is set to 0.05 by default. Under **LGN name**, type the _gene symbol_ of you Hs gene.
For example:
1. Type MFSD2A in the **LGN name** box.
2.  You are then redirected to the output page where the first result contains the expansion list of MFSD2A. Click on its **pcim_id** (193111) and the download will start automatically. 
3. The expansion list will appear in your **Donwnload** folder as a zip compressed folder (_193111_Hs.zip_) . Unzip the folder to get access to its content: the _.interactions_ file (_193111_Hs.interactions_) is the output of NES2RA, while the _.expansion_ file (_193111_Hs.expansion_) is the actual expansion list, which is not in its working form, you have to annotate it first.
4. Move the _.interactions_ file to the **MAIN** folder of this project, open a new terminal panel in this folder and type the following command:

  ```
  % php anno-hsf5.php  file.interactions
  ```
5. You expansion list is now available in **MAIN** in csv format (_193111_Hs_p1@<!-- -->MFSD2A.csv_) and you can provide it as input to the data processing procedure.


## Transcriptomic dataset preparation for the data processing procedure

### Homo sapiens

Right now, the file _fantom_mat.csv_ is a placeholder for the actual FANTOM-full transcriptomic dataset, a gene@home version of the FANTOM5 transcriptomic dataset. This file should be replaced and renamed accordingly in order for the data processing procedure to work. The FANTOM-full transcriptomic dataset can be downoloaded from: [Human OneGenE download page](https://gene.disi.unitn.it/test/download/bc/hgnc_data_mat.csv.gz). The file will need to be extracted, renamed (_fantom_mat.csv_) and placed in the **MAIN** folder.

### Vitis vinifera

Right now, the file _vespucci_mat.csv_ is a placeholder for the actual VESPUCCI transcriptomic dataset. This file should be replaced and renamed accordingly in order for the data processing procedure to work. The VESPUCCI transcriptomic dataset can be downoloaded from: [Vitis OneGenE downolad page](http://ibdm.disi.unitn.it/download-2/). The file will need to be extracted, renamed (_vespucci_mat.csv_) and placed in the **MAIN** folder.


## Input submission to the data processing procedure

To run the data processing procedure, make sure you have a terminal panel open in the **MAIN** folder and type the following command:

  ```
  % Rscript install_packages.R
  ```

In this way, all the necessary R pacakges will be installed, if not present.

Next, you can type the actual command that runs the data processing procedure:

  ```
  % Rscript --vanilla data_processing_procedure.R explist.csv organism_type n
  ```
The arguments you can provide are:

- **explist.csv**, which corresponds to the annotated expansion list of the gene under investigation. For Vv gene _VIT_04s0008g06000_ it is _54651_Vv-VIT_04s0008g06000.exp.csv_ and for Hs gene _MFSD2A_, it is _193111_Hs_p1@<!-- -->MFSD2A.csv_ (as explained in [input-preparation](https://github.com/Camilla9347/JBI-Special-Issue-BITS-2022/edit/main/README.md#input-preparation-for-the-data-processing-procedure));
- **organism_type**, which corresponds to the organism to which the gene under investigation belongs;
- **n**, which can be the first **n** genes you select from the expansion list (make sure that **n** is not greater than the expansion list length) or the relative frequency threshold according to which you can cut the expansion list, by selecting only the candidate genes with relative frequency >= **n** (0 <**n** <= 1).

Here are some examples:

#### Vitis vinifera

 ```
  % Rscript --vanilla data_processing_procedure.R 54651_Vv-VIT_04s0008g06000.exp.csv Vv 0.7

  ```
  
  ```
  % Rscript --vanilla data_processing_procedure.R 54651_Vv-VIT_04s0008g06000.exp.csv Vv 150

  ```

#### Homo Sapiens

  ```
  % Rscript --vanilla data_processing_procedure.R 193111_Hs_p1@MFSD2A.csv Hs 0.5

  ```
  
  ```
  % Rscript --vanilla data_processing_procedure.R 193111_Hs_p1@MFSD2A.csv Hs 200

  ```
  
**Note:** The expansion lists of NFKB1 and TNF, used in the biological validation of the paper are made available for testing.

## Output visualisation in Cytoscape

The data processing procedure has two output files:

 - a list of edges (_gene_edges.csv_), which represents the interactions retrieved by pc_parallel() between the surviving input gene nodes, divided into _source_ and _target_, and the direction of their interaction, --- if undirected or  --> if directed. Also the pearson correlation (_cor_) computed between the input genes, as the zero-order conditional independence test, is provided, along with its sign (_cor_sign_)
 
##### Homo sapiens example _p1@<!-- -->MFSD2A_edges.csv_

| spurce | interaction | target | cor | cor_sign |
| :---   |     :---  |   :--- |  :---  |:---    | 
| T178190| --- | T009518 | 0.455 | + |
| T032201 | --> | T054717 | 0.699 | + |

 - a list of nodes (_gene_nodes.csv_), which represents the input gene nodes that survived after pc_parallel() application and for which an interaction was found in the output graph. Additional information, extracted from human and grapevine annotation files, is added for the biological interpretation of the results.

These two files are contained in the **Vv** folder inside the **MAIN** folder, if a grapevine expansion list was chosen as input to the data processing procedure, or they are contained in the **Hs** folder inside the **MAIN** folder, if a human expansion list was chosen as input.

##### Homo sapiens example _p1@<!-- -->MFSD2A_nodes.csv_

| ID | association_with_transcript | entrezgene_id | hgnc_id | uniprot_id | description | rank | Frel | type |
| :---    |     :---  |   :--- |  :---  |:---    | :---                      | :--- | :---  |  :---                     |  
| T178190 | p4@CDK6   | CDK6   | 1777   | Q00534 | cyclin dependent kinase 6 | 14   | 0.996 | gene with protein product |


To visualize the pc_parallel() ouptut graph on Cytoscape, do the following steps:

1. Open Cytoscape and allow the app to accept incoming network connections;
2. Select the _network_ icon from the main horizontal toolbar, which stands for _Import Network from File System_ (or from File -> Import -> Network from file) and select the _gene_edges.csv_ file from the **Vv** folder or **Hs** folder inside **MAIN**;
3. Click the OK button in the _Import Network from Table_ panel and wait for the network to load;
4. Select the _table_ icon from the main horizontal toolbar, which stands for _Import Table from File_ (or from File -> Import -> Table from file) and select the _gene_nodes.csv_ file from the **Vv** folder or **Hs** folder inside **MAIN**;
5. Click the OK button in the _Import Columns from Table_ panel and wait for the table to load;
6. The **Vv** and **Hs** folders contain respectively a _Vv_style.xml_ and a _Hs_style.xml_ that can be uploaded in Cytoscape to customize the network appearance (From File -> Import -> Styles from file...). This feature is managed by the _Style_ panel (under _Network_ in the main vertical toolbar), from which you can selected the uploaded style and visualize the network in a more human-friendly and enriched way.
