# Biospark Challenge 2022/23

## Data download
Taken from here: https://downloads.thebiogrid.org/File/BioGRID-ORCS/Release-Archive/BIOGRID-ORCS-1.1.13/BIOGRID-ORCS-ALL-homo_sapiens-1.1.13.screens.tar.gz

Information about the file format: https://wiki.thebiogrid.org/doku.php/orcs:downloads

```sh
mkdir data/
cd data/

# download data and unzip
wget https://downloads.thebiogrid.org/Download/BioGRID-ORCS/Release-Archive/BIOGRID-ORCS-1.1.13/BIOGRID-ORCS-ALL-homo_sapiens-1.1.13.screens.tar.gz
tar -xzvf BIOGRID-ORCS-ALL-homo_sapiens-1.1.13.screens.tar.gz

# remove unzipped file afterwards
rm BIOGRID-ORCS-ALL-homo_sapiens-1.1.13.screens.tar.gz
```

By default, these raw data files are not committed onto GitHub because of their large size. However, present in this repository is an Excel version of the Index file, which is easier to browse.

Also included in `data/` is `metadata_celltype.csv`, which is a (draft) mapping of the cell lines (`CELL_TYPE`) into broader cell types (`CELL_TYPE_BROAD`).

Manually download `homo_sapiens_gene_id.txt` from NCBI https://www.ncbi.nlm.nih.gov/gene/?term=Homo+sapiens.

## Environment setup
### Conda (Python)
Download Miniconda or Anaconda and use the `conda_env.yml` to create the environment. I've included some commonly used data science and plotting related packages for the time being.

```sh
conda install -c conda-forge mamba
mamba env create -f conda_env.yml
```

### R
Using the Terminal, `cd` into the root of this repository and start R. Then, run the following commands to install the dependency manager `renv` and use these commands to recreate the enviroment:
```r
install.packages("renv")

# then restart R
renv::activate()
renv::restore()
```

To use this environment in the RStudio IDE, use `File` > `Open Project` to choose the root of this repository to start a new project. You will see a `.proj` file created.

More information about `renv`: https://rstudio.github.io/renv/articles/renv.html

### GitHub
We should work with our independent branches and not send any direct commits to the `master` branch. The latter should only updated via pull requests.

I recommend putting in pull requests often for small changes, rather than leaving everything as a huge chunk.

Good tutorials about GitHub branching and collaboration.
* Branching: https://www.youtube.com/watch?v=QV0kVNvkMxc
* GitHub collaboration: https://www.youtube.com/watch?v=MnUd31TvBoU

