# ihouaga_adgg_spatial

Files for ihouaga_adgg_spatial project where we want to evaluate the benefit of spatial modelling on the accuracy of genomic prediction in African small holder dairy farms

data
  * data folder
  * see README.txt for details
  
To get a copy of the data from Datastore do:
  * On Laptop (Apple/Mac)
    ```
    cd to_your_working_directory
    git clone https://github.com/HighlanderLab/ihouaga_adgg_spatial.git
    cd ihouaga_adgg_spatial
    mkdir data
    # load Datastore as mentioned in https://github.com/HighlanderLab/HighlanderLab_Handbook/blob/main/Data/Workspaces.md#Datastore
    cp /Volumes/HighlanderLab/ADGG_spatial/* /exports/cmvm/eddie/eb/groups/HighlanderLab/your_uun/ihouaga_adgg_spatial/data/.
    ```
  * On Eddie
    ```
    ssh your_uun@eddie.ecdf.ed.ac.uk
    cd /exports/cmvm/eddie/eb/groups/HighlanderLab/your_uun/
    git clone https://github.com/HighlanderLab/ihouaga_adgg_spatial.git
    cd ihouaga_adgg_spatial
    mkdir data
    qlogin -q staging
    cp /exports/cmvm/datastore/eb/groups/HighlanderLab/ADGG_spatial/* /exports/cmvm/eddie/eb/groups/HighlanderLab/your_uun/ihouaga_adgg_spatial/data/.
    ```

scripts
  * all the scripts folder
  * see README.txt for details

results
  * all intermediate and final results
  * see README.txt for details