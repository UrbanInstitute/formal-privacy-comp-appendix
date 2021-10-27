#!/bin/bash

# install script for cloud syntheses - run from the RStudio Docker Container Terminal Pane

# install R packages
sudo su - -c "R -e \"install.packages('rmutil', repos='http://cran.rstudio.com/')\""
sudo su - -c "R -e \"install.packages('LaplacesDemon', repos='http://cran.rstudio.com/')\""
sudo su - -c "R -e \"install.packages('reticulate', repos='http://cran.rstudio.com/')\""
sudo su - -c "R -e \"install.packages('furrr', repos='http://cran.rstudio.com/')\""
sudo su - -c "R -e \"install.packages('here', repos='http://cran.rstudio.com/')\""
sudo su - -c "R -e \"install.packages('tictoc', repos='http://cran.rstudio.com/')\""
sudo su - -c "R -e \"install.packages('ipumsr', repos='http://cran.rstudio.com/')\""
sudo su - -c "R -e \"install.packages('srvyr', repos='http://cran.rstudio.com/')\""
sudo su - -c "R -e \"install.packages('synthpop', repos='http://cran.rstudio.com/')\""
sudo su - -c "R -e \"install.packages('caret', repos='http://cran.rstudio.com/')\""
sudo su - -c "R -e \"install.packages('moments', repos='http://cran.rstudio.com/')\""
sudo su - -c "R -e \"install.packages('randomForest', repos='http://cran.rstudio.com/')\""
sudo su - -c "R -e \"install.packages('rpart.LAD', repos='http://cran.rstudio.com/')\""
sudo su - -c "R -e \"install.packages('reticulate', repos='http://cran.rstudio.com/')\""
sudo su - -c "R -e \"install.packages('yardstick', repos='http://cran.rstudio.com/')\""

# download Tax-Calculator from github
cd $HOME
git clone https://github.com/UI-Research/Tax-Calculator.git

# install anaconda
wget https://repo.anaconda.com/archive/Anaconda3-2019.10-Linux-x86_64.sh
bash Anaconda3-2019.10-Linux-x86_64.sh -b
source $HOME/anaconda3/bin/activate
conda init
source ~/.bashrc
cd $HOME

# install tax conda env and taxcalc
cd Tax-Calculator
conda env create -f environment.yml
conda activate taxcalc-dev
python setup.py install
cd $HOME

# install privacy conda env
cd formal-privacy-comp
conda env create -f environment.yml
cd $HOME
