M2R-ParallelQuicksort
=====================

This project is done originally as a homework for the [Scientific Methodology and Performance Evaluation](https://github.com/alegrand/SMPE) course taken in my last year of master. The goal of this project is to teach the design of sound experiments and the performance evaluation of a system (or algorithm) and then visualize the results in a clear way.

This also should be done using the best practice of [Reproducible Research](https://en.wikipedia.org/wiki/Reproducibility#Reproducible_research) as to document and specify all the details of the working/testing environment in a way that anyone can redo the same experiments and get the same results.

The motivation of the course become obvious during the project defense, when the students present their work and experiments showing (sometimes, completely) different results! Even though, the studied problem is considered simple.


## Software requirements

* C++ compiler: to compile and run the system under evaluation.
* Python: to run the experiments.
* R with the following packages:
	+ DoE.base: to design the experiments
	+ plyr and reshape2: to analyze the data
	+ ggplot2: to visualize the information
	+ rmarkdown: to produce nice reports
* (optional) pandoc with at least version 1.12.3

### On ubuntu

To install R
```
sudo apt-get install r-base
```
After that type `R` in the terminal then enter to following command to install the required libraries.
```
install.packages(c("DoE.base","ggplot2","plyr","reshape2","rmarkdown"))
```
To install pandoc
```
sudo apt-get install pandoc
```


## Run everything

A bash script is included to facilitate the whole process. It compiles the source code, generates the Design of Experiments, runs the experiments, analyzes the data and generates the final report.

Just go to the project folder in you terminal and type
```
sh build.sh
```
then go and take a cup of tea :)
