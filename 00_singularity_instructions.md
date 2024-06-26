# Instructions for setting up the container 
## Verify singularity is installed
First, ensure you have the singularity program installed or loaded if you are working in an HPCC environment. If you wish to run the container on your local machine but don't have singularity installed then please refer to 
the link below for instructions on how to setup singularity on your local machine. The only necessary sections you need to follow are the following.

* Install system dependencies
* Install Go
* Clone the repo
* Compiling SingularityCE

Note that the last command in the section Compiling SingularityCE doesn’t need to be run. This command is for if you wish to modify where singularity is installed.

Singularity install instructions: https://github.com/sylabs/singularity/blob/main/INSTALL.md
<br></br>
## Building the container
In order to build the container, you will need a **definition file** (.def), as this file provides the instructions from which singularity will follow to build the container. It is **suggested** that you build the container with **sudo** privileges. If you don't have sudo privileges on the machine you are working with, you can still build the container; however, the singularity documentation suggests this method may not always work as intended. There are two types of containers you can build: **sandbox images** and **singularity image files**. This section will cover how to build both, starting with singularity image files.
<br></br>
A singularity image file (.SIF for short) is the standard container format most widely used and accessible (some HPCC environments may restrict sandbox images). A SIF container is compressed as much as possible by singularity and is immutable, meaning once the container has been built, root files can't be edited. In order to build the container with sudo privileges, open a terminal from which the definition file is located and type the following command: 

``` console
sudo singularity build <container-name>.sif pipeline.def
```

In order to build the container with non sudo privileges, open a terminal from which the definition file is located and type the following command:

``` console
singularity build --fakeroot <container-name>.sif pipeline.def
```

Refer to the Singularity documentation for more information not covered here.
SingularityCE documentation: https://sylabs.io/docs/

<br></br>
## Running the container
There are a couple of things to make note of when running the container. By default, singularity will **mount your home directory** to the container, meaning that any changes you make to files or directories inside the container will **persist 
even after the container has been terminated**. Host environment variables will also be passed into the container and may sometimes work as defined or not, depending on whether the path a environment variable points to exist within the container.

In some cases, this is great, as directories like scratch (a common directory when working in an HPCC environment) will be available in the container; however, if you have any of the following programs **loaded or installed** then you must execute the container with some special flags passed or if you're in an HPCC environment unload them.

* Conda
* R
* Blast
* E Utilities 

Having either of these programs installed or loaded without passing the flags in the following command will cause issues inside the container as the container won't know which programs to reference and environment variables already written inside the container could be overwritten. These flags ensure such issues don't occur as they provide even more isolation from the host environment. In general,  if you're running the container on your **local machine** this is how you should be running the container.

``` console
singularity shell --cleanenv --no-home <container-name>.sif
```

If you're running the container in an **HPCC environment** then you can just unload all the programs mentioned, or if you don't have **any of these programs** installed on your local machine, then you can run the container using the following command:
``` console
singularity shell <container-name>.sif
```

<br></br>
## Binding external software or databases to the container
Since this pipeline requires the use of **usearch 64 bit**, in order to access this program inside of the container, it will need to be **binded to a directory inside the container**. The directory you should bind usearch to is the /bind/ directory. You could
also bind other executable programs to this directory as well. 

In order to bind usearch to the container, you will need to know the **path** to your usearch executable. If you have a local copy of usearch then simply store it in a directory and pass in
the relative path to usearch to the container. 

``` console
singularity shell --bind path/to/usearch:/bind/ <container-name>.sif
```
Similarly, if you wish to access a database such as the NCBI, Uniprot, Pfam, etc, then you will need to bind the path to whatever database you wish to access to the /database/ directory inside of the container as follows:
``` console
singularity shell --bind path/to/usearch:/bind/ --bind path/to/database/:/database/ <container-name>.sif
```

## Developing within a container
Typically, when developing software specifically for the container, you'll want to be working with sandbox containers. These are containers that are writeable (basically just a writeable directory on your system), meaning users can make live changes to the environment within the container. You **will need sudo privileges** in order to even make changes to the sandbox container, as without these privileges you won't have the privileges necessary inside of the container to make changes to the environment. You can build the container without these privileges; however, like with any container, it is suggested that you build the container with these privileges.

Sandbox containers can be built without adding a file extension, as they are just directories on your system. However, if you wish to build the container from within the repository directory, then I recommend adding the **.simg** extension to this container, as any file or directory ending in this extension will be ignored when making changes to the repository. In order to build the container with sudo privileges, open a terminal from which the definition file is located and type the following command: 

``` console
sudo singularity build --sandbox <container-name>.simg pipeline.def
```

In order to build the container with non sudo privileges, open a terminal from which the definition file is located and type the following command:

``` console
singularity build --fakeroot --sandbox <container-name>.simg pipeline.def
```

Interacting with the sandbox container can be done similarly to how you would interact with a regular container, with the following exceptions: you need the --writeable flag, and you need sudo privileges to make any changes to the container's environment. You don't need the --cleanenv or --no-home flags when interacting with the container; however, it is best to pass these in either way for the same reasons stated in the previous section. To run the container, enter the following command into a terminal from where the container is located:

``` console
singularity shell --cleanenv --no-home --writeable <container-name>.sif
```

<br></br>
## Verifying that the container is working properly
First, check to make sure the /bind/ and /database/ directories did in fact bind the directories passed in. Next, check the output of the conda environment list by pasting in the following command:
``` console
conda env list
```
Expected output:
``` console
# conda environments:
#
base                     /opt/miniconda
qiime1                *  /opt/miniconda/envs/qiime1
```
If there are only **two** conda environments and the active conda environment is **qiime1**, then conda is working as intended. After this, ensure the version of R running is 4.3.3, and the output of checking the version is outputting appropriately. 
``` console
R --version
```
Expected output:
``` console
R version 4.3.3 (2024-02-29) -- "Angel Food Cake"
Copyright (C) 2024 The R Foundation for Statistical Computing
Platform: x86_64-pc-linux-gnu (64-bit)

R is free software and comes with ABSOLUTELY NO WARRANTY.
You are welcome to redistribute it under the terms of the
GNU General Public License versions 2 or 3.
For more information about these matters see
https://www.gnu.org/licenses/.
```
