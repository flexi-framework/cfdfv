# Installation on Ubuntu 20.04

1. Install and set up neccessary software:
```
$ sudo apt update && sudo apt upgrade
$ sudo apt install git gfortran cmake 
$ wget "https://www.paraview.org/paraview-downloads/download.php?submit=Download&version=v5.8&type=binary&os=Linux&downloadFile=ParaView-5.8.0-MPI-Linux-Python3.7-64bit.tar.gz" -O ParaView-5.8.0-MPI-Linux-Python3.7-64bit.tar.gz
$ tar xzvf ParaView-5.8.0-MPI-Linux-Python3.7-64bit.tar.gz
$ sudo mv ParaView-5.8.0-MPI-Linux-Python3.7-64bit /opt/ParaView
$ echo "export PATH='\$PATH:/opt/ParaView/bin'" >> .bashrc
```

2. Clone `cfdfv` repository:
```
$ git clone https://github.com/flexi-framework/cfdfv.git
```

3. Compile the CGNS library and the code:
```
$ cd cfdfv
$ make shared && make
```
