#from here this has openfpm inside

#FROM danielweilandt/openbread_docker:latest
FROM python:3.6-stretch


USER root

RUN apt-get update && apt-get install -y --no-install-recommends \
        libxml2-dev     \
        libxslt1-dev    \
        libopenblas-dev \
        liblapack-dev   \
		less			\
		build-essential \
		gfortran        \
		fort77          \
		wget            \
		cmake           \
        libflint-2.5.2  \
        libflint-dev    \
        libgmp-dev      \
		yasm            \
		xvfb            \
		xauth           \
		ffmpeg          \
        firefox-esr

# Clean up lists
RUN rm -rf /var/lib/apt/lists/*

ENV USER open
ENV HOME /home/$USER

RUN useradd -ms "/bin/bash" "$USER"
USER $USER
WORKDIR $HOME


USER root

#RUN ln -s /usr/bin/python3 python
#RUN pip3 install --user --upgrade pip
#RUN mkdir /src
#COPY src/ /src/


RUN python3 -m pip install --upgrade pip
RUN python3 -m pip install --upgrade pipenv


#install hdf package to store data as hdf files
#RUN pip3 install 'HDF5'
# Take care of the solvers and HDF5
COPY ./HDF5Files /HDF5Files
COPY ./solvers /solvers
COPY ./utils /utils
RUN chmod u+x /utils/*.sh
# Install CPLEX
RUN /utils/install_cplex.sh
#install gurobi 800
COPY ./utils/gurobi.lic* ./
RUN /utils/install_gurobi.sh
#activate required solvers
#install hdf5
RUN /utils/install_hdf5.sh

#install openbread
#first already put this id_rsa your public ssh key here to clone from c4science openbread
RUN mkdir /root/.ssh/
#RUN echo "${SSH_PRIVATE_KEY}" > /root/.ssh/id_rsa
# make sure your domain is accepted
RUN touch /root/.ssh/known_hosts
RUN ssh-keyscan c4science.ch >> /root/.ssh/known_hosts
#this id rsa you need to generate it and also paste it to c4science public ssh and put private file here so that it can login with your user name!!


RUN rm -rf /solvers

# these for hdf for tables to find hdf
ENV HDF5_DIR='/usr/local/hdf5'
ENV LD_LIBRARY_PATH='/usr/local/hdf5/lib'

ENV LD_LIBRARY_PATH=$LD_LIBRARY_PATH:"/usr/local/lib"
#this lets see if adter build run will be super fast or not if yes leave it in built

COPY requirements.txt .
RUN pip install -r requirements.txt
#when from openbread docker
#RUN pip install --user -r requirements.txt

#USER $USER
RUN mkdir /open
#RUN mkdir ./work

#for colorful
COPY .bashrc $HOME
RUN chown "$USER" "$HOME/.bashrc"
#RUN chown -R "$USER" "$HOME/.cache

#permission is root normally change it #somehow $USER doesnt work need to five string!
#RUN chown -R open:open /home/open/.local

#put this at last as it has export LD_Library_Path
RUN /utils/activate_gurobi.sh

#Finalizing installation
#USER user
USER $USER
RUN mkdir ./work
WORKDIR ./work


ENTRYPOINT ["/bin/bash", "-c",  "pip install --user -e /open  &&  $0 $*"]

CMD /bin/bash











