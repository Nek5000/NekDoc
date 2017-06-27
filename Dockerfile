FROM ubuntu:latest

# SET WORKING DIRECTORY TO /home
WORKDIR /home

# COPY CURRENT DIRECTORY CONTENTS TO /home
ADD . /home

# INSTALL NECESSARY DEPENDENCIES
RUN apt-get -y update
RUN apt-get install -y texlive-full
RUN apt-get install -y python-pip
RUN pip install sphinx

# MAKE PORT 80 AVAILABLE
EXPOSE 80

# RUN SPHINX MAKE WHEN CONTAINER LAUNCHES
CMD ["./docker_test"]
