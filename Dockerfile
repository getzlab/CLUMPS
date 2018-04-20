FROM centos

## YUM INSTALLS
RUN yum -y install gcc && yum clean all
RUN yum -y install nano && yum clean all
RUN yum -y install java-1.8.0 && yum clean all
RUN yum -y install epel-release && yum clean all
RUN yum -y install python-devel && yum clean all
RUN yum -y install python-pip && yum clean all
RUN yum -y install python2-jpype && yum clean all
RUN yum -y install python-biopython && yum clean all
RUN yum -y install wget && yum clean all

#RUN rpm -i http://download.fedoraproject.org/pub/epel/6/i386/epel-release-6-5.noarch.rpm
RUN yum --enablerepo=epel install pymol && yum clean all


## PYTHON LIB INSTALLS
RUN pip install numpy
RUN pip install scipy
RUN pip install twobitreader
RUN pip install prody
RUN pip install lxml
RUN pip install statsmodels
 
## COPY CLUMPS/EMPRINT CODE
RUN mkdir /sw
COPY sw/ /sw/

