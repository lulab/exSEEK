# refer to biocontainer
FROM dongzhuoer/ubuntu-cn:rolling
# may need if you `FROM ubuntu:rolling` -----------
# # install locales
# RUN apt update && apt -y upgrade && apt -y install locales && rm -r /var/lib/apt/lists/ \
#     && echo "en_US.UTF-8 UTF-8" > /etc/locale.gen && locale-gen

# ENV LANG=en_US.UTF-8 LANGUAGE=en_US:en LC_ALL=en_US.UTF-8 

# # enable extended glob: http://mywiki.wooledge.org/glob#extglob
# RUN echo "\n#enable wildcard\nshopt -s extglob" >> /etc/bash.bashrc

# Set a default user ---------------
RUN useradd exseek -m -s /bin/bash
WORKDIR /home/exseek

# install R ----------
# tzdata is a dependency of r-base-dev, so it's the same as we just install the latter
RUN apt update \
 && DEBIAN_FRONTEND=noninteractive apt -y install tzdata \
 && apt -y install add-apt-key && apt-key adv --keyserver hkp://keyserver.ubuntu.com:80 --recv-keys E084DAB9 \
 && apt -y install software-properties-common && apt-add-repository -y "deb https://mirrors4.tuna.tsinghua.edu.cn/CRAN/bin/linux/ubuntu bionic-cran35/" \
 && apt update && apt -y install r-base-dev && rm -r /usr/local/lib/R/site-library/ \
 && apt -y purge add-apt-key software-properties-common && apt -y autoremove && rm -r /var/lib/apt/lists/

RUN echo 'R_LIBS_USER="~/.local/lib/R"' > /usr/lib/R/etc/Renviron.site 

RUN echo "options(Ncpus = 8L)" >> /usr/lib/R/etc/Rprofile.site 
RUN echo "options(BioC_mirror = 'https://mirrors4.tuna.tsinghua.edu.cn/bioconductor')" >> /usr/lib/R/etc/Rprofile.site 
RUN echo "options(repos = c('CRAN' = 'https://mirrors4.tuna.tsinghua.edu.cn/CRAN'))" >> /usr/lib/R/etc/Rprofile.site 



RUN mkdir -p /home/exseek/.local/lib/R && chown -R exseek:exseek /home/exseek

 