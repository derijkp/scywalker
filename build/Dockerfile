FROM ubuntu:20.04

RUN ln -snf /usr/share/zoneinfo/$CONTAINER_TIMEZONE /etc/localtime && echo $CONTAINER_TIMEZONE > /etc/timezone

RUN apt-get update \
    && apt-get upgrade -y \
    && apt-get install -y \
      wget \
      build-essential \
      curl libsm6 libgtk2.0-0 libegl1 libxt6 libxtst6 libxext6 libgl1-mesa-glx less \
      wget \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /usr/local/bin

RUN wget https://github.com/derijkp/scywalker/releases/download/0.110.0/scywalker-0.110.0-linux-x86_64.tar.gz \
  && tar xvzf scywalker-0.110.0-linux-x86_64.tar.gz \
  && rm -f scywalker-0.110.0-linux-x86_64.tar.gz
RUN chmod go+rX /usr/local/bin/scywalker-0.110.0-linux-x86_64

RUN mkdir -p /home/user
RUN chmod ugo+rwX /home/user
ENV HOME=/home/user

RUN ln -sf scywalker-0.110.0-linux-x86_64/scywalker_makerefdir scywalker_makerefdir

RUN ln -sf scywalker-0.110.0-linux-x86_64/sw sw

RUN ln -sf scywalker-0.110.0-linux-x86_64/cg cg

CMD ["scywalker", "help"]
