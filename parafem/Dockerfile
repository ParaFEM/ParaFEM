FROM debian

WORKDIR /app
RUN apt update && apt install -y make mpich curl

COPY  . /app
RUN /app/make-parafem MACHINE=linuxdesktop

RUN cp /app/bin/* /usr/local/bin

CMD ["sh"]

