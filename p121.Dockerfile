FROM debian

WORKDIR /app
RUN apt update && apt install -y make mpich curl && \
    curl -L https://github.com/ParaFEM/parafem/releases/download/5.0.3/parafem.5.0.3.tar.gz -o parafem.5.0.3.tar.gz && \
    tar -xf parafem.5.0.3.tar.gz

WORKDIR /app/parafem/
RUN /app/parafem/make-parafem MACHINE=linuxdesktop

RUN cp /app/parafem/bin/* /usr/local/bin
WORKDIR /
COPY entrypoint.sh entrypoint.sh
WORKDIR /app/test121
RUN cp  /app/parafem/examples/5th_ed/p121/demo/*.bnd /app/test121
RUN cp  /app/parafem/examples/5th_ed/p121/demo/*.lds /app/test121
RUN cp  /app/parafem/examples/5th_ed/p121/demo/*.dat  /app/test121
RUN cp  /app/parafem/examples/5th_ed/p121/demo/*.d  /app/test121
ENTRYPOINT ["/entrypoint.sh"]
CMD ["sh"]

