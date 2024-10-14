FROM julia:1.11

COPY . /src

ENV JULIA_DEPOT_PATH=/usr/local/share/julia

RUN julia -e 'using Pkg; Pkg.add(path="/src"); Pkg.precompile()' \
 && rm -rf /src \
# smoke test
 && julia -e 'using Fronts'
