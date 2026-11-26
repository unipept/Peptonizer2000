FROM rust:1-bookworm AS build
RUN apt-get update && apt-get install -y --no-install-recommends python3 \
    && rm -rf /var/lib/apt/lists/*
WORKDIR /src
COPY Cargo.toml Cargo.lock ./
COPY peptonizer_rust ./peptonizer_rust
COPY peptonizer_analysis ./peptonizer_analysis
RUN cargo build --release --locked -p peptonizer_analysis

FROM debian:bookworm-slim
ARG GIT_COMMIT=unknown
LABEL org.opencontainers.image.source="https://github.com/unipept/Peptonizer2000" \
      org.opencontainers.image.revision=$GIT_COMMIT
RUN apt-get update && apt-get install -y --no-install-recommends ca-certificates \
    && rm -rf /var/lib/apt/lists/*
COPY --from=build /src/target/release/taxonomic_analysis \
                  /src/target/release/protein_inference \
                  /src/target/release/functional_analysis \
                  /usr/local/bin/
RUN useradd --create-home peptonizer
USER peptonizer
WORKDIR /data
