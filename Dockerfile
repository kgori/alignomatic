# Build stage
FROM rust:bookworm AS builder

# Install build dependencies
# alignomatic needs clang/llvm for bindgen, and zlib for htslib/bio
RUN apt-get update && apt-get install -y \
    clang \
    llvm \
    libclang-dev \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /usr/src/alignomatic

# Create a blank project to cache dependencies
RUN cargo init

# Copy manifests
COPY Cargo.toml Cargo.lock ./

# Build dependencies (this layer is cached if manifests don't change)
RUN cargo build --release && \
    rm src/*.rs && \
    rm ./target/release/deps/alignomatic*

# Copy actual source code
COPY src ./src
# Copy other necessary files if any (e.g., build.rs, C libraries)
# COPY build.rs ./ 

# Build the actual application
RUN cargo build --release

# Runtime stage
FROM debian:bookworm-slim

# Install runtime dependencies
RUN apt-get update && apt-get install -y \
    zlib1g \
    libbz2-1.0 \
    liblzma5 \
    samtools \
    ca-certificates \
    && rm -rf /var/lib/apt/lists/*

# Copy the binary from the builder
COPY --from=builder /usr/src/alignomatic/target/release/alignomatic /usr/local/bin/alignomatic

# Set the entrypoint
ENTRYPOINT ["alignomatic"]
