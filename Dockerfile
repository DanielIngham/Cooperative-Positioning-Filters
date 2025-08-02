FROM ubuntu:24.04

# Install system dependencies
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        build-essential \
        cmake \
        gnuplot \
        libeigen3-dev \
        && \
    apt-get clean && rm -rf /var/lib/apt/lists/*

# Set working directory
WORKDIR /app

# Copy source code into container
COPY . .

# Create a separate build directory
RUN rm -rf ./build && mkdir -p build && cd build && \
    cmake .. && \
    cmake --build .

# Set working directory to run built app
WORKDIR /app/build

# Set the default command to run your application
# Replace `your_executable` with the actual binary name
CMD ["bash"]
