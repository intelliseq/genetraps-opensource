FROM    ubuntu:18.04

# Install python3
RUN     apt-get update
RUN     apt-get install -y python3
RUN     apt-get install -y curl

# Copy html
ADD src/ /src
WORKDIR /src

# Run http server on port 8080
CMD ["python3", "-m", "http.server", "8087"]
