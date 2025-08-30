# Dockerfile

# Use a Conda-based Python image as our foundation
FROM continuumio/miniconda3

# Set the working directory inside the container
WORKDIR /app

# Copy our dependency file
COPY requirements.txt .

# Install RDKit using Conda, which is the most reliable method
RUN conda install -c conda-forge rdkit

# Install the rest of the dependencies using pip
RUN pip install -r requirements.txt

# Copy the rest of our application code into the container
COPY . .

# Tell Render what command to run when the container starts
CMD ["gunicorn", "--workers", "1", "--bind", "0.0.0.0:10000", "app:app"]
