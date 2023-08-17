# Use an official Python runtime as a parent image
FROM python:3.9-slim-buster

# Install PostgreSQL development package
RUN apt-get update && apt-get install -y libpq-dev

# Install the requirements
COPY requirements.txt requirements.txt
RUN pip install --trusted-host pypi.python.org -r requirements.txt

# Expose the port the app runs on
EXPOSE 80

# Command to run the application
CMD ["python", "upload/app.py"]
