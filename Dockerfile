# docker build -t toxindex .
# docker run -p 6513:6513 toxindex
FROM python:3.10-slim-buster

# Install PostgreSQL development package
RUN apt-get update && apt-get install -y libpq-dev

# Install required packages
COPY requirements.txt /requirements.txt
RUN pip install --trusted-host pypi.python.org -r /requirements.txt

# Make port 6513 available to the world outside this container
EXPOSE 6513

# Set the FLASK_APP environment variable
ENV FLASK_APP=webserver.app
ENV FLASK_ENV=development
ENV ROOT_URL=http://localhost:6513

# Command to run the application
# CMD ["flask", "run", "--host=0.0.0.0", "--port=6513"]
CMD ["flask", "run", "--debug", "--host=0.0.0.0", "--port=6513"]
