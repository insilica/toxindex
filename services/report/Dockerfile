FROM python:3.10-slim-buster

# Install PostgreSQL development package
RUN apt-get update && apt-get install -y libpq-dev

# Install the requirements
COPY requirements.txt requirements.txt
RUN pip install --trusted-host pypi.python.org -r requirements.txt

# Expose the port the app runs on
EXPOSE 6515

# Command to run the application
ENV FLASK_APP=report.app
ENV FLASK_ENV=development
ENV ROOT_URL=http://localhost:6515

CMD ["flask", "run", "--debug", "--host=0.0.0.0", "--port=6515"]
