-- Create a custom ENUM type for the report status

CREATE TABLE report_object (
    report_object_id SERIAL PRIMARY KEY,  -- Unique identifier for the report
    report_id INTEGER NOT NULL REFERENCES reports(report_id), -- Reference to the associated report
    s3_reference VARCHAR(255) NOT NULL,  -- S3 reference or link to the report file (e.g., 's3://toxindex/report1234.pdf')
    mimetype VARCHAR(255) NOT NULL -- S3 mimetype of the report file (e.g., 'application/pdf')    
);
