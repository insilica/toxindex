-- Create a custom ENUM type for the report status
CREATE TYPE report_status AS ENUM ('PENDING', 'GENERATED', 'FAILED');

CREATE TABLE reports (
    report_id SERIAL PRIMARY KEY,  -- Unique identifier for the report
    s3_reference VARCHAR(255) NOT NULL,  -- S3 reference or link to the report file (e.g., 's3://toxindex/report1234.pdf')
    user_id INTEGER NOT NULL, -- User ID associated with the report (you may want to set up a foreign key relationship here)
    title VARCHAR(255), -- Title or name of the report
    description TEXT, -- Description or summary of the report
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP, -- Timestamp when the report was created
    updated_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP, -- Timestamp when the report was last updated
    status report_status DEFAULT 'PENDING' -- Status of the report generation (e.g., pending, generated, failed)
);


-- Function to update the 'updated_at' column to the current timestamp
CREATE OR REPLACE FUNCTION update_updated_at_column()
RETURNS TRIGGER AS $$
BEGIN
   NEW.updated_at = now(); -- Set the 'updated_at' field to the current time
   RETURN NEW;
END;
$$ language 'plpgsql';

-- Trigger to run the 'update_updated_at_column' function before each update on the 'reports' table
CREATE TRIGGER trigger_name
BEFORE UPDATE ON reports
FOR EACH ROW EXECUTE FUNCTION update_updated_at_column();
