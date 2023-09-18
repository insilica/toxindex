CREATE TABLE project_reports (
    project_report_id SERIAL PRIMARY KEY,  -- Unique identifier for the project report association
    project_id INTEGER NOT NULL REFERENCES projects(project_id), -- Reference to the associated project
    report_id INTEGER NOT NULL REFERENCES reports(report_id), -- Reference to the associated report
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP, -- Timestamp when the association was created
    updated_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP -- Timestamp when the association was last updated
);

-- Trigger to update the updated_at field
CREATE OR REPLACE FUNCTION update_updated_at_column()
RETURNS TRIGGER AS $$
BEGIN
   NEW.updated_at = now();
   RETURN NEW;
END;
$$ language 'plpgsql';

CREATE TRIGGER update_project_reports_updated_at
BEFORE UPDATE ON project_reports
FOR EACH ROW EXECUTE FUNCTION update_updated_at_column();

-- Example indexes to optimize querying by project or report
CREATE INDEX idx_project_reports_project_id ON project_reports(project_id);
CREATE INDEX idx_project_reports_report_id ON project_reports(report_id);
