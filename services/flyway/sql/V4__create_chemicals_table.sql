CREATE TABLE chemicals (
    chemical_id SERIAL PRIMARY KEY,
    project_id INTEGER REFERENCES projects(project_id),
    name VARCHAR(255) NOT NULL,
    formula VARCHAR(255),
    known_properties TEXT, -- consider defining a more structured format
    predicted_properties TEXT, -- consider defining a more structured format
    upload_date TIMESTAMP DEFAULT current_timestamp
);
