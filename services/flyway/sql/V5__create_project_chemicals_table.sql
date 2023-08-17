CREATE TABLE source_chemicals (
    source_chemical_id SERIAL PRIMARY KEY,
    source_id INTEGER REFERENCES projects(project_id),
    chemical_id INTEGER REFERENCES chemicals(chemical_id)
);
