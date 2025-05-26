ALTER TABLE tasks ADD COLUMN environment_id INTEGER REFERENCES environments(environment_id);
