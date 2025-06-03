ALTER TABLE tasks DROP CONSTRAINT tasks_environment_id_fkey;
ALTER TABLE tasks
ADD CONSTRAINT tasks_environment_id_fkey
FOREIGN KEY (environment_id) REFERENCES environments(environment_id) ON DELETE CASCADE; 