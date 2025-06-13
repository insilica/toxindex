-- V2__cascade_delete_environments.sql
-- Update foreign keys to cascade deletes from environments

-- Drop old foreign key constraints if they exist
ALTER TABLE tasks DROP CONSTRAINT IF EXISTS tasks_environment_id_fkey;
ALTER TABLE files DROP CONSTRAINT IF EXISTS files_environment_id_fkey;

-- Add new foreign key constraints with ON DELETE CASCADE
ALTER TABLE tasks
  ADD CONSTRAINT tasks_environment_id_fkey
  FOREIGN KEY (environment_id) REFERENCES environments(environment_id) ON DELETE CASCADE;

ALTER TABLE files
  ADD CONSTRAINT files_environment_id_fkey
  FOREIGN KEY (environment_id) REFERENCES environments(environment_id) ON DELETE CASCADE; 