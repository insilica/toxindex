-- Add archived and last_accessed columns to tasks
ALTER TABLE tasks ADD COLUMN archived BOOLEAN DEFAULT FALSE;
ALTER TABLE tasks ADD COLUMN last_accessed TIMESTAMP; 